#include "correct_errors.hpp"



//----------------------------------------------------------------------
// determine_bloom_filter_parameters
//----------------------------------------------------------------------
void C_correct_errors::determine_bloom_filter_parameters(const C_arg& c_inst_args, C_time& c_inst_time) {
   if (rank_node == 0) {
      std::cout << "Determining Bloom filter parameters" << std::endl;
   }

   // find optimal parameters of for this Bloom filter
   // initialize variables
   double min_bit_vector_width_element = std::numeric_limits<double>::infinity();
   double min_num_hash_func = 1.0;
   double current_bit_vector_width_element = 0.0;

   // find the number of hash functions for minimum bit-vector size
   for (double j = 1.0; j < 1000.0; ++j) {
      if ((current_bit_vector_width_element = ((- j * num_unique_solid_kmers) / std::log(1.0 - std::pow(c_inst_args.target_false_positive_prob, 1.0 / j)))) < min_bit_vector_width_element) {
         min_bit_vector_width_element = current_bit_vector_width_element;
         min_num_hash_func = j;
      }
   }

   // find the tradeoff point between the bit vector size and the number of hash functions
   for (double j = 1.0; j < min_num_hash_func; ++j) {
      if ((current_bit_vector_width_element = ((- j * num_unique_solid_kmers) / std::log(1.0 - std::pow(c_inst_args.target_false_positive_prob, 1.0 / j)))) < (min_bit_vector_width_element * BIT_VEC_INC_RATIO)) {
         min_bit_vector_width_element = current_bit_vector_width_element;
         min_num_hash_func = j;
         break;
      }
   }

   // determine the size of the bloom filter
   num_hash_func         = static_cast<unsigned int>(min_num_hash_func);
   bit_vector_width      = static_cast<bloom_type>(min_bit_vector_width_element);
   bit_vector_width      += (((bit_vector_width % BITS_PER_CHAR) != 0) ? (BITS_PER_CHAR - (bit_vector_width % BITS_PER_CHAR)) : 0);
   bit_vector_width_byte = bit_vector_width / BITS_PER_CHAR;

   bloom_type bit_vector_size_mb(bit_vector_width_byte / (1024 * 1024));
   if (bit_vector_size_mb < 1) {
      bit_vector_size_mb = 1;
   }

   if (rank_node == 0) {
      std::cout << "     Bloom filter" << std::endl;
      std::cout << "     Number of keys           : " << num_unique_solid_kmers << std::endl;
      std::cout << "     Bit-vector size          : " << bit_vector_size_mb << " MB" << std::endl;
      std::cout << "     Number of hash functions : " << num_hash_func << std::endl;
      std::cout << "     Determining Bloom filter parameters: done" << std::endl;
      std::cout << std::endl;
   }
}



//----------------------------------------------------------------------
// program_kmers
//----------------------------------------------------------------------
void C_correct_errors::program_kmers(const C_arg& c_inst_args, C_time& c_inst_time) {
   time_t rawtime;
   time(&rawtime);
   c_inst_time.start_program_kmers_into_bloom_filter = asctime(localtime(&rawtime));

   if (rank_node == 0) {
      std::cout << "Programming k-mers to the Bloom filter" << std::endl;
   }

   CKMCFile kmer_database;

   //--------------------------------------------------
   // open the database
   //--------------------------------------------------
   // set rank_node_text
   std::stringstream rank_node_stream;
   rank_node_stream << std::setw(5) << std::setfill('0') << rank_node;
   rank_node_text = rank_node_stream.str();

   kmc_prefix = c_inst_args.kmc_prefix + "." + rank_node_text;

   // open the database
   if (!kmer_database.OpenForListing(kmc_prefix)) {
      std::cout << std::endl << "ERROR: Cannot open " << kmc_prefix << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 200);
   }

   //--------------------------------------------------
   // generate a temporary bit vector
   //--------------------------------------------------
   // generate unique hash seeds
   std::vector<unsigned int> hash_seed;
   generate_hash_seed(c_inst_args.random_seed, hash_seed);

   // generate a bloom filter
   std::vector<unsigned char> bit_vector(bit_vector_width_byte, 0);

   //--------------------------------------------------
   // iterate the database
   //--------------------------------------------------
   CKmerAPI kmer((uint32)kmer_length);

   float num_occurrences;

   std::vector<CKmerAPI> kmer_vector;

   std::size_t kmer_vector_index(0);

   std::string kmer_string;

   bloom_type original_index;
   bloom_type hash [2];

   unsigned short int bit_index;

   kmer_vector.resize(KMER_BLOCK_SIZE);

   // iterate k-mers
   while (kmer_database.ReadNextKmer(kmer, num_occurrences)) {
      if (num_occurrences >= kmer_occurrence_threshold) {
         // add the k-mer to the vector
         kmer_vector[kmer_vector_index] = kmer;

         kmer_vector_index++;

         // kmer_vector is full
         if (kmer_vector_index == KMER_BLOCK_SIZE) {
            //--------------------------------------------------
            // correct reads in a block
            //--------------------------------------------------
            #pragma omp parallel shared(bit_vector, hash_seed, kmer_vector) private(kmer_string, original_index, hash, bit_index)
            {
               // iterate reads
               #pragma omp for schedule(dynamic, OPENMP_CHUNK_SIZE)
               for (std::size_t it_kmer = 0; it_kmer < KMER_BLOCK_SIZE; it_kmer++) {
                  // convert the k-mer in the db to a string
                  (kmer_vector[it_kmer]).to_string(kmer_string);
                  kmer_string = kmer_string.substr(0, kmer_length);

                  // program the k-mer into the bloom filter
                  for (unsigned short int it_hash_func = 0; it_hash_func < num_hash_func; it_hash_func++) {
                     MurmurHash3_x64_128(kmer_string.c_str(), kmer_length, hash_seed[it_hash_func], hash);

                     original_index = hash[0] % bit_vector_width;
                     bit_index      = original_index % BITS_PER_CHAR;

                     #pragma omp atomic
                     bit_vector[original_index / BITS_PER_CHAR] |= BIT_MASK[bit_index];
                  }
               }
            }

            kmer_vector_index = 0;
         }
      }
   }

   // program remaining k-mers
   if (kmer_vector_index > 0) {
      #pragma omp parallel shared(bit_vector, hash_seed, kmer_vector, kmer_vector_index) private(kmer_string, original_index, hash, bit_index)
      {
         // iterate reads
         #pragma omp for schedule(dynamic, OPENMP_CHUNK_SIZE)
         for (std::size_t it_kmer = 0; it_kmer < kmer_vector_index; it_kmer++) {
            // convert the k-mer in the db to a string
            kmer_vector[it_kmer].to_string(kmer_string);
            kmer_string = kmer_string.substr(0, kmer_length);

            // program the k-mer into the bloom filter
            for (unsigned short int it_hash_func = 0; it_hash_func < num_hash_func; it_hash_func++) {
               MurmurHash3_x64_128(kmer_string.c_str(), kmer_length, hash_seed[it_hash_func], hash);

               original_index = hash[0] % bit_vector_width;
               bit_index      = original_index % BITS_PER_CHAR;

               #pragma omp atomic
               bit_vector[original_index / BITS_PER_CHAR] |= BIT_MASK[bit_index];
            }
         }
      }
   }

   kmer_vector.clear();

   //--------------------------------------------------
   // do bit-wise or for the bloom filter data
   //--------------------------------------------------
   std::size_t receive_buffer_size_byte;
   std::size_t num_iterations;
   std::size_t buffer_size_residue;

   // determine the size of the receive buffer
   if (bit_vector_width_byte >= BLOOM_RCV_BUFFER_SIZE) {
      receive_buffer_size_byte = BLOOM_RCV_BUFFER_SIZE;
   }
   else {
      receive_buffer_size_byte = bit_vector_width_byte;
   }

   num_iterations      = bit_vector_width_byte / receive_buffer_size_byte;
   buffer_size_residue = bit_vector_width_byte % receive_buffer_size_byte;

   unsigned char* receive_buffer(new unsigned char[static_cast<std::size_t>(receive_buffer_size_byte)]);

   //MPI_Barrier(comm_node);

   // bit-wise or
   for (std::size_t it_iter = 0; it_iter < num_iterations; it_iter++) {
      MPI_Allreduce(&bit_vector[it_iter * receive_buffer_size_byte], receive_buffer, receive_buffer_size_byte, MPI_UNSIGNED_CHAR, MPI_BOR, comm_node);
      std::memcpy(&bit_vector[it_iter * receive_buffer_size_byte], receive_buffer, receive_buffer_size_byte);
   }

   MPI_Allreduce(&bit_vector[num_iterations * receive_buffer_size_byte], receive_buffer, buffer_size_residue, MPI_UNSIGNED_CHAR, MPI_BOR, comm_node);
   std::memcpy(&bit_vector[num_iterations * receive_buffer_size_byte], receive_buffer, buffer_size_residue);

   boost::filesystem::path path_kmc_pre(kmc_prefix + ".kmc_pre");
   boost::filesystem::path path_kmc_suf(kmc_prefix + ".kmc_suf");
   boost::filesystem::remove(path_kmc_pre);
   boost::filesystem::remove(path_kmc_suf);

   //--------------------------------------------------
   // dump the bloom filter data
   //--------------------------------------------------
   if (rank_node == 0) {
      // open the bloom filter dump file
      std::ofstream f_bf_dump_size;
      std::ofstream f_bf_dump_data;
      f_bf_dump_size.open(c_inst_args.bf_size_file_name.c_str());
      f_bf_dump_data.open(c_inst_args.bf_data_file_name.c_str(), std::ios::binary);

      // check the bloom filter dump file
      if (f_bf_dump_size.is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.bf_size_file_name << std::endl << std::endl;
         MPI_Abort(MPI_COMM_WORLD, 201);
      }

      if (f_bf_dump_data.is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.bf_data_file_name << std::endl << std::endl;
         MPI_Abort(MPI_COMM_WORLD, 202);
      }

      // write the size file
      // order: <bit vector size in bits> <bit vector size in bytes> <# of hash functions> <k-mer occurrence threshold> <hash function random seed> <# of unique solid k-mers>
      f_bf_dump_size << bit_vector_width << " ";
      f_bf_dump_size << bit_vector_width_byte << " ";
      f_bf_dump_size << num_hash_func << " ";
      f_bf_dump_size << kmer_length << " ";
      f_bf_dump_size << c_inst_args.kmer_occurrence_threshold << " ";
      f_bf_dump_size << c_inst_args.random_seed << " ";
      f_bf_dump_size << num_unique_solid_kmers << std::endl;

      // dump the bit vector
      f_bf_dump_data.write((const char*)&bit_vector[0], bit_vector_width_byte);

      f_bf_dump_data.close();
      f_bf_dump_size.close();
   }

   // purge allocated memory
   bit_vector.clear();

   delete[] receive_buffer;

   time(&rawtime);
   c_inst_time.end_program_kmers_into_bloom_filter = asctime(localtime(&rawtime));

   if (rank_node == 0) {
      std::cout << "     Programming k-mers to the Bloom filter: done" << std::endl;
      std::cout << std::endl;
   }
}



//----------------------------------------------------------------------
// read_bloom_filter_parameters
//----------------------------------------------------------------------
void C_correct_errors::read_bloom_filter_parameters(const C_arg& c_inst_args, C_time& c_inst_time) {
   // open the bloom filter dump file
   std::ifstream f_bf_dump_size;
   f_bf_dump_size.open(c_inst_args.loaded_bf_size_file_name.c_str());

   // check the bloom filter dump file
   if (f_bf_dump_size.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.loaded_bf_size_file_name << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 203);
   }

   // read the size file
   // order: <bit vector size> <# of hash functions> <k-mer occurrence threshold> <# of unique solid k-mers>
   std::size_t kmer_length_read;

   f_bf_dump_size
                  >> bit_vector_width
                  >> bit_vector_width_byte
                  >> num_hash_func
                  >> kmer_length_read
                  >> kmer_occurrence_threshold
                  >> random_seed
                  >> num_unique_solid_kmers;

   if (f_bf_dump_size.good() == false) {
      std::cout << std::endl << "ERROR: Wrong number of items in " << c_inst_args.loaded_bf_size_file_name << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 204);
   }

   // compare the k-mer length between the loaded data and the command line input
   if (kmer_length_read != kmer_length) {
      std::cout << std::endl << "ERROR: k-mer length in the loaded data: " << kmer_length << " vs in the command line: " << kmer_length << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 205);
   }

   // set rank_node_text
   // this will be needed when errors are corrected
   std::stringstream rank_node_stream;
   rank_node_stream << std::setw(5) << std::setfill('0') << rank_node;
   rank_node_text = rank_node_stream.str();

   f_bf_dump_size.close();
}



//----------------------------------------------------------------------
// remove_new_line
//----------------------------------------------------------------------
std::string C_correct_errors::remove_new_line(std::string in_string) {
   in_string.erase(std::remove(in_string.begin(), in_string.end(), '\n'), in_string.end());
   return in_string;
}



//----------------------------------------------------------------------
// summarize_outputs
//----------------------------------------------------------------------
void C_correct_errors::summarize_outputs(const C_arg& c_inst_args, C_time& c_inst_time) {
   std::string last_ending_time;

   //--------------------------------------------------
   // stdout
   //--------------------------------------------------
   std::cout << "Running Time" << std::endl;
   std::cout << "     Checking reads" << std::endl;
   std::cout << "          Start: " << remove_new_line(c_inst_time.start_parse_args) << std::endl;
   std::cout << "          End  : " << remove_new_line(c_inst_time.end_check_read_file) << std::endl;

   if (c_inst_args.load_bf == false) {
      std::cout << "     Counting the number of unique solid k-mers" << std::endl;
      std::cout << "          Start: " << remove_new_line(c_inst_time.end_check_read_file) << std::endl;
      std::cout << "          End  : " << remove_new_line(c_inst_time.end_count_kmers) << std::endl;

      std::cout << "     Programming k-mers into the Bloom filter" << std::endl;
      std::cout << "          Start: " << remove_new_line(c_inst_time.end_count_kmers) << std::endl;
      std::cout << "          End  : " << remove_new_line(c_inst_time.end_program_kmers_into_bloom_filter) << std::endl;
      last_ending_time = c_inst_time.end_program_kmers_into_bloom_filter;
   }
   else {
      last_ending_time = c_inst_time.end_check_read_file;
   }

   std::cout << "     Correcting errors in reads" << std::endl;
   std::cout << "          Start: " << remove_new_line(last_ending_time) << std::endl;
   std::cout << "          End  : " << remove_new_line(c_inst_time.end_correct_errors_in_reads) << std::endl;

   std::cout << "     Merging output files" << std::endl;
   std::cout << "          Start: " << remove_new_line(c_inst_time.end_correct_errors_in_reads) << std::endl;
   std::cout << "          End  : " << remove_new_line(c_inst_time.end_merge_output_files) << std::endl;

   std::cout << std::endl;
   std::cout << "BLESS is successfully completed" << std::endl << std::endl;
}



//----------------------------------------------------------------------
// correct_errors_in_reads
//----------------------------------------------------------------------
void C_correct_errors::correct_errors_in_reads(const C_arg& c_inst_args, C_time& c_inst_time) {
   time_t rawtime;
   time(&rawtime);
   c_inst_time.start_correct_errors_in_reads = asctime(localtime(&rawtime));

   if (rank_node == 0) {
      std::cout << "Correcting errors in reads" << std::endl;
   }

   //----------------------------------------------------------------------
   // check the rank of this core and the number of cores
   //----------------------------------------------------------------------
   if (rank_node < 0) {
      std::cout << std::endl << "ERROR: The rank of this core is " << rank_node << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 206);
   }

   if (size_node < 0) {
      std::cout << std::endl << "ERROR: The number of cores is " << size_node << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 207);
   }

   //----------------------------------------------------------------------
   // correct reads
   //----------------------------------------------------------------------
   read_block_size = OPENMP_CHUNK_SIZE * READ_BLOCK_SIZE_RATIO * c_inst_args.smpthread;

   if (c_inst_args.paired_read == true) {
      if (c_inst_args.gzipped_input_read) {
         correct_errors_in_reads_paired_fastq_gzipped(c_inst_args);
      }
      else {
         correct_errors_in_reads_paired_fastq_unzipped(c_inst_args);
      }
   }
   else {
      if (c_inst_args.gzipped_input_read) {
         correct_errors_in_reads_single_fastq_gzipped(c_inst_args);
      }
      else {
         correct_errors_in_reads_single_fastq_unzipped(c_inst_args);
      }
   }

   //MPI_Barrier(comm_node);

   //----------------------------------------------------------------------
   // reduce all the numbers
   //----------------------------------------------------------------------
   std::size_t reduced_num_corrected_errors;
   std::size_t reduced_num_trimmed_bases;
   std::size_t reduced_num_corrected_reads;

   MPI_Reduce(&num_corrected_errors, &reduced_num_corrected_errors, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm_node);
   MPI_Reduce(&num_trimmed_bases,    &reduced_num_trimmed_bases,    1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm_node);
   MPI_Reduce(&num_corrected_reads,  &reduced_num_corrected_reads,  1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm_node);

   if (rank_node == 0) {
      std::cout << "     Number of corrected errors: " << reduced_num_corrected_errors << std::endl;
      std::cout << "     Number of trimmed bases   : " << reduced_num_trimmed_bases << std::endl;
      std::cout << "     Number of corrected reads : " << reduced_num_corrected_reads << std::endl;
      std::cout << "     Correcting errors in reads: done" << std::endl;
      std::cout << std::endl;
   }

   time(&rawtime);
   c_inst_time.end_correct_errors_in_reads = asctime(localtime(&rawtime));

   //----------------------------------------------------------------------
   // merge output files
   //----------------------------------------------------------------------
   time(&rawtime);
   c_inst_time.start_merge_output_files = asctime(localtime(&rawtime));

   if (rank_node == 0) {
      std::cout << "Merging output files" << std::endl;
   }

   //----------------------------------------------------------------------
   // multiple nodes
   //----------------------------------------------------------------------
   if (size_node > 1) {
      //----------------------------------------------------------------------
      // paired-end output file
      //----------------------------------------------------------------------
      if (c_inst_args.paired_read == true) {
         // gzip outputs
         if (c_inst_args.gzipped_output_read) {
            // only rank 0 does this
            // compressions cannot be parallelized
            if (rank_node == 0) {
               // variables
               std::string node_text;
               std::string tmp_file;

               std::stringstream rank_node_stream;

               std::ifstream f_in;

               std::ofstream f_out;

               boost::iostreams::filtering_streambuf<boost::iostreams::input> f_in_filter;

               //
               // forward file
               //
               // open output files
               f_out.open(c_inst_args.corrected_read_file_name1.c_str(), std::ofstream::binary);
               if (!f_out.is_open()) {
                  std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.corrected_read_file_name1 << std::endl << std::endl;
                  MPI_Abort(MPI_COMM_WORLD, 208);
               }

               for (int it_rank = 0; it_rank < size_node; it_rank++) {
                  // open a partial output file
                  rank_node_stream.str(std::string());
                  rank_node_stream.clear();
                  rank_node_stream << std::setw(5) << std::setfill('0') << it_rank;
                  node_text = rank_node_stream.str();
                  tmp_file  = c_inst_args.corrected_read_file_name1 + '.' + node_text;

                  f_in.open(tmp_file.c_str());
                  if (!f_out.is_open()) {
                     std::cout << std::endl << "ERROR: Cannot open " << tmp_file << std::endl << std::endl;
                     MPI_Abort(MPI_COMM_WORLD, 208);
                  }

                  f_in_filter.reset();
                  f_in_filter.push(boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(boost::iostreams::gzip::best_speed)));
                  f_in_filter.push(f_in);

                  // write compressed temporarily filed to the output
                  boost::iostreams::copy(f_in_filter, f_out);

                  f_in.close();

                  // remove the temporary file
                  boost::filesystem::path path_tmp1(tmp_file);
                  boost::filesystem::remove(path_tmp1);
               }

               //
               // reverse file
               //
               // open output files
               f_out.open(c_inst_args.corrected_read_file_name1.c_str(), std::ofstream::binary);
               if (!f_out.is_open()) {
                  std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.corrected_read_file_name1 << std::endl << std::endl;
                  MPI_Abort(MPI_COMM_WORLD, 208);
               }

               for (int it_rank = 0; it_rank < size_node; it_rank++) {
                  // open a partial output file
                  rank_node_stream.str(std::string());
                  rank_node_stream.clear();
                  rank_node_stream << std::setw(5) << std::setfill('0') << it_rank;
                  node_text = rank_node_stream.str();
                  tmp_file  = c_inst_args.corrected_read_file_name2 + '.' + node_text;

                  f_in.open(tmp_file.c_str());
                  if (!f_out.is_open()) {
                     std::cout << std::endl << "ERROR: Cannot open " << tmp_file << std::endl << std::endl;
                     MPI_Abort(MPI_COMM_WORLD, 208);
                  }

                  f_in_filter.reset();
                  f_in_filter.push(boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(boost::iostreams::gzip::best_speed)));
                  f_in_filter.push(f_in);

                  // write compressed temporarily filed to the output
                  boost::iostreams::copy(f_in_filter, f_out);

                  f_in.close();

                  // remove the temporary file
                  boost::filesystem::path path_tmp2(tmp_file);
                  boost::filesystem::remove(path_tmp2);
               }

               f_out.close();
            }
         }
         // do not gzip outputs
         else {
            // variables
            MPI_File f_out;

            MPI_Offset write_offset;

            MPI_Status status;

            int err;

            std::string node_text;
            std::string tmp_file;

            std::size_t total_size;
            std::size_t num_buffers;
            std::size_t residue;

            boost::iostreams::mapped_file_source mmap;

            std::stringstream rank_node_stream;

            //
            // forward file
            //
            // open output files
            err = MPI_File_open(comm_node, c_inst_args.corrected_read_file_name1.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &f_out);
            if (err != MPI_SUCCESS) {
               std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.corrected_read_file_name1 << std::endl << std::endl;
               MPI_Abort(MPI_COMM_WORLD, 208);
            }

            // calculate the output offset
            write_offset = 0;

            for (int it_nodes = 0; it_nodes < rank_node; it_nodes++) {
               // set temporary file names
               rank_node_stream.str(std::string());
               rank_node_stream.clear();
               rank_node_stream << std::setw(5) << std::setfill('0') << it_nodes;
               node_text = rank_node_stream.str();
               tmp_file  = c_inst_args.corrected_read_file_name1 + '.' + node_text;

               // add the size of each file
               mmap.open(tmp_file.c_str());
               if (!mmap.is_open()) {
                  std::cout << std::endl << "ERROR: Cannot open " << tmp_file << std::endl << std::endl;
               }
               write_offset += mmap.size();
               mmap.close();
            }

            // open a partial output file
            rank_node_stream.str(std::string());
            rank_node_stream.clear();
            rank_node_stream << std::setw(5) << std::setfill('0') << rank_node;
            node_text = rank_node_stream.str();
            tmp_file  = c_inst_args.corrected_read_file_name1 + '.' + node_text;

            // find out the total size of the file
            mmap.open(tmp_file.c_str());
            if (!mmap.is_open()) {
               std::cout << std::endl << "ERROR: Cannot open " << tmp_file << std::endl << std::endl;
            }
            total_size = mmap.size();
            mmap.close();

            num_buffers = total_size / MMAP_FILE_SIZE;
            residue     = total_size % MMAP_FILE_SIZE;

            // iterate the file
            for (std::size_t it_buffer = 0; it_buffer < num_buffers; it_buffer++) {
               // open(<file name>, <length>, <offset>)
               mmap.open(tmp_file.c_str(), MMAP_FILE_SIZE, it_buffer * MMAP_FILE_SIZE);
               if (!mmap.is_open()) {
                  std::cout << std::endl << "ERROR: Cannot open " << tmp_file << std::endl << std::endl;
               }
               MPI_File_write_at(f_out, write_offset, mmap.data(), MMAP_FILE_SIZE, MPI_CHAR, &status);
               mmap.close();

               write_offset += MMAP_FILE_SIZE;
            }

            // residue lines
            // open(<file name>, <length>, <offset>)
            mmap.open(tmp_file.c_str(), residue, num_buffers * MMAP_FILE_SIZE);
            if (!mmap.is_open()) {
               std::cout << std::endl << "ERROR: Cannot open " << tmp_file << std::endl << std::endl;
            }
            MPI_File_write_at(f_out, write_offset, mmap.data(), residue, MPI_CHAR, &status);
            mmap.close();

            MPI_File_close(&f_out);

            // remove the temporary file
            boost::filesystem::path path_tmp1(tmp_file);
            boost::filesystem::remove(path_tmp1);

            //
            // reverse file
            //
            // open output files
            err = MPI_File_open(comm_node, c_inst_args.corrected_read_file_name2.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &f_out);
            if (err != MPI_SUCCESS) {
               std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.corrected_read_file_name2 << std::endl << std::endl;
               MPI_Abort(MPI_COMM_WORLD, 209);
            }

            // calculate the output offset
            write_offset = 0;

            for (int it_nodes = 0; it_nodes < rank_node; it_nodes++) {
               // set temporary file names
               rank_node_stream.str(std::string());
               rank_node_stream.clear();
               rank_node_stream << std::setw(5) << std::setfill('0') << it_nodes;
               node_text = rank_node_stream.str();
               tmp_file  = c_inst_args.corrected_read_file_name2 + '.' + node_text;

               // add the size of each file
               mmap.open(tmp_file.c_str());
               if (!mmap.is_open()) {
                  std::cout << std::endl << "ERROR: Cannot open " << tmp_file << std::endl << std::endl;
               }
               write_offset += mmap.size();
               mmap.close();
            }

            // open a partial output file
            rank_node_stream.str(std::string());
            rank_node_stream.clear();
            rank_node_stream << std::setw(5) << std::setfill('0') << rank_node;
            node_text = rank_node_stream.str();
            tmp_file  = c_inst_args.corrected_read_file_name2 + '.' + node_text;

            // find out the total size of the file
            mmap.open(tmp_file.c_str());
            if (!mmap.is_open()) {
               std::cout << std::endl << "ERROR: Cannot open " << tmp_file << std::endl << std::endl;
            }
            total_size = mmap.size();
            mmap.close();

            num_buffers = total_size / MMAP_FILE_SIZE;
            residue     = total_size % MMAP_FILE_SIZE;

            // iterate the file
            for (std::size_t it_buffer = 0; it_buffer < num_buffers; it_buffer++) {
               // open(<file name>, <length>, <offset>)
               mmap.open(tmp_file.c_str(), MMAP_FILE_SIZE, it_buffer * MMAP_FILE_SIZE);
               if (!mmap.is_open()) {
                  std::cout << std::endl << "ERROR: Cannot open " << tmp_file << std::endl << std::endl;
               }
               MPI_File_write_at(f_out, write_offset, mmap.data(), MMAP_FILE_SIZE, MPI_CHAR, &status);
               mmap.close();

               write_offset += MMAP_FILE_SIZE;
            }

            // residue lines
            // open(<file name>, <length>, <offset>)
            mmap.open(tmp_file.c_str(), residue, num_buffers * MMAP_FILE_SIZE);
            if (!mmap.is_open()) {
               std::cout << std::endl << "ERROR: Cannot open " << tmp_file << std::endl << std::endl;
            }
            MPI_File_write_at(f_out, write_offset, mmap.data(), residue, MPI_CHAR, &status);
            mmap.close();

            MPI_File_close(&f_out);

            // remove the temporary file
            boost::filesystem::path path_tmp2(tmp_file);
            boost::filesystem::remove(path_tmp2);
         }
      }
      //----------------------------------------------------------------------
      // single-end output file
      //----------------------------------------------------------------------
      else {
         // gzip outputs
         if (c_inst_args.gzipped_output_read) {
            // only rank 0 does this
            // compressions cannot be parallelized
            if (rank_node == 0) {
               // variables
               std::string node_text;
               std::string tmp_file;

               std::stringstream rank_node_stream;

               std::ifstream f_in;

               std::ofstream f_out;

               boost::iostreams::filtering_streambuf<boost::iostreams::input> f_in_filter;

               // open output files
               f_out.open(c_inst_args.corrected_read_file_name.c_str(), std::ofstream::binary);
               if (!f_out.is_open()) {
                  std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.corrected_read_file_name << std::endl << std::endl;
                  MPI_Abort(MPI_COMM_WORLD, 208);
               }

               for (int it_rank = 0; it_rank < size_node; it_rank++) {
                  // open a partial output file
                  rank_node_stream.str(std::string());
                  rank_node_stream.clear();
                  rank_node_stream << std::setw(5) << std::setfill('0') << it_rank;
                  node_text = rank_node_stream.str();
                  tmp_file  = c_inst_args.corrected_read_file_name + '.' + node_text;

                  f_in.open(tmp_file.c_str());
                  if (!f_out.is_open()) {
                     std::cout << std::endl << "ERROR: Cannot open " << tmp_file << std::endl << std::endl;
                     MPI_Abort(MPI_COMM_WORLD, 208);
                  }

                  f_in_filter.reset();
                  f_in_filter.push(boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(boost::iostreams::gzip::best_speed)));
                  f_in_filter.push(f_in);

                  // write compressed temporarily filed to the output
                  boost::iostreams::copy(f_in_filter, f_out);

                  // remove the temporary file
                  boost::filesystem::path path_tmp(tmp_file);
                  boost::filesystem::remove(path_tmp);
               }
            }
         }
         // do not gzip outputs
         else {
            // open output files
            MPI_File f_out;

            int err;
            err = MPI_File_open(comm_node, c_inst_args.corrected_read_file_name.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &f_out);
            if (err != MPI_SUCCESS) {
               std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.corrected_read_file_name << std::endl << std::endl;
               MPI_Abort(MPI_COMM_WORLD, 210);
            }

            // calculate the output offset
            MPI_Offset write_offset(0);

            boost::iostreams::mapped_file_source mmap;

            std::stringstream rank_node_stream;

            std::string node_text;
            std::string tmp_file;

            for (int it_nodes = 0; it_nodes < rank_node; it_nodes++) {
               // set temporary file names
               rank_node_stream.str(std::string());
               rank_node_stream.clear();
               rank_node_stream << std::setw(5) << std::setfill('0') << it_nodes;
               node_text = rank_node_stream.str();
               tmp_file  = c_inst_args.corrected_read_file_name + '.' + node_text;

               // add the size of each file
               mmap.open(tmp_file.c_str());
               if (!mmap.is_open()) {
                  std::cout << std::endl << "ERROR: Cannot open " << tmp_file << std::endl << std::endl;
               }
               write_offset += mmap.size();
               mmap.close();
            }

            // open a partial output file
            MPI_Status status;

            rank_node_stream.str(std::string());
            rank_node_stream.clear();
            rank_node_stream << std::setw(5) << std::setfill('0') << rank_node;
            node_text = rank_node_stream.str();
            tmp_file  = c_inst_args.corrected_read_file_name + '.' + node_text;


            // find out the total size of the file
            mmap.open(tmp_file.c_str());
            if (!mmap.is_open()) {
               std::cout << std::endl << "ERROR: Cannot open " << tmp_file << std::endl << std::endl;
            }
            std::size_t total_size(mmap.size());
            mmap.close();

            std::size_t num_buffers(total_size / MMAP_FILE_SIZE);
            std::size_t residue(total_size % MMAP_FILE_SIZE);

            // iterate the file
            for (std::size_t it_buffer = 0; it_buffer < num_buffers; it_buffer++) {
               // open(<file name>, <length>, <offset>)
               mmap.open(tmp_file.c_str(), MMAP_FILE_SIZE, it_buffer * MMAP_FILE_SIZE);
               if (!mmap.is_open()) {
                  std::cout << std::endl << "ERROR: Cannot open " << tmp_file << std::endl << std::endl;
               }
               MPI_File_write_at(f_out, write_offset, mmap.data(), MMAP_FILE_SIZE, MPI_CHAR, &status);
               mmap.close();

               write_offset += MMAP_FILE_SIZE;
            }

            // residue lines
            // open(<file name>, <length>, <offset>)
            mmap.open(tmp_file.c_str(), residue, num_buffers * MMAP_FILE_SIZE);
            if (!mmap.is_open()) {
               std::cout << std::endl << "ERROR: Cannot open " << tmp_file << std::endl << std::endl;
            }
            MPI_File_write_at(f_out, write_offset, mmap.data(), residue, MPI_CHAR, &status);
            mmap.close();

            MPI_File_close(&f_out);

            // remove the temporary file
            boost::filesystem::path path_tmp(tmp_file);
            boost::filesystem::remove(path_tmp);
         }
      } // single-end output file
   }
   //----------------------------------------------------------------------
   // only one node
   //----------------------------------------------------------------------
   else {
      //----------------------------------------------------------------------
      // paired-end output file
      //----------------------------------------------------------------------
      if (c_inst_args.paired_read == true) {
         // gzip outputs
         if (c_inst_args.gzipped_output_read) {
            std::string file_name;

            std::ifstream f_in;

            std::ofstream f_out;

            boost::iostreams::filtering_streambuf<boost::iostreams::input> f_in_filter;

            // forward
            file_name = c_inst_args.corrected_read_file_name1 + '.' + rank_node_text;

            f_in.open(file_name.c_str());
            if (!f_in.is_open()) {
               std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.corrected_read_file_name1 << std::endl << std::endl;
               MPI_Abort(MPI_COMM_WORLD, 208);
            }

            f_out.open(c_inst_args.corrected_read_file_name1.c_str(), std::ofstream::binary);
            if (!f_out.is_open()) {
               std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.corrected_read_file_name1 << std::endl << std::endl;
               MPI_Abort(MPI_COMM_WORLD, 208);
            }

            f_in_filter.reset();
            f_in_filter.push(boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(boost::iostreams::gzip::best_speed)));
            f_in_filter.push(f_in);
            boost::iostreams::copy(f_in_filter, f_out);

            f_in.close();
            f_out.close();

            boost::filesystem::path path_source1(c_inst_args.corrected_read_file_name1 + '.' + rank_node_text);
            boost::filesystem::remove(path_source1);

            // reverse
            file_name = c_inst_args.corrected_read_file_name2 + '.' + rank_node_text;

            f_in.open(file_name.c_str());
            if (!f_in.is_open()) {
               std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.corrected_read_file_name2 << std::endl << std::endl;
               MPI_Abort(MPI_COMM_WORLD, 208);
            }

            f_out.open(c_inst_args.corrected_read_file_name2.c_str(), std::ofstream::binary);
            if (!f_out.is_open()) {
               std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.corrected_read_file_name2 << std::endl << std::endl;
               MPI_Abort(MPI_COMM_WORLD, 208);
            }

            f_in_filter.reset();
            f_in_filter.push(boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(boost::iostreams::gzip::best_speed)));
            f_in_filter.push(f_in);
            boost::iostreams::copy(f_in_filter, f_out);

            f_in.close();
            f_out.close();

            boost::filesystem::path path_source2(c_inst_args.corrected_read_file_name2 + '.' + rank_node_text);
            boost::filesystem::remove(path_source2);
         }
         // do not gzip outputs
         else {
            boost::filesystem::path path_source1(c_inst_args.corrected_read_file_name1 + '.' + rank_node_text);
            boost::filesystem::path path_source2(c_inst_args.corrected_read_file_name2 + '.' + rank_node_text);
            boost::filesystem::path path_destination1(c_inst_args.corrected_read_file_name1);
            boost::filesystem::path path_destination2(c_inst_args.corrected_read_file_name2);

            if (boost::filesystem::exists(path_destination1)) {
               boost::filesystem::remove(path_destination1);
            }
            if (boost::filesystem::exists(path_destination2)) {
               boost::filesystem::remove(path_destination2);
            }

            boost::filesystem::rename(path_source1, path_destination1);
            boost::filesystem::rename(path_source2, path_destination2);
         }
      }
      //----------------------------------------------------------------------
      // single-end output file
      //----------------------------------------------------------------------
      else {
         // gzip outputs
         if (c_inst_args.gzipped_output_read) {
            std::string file_name;

            std::ifstream f_in;

            std::ofstream f_out;

            boost::iostreams::filtering_streambuf<boost::iostreams::input> f_in_filter;

            file_name = c_inst_args.corrected_read_file_name + '.' + rank_node_text;

            f_in.open(file_name.c_str());
            if (!f_in.is_open()) {
               std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.corrected_read_file_name << std::endl << std::endl;
               MPI_Abort(MPI_COMM_WORLD, 208);
            }

            f_out.open(c_inst_args.corrected_read_file_name.c_str(), std::ofstream::binary);
            if (!f_out.is_open()) {
               std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.corrected_read_file_name << std::endl << std::endl;
               MPI_Abort(MPI_COMM_WORLD, 208);
            }

            f_in_filter.reset();
            f_in_filter.push(boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(boost::iostreams::gzip::best_speed)));
            f_in_filter.push(f_in);
            boost::iostreams::copy(f_in_filter, f_out);

            f_in.close();
            f_out.close();

            boost::filesystem::path path_source(c_inst_args.corrected_read_file_name + '.' + rank_node_text);
            boost::filesystem::remove(path_source);
         }
         // do not gzip outputs
         else {
            boost::filesystem::path path_source(c_inst_args.corrected_read_file_name + '.' + rank_node_text);
            boost::filesystem::path path_destination(c_inst_args.corrected_read_file_name);

            if (boost::filesystem::exists(path_destination)) {
               boost::filesystem::remove(path_destination);
            }

            boost::filesystem::rename(path_source, path_destination);
         }
      }
   }

   if (rank_node == 0) {
      std::cout << "     Merging output files: done" << std::endl;
      std::cout << std::endl;
   }

   //MPI_Barrier(comm_node);

   time(&rawtime);
   c_inst_time.end_merge_output_files = asctime(localtime(&rawtime));
}



//----------------------------------------------------------------------
// correct_errors_in_reads_single_fastq_unzipped
//----------------------------------------------------------------------
void C_correct_errors::correct_errors_in_reads_single_fastq_unzipped(const C_arg& c_inst_args) {
   //----------------------------------------------------------------------
   // load bloom filter data
   //----------------------------------------------------------------------
   // generate unique hash seeds
   std::vector<unsigned int> hash_seed;
   // no bloom filter data load
   // use a command line random seed
   if (c_inst_args.load_bf == false) {
      generate_hash_seed(c_inst_args.random_seed, hash_seed);
   }
   // load bloom filter data
   // use a loaded random seed
   else {
      generate_hash_seed(random_seed, hash_seed);
   }

   // generate a bloom filter
   unsigned char* bit_vector(new unsigned char[static_cast<std::size_t>(bit_vector_width_byte)]);
   std::fill_n(bit_vector, bit_vector_width_byte, 0x00);

   // open the bloom filter data file
   std::string bf_data_file_name;
   std::string bf_size_file_name;
   if (c_inst_args.load_bf == false) {
      bf_data_file_name = c_inst_args.bf_data_file_name;
      bf_size_file_name = c_inst_args.bf_size_file_name;
   }
   else {
      bf_data_file_name = c_inst_args.loaded_bf_data_file_name;
      bf_size_file_name = c_inst_args.loaded_bf_size_file_name;
   }

   std::ifstream f_bf_dump_data;
   f_bf_dump_data.open(bf_data_file_name.c_str(), std::ios::binary);
   if (f_bf_dump_data.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << bf_data_file_name << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 211);
   }

   // load bit vector data
   f_bf_dump_data.read(reinterpret_cast<char*>(bit_vector), bit_vector_width_byte);

   f_bf_dump_data.close();

   if (rank_node == 0) {
      // calculate the size of the bloom filter in megabyte
      bloom_type bit_vector_size_mb(bit_vector_width_byte / (1024 * 1024));

      if (c_inst_args.load_bf == true) {
         std::cout << "     Bloom filter data file    : " << c_inst_args.bf_data_file_name << std::endl;
         std::cout << "     Bloom filter size file    : " << c_inst_args.bf_size_file_name << std::endl;
         std::cout << "     Number of keys            : " << num_unique_solid_kmers << std::endl;
         std::cout << "     Bit-vector size           : " << bit_vector_size_mb << " MB" << std::endl;
         std::cout << "     Number of hash functions  : " << num_hash_func << std::endl;
         std::cout << "     k-mer threshold           : " << kmer_occurrence_threshold << std::endl;
      }
   }

   //----------------------------------------------------------------------
   // correct errors
   //----------------------------------------------------------------------
   //
   // variables
   //
   boost::iostreams::mapped_file_source f_read;

   std::vector<std::string> read_vector;

   std::string sequence_modification;
   std::string trimmed_seq;
   std::string trimmed_qs;

   std::size_t read_length;
   std::size_t min_check_length;
   std::size_t max_allowed_ns;
   std::size_t max_trimmed_bases;
   std::size_t num_corrected_errors_local;
   std::size_t trim_5_end;
   std::size_t trim_3_end;
   std::size_t num_corrected_errors_tmp;
   std::size_t num_corrected_reads_tmp;
   std::size_t num_trimmed_bases_tmp;
   std::size_t read_vector_index;
   std::size_t current_read_index;
   std::size_t current_read_index_write;

   bool keep_going;
   bool all_reads;
   bool too_many_errors;

   short int line_order;

   long long int processed_bytes;
   long long int alignment_offset;
   long long int remaining_bytes;
   long long int prev_remaining_bytes;
   long long int this_iteration_bytes;

   const char* pt_current;
   const char* pt_last;
   const char* pt_current_read_start;
   const char* pt_current_line_start;

   std::size_t num_reads_local;

   std::regex regex_non_acgtn("[^ACGTN]");

   read_vector.resize(read_block_size * 3);

   num_corrected_errors_tmp = 0;
   num_corrected_reads_tmp  = 0;
   num_trimmed_bases_tmp    = 0;

   //
   // file specific
   //
   // open a temporary output file
   std::string corrected_read_file_name(c_inst_args.corrected_read_file_name + '.' + rank_node_text);

   FILE* f_corrected_read(fopen(corrected_read_file_name.c_str(), "w"));
   if (f_corrected_read == NULL) {
      std::cout << std::endl << "ERROR: Cannot open " << corrected_read_file_name << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 212);
   }

   // initialize variables
   num_reads_local          = 0;
   keep_going               = true;
   all_reads                = false;
   processed_bytes          = (starting_point_vector[rank_node] / read_file_unit_size_byte) * read_file_unit_size_byte;
   alignment_offset         = starting_point_vector[rank_node] % read_file_unit_size_byte;
   read_vector_index        = 0;
   current_read_index_write = 0;
   remaining_bytes          = read_file_size_byte - starting_point_vector[rank_node] + 1;
   prev_remaining_bytes     = read_file_size_byte - starting_point_vector[rank_node] + 1;

   // check whether the input file is empty
   if (remaining_bytes <= 0) {
      std::cout << std::endl << "ERROR: Remaining bytes are " << remaining_bytes << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 213);
   }

   // read loop
   while (keep_going) {
      // calculate the number of bytes that will be processed in this iteration
      if (remaining_bytes <= MMAP_FILE_SIZE) {
         this_iteration_bytes = remaining_bytes;
      }
      else {
         this_iteration_bytes = MMAP_FILE_SIZE;
      }

      // open the file
      // open(<file name>, <length>, <offset>)
      f_read.open(c_inst_args.read_file_name.c_str(), this_iteration_bytes, processed_bytes);

      if (!f_read.is_open()) {
         std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.read_file_name << std::endl << std::endl;
         MPI_Abort(MPI_COMM_WORLD, 214);
      }

      // initialize variables
      pt_current            = f_read.data() + alignment_offset;
      pt_last               = f_read.data() + f_read.size();
      pt_current_read_start = pt_current;
      pt_current_line_start = pt_current;

      processed_bytes += alignment_offset;

      line_order = 0;

      while (pt_current && (pt_current != pt_last) && (!all_reads)) {
         // find a new line
         if ((pt_current = static_cast<const char*>(memchr(pt_current, '\n', pt_last - pt_current)))) {
            // header
            if (line_order == 0) {
               current_read_index_write = read_vector_index * 3;

               // include '\n'
               read_vector[current_read_index_write].assign(pt_current_line_start, pt_current - pt_current_line_start + 1);

               line_order++;
            }
            // sequence
            else if (line_order == 1) {
               read_vector[current_read_index_write + 1].assign(pt_current_line_start, pt_current - pt_current_line_start);

               line_order++;
            }
            // connector
            else if (line_order == 2) {
               line_order++;
            }
            // quality score
            else {
               read_vector[current_read_index_write + 2].assign(pt_current_line_start, pt_current - pt_current_line_start);

               line_order = 0;
               processed_bytes += (pt_current - pt_current_read_start + 1);
               pt_current_read_start = pt_current + 1;
               read_vector_index++;

               num_reads_local++;

               // all reads have been processed
               // stop the process
               if (num_reads_local == num_reads_vector[rank_node]) {
                  all_reads  = true;
                  keep_going = false;
               }

               // read_vector is full
               if (read_vector_index == read_block_size) {
                  //--------------------------------------------------
                  // correct reads in a block
                  //--------------------------------------------------
                  #pragma omp parallel shared(bit_vector, hash_seed, read_vector) private(current_read_index, read_length, min_check_length, max_allowed_ns, max_trimmed_bases, num_corrected_errors_local, trim_5_end, trim_3_end, sequence_modification, too_many_errors) reduction(+: num_corrected_errors_tmp, num_corrected_reads_tmp, num_trimmed_bases_tmp)
                  {
                     // iterate reads
                     #pragma omp for schedule(dynamic, OPENMP_CHUNK_SIZE)
                     for (std::size_t it_read = 0; it_read < read_block_size; it_read++) {
                        // calculate the current index
                        current_read_index = it_read * 3;

                        // change sequences to upper case
                        std::transform(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), read_vector[current_read_index + 1].begin(), static_cast<int(*)(int)>(std::toupper));

                        // substitute non-standard characters with Ns
                        read_vector[current_read_index + 1] = std::regex_replace(read_vector[current_read_index + 1], regex_non_acgtn, "N");

                        // set various thresholds
                        read_length      = read_vector[current_read_index + 1].length();
                        min_check_length = read_length * CHECK_RANGE_RATIO;
                        max_allowed_ns   = read_length * MAX_N_RATIO;

                        // forward
                        // too short read: no trimming
                        if (read_length <= MIN_BASES_AFTER_TRIMMING) {
                           max_trimmed_bases = 0;
                        }
                        else {
                           max_trimmed_bases = std::min((std::size_t)(read_length * MAX_TRIMMING_RATE), read_length - MIN_BASES_AFTER_TRIMMING);
                        }

                        // check the number of Ns in the read
                        if (((std::size_t)std::count(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), 'N') <= max_allowed_ns) &&
                            (read_length >= kmer_length)) {
                           // substitute Ns other characters
                           std::replace(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), 'N', SUBST_CHAR);

                           // initialize variables
                           num_corrected_errors_local = 0;
                           trim_5_end                 = 0;
                           trim_3_end                 = 0;

                           // storage for the modification of the reads
                           // # of entries = read length
                           sequence_modification.assign(read_length, '0');

                           //----------------------------------------------------------------------
                           // correct errors in a read
                           //----------------------------------------------------------------------
                           // forward read
                           // errors cannot be corrected if k is equal to read length
                           if (read_length > kmer_length) {
                              correct_errors_in_a_read_fastq(
                                                             read_vector[current_read_index + 1],
                                                             sequence_modification,
                                                             read_vector[current_read_index + 2],
                                                             trim_5_end,
                                                             trim_3_end,
                                                             read_length,
                                                             max_trimmed_bases,
                                                             min_check_length,
                                                             num_corrected_errors_local,
                                                             bit_vector,
                                                             hash_seed
                                                            );
                           }

                           // no trim
                           if (c_inst_args.notrim == true) {
                              trim_5_end = 0;
                              trim_3_end = 0;
                           }
                           // adjust the number of trimmed bases
                           else {
                              if ((trim_5_end + trim_3_end) > max_trimmed_bases) {
                                 if (trim_3_end <= max_trimmed_bases) {
                                    trim_5_end = 0;
                                 }
                                 else if (trim_5_end <= max_trimmed_bases) {
                                    trim_3_end = 0;
                                 }
                                 else {
                                    trim_5_end = 0;
                                    trim_3_end = max_trimmed_bases;
                                 }
                              }
                           }

                           num_trimmed_bases_tmp += (trim_5_end + trim_3_end);

                           // update num_corrected_reads
                           too_many_errors = false;
                           if (num_corrected_errors_local > (read_vector[current_read_index + 1].length() * MAX_ERROR_RATE)) {
                              too_many_errors = true;
                           }
                           else if (num_corrected_errors_local > 0) {
                              num_corrected_errors_tmp += num_corrected_errors_local;

                              num_corrected_reads_tmp++;
                           }
                           else if (c_inst_args.notrim == false) {
                              if ((trim_5_end > 0) || (trim_3_end > 0)) {
                                 num_corrected_reads_tmp++;
                              }
                           }

                           // make a corrected read
                           if (too_many_errors == false) {
                              // apply modifications to the read
                              for (unsigned int it_base = trim_5_end; it_base < (read_vector[current_read_index + 1].length() - trim_3_end); it_base++) {
                                 if (sequence_modification[it_base] != '0') {
                                    read_vector[current_read_index + 1][it_base] = sequence_modification[it_base];
                                 }
                              }
                           }

                           // make a trimmed read
                           if ((trim_5_end + trim_3_end) > 0) {
                              read_vector[current_read_index + 1] = read_vector[current_read_index + 1].substr(trim_5_end, read_vector[current_read_index + 1].length() - trim_5_end - trim_3_end) + '\n';
                              read_vector[current_read_index + 2] = read_vector[current_read_index + 2].substr(trim_5_end, read_vector[current_read_index + 2].length() - trim_5_end - trim_3_end) + '\n';
                           }
                           else {
                              read_vector[current_read_index + 1] = read_vector[current_read_index + 1] + '\n';
                              read_vector[current_read_index + 2] = read_vector[current_read_index + 2] + '\n';
                           }
                        }
                        // too many Ns
                        // write reads without modification
                        else {
                           // sequence
                           read_vector[current_read_index + 1] += '\n';

                           // quality score
                           read_vector[current_read_index + 2] += '\n';
                        }
                     // it_read
                     }
                  // omp parallel
                  }

                  // write corrected_reads
                  for (std::size_t it_write = 0; it_write < read_block_size; it_write++) {
                     current_read_index = it_write * 3;

                     // header
                     fwrite(read_vector[current_read_index].c_str(), 1, read_vector[current_read_index].length(), f_corrected_read);
                     // sequence
                     fwrite(read_vector[current_read_index + 1].c_str(), 1, read_vector[current_read_index + 1].length(), f_corrected_read);
                     // connector
                     fwrite("+\n", 1, 2, f_corrected_read);
                     // quality score
                     fwrite(read_vector[current_read_index + 2].c_str(), 1, read_vector[current_read_index + 2].length(), f_corrected_read);
                  }

                  read_vector_index = 0;
               }
            }

            pt_current++;
            pt_current_line_start = pt_current;
         }
      }

      remaining_bytes  = read_file_size_byte - processed_bytes;

      if (remaining_bytes <= 0) {
         keep_going = false;
      }
      else if (remaining_bytes == prev_remaining_bytes) {
         std::cout << std::endl << "ERROR: The number of lines in " << c_inst_args.read_file_name << " is wrong" << std::endl << std::endl;
         MPI_Abort(MPI_COMM_WORLD, 215);
      }
      else {
         // align processed bytes to the multiple of unit_alignment_offset
         alignment_offset     = processed_bytes % read_file_unit_size_byte;
         processed_bytes     -= alignment_offset;
         prev_remaining_bytes = remaining_bytes;
         remaining_bytes     += alignment_offset;
      }

      f_read.close();
   }

   // correct errors in remaining reads
   if (read_vector_index > 0) {
      #pragma omp parallel shared(bit_vector, hash_seed, read_vector) private(current_read_index, read_length, min_check_length, max_allowed_ns, max_trimmed_bases, num_corrected_errors_local, trim_5_end, trim_3_end, sequence_modification, too_many_errors) reduction(+: num_corrected_errors_tmp, num_corrected_reads_tmp, num_trimmed_bases_tmp)
      {
         // iterate reads
         #pragma omp for schedule(dynamic, OPENMP_CHUNK_SIZE)
         for (std::size_t it_read = 0; it_read < read_vector_index; it_read++) {
            // calculate the current index
            current_read_index = it_read * 3;

            // change sequences to upper case
            std::transform(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), read_vector[current_read_index + 1].begin(), static_cast<int(*)(int)>(std::toupper));

            // substitute non-standard characters with Ns
            read_vector[current_read_index + 1] = std::regex_replace(read_vector[current_read_index + 1], regex_non_acgtn, "N");

            // set various thresholds
            read_length      = read_vector[current_read_index + 1].length();
            min_check_length = read_length * CHECK_RANGE_RATIO;
            max_allowed_ns   = read_length * MAX_N_RATIO;

            // forward
            // too short read: no trimming
            if (read_length <= MIN_BASES_AFTER_TRIMMING) {
               max_trimmed_bases = 0;
            }
            else {
               max_trimmed_bases = std::min((std::size_t)(read_length * MAX_TRIMMING_RATE), read_length - MIN_BASES_AFTER_TRIMMING);
            }

            // check the number of Ns in the read
            if (((std::size_t)std::count(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), 'N') <= max_allowed_ns) &&
                (read_length >= kmer_length)) {
               // substitute Ns other characters
               std::replace(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), 'N', SUBST_CHAR);

               // initialize variables
               num_corrected_errors_local = 0;
               trim_5_end                 = 0;
               trim_3_end                 = 0;

               // storage for the modification of the reads
               // # of entries = read length
               sequence_modification.assign(read_length, '0');

               //----------------------------------------------------------------------
               // correct errors in a read
               //----------------------------------------------------------------------
               // forward read
               // errors cannot be corrected if k is equal to read length
               if (read_length > kmer_length) {
                  correct_errors_in_a_read_fastq(
                                                 read_vector[current_read_index + 1],
                                                 sequence_modification,
                                                 read_vector[current_read_index + 2],
                                                 trim_5_end,
                                                 trim_3_end,
                                                 read_length,
                                                 max_trimmed_bases,
                                                 min_check_length,
                                                 num_corrected_errors_local,
                                                 bit_vector,
                                                 hash_seed
                                                );
               }

               // no trim
               if (c_inst_args.notrim == true) {
                  trim_5_end = 0;
                  trim_3_end = 0;
               }
               // adjust the number of trimmed bases
               else {
                  if ((trim_5_end + trim_3_end) > max_trimmed_bases) {
                     if (trim_3_end <= max_trimmed_bases) {
                        trim_5_end = 0;
                     }
                     else if (trim_5_end <= max_trimmed_bases) {
                        trim_3_end = 0;
                     }
                     else {
                        trim_5_end = 0;
                        trim_3_end = max_trimmed_bases;
                     }
                  }
               }

               num_trimmed_bases_tmp += (trim_5_end + trim_3_end);

               // update num_corrected_reads
               too_many_errors = false;
               if (num_corrected_errors_local > (read_vector[current_read_index + 1].length() * MAX_ERROR_RATE)) {
                  too_many_errors = true;
               }
               else if (num_corrected_errors_local > 0) {
                  num_corrected_errors_tmp += num_corrected_errors_local;

                  num_corrected_reads_tmp++;
               }
               else if (c_inst_args.notrim == false) {
                  if ((trim_5_end > 0) || (trim_3_end > 0)) {
                     num_corrected_reads_tmp++;
                  }
               }

               // make a corrected read
               if (too_many_errors == false) {
                  // apply modifications to the read
                  for (unsigned int it_base = trim_5_end; it_base < (read_vector[current_read_index + 1].length() - trim_3_end); it_base++) {
                     if (sequence_modification[it_base] != '0') {
                        read_vector[current_read_index + 1][it_base] = sequence_modification[it_base];
                     }
                  }
               }

               // make a trimmed read
               if ((trim_5_end + trim_3_end) > 0) {
                  read_vector[current_read_index + 1] = read_vector[current_read_index + 1].substr(trim_5_end, read_vector[current_read_index + 1].length() - trim_5_end - trim_3_end) + '\n';
                  read_vector[current_read_index + 2] = read_vector[current_read_index + 2].substr(trim_5_end, read_vector[current_read_index + 2].length() - trim_5_end - trim_3_end) + '\n';
               }
               else {
                  read_vector[current_read_index + 1] = read_vector[current_read_index + 1] + '\n';
                  read_vector[current_read_index + 2] = read_vector[current_read_index + 2] + '\n';
               }
            }
            // too many Ns
            // write reads without modification
            else {
               // sequence
               read_vector[current_read_index + 1] += '\n';

               // quality score
               read_vector[current_read_index + 2] += '\n';
            }
         // it_read
         }

      // omp parallel
      }

      // write corrected_reads
      for (std::size_t it_read = 0; it_read < read_vector_index; it_read++) {
         current_read_index = it_read * 3;

         // header
         fwrite(read_vector[current_read_index].c_str(), 1, read_vector[current_read_index].length(), f_corrected_read);
         // sequence
         fwrite(read_vector[current_read_index + 1].c_str(), 1, read_vector[current_read_index + 1].length(), f_corrected_read);
         // connector
         fwrite("+\n", 1, 2, f_corrected_read);
         // quality score
         fwrite(read_vector[current_read_index + 2].c_str(), 1, read_vector[current_read_index + 2].length(), f_corrected_read);
      }
   }

   fclose(f_corrected_read);

   read_vector.clear();

   num_corrected_errors = num_corrected_errors_tmp;
   num_corrected_reads  = num_corrected_reads_tmp;
   num_trimmed_bases    = num_trimmed_bases_tmp;
}



//----------------------------------------------------------------------
// correct_errors_in_reads_paired_fastq_unzipped
//----------------------------------------------------------------------
void C_correct_errors::correct_errors_in_reads_paired_fastq_unzipped(const C_arg& c_inst_args) {
   //----------------------------------------------------------------------
   // load bloom filter data
   //----------------------------------------------------------------------
   // generate unique hash seeds
   std::vector<unsigned int> hash_seed;
   // no bloom filter data load
   // use a command line random seed
   if (c_inst_args.load_bf == false) {
      generate_hash_seed(c_inst_args.random_seed, hash_seed);
   }
   // load bloom filter data
   // use a loaded random seed
   else {
      generate_hash_seed(random_seed, hash_seed);
   }

   // generate a bloom filter
   unsigned char* bit_vector(new unsigned char[static_cast<std::size_t>(bit_vector_width_byte)]);
   std::fill_n(bit_vector, bit_vector_width_byte, 0x00);

   // open the bloom filter data file
   std::string bf_data_file_name;
   std::string bf_size_file_name;
   if (c_inst_args.load_bf == false) {
      bf_data_file_name = c_inst_args.bf_data_file_name;
      bf_size_file_name = c_inst_args.bf_size_file_name;
   }
   else {
      bf_data_file_name = c_inst_args.loaded_bf_data_file_name;
      bf_size_file_name = c_inst_args.loaded_bf_size_file_name;
   }

   std::ifstream f_bf_dump_data;
   f_bf_dump_data.open(bf_data_file_name.c_str(), std::ios::binary);
   if (f_bf_dump_data.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << bf_data_file_name << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 216);
   }

   // load bit vector data
   f_bf_dump_data.read(reinterpret_cast<char*>(bit_vector), bit_vector_width_byte);

   f_bf_dump_data.close();

   if (rank_node == 0) {
      // calculate the size of the bloom filter in megabyte
      bloom_type bit_vector_size_mb(bit_vector_width_byte / (1024 * 1024));

      if (c_inst_args.load_bf == true) {
         std::cout << "     Bloom filter data file    : " << c_inst_args.bf_data_file_name << std::endl;
         std::cout << "     Bloom filter size file    : " << c_inst_args.bf_size_file_name << std::endl;
         std::cout << "     Number of keys            : " << num_unique_solid_kmers << std::endl;
         std::cout << "     Bit-vector size           : " << bit_vector_size_mb << " MB" << std::endl;
         std::cout << "     Number of hash functions  : " << num_hash_func << std::endl;
         std::cout << "     k-mer threshold           : " << kmer_occurrence_threshold << std::endl;
      }
   }

   //--------------------------------------------------
   // correct errors
   //--------------------------------------------------
   //
   // variables
   //
   boost::iostreams::mapped_file_source f_read;

   std::vector<std::string> read_vector;

   std::string sequence_modification;
   std::string trimmed_seq;
   std::string trimmed_qs;

   std::size_t read_length;
   std::size_t min_check_length;
   std::size_t max_allowed_ns;
   std::size_t max_trimmed_bases;
   std::size_t num_corrected_errors_local;
   std::size_t trim_5_end;
   std::size_t trim_3_end;
   std::size_t num_corrected_errors_tmp;
   std::size_t num_corrected_reads_tmp;
   std::size_t num_trimmed_bases_tmp;
   //std::size_t num_read_blocks;
   std::size_t read_vector_index;
   std::size_t current_read_index;
   std::size_t current_read_index_write;

   bool keep_going;
   bool all_reads;
   bool too_many_errors;

   short int line_order;

   long long int processed_bytes;
   long long int alignment_offset;
   long long int remaining_bytes;
   long long int prev_remaining_bytes;
   long long int this_iteration_bytes;

   const char* pt_current;
   const char* pt_last;
   const char* pt_current_read_start;
   const char* pt_current_line_start;

   std::size_t num_reads_local;

   std::regex regex_non_acgtn("[^ACGTN]");

   read_vector.resize(read_block_size * 3);

   num_corrected_errors_tmp = 0;
   num_corrected_reads_tmp  = 0;
   num_trimmed_bases_tmp    = 0;

   //
   // forward
   //
   // open a temporary output file
   std::string corrected_read_file_name1(c_inst_args.corrected_read_file_name1 + '.' + rank_node_text);

   FILE* f_corrected_read1(fopen(corrected_read_file_name1.c_str(), "w"));
   if (f_corrected_read1 == NULL) {
      std::cout << std::endl << "ERROR: Cannot open " << corrected_read_file_name1 << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 217);
   }

   // initialize variables
   num_reads_local          = 0;
   keep_going               = true;
   all_reads                = false;
   processed_bytes          = (starting_point_vector1[rank_node] / read_file_unit_size_byte1) * read_file_unit_size_byte1;
   alignment_offset         = starting_point_vector1[rank_node] % read_file_unit_size_byte1;
   read_vector_index        = 0;
   current_read_index_write = 0;
   remaining_bytes          = read_file_size_byte1 - starting_point_vector1[rank_node] + 1;
   prev_remaining_bytes     = read_file_size_byte1 - starting_point_vector1[rank_node] + 1;

   // check whether the input file is empty
   if (remaining_bytes <= 0) {
      std::cout << std::endl << "ERROR: Remaining bytes are " << remaining_bytes << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 218);
   }

   // read loop
   while (keep_going) {
      // calculate the number of bytes that will be processed in this iteration
      if (remaining_bytes <= MMAP_FILE_SIZE) {
         this_iteration_bytes = remaining_bytes;
      }
      else {
         this_iteration_bytes = MMAP_FILE_SIZE;
      }

      // open the file
      // open(<file name>, <length>, <offset>)
      f_read.open(c_inst_args.read_file_name1.c_str(), this_iteration_bytes, processed_bytes);

      if (!f_read.is_open()) {
         std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.read_file_name1 << std::endl << std::endl;
         MPI_Abort(MPI_COMM_WORLD, 219);
      }

      // initialize variables
      pt_current            = f_read.data() + alignment_offset;
      pt_last               = f_read.data() + f_read.size();
      pt_current_read_start = pt_current;
      pt_current_line_start = pt_current;

      processed_bytes += alignment_offset;

      line_order = 0;

      while (pt_current && (pt_current != pt_last) && (!all_reads)) {
         // find a new line
         if ((pt_current = static_cast<const char*>(memchr(pt_current, '\n', pt_last - pt_current)))) {
            // header
            if (line_order == 0) {
               current_read_index_write = read_vector_index * 3;

               // include '\n'
               read_vector[current_read_index_write].assign(pt_current_line_start, pt_current - pt_current_line_start + 1);

               line_order++;
            }
            // sequence
            else if (line_order == 1) {
               read_vector[current_read_index_write + 1].assign(pt_current_line_start, pt_current - pt_current_line_start);

               line_order++;
            }
            // connector
            else if (line_order == 2) {
               line_order++;
            }
            // quality score
            else {
               read_vector[current_read_index_write + 2].assign(pt_current_line_start, pt_current - pt_current_line_start);

               line_order = 0;
               processed_bytes += (pt_current - pt_current_read_start + 1);
               pt_current_read_start = pt_current + 1;
               read_vector_index++;

               num_reads_local++;

               // all reads have been processed
               // stop the process
               if (num_reads_local == num_reads_vector1[rank_node]) {
                  all_reads  = true;
                  keep_going = false;
               }

               // read_vector is full
               if (read_vector_index == read_block_size) {
                  //--------------------------------------------------
                  // correct reads in a block
                  //--------------------------------------------------
                  #pragma omp parallel shared(bit_vector, hash_seed, read_vector) private(current_read_index, read_length, min_check_length, max_allowed_ns, max_trimmed_bases, num_corrected_errors_local, trim_5_end, trim_3_end, sequence_modification, too_many_errors) reduction(+: num_corrected_errors_tmp, num_corrected_reads_tmp, num_trimmed_bases_tmp)
                  {
                     // iterate reads
                     #pragma omp for schedule(dynamic, OPENMP_CHUNK_SIZE)
                     for (std::size_t it_read = 0; it_read < read_block_size; it_read++) {
                        // calculate the current index
                        current_read_index = it_read * 3;

                        // change sequences to upper case
                        std::transform(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), read_vector[current_read_index + 1].begin(), static_cast<int(*)(int)>(std::toupper));

                        // substitute non-standard characters with Ns
                        read_vector[current_read_index + 1] = std::regex_replace(read_vector[current_read_index + 1], regex_non_acgtn, "N");

                        // set various thresholds
                        read_length      = read_vector[current_read_index + 1].length();
                        min_check_length = read_length * CHECK_RANGE_RATIO;
                        max_allowed_ns   = read_length * MAX_N_RATIO;

                        // forward
                        // too short read: no trimming
                        if (read_length <= MIN_BASES_AFTER_TRIMMING) {
                           max_trimmed_bases = 0;
                        }
                        else {
                           max_trimmed_bases = std::min((std::size_t)(read_length * MAX_TRIMMING_RATE), read_length - MIN_BASES_AFTER_TRIMMING);
                        }

                        // check the number of Ns in the read
                        if (((std::size_t)std::count(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), 'N') <= max_allowed_ns) &&
                            (read_length >= kmer_length)) {
                           // substitute Ns other characters
                           std::replace(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), 'N', SUBST_CHAR);

                           // initialize variables
                           num_corrected_errors_local = 0;
                           trim_5_end                 = 0;
                           trim_3_end                 = 0;

                           // storage for the modification of the reads
                           // # of entries = read length
                           sequence_modification.assign(read_length, '0');

                           //----------------------------------------------------------------------
                           // correct errors in a read
                           //----------------------------------------------------------------------
                           // forward read
                           // errors cannot be corrected if k is equal to read length
                           if (read_length > kmer_length) {
                              correct_errors_in_a_read_fastq(
                                                             read_vector[current_read_index + 1],
                                                             sequence_modification,
                                                             read_vector[current_read_index + 2],
                                                             trim_5_end,
                                                             trim_3_end,
                                                             read_length,
                                                             max_trimmed_bases,
                                                             min_check_length,
                                                             num_corrected_errors_local,
                                                             bit_vector,
                                                             hash_seed
                                                            );
                           }

                           // no trim
                           if (c_inst_args.notrim == true) {
                              trim_5_end = 0;
                              trim_3_end = 0;
                           }
                           // adjust the number of trimmed bases
                           else {
                              if ((trim_5_end + trim_3_end) > max_trimmed_bases) {
                                 if (trim_3_end <= max_trimmed_bases) {
                                    trim_5_end = 0;
                                 }
                                 else if (trim_5_end <= max_trimmed_bases) {
                                    trim_3_end = 0;
                                 }
                                 else {
                                    trim_5_end = 0;
                                    trim_3_end = max_trimmed_bases;
                                 }
                              }
                           }

                           num_trimmed_bases_tmp += (trim_5_end + trim_3_end);

                           // update num_corrected_reads
                           too_many_errors = false;
                           if (num_corrected_errors_local > (read_vector[current_read_index + 1].length() * MAX_ERROR_RATE)) {
                              too_many_errors = true;
                           }
                           else if (num_corrected_errors_local > 0) {
                              num_corrected_errors_tmp += num_corrected_errors_local;

                              num_corrected_reads_tmp++;
                           }
                           else if (c_inst_args.notrim == false) {
                              if ((trim_5_end > 0) || (trim_3_end > 0)) {
                                 num_corrected_reads_tmp++;
                              }
                           }

                           // make a corrected read
                           if (too_many_errors == false) {
                              // apply modifications to the read
                              for (unsigned int it_base = trim_5_end; it_base < (read_vector[current_read_index + 1].length() - trim_3_end); it_base++) {
                                 if (sequence_modification[it_base] != '0') {
                                    read_vector[current_read_index + 1][it_base] = sequence_modification[it_base];
                                 }
                              }
                           }

                           // make a trimmed read
                           if ((trim_5_end + trim_3_end) > 0) {
                              read_vector[current_read_index + 1] = read_vector[current_read_index + 1].substr(trim_5_end, read_vector[current_read_index + 1].length() - trim_5_end - trim_3_end) + '\n';
                              read_vector[current_read_index + 2] = read_vector[current_read_index + 2].substr(trim_5_end, read_vector[current_read_index + 2].length() - trim_5_end - trim_3_end) + '\n';
                           }
                           else {
                              read_vector[current_read_index + 1] = read_vector[current_read_index + 1] + '\n';
                              read_vector[current_read_index + 2] = read_vector[current_read_index + 2] + '\n';
                           }
                        }
                        // too many Ns
                        // write reads without modification
                        else {
                           // sequence
                           read_vector[current_read_index + 1] += '\n';

                           // quality score
                           read_vector[current_read_index + 2] += '\n';
                        }
                     // it_read
                     }
                  // omp parallel
                  }

                  // write corrected_reads
                  for (std::size_t it_write = 0; it_write < read_block_size; it_write++) {
                     current_read_index = it_write * 3;

                     // header
                     fwrite(read_vector[current_read_index].c_str(), 1, read_vector[current_read_index].length(), f_corrected_read1);
                     // sequence
                     fwrite(read_vector[current_read_index + 1].c_str(), 1, read_vector[current_read_index + 1].length(), f_corrected_read1);
                     // connector
                     fwrite("+\n", 1, 2, f_corrected_read1);
                     // quality score
                     fwrite(read_vector[current_read_index + 2].c_str(), 1, read_vector[current_read_index + 2].length(), f_corrected_read1);
                  }

                  read_vector_index = 0;
               }
            }

            pt_current++;
            pt_current_line_start = pt_current;
         }
      }

      remaining_bytes  = read_file_size_byte1 - processed_bytes;

      if (remaining_bytes <= 0) {
         keep_going = false;
      }
      else if (remaining_bytes == prev_remaining_bytes) {
         std::cout << std::endl << "ERROR: The number of lines in " << c_inst_args.read_file_name1 << " is wrong" << std::endl << std::endl;
         MPI_Abort(MPI_COMM_WORLD, 220);
      }
      else {
         // align processed bytes to the multiple of unit_alignment_offset
         alignment_offset     = processed_bytes % read_file_unit_size_byte1;
         processed_bytes     -= alignment_offset;
         prev_remaining_bytes = remaining_bytes;
         remaining_bytes     += alignment_offset;
      }

      f_read.close();
   }

   // correct errors in remaining reads
   if (read_vector_index > 0) {
      #pragma omp parallel shared(bit_vector, hash_seed, read_vector) private(current_read_index, read_length, min_check_length, max_allowed_ns, max_trimmed_bases, num_corrected_errors_local, trim_5_end, trim_3_end, sequence_modification, too_many_errors) reduction(+: num_corrected_errors_tmp, num_corrected_reads_tmp, num_trimmed_bases_tmp)
      {
         // iterate reads
         #pragma omp for schedule(dynamic, OPENMP_CHUNK_SIZE)
         for (std::size_t it_read = 0; it_read < read_vector_index; it_read++) {
            // calculate the current index
            current_read_index = it_read * 3;

            // change sequences to upper case
            std::transform(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), read_vector[current_read_index + 1].begin(), static_cast<int(*)(int)>(std::toupper));

            // substitute non-standard characters with Ns
            read_vector[current_read_index + 1] = std::regex_replace(read_vector[current_read_index + 1], regex_non_acgtn, "N");

            // set various thresholds
            read_length      = read_vector[current_read_index + 1].length();
            min_check_length = read_length * CHECK_RANGE_RATIO;
            max_allowed_ns   = read_length * MAX_N_RATIO;

            // forward
            // too short read: no trimming
            if (read_length <= MIN_BASES_AFTER_TRIMMING) {
               max_trimmed_bases = 0;
            }
            else {
               max_trimmed_bases = std::min((std::size_t)(read_length * MAX_TRIMMING_RATE), read_length - MIN_BASES_AFTER_TRIMMING);
            }

            // check the number of Ns in the read
            if (((std::size_t)std::count(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), 'N') <= max_allowed_ns) &&
                (read_length >= kmer_length)) {
               // substitute Ns other characters
               std::replace(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), 'N', SUBST_CHAR);

               // initialize variables
               num_corrected_errors_local = 0;
               trim_5_end                 = 0;
               trim_3_end                 = 0;

               // storage for the modification of the reads
               // # of entries = read length
               sequence_modification.assign(read_length, '0');

               //----------------------------------------------------------------------
               // correct errors in a read
               //----------------------------------------------------------------------
               // forward read
               // errors cannot be corrected if k is equal to read length
               if (read_length > kmer_length) {
                  correct_errors_in_a_read_fastq(
                                                 read_vector[current_read_index + 1],
                                                 sequence_modification,
                                                 read_vector[current_read_index + 2],
                                                 trim_5_end,
                                                 trim_3_end,
                                                 read_length,
                                                 max_trimmed_bases,
                                                 min_check_length,
                                                 num_corrected_errors_local,
                                                 bit_vector,
                                                 hash_seed
                                                );
               }

               // no trim
               if (c_inst_args.notrim == true) {
                  trim_5_end = 0;
                  trim_3_end = 0;
               }
               // adjust the number of trimmed bases
               else {
                  if ((trim_5_end + trim_3_end) > max_trimmed_bases) {
                     if (trim_3_end <= max_trimmed_bases) {
                        trim_5_end = 0;
                     }
                     else if (trim_5_end <= max_trimmed_bases) {
                        trim_3_end = 0;
                     }
                     else {
                        trim_5_end = 0;
                        trim_3_end = max_trimmed_bases;
                     }
                  }
               }

               num_trimmed_bases_tmp += (trim_5_end + trim_3_end);

               // update num_corrected_reads
               too_many_errors = false;
               if (num_corrected_errors_local > (read_vector[current_read_index + 1].length() * MAX_ERROR_RATE)) {
                  too_many_errors = true;
               }
               else if (num_corrected_errors_local > 0) {
                  num_corrected_errors_tmp += num_corrected_errors_local;

                  num_corrected_reads_tmp++;
               }
               else if (c_inst_args.notrim == false) {
                  if ((trim_5_end > 0) || (trim_3_end > 0)) {
                     num_corrected_reads_tmp++;
                  }
               }

               // make a corrected read
               if (too_many_errors == false) {
                  // apply modifications to the read
                  for (unsigned int it_base = trim_5_end; it_base < (read_vector[current_read_index + 1].length() - trim_3_end); it_base++) {
                     if (sequence_modification[it_base] != '0') {
                        read_vector[current_read_index + 1][it_base] = sequence_modification[it_base];
                     }
                  }
               }

               // make a trimmed read
               if ((trim_5_end + trim_3_end) > 0) {
                  read_vector[current_read_index + 1] = read_vector[current_read_index + 1].substr(trim_5_end, read_vector[current_read_index + 1].length() - trim_5_end - trim_3_end) + '\n';
                  read_vector[current_read_index + 2] = read_vector[current_read_index + 2].substr(trim_5_end, read_vector[current_read_index + 2].length() - trim_5_end - trim_3_end) + '\n';
               }
               else {
                  read_vector[current_read_index + 1] = read_vector[current_read_index + 1] + '\n';
                  read_vector[current_read_index + 2] = read_vector[current_read_index + 2] + '\n';
               }
            }
            // too many Ns
            // write reads without modification
            else {
               // sequence
               read_vector[current_read_index + 1] += '\n';

               // quality score
               read_vector[current_read_index + 2] += '\n';
            }
         // it_read
         }

      // omp parallel
      }

      // write corrected_reads
      for (std::size_t it_read = 0; it_read < read_vector_index; it_read++) {
         current_read_index = it_read * 3;

         // header
         fwrite(read_vector[current_read_index].c_str(), 1, read_vector[current_read_index].length(), f_corrected_read1);
         // sequence
         fwrite(read_vector[current_read_index + 1].c_str(), 1, read_vector[current_read_index + 1].length(), f_corrected_read1);
         // connector
         fwrite("+\n", 1, 2, f_corrected_read1);
         // quality score
         fwrite(read_vector[current_read_index + 2].c_str(), 1, read_vector[current_read_index + 2].length(), f_corrected_read1);
      }
   }

   fclose(f_corrected_read1);

   //
   // reverse
   //
   // open a temporary output file
   std::string corrected_read_file_name2(c_inst_args.corrected_read_file_name2 + '.' + rank_node_text);

   FILE* f_corrected_read2(fopen(corrected_read_file_name2.c_str(), "w"));
   if (f_corrected_read2 == NULL) {
      std::cout << std::endl << "ERROR: Cannot open " << corrected_read_file_name2 << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 221);
   }

   // initialize variables
   num_reads_local          = 0;
   keep_going               = true;
   all_reads                = false;
   processed_bytes          = (starting_point_vector2[rank_node] / read_file_unit_size_byte2) * read_file_unit_size_byte2;
   alignment_offset         = starting_point_vector2[rank_node] % read_file_unit_size_byte2;
   read_vector_index        = 0;
   current_read_index_write = 0;
   remaining_bytes          = read_file_size_byte2 - starting_point_vector2[rank_node] + 1;
   prev_remaining_bytes     = read_file_size_byte2 - starting_point_vector2[rank_node] + 1;

   // check whether the input file is empty
   if (remaining_bytes <= 0) {
      std::cout << std::endl << "ERROR: Remaining bytes are " << remaining_bytes << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 222);
   }

   // read loop
   while (keep_going) {
      // calculate the number of bytes that will be processed in this iteration
      if (remaining_bytes <= MMAP_FILE_SIZE) {
         this_iteration_bytes = remaining_bytes;
      }
      else {
         this_iteration_bytes = MMAP_FILE_SIZE;
      }

      // open the file
      //open(<file name>, <length>, <offset>)
      f_read.open(c_inst_args.read_file_name2.c_str(), this_iteration_bytes, processed_bytes);

      if (!f_read.is_open()) {
         std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.read_file_name2 << std::endl << std::endl;
         MPI_Abort(MPI_COMM_WORLD, 223);
      }

      // initialize variables
      pt_current            = f_read.data() + alignment_offset;
      pt_last               = f_read.data() + f_read.size();
      pt_current_read_start = pt_current;
      pt_current_line_start = pt_current;

      processed_bytes += alignment_offset;

      line_order = 0;

      while (pt_current && (pt_current != pt_last) && (!all_reads)) {
         // find a new line
         if ((pt_current = static_cast<const char*>(memchr(pt_current, '\n', pt_last - pt_current)))) {
            // header
            if (line_order == 0) {
               current_read_index_write = read_vector_index * 3;

               // include '\n'
               read_vector[current_read_index_write].assign(pt_current_line_start, pt_current - pt_current_line_start + 1);

               line_order++;
            }
            // sequence
            else if (line_order == 1) {
               read_vector[current_read_index_write + 1].assign(pt_current_line_start, pt_current - pt_current_line_start);

               line_order++;
            }
            // connector
            else if (line_order == 2) {
               line_order++;
            }
            // quality score
            else {
               read_vector[current_read_index_write + 2].assign(pt_current_line_start, pt_current - pt_current_line_start);

               line_order = 0;
               processed_bytes += (pt_current - pt_current_read_start + 1);
               pt_current_read_start = pt_current + 1;
               read_vector_index++;

               num_reads_local++;

               // all reads have been processed
               // stop the process
               if (num_reads_local == num_reads_vector2[rank_node]) {
                  all_reads  = true;
                  keep_going = false;
               }

               // read_vector is full
               if (read_vector_index == read_block_size) {
                  //--------------------------------------------------
                  // correct reads in a block
                  //--------------------------------------------------
                  #pragma omp parallel shared(bit_vector, hash_seed, read_vector) private(current_read_index, read_length, min_check_length, max_allowed_ns, max_trimmed_bases, num_corrected_errors_local, trim_5_end, trim_3_end, sequence_modification, too_many_errors) reduction(+: num_corrected_errors_tmp, num_corrected_reads_tmp, num_trimmed_bases_tmp)
                  {
                     // iterate reads
                     #pragma omp for schedule(dynamic, OPENMP_CHUNK_SIZE)
                     for (std::size_t it_read = 0; it_read < read_block_size; it_read++) {
                        // calculate the current index
                        current_read_index = it_read * 3;

                        // change sequences to upper case
                        std::transform(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), read_vector[current_read_index + 1].begin(), static_cast<int(*)(int)>(std::toupper));

                        // substitute non-standard characters with Ns
                        read_vector[current_read_index + 1] = std::regex_replace(read_vector[current_read_index + 1], regex_non_acgtn, "N");

                        // set various thresholds
                        read_length      = read_vector[current_read_index + 1].length();
                        min_check_length = read_length * CHECK_RANGE_RATIO;
                        max_allowed_ns   = read_length * MAX_N_RATIO;

                        // forward
                        // too short read: no trimming
                        if (read_length <= MIN_BASES_AFTER_TRIMMING) {
                           max_trimmed_bases = 0;
                        }
                        else {
                           max_trimmed_bases = std::min((std::size_t)(read_length * MAX_TRIMMING_RATE), read_length - MIN_BASES_AFTER_TRIMMING);
                        }

                        // check the number of Ns in the read
                        if (((std::size_t)std::count(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), 'N') <= max_allowed_ns) &&
                            (read_length >= kmer_length)) {
                           // substitute Ns other characters
                           std::replace(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), 'N', SUBST_CHAR);

                           // initialize variables
                           num_corrected_errors_local = 0;
                           trim_5_end                 = 0;
                           trim_3_end                 = 0;

                           // storage for the modification of the reads
                           // # of entries = read length
                           sequence_modification.assign(read_length, '0');

                           //----------------------------------------------------------------------
                           // correct errors in a read
                           //----------------------------------------------------------------------
                           // forward read
                           // errors cannot be corrected if k is equal to read length
                           if (read_length > kmer_length) {
                              correct_errors_in_a_read_fastq(
                                                             read_vector[current_read_index + 1],
                                                             sequence_modification,
                                                             read_vector[current_read_index + 2],
                                                             trim_5_end,
                                                             trim_3_end,
                                                             read_length,
                                                             max_trimmed_bases,
                                                             min_check_length,
                                                             num_corrected_errors_local,
                                                             bit_vector,
                                                             hash_seed
                                                            );
                           }

                           // no trim
                           if (c_inst_args.notrim == true) {
                              trim_5_end = 0;
                              trim_3_end = 0;
                           }
                           // adjust the number of trimmed bases
                           else {
                              if ((trim_5_end + trim_3_end) > max_trimmed_bases) {
                                 if (trim_3_end <= max_trimmed_bases) {
                                    trim_5_end = 0;
                                 }
                                 else if (trim_5_end <= max_trimmed_bases) {
                                    trim_3_end = 0;
                                 }
                                 else {
                                    trim_5_end = 0;
                                    trim_3_end = max_trimmed_bases;
                                 }
                              }
                           }

                           num_trimmed_bases_tmp += (trim_5_end + trim_3_end);

                           // update num_corrected_reads
                           too_many_errors = false;
                           if (num_corrected_errors_local > (read_vector[current_read_index + 1].length() * MAX_ERROR_RATE)) {
                              too_many_errors = true;
                           }
                           else if (num_corrected_errors_local > 0) {
                              num_corrected_errors_tmp += num_corrected_errors_local;

                              num_corrected_reads_tmp++;
                           }
                           else if (c_inst_args.notrim == false) {
                              if ((trim_5_end > 0) || (trim_3_end > 0)) {
                                 num_corrected_reads_tmp++;
                              }
                           }

                           // make a corrected read
                           if (too_many_errors == false) {
                              // apply modifications to the read
                              for (unsigned int it_base = trim_5_end; it_base < (read_vector[current_read_index + 1].length() - trim_3_end); it_base++) {
                                 if (sequence_modification[it_base] != '0') {
                                    read_vector[current_read_index + 1][it_base] = sequence_modification[it_base];
                                 }
                              }
                           }

                           // make a trimmed read
                           if ((trim_5_end + trim_3_end) > 0) {
                              read_vector[current_read_index + 1] = read_vector[current_read_index + 1].substr(trim_5_end, read_vector[current_read_index + 1].length() - trim_5_end - trim_3_end) + '\n';
                              read_vector[current_read_index + 2] = read_vector[current_read_index + 2].substr(trim_5_end, read_vector[current_read_index + 2].length() - trim_5_end - trim_3_end) + '\n';
                           }
                           else {
                              read_vector[current_read_index + 1] = read_vector[current_read_index + 1] + '\n';
                              read_vector[current_read_index + 2] = read_vector[current_read_index + 2] + '\n';
                           }
                        }
                        // too many Ns
                        // write reads without modification
                        else {
                           // sequence
                           read_vector[current_read_index + 1] += '\n';

                           // quality score
                           read_vector[current_read_index + 2] += '\n';
                        }
                     // it_read
                     }
                  // omp parallel
                  }

                  // write corrected_reads
                  for (std::size_t it_write = 0; it_write < read_block_size; it_write++) {
                     current_read_index = it_write * 3;

                     // header
                     fwrite(read_vector[current_read_index].c_str(), 1, read_vector[current_read_index].length(), f_corrected_read2);
                     // sequence
                     fwrite(read_vector[current_read_index + 1].c_str(), 1, read_vector[current_read_index + 1].length(), f_corrected_read2);
                     // connector
                     fwrite("+\n", 1, 2, f_corrected_read2);
                     // quality score
                     fwrite(read_vector[current_read_index + 2].c_str(), 1, read_vector[current_read_index + 2].length(), f_corrected_read2);
                  }

                  read_vector_index = 0;
               }
            }

            pt_current++;
            pt_current_line_start = pt_current;
         }
      }

      remaining_bytes  = read_file_size_byte2 - processed_bytes;

      if (remaining_bytes <= 0) {
         keep_going = false;
      }
      else if (remaining_bytes == prev_remaining_bytes) {
         std::cout << std::endl << "ERROR: The number of lines in " << c_inst_args.read_file_name2 << " is wrong" << std::endl << std::endl;
         MPI_Abort(MPI_COMM_WORLD, 224);
      }
      else {
         // align processed bytes to the multiple of unit_alignment_offset
         alignment_offset     = processed_bytes % read_file_unit_size_byte2;
         processed_bytes     -= alignment_offset;
         prev_remaining_bytes = remaining_bytes;
         remaining_bytes     += alignment_offset;
      }

      f_read.close();
   }

   // correct errors in remaining reads
   if (read_vector_index > 0) {
      #pragma omp parallel shared(bit_vector, hash_seed, read_vector) private(current_read_index, read_length, min_check_length, max_allowed_ns, max_trimmed_bases, num_corrected_errors_local, trim_5_end, trim_3_end, sequence_modification, too_many_errors) reduction(+: num_corrected_errors_tmp, num_corrected_reads_tmp, num_trimmed_bases_tmp)
      {
         // iterate reads
         #pragma omp for schedule(dynamic, OPENMP_CHUNK_SIZE)
         for (std::size_t it_read = 0; it_read < read_vector_index; it_read++) {
            // calculate the current index
            current_read_index = it_read * 3;

            // change sequences to upper case
            std::transform(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), read_vector[current_read_index + 1].begin(), static_cast<int(*)(int)>(std::toupper));

            // substitute non-standard characters with Ns
            read_vector[current_read_index + 1] = std::regex_replace(read_vector[current_read_index + 1], regex_non_acgtn, "N");

            // set various thresholds
            read_length      = read_vector[current_read_index + 1].length();
            min_check_length = read_length * CHECK_RANGE_RATIO;
            max_allowed_ns   = read_length * MAX_N_RATIO;

            // forward
            // too short read: no trimming
            if (read_length <= MIN_BASES_AFTER_TRIMMING) {
               max_trimmed_bases = 0;
            }
            else {
               max_trimmed_bases = std::min((std::size_t)(read_length * MAX_TRIMMING_RATE), read_length - MIN_BASES_AFTER_TRIMMING);
            }

            // check the number of Ns in the read
            if (((std::size_t)std::count(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), 'N') <= max_allowed_ns) &&
                (read_length >= kmer_length)) {
               // substitute Ns other characters
               std::replace(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), 'N', SUBST_CHAR);

               // initialize variables
               num_corrected_errors_local = 0;
               trim_5_end                 = 0;
               trim_3_end                 = 0;

               // storage for the modification of the reads
               // # of entries = read length
               sequence_modification.assign(read_length, '0');

               //----------------------------------------------------------------------
               // correct errors in a read
               //----------------------------------------------------------------------
               // forward read
               // errors cannot be corrected if k is equal to read length
               if (read_length > kmer_length) {
                  correct_errors_in_a_read_fastq(
                                                 read_vector[current_read_index + 1],
                                                 sequence_modification,
                                                 read_vector[current_read_index + 2],
                                                 trim_5_end,
                                                 trim_3_end,
                                                 read_length,
                                                 max_trimmed_bases,
                                                 min_check_length,
                                                 num_corrected_errors_local,
                                                 bit_vector,
                                                 hash_seed
                                                );
               }

               // no trim
               if (c_inst_args.notrim == true) {
                  trim_5_end = 0;
                  trim_3_end = 0;
               }
               // adjust the number of trimmed bases
               else {
                  if ((trim_5_end + trim_3_end) > max_trimmed_bases) {
                     if (trim_3_end <= max_trimmed_bases) {
                        trim_5_end = 0;
                     }
                     else if (trim_5_end <= max_trimmed_bases) {
                        trim_3_end = 0;
                     }
                     else {
                        trim_5_end = 0;
                        trim_3_end = max_trimmed_bases;
                     }
                  }
               }

               num_trimmed_bases_tmp += (trim_5_end + trim_3_end);

               // update num_corrected_reads
               too_many_errors = false;
               if (num_corrected_errors_local > (read_vector[current_read_index + 1].length() * MAX_ERROR_RATE)) {
                  too_many_errors = true;
               }
               else if (num_corrected_errors_local > 0) {
                  num_corrected_errors_tmp += num_corrected_errors_local;

                  num_corrected_reads_tmp++;
               }
               else if (c_inst_args.notrim == false) {
                  if ((trim_5_end > 0) || (trim_3_end > 0)) {
                     num_corrected_reads_tmp++;
                  }
               }

               // make a corrected read
               if (too_many_errors == false) {
                  // apply modifications to the read
                  for (unsigned int it_base = trim_5_end; it_base < (read_vector[current_read_index + 1].length() - trim_3_end); it_base++) {
                     if (sequence_modification[it_base] != '0') {
                        read_vector[current_read_index + 1][it_base] = sequence_modification[it_base];
                     }
                  }
               }

               // make a trimmed quality score
               if ((trim_5_end + trim_3_end) > 0) {
                  read_vector[current_read_index + 1] = read_vector[current_read_index + 1].substr(trim_5_end, read_vector[current_read_index + 1].length() - trim_5_end - trim_3_end) + '\n';
                  read_vector[current_read_index + 2] = read_vector[current_read_index + 2].substr(trim_5_end, read_vector[current_read_index + 2].length() - trim_5_end - trim_3_end) + '\n';
               }
               else {
                  read_vector[current_read_index + 1] = read_vector[current_read_index + 1] + '\n';
                  read_vector[current_read_index + 2] = read_vector[current_read_index + 2] + '\n';
               }
            }
            // too many Ns
            // write reads without modification
            else {
               // sequence
               read_vector[current_read_index + 1] += '\n';

               // quality score
               read_vector[current_read_index + 2] += '\n';
            }
         // it_read
         }

      // omp parallel
      }

      // write corrected_reads
      for (std::size_t it_read = 0; it_read < read_vector_index; it_read++) {
         current_read_index = it_read * 3;

         // header
         fwrite(read_vector[current_read_index].c_str(), 1, read_vector[current_read_index].length(), f_corrected_read2);
         // sequence
         fwrite(read_vector[current_read_index + 1].c_str(), 1, read_vector[current_read_index + 1].length(), f_corrected_read2);
         // connector
         fwrite("+\n", 1, 2, f_corrected_read2);
         // quality score
         fwrite(read_vector[current_read_index + 2].c_str(), 1, read_vector[current_read_index + 2].length(), f_corrected_read2);
      }
   }

   fclose(f_corrected_read2);

   read_vector.clear();

   num_corrected_errors = num_corrected_errors_tmp;
   num_corrected_reads  = num_corrected_reads_tmp;
   num_trimmed_bases    = num_trimmed_bases_tmp;
}



//----------------------------------------------------------------------
// correct_errors_in_a_read_fastq
//----------------------------------------------------------------------
inline void C_correct_errors::correct_errors_in_a_read_fastq(const std::string& sequence, std::string& sequence_modification, const std::string& quality_score, std::size_t& trim_5_end, std::size_t& trim_3_end, const std::size_t& read_length, const std::size_t& max_trimmed_bases, const std::size_t& min_check_length, std::size_t& num_corrected_errors_local, const unsigned char* bit_vector, const std::vector<unsigned int>& hash_seed) {
   //--------------------------------------------------
   // STEP 0-0: find solid k-mers in this read
   //--------------------------------------------------
   // variables
   std::size_t num_kmers(read_length - kmer_length + 1);

   std::vector< std::pair<std::size_t, std::size_t> > solid_regions;
   std::vector< std::pair<std::size_t, std::size_t> > solid_regions_org;

   bool is_solid_kmer_prev(false);

   // find solid regions
   std::pair<std::size_t, std::size_t> new_solid_region;
   for (std::size_t it_kmer = 0; it_kmer < num_kmers; it_kmer++) {
      std::string current_kmer(sequence.substr(it_kmer, kmer_length));

      // k-mer is solid
      if (query_text(current_kmer, bit_vector, hash_seed) == true) {
         // start point of a solid region
         if (is_solid_kmer_prev == false) {
            new_solid_region.first = it_kmer;
            is_solid_kmer_prev = true;
         }
      }
      else {
         // end point of a solid region
         if (is_solid_kmer_prev == true) {
            new_solid_region.second = it_kmer - 1;

            if (new_solid_region.second < new_solid_region.first) {
               std::cout << std::endl << "ERROR: The second index is smaller than the first" << std::endl << std::endl;
               MPI_Abort(MPI_COMM_WORLD, 225);
            }
            solid_regions.push_back(new_solid_region);
            solid_regions_org.push_back(new_solid_region);

            is_solid_kmer_prev = false;
         }
      }
   }

   // last solid region
   if (is_solid_kmer_prev == true) {
      new_solid_region.second = num_kmers - 1;
      solid_regions.push_back(new_solid_region);
      solid_regions_org.push_back(new_solid_region);
   }

   //--------------------------------------------------
   // STEP 0-1: adjust solid regioins using quality scores
   //--------------------------------------------------
   // solid_regions_org: indexes that are not modified
   // when the indices are reached, all kinds of modifications (A/C/G/T) should be made and checked
   // at least one solid region
   if (solid_regions.size() > 0) {
      // exceptional case: only one solid island that covers the entire read
      if ((solid_regions.size() != 1) || (solid_regions[0].first != 0) || (solid_regions[0].second != (read_length - kmer_length))) {
         // at least two solid k-mer islands
         if (solid_regions.size() > 1) {
            // check the distance between every two solid k-mer islands
            bool flag_short_distance(false);

            for (std::size_t it_sr = 0; it_sr < solid_regions.size() - 1; it_sr++) {
               if ((solid_regions[it_sr + 1].first - solid_regions[it_sr].second) < kmer_length) {
                  flag_short_distance = true;
               }
            }

            if (flag_short_distance == true) {
               std::vector< std::pair<std::size_t, std::size_t> > solid_regions_tmp;
               std::vector< std::pair<std::size_t, std::size_t> > solid_regions_org_tmp;

               // each solid island
               for (std::size_t it_sr = 0; it_sr < solid_regions.size(); it_sr++) {
                  // each base in the solid island (0-base)
                  std::size_t num_low_quality_base(0);
                  // set an initial value to avoid a compilation warning
                  std::size_t index_prev_low_quality_base(0);
                  for (unsigned int it_base = solid_regions[it_sr].first; it_base < solid_regions[it_sr].second + kmer_length; it_base++) {
                     // a current base has a low quality score
                     if (((unsigned short int)quality_score[it_base] - quality_score_offset) < quality_score_cutoff) {
                        num_low_quality_base++;

                        // first low quality base
                        if (num_low_quality_base == 1) {
                           // the low quality base is not in the first k-mer of the solid island
                           if (it_base >= (solid_regions[it_sr].first + kmer_length)) {
                              // add the left most high quality region to a temporary vector
                              std::pair<std::size_t, std::size_t> new_solid_region;

                              new_solid_region.first  = solid_regions[it_sr].first;
                              new_solid_region.second = it_base - kmer_length;
                              solid_regions_tmp.push_back(new_solid_region);

                              solid_regions_org_tmp.push_back(solid_regions_org[it_sr]);
                           }
                        }
                        // not first low quality base
                        else {
                           if ((it_base - index_prev_low_quality_base) > kmer_length) {
                              std::pair<std::size_t, std::size_t> new_solid_region;

                              new_solid_region.first = index_prev_low_quality_base + 1;
                              new_solid_region.second = it_base - kmer_length;
                              solid_regions_tmp.push_back(new_solid_region);

                              solid_regions_org_tmp.push_back(solid_regions_org[it_sr]);
                           }
                        }

                        index_prev_low_quality_base = it_base;
                     }
                  }

                  // process the bases to the right of the rightmost low quality base
                  if (num_low_quality_base > 0) {
                     if (solid_regions[it_sr].second >= (index_prev_low_quality_base + kmer_length)) {
                        std::pair<std::size_t, std::size_t> new_solid_region;

                        new_solid_region.first  = index_prev_low_quality_base + kmer_length;
                        new_solid_region.second = solid_regions[it_sr].second;
                        solid_regions_tmp.push_back(new_solid_region);

                        solid_regions_org_tmp.push_back(solid_regions_org[it_sr]);
                     }
                  }
                  // no low quality base
                  // add the current solid island
                  else {
                     solid_regions_tmp.push_back(solid_regions[it_sr]);

                     solid_regions_org_tmp.push_back(solid_regions_org[it_sr]);
                  }
               }

               solid_regions     = solid_regions_tmp;
               solid_regions_org = solid_regions_org_tmp;
            }
         }
         // only one solid k-mer island
         else if (solid_regions.size() == 1) {
            std::vector< std::pair<std::size_t, std::size_t> > solid_regions_tmp;
            std::vector< std::pair<std::size_t, std::size_t> > solid_regions_org_tmp;

            std::size_t num_low_quality_base(0);
            std::size_t prev_low_quality_index(0);

            // each base in the solid island (0-base)
            for (unsigned int it_base = solid_regions[0].first; it_base < solid_regions[0].second + kmer_length; it_base++) {
               // a current base has a low quality score
               if (((unsigned short int)quality_score[it_base] - quality_score_offset) < quality_score_cutoff) {
                  num_low_quality_base++;

                  // first low quality base
                  if (num_low_quality_base == 1) {
                     if ((it_base - solid_regions[0].first) >= (kmer_length + MIN_SOLID_LENGTH - 1)) {
                        std::pair<std::size_t, std::size_t> new_solid_region;

                        new_solid_region.first = solid_regions[0].first;
                        new_solid_region.second = it_base - kmer_length;
                        solid_regions_tmp.push_back(new_solid_region);

                        solid_regions_org_tmp.push_back(solid_regions_org[0]);
                     }

                     prev_low_quality_index = it_base;
                  }
                  // not first low quality base
                  else {
                     if ((it_base - prev_low_quality_index) >= (kmer_length + MIN_SOLID_LENGTH)) {
                        std::pair<std::size_t, std::size_t> new_solid_region;

                        new_solid_region.first = prev_low_quality_index + 1;
                        new_solid_region.second = it_base - kmer_length;
                        solid_regions_tmp.push_back(new_solid_region);

                        solid_regions_org_tmp.push_back(solid_regions_org[0]);
                     }

                     prev_low_quality_index = it_base;
                  }
               }
            }

            // the above is done only when this procedure does not remove the only solid island
            if (solid_regions_tmp.size() > 0) {
               solid_regions     = solid_regions_tmp;
               solid_regions_org = solid_regions_org_tmp;
            }
         }
      }
   }

   //--------------------------------------------------
   // STEP 0-2: remove short solid regions
   //--------------------------------------------------
   if (solid_regions.size() > 0) {
      std::vector< std::pair<std::size_t, std::size_t> > solid_regions_tmp;
      std::vector< std::pair<std::size_t, std::size_t> > solid_regions_org_tmp;

      for (std::size_t it_region = 0; it_region < solid_regions.size(); it_region++) {
         if ((solid_regions[it_region].second - solid_regions[it_region].first + 1) >= MIN_SOLID_LENGTH) {
            solid_regions_tmp.push_back(solid_regions[it_region]);
            solid_regions_org_tmp.push_back(solid_regions_org[it_region]);
         }
      }

      solid_regions     = solid_regions_tmp;
      solid_regions_org = solid_regions_org_tmp;
   }

   //--------------------------------------------------
   // STEP 0-3: remove short non-solid regions
   //--------------------------------------------------
   if (solid_regions.size() > 0) {
      std::vector< std::pair<std::size_t, std::size_t> > solid_regions_tmp;
      std::vector< std::pair<std::size_t, std::size_t> > solid_regions_org_tmp;

      solid_regions_tmp.push_back(solid_regions[0]);
      solid_regions_org_tmp.push_back(solid_regions_org[0]);

      if (solid_regions.size() > 1) {
         for (std::size_t it_region = 1; it_region < solid_regions.size(); it_region++) {
            if ((solid_regions[it_region].first - solid_regions[it_region - 1].second - 1) < MIN_NON_SOLID_LENGTH) {
               solid_regions_tmp[solid_regions_tmp.size() - 1].second = solid_regions[it_region].second;
            }
            else {
               solid_regions_tmp.push_back(solid_regions[it_region]);
               solid_regions_org_tmp.push_back(solid_regions_org[it_region]);
            }
         }
      }
      solid_regions     = solid_regions_tmp;
      solid_regions_org = solid_regions_org_tmp;
   }

   //--------------------------------------------------
   // STEP 0-4: reduce the size of solid regions
   //--------------------------------------------------
   if (solid_regions.size() > 1) {
      for (std::size_t it_region = 1; it_region < solid_regions.size(); it_region++) {
         // (length of a non-solid region < kmer_length) && (length of a non-solid region >= kmer_length - FP_SUSPECT_LENGTH(default: 1))
         if (((solid_regions[it_region].first - solid_regions[it_region - 1].second - 1) < kmer_length) &&
             ((solid_regions[it_region].first - solid_regions[it_region - 1].second - 1) >= kmer_length - FP_SUSPECT_LENGTH)) {
            // length of the right solid region > FP_SUSPECT_LENGTH(default: 1)
            if ((solid_regions[it_region].second - solid_regions[it_region].first + 1) > FP_SUSPECT_LENGTH) {
               solid_regions[it_region].first += FP_SUSPECT_LENGTH;
            }

            // length of the left solid region > FP_SUSPECT_LENGTH(default: 1)
            if ((solid_regions[it_region - 1].second - solid_regions[it_region - 1].first + 1) > FP_SUSPECT_LENGTH) {
               solid_regions[it_region - 1].second -= FP_SUSPECT_LENGTH;
            }
         }
      }
   }

   //--------------------------------------------------
   // STEP 0-5: remove a solid region that makes a non-solid reiong shorter than k
   //--------------------------------------------------
   if (solid_regions.size() == 2) {
      // the first solid region starts from the first k-mer
      if (solid_regions[0].first == 0) {
         // the distance between two regions is shorter than k
         if ((solid_regions[1].first - solid_regions[0].second) < (kmer_length + 1)) {
            // remove the second solid region
            solid_regions.erase(solid_regions.begin() + 1);
         }
      }
      // the second solid region ends in the last k-mer
      else if (solid_regions[1].second == (sequence.length() - kmer_length)) {
         // the distance between two regions is shorter than k
         if ((solid_regions[1].first - solid_regions[0].second) < (kmer_length + 1)) {
            // the length of the second solid region is >= 10% of the sequence length
            if ((solid_regions[1].second - solid_regions[1].first + 1) >= (sequence.length() * 0.1)) {
               // the length of the first solid region is < 10% of the sequence length
               if ((solid_regions[0].second - solid_regions[0].first + 1) < (sequence.length() * 0.1)) {
                  // remove the second solid region
                  solid_regions.erase(solid_regions.begin());
               }
            }
         }
      }
   }

   //--------------------------------------------------
   // STEP 0-6: check the quality scores of right side of each solid k-mer region
   //--------------------------------------------------
   // at least one solid region
   if (solid_regions.size() > 0) {
      // 1 - (n - 1) solid region
      for (std::size_t it_sr = 0; it_sr < (solid_regions.size() - 1); it_sr++) {
         // sufficient solid regions length
         if ((solid_regions[it_sr].second - solid_regions[it_sr].first) > SOLID_REGION_ADJUST_RANGE) {
            for (std::size_t it_adjust = solid_regions[it_sr].second; it_adjust > (solid_regions[it_sr].second - SOLID_REGION_ADJUST_RANGE); it_adjust--) {
               // low quality score
               if ((((unsigned short int)quality_score[it_adjust + kmer_length - 1] - quality_score_offset) < quality_score_cutoff) ||
                   (((unsigned short int)quality_score[it_adjust] - quality_score_offset) < quality_score_cutoff)
                  ) {
                  solid_regions[it_sr].second = it_adjust - 1;
                  break;
               }
            }
         }
      }

      // last solid region
      std::size_t index_solid_region(solid_regions.size() - 1);

      // non-solid k-mers exist at the 3-prime end
      if (solid_regions[index_solid_region].second < (sequence.length() - kmer_length)) {
         // sufficient solid regions length
         if ((solid_regions[index_solid_region].second - solid_regions[index_solid_region].first) > SOLID_REGION_ADJUST_RANGE) {
            for (std::size_t it_adjust = solid_regions[index_solid_region].second; it_adjust > (solid_regions[index_solid_region].second - SOLID_REGION_ADJUST_RANGE); it_adjust--) {
               // low quality score
               if (((unsigned short int)quality_score[it_adjust + kmer_length - 1] - quality_score_offset) < quality_score_cutoff) {
                  solid_regions[index_solid_region].second = it_adjust - 1;
                  break;
               }
            }
         }
      }

      // non-solid k-mers exist at the 5-prime end
      if (solid_regions[0].first > 0) {
         // sufficient solid regions length
         if ((solid_regions[0].second - solid_regions[0].first) > SOLID_REGION_ADJUST_RANGE) {
            for (std::size_t it_adjust = solid_regions[0].first; it_adjust < (solid_regions[0].first + SOLID_REGION_ADJUST_RANGE); it_adjust++) {
               // low quality score
               if (((unsigned short int)quality_score[it_adjust + kmer_length - 1] - quality_score_offset) < quality_score_cutoff) {
                  solid_regions[0].first = it_adjust + 1;
                  break;
               }
            }
         }
      }
   }

   //--------------------------------------------------
   // STEP 0-7: check whether a non-solid region < k still exists
   //--------------------------------------------------
   bool short_non_solid_region(false);
   if (solid_regions.size() > 1) {
      for (std::size_t it_sr = 1; it_sr < (solid_regions.size() - 1); it_sr++) {
         if ((solid_regions[it_sr].first - solid_regions[it_sr - 1].second) <= kmer_length) {
            short_non_solid_region = true;
            break;
         }
      }
   }

   //--------------------------------------------------
   // correct errors
   //--------------------------------------------------
   std::string sequence_modified(sequence);

   if ((solid_regions.size() > 0) && (short_non_solid_region == false)) {
      //--------------------------------------------------
      // STEP 1-1: Correct errors between solid regions
      //--------------------------------------------------
      if (solid_regions.size() > 1) {
         // for each solid region
         for (std::size_t it_region = 1; it_region < solid_regions.size(); it_region++) {
            if ((((solid_regions[it_region].first - 1) - (solid_regions[it_region - 1].second + 1)) + 1) >= kmer_length) {
               // the bases that may be modified: from (solid_regions[it_region - 1].second + kmer_length) to (solid_regions[it_region].first - 1)
               // they should not be in trimmed regions
               if (((solid_regions[it_region - 1].second + kmer_length) < (read_length - trim_3_end)) && ((solid_regions[it_region].first) > trim_5_end)) {
                  correct_errors_between_solid_regions(
                                                       sequence,
                                                       sequence_modified,
                                                       quality_score,
                                                       solid_regions[it_region - 1].first,
                                                       solid_regions[it_region - 1].second + 1,
                                                       solid_regions[it_region].first - 1,
                                                       solid_regions[it_region].second,
                                                       solid_regions_org[it_region - 1].second - 1,
                                                       solid_regions_org[it_region].first - 1,
                                                       sequence_modification,
                                                       trim_5_end,
                                                       trim_3_end,
                                                       solid_regions.size(),
                                                       read_length,
                                                       max_trimmed_bases,
                                                       min_check_length,
                                                       num_corrected_errors_local,
                                                       bit_vector,
                                                       hash_seed
                                                      );
               }
            }
         }
      }

      //--------------------------------------------------
      // STEP 1-2: Correct errors in the 5' end
      //--------------------------------------------------
      // number of solid regions is >= 1
      if (solid_regions.size() >= 1) {
         // the first solid region does not start from the 0-th k-mer in a read
         if (solid_regions[0].first > 0) {
            // the bases that may be modified: from 0 to (solid_regions[0].first - 1)
            // they should not be in trimmed regions
            if (solid_regions[0].first > trim_5_end) {
               correct_errors_5_prime_end(
                                          sequence,
                                          sequence_modified,
                                          quality_score,
                                          solid_regions[0].first - 1,
                                          sequence_modification,
                                          trim_5_end,
                                          trim_3_end,
                                          solid_regions_org[0].first - 1,
                                          read_length,
                                          max_trimmed_bases,
                                          min_check_length,
                                          num_corrected_errors_local,
                                          bit_vector,
                                          hash_seed
                                         );
            }
         }
      }

      //--------------------------------------------------
      // STEP 1-3: Correct errors in the 3' end
      //--------------------------------------------------
      // number of solid regions is >= 1
      if (solid_regions.size() >= 1) {
         // the last solid region does not end in the last k-mer in a read
         if (solid_regions[solid_regions.size() - 1].second < (sequence.length() - kmer_length)) {
            // the bases that may be modified: from (solid_regions[solid_regions.size() - 1].second + kmer_length) to (read_length - 1)
            // they should not be in trimmed regions
            if ((solid_regions[solid_regions.size() - 1].second + kmer_length) < (read_length - trim_3_end)) {
               correct_errors_3_prime_end(
                                          sequence,
                                          sequence_modified,
                                          quality_score,
                                          solid_regions[solid_regions.size() - 1].second + 1,
                                          sequence_modification,
                                          trim_5_end,
                                          trim_3_end,
                                          solid_regions_org[solid_regions.size() - 1].second + 1,
                                          read_length,
                                          max_trimmed_bases,
                                          min_check_length,
                                          num_corrected_errors_local,
                                          bit_vector,
                                          hash_seed
                                         );
            }
         }
      }
   }
   //--------------------------------------------------
   // no solid region or short weak regions
   //--------------------------------------------------
   else {
      //--------------------------------------------------
      // STEP 2-1: Correct errors in the first k-mer
      //--------------------------------------------------
      // find potentially wrong bases
      std::vector<C_candidate_path> candidate_path_vector_tmp;

      correct_errors_first_kmer(
                                sequence_modified,
                                quality_score,
                                sequence_modification,
                                candidate_path_vector_tmp,
                                bit_vector,
                                hash_seed
                               );

      // candidiate_path_vector_tmp: differently modified versions of the first k-mer

      // filter some candidates by extending the first k-mer to the left
      std::vector<C_candidate_path> candidate_path_vector_tmp_tmp;

      if (candidate_path_vector_tmp.size() > 0) {
         // each path
         for (std::size_t it_candidates = 0; it_candidates < candidate_path_vector_tmp.size(); it_candidates++) {
            // no modified base
            // if the original first k-mer is added to the candidate vector
            if (candidate_path_vector_tmp[it_candidates].modified_bases.size() == 0) {
               candidate_path_vector_tmp_tmp.push_back(candidate_path_vector_tmp[it_candidates]);
            }
            // there exist modified bases
            else {
               bool really_modified(false);
               for (unsigned int it_mod_base = 0; it_mod_base < candidate_path_vector_tmp[it_candidates].modified_bases.size(); it_mod_base++) {
                  if (sequence[candidate_path_vector_tmp[it_candidates].modified_bases[it_mod_base].first] != candidate_path_vector_tmp[it_candidates].modified_bases[it_mod_base].second) {
                     really_modified = true;
                  }
               }

               if (really_modified == true) {
                  // check the index of the first modified base
                  // extension is needed
                  if (candidate_path_vector_tmp[it_candidates].modified_bases[0].first < (kmer_length - 1)) {
                     bool extension_success(false);
                     solid_first_kmer(
                                      candidate_path_vector_tmp[it_candidates],
                                      sequence_modified,
                                      extension_success,
                                      bit_vector,
                                      hash_seed
                                     );

                     if (extension_success == true) {
                        candidate_path_vector_tmp_tmp.push_back(candidate_path_vector_tmp[it_candidates]);
                     }
                  }
                  // extension is not needed
                  else {
                     candidate_path_vector_tmp_tmp.push_back(candidate_path_vector_tmp[it_candidates]);
                  }
               }
            }
         }
      }

      // candidiate_path_vector_tmp_tmp: solid k-mers in candidate_path_vector_tmp

      // candidates in candidiate_path_vector_tmp_tmp are moved to candidate_path_vector_tmp
      candidate_path_vector_tmp = candidate_path_vector_tmp_tmp;
      candidate_path_vector_tmp_tmp.clear();

      //--------------------------------------------------
      // STEP 2-2: extend candidate paths to the right
      //--------------------------------------------------
      if (candidate_path_vector_tmp.size() > 0) {
         // each path
         run_exploration = true;

         for (std::size_t it_candidates = 0; it_candidates < candidate_path_vector_tmp.size(); it_candidates++) {
            if (run_exploration == true) {
               extend_first_kmer_to_right(
                                          sequence_modified,
                                          quality_score,
                                          candidate_path_vector_tmp[it_candidates],
                                          candidate_path_vector_tmp_tmp,
                                          read_length,
                                          bit_vector,
                                          hash_seed
                                         );
            }
         }

         // complete exploration was not done because there are too many candidata paths
         // remove all the paths in candidate_path_vector_tmp
         if (run_exploration == false) {
            candidate_path_vector_tmp_tmp.clear();
         }
      }

      // candidiate_path_vector_tmp_tmp: successfully right extended candidates

      // remain only really modified paths
      candidate_path_vector_tmp.clear();
      // each path
      for (std::size_t it_candidate = 0; it_candidate < candidate_path_vector_tmp_tmp.size(); it_candidate++) {
         // each modification
         bool really_modified(false);
         for (unsigned int it_mod_base = 0; it_mod_base < candidate_path_vector_tmp_tmp[it_candidate].modified_bases.size(); it_mod_base++) {
            if (sequence[candidate_path_vector_tmp_tmp[it_candidate].modified_bases[it_mod_base].first] != candidate_path_vector_tmp_tmp[it_candidate].modified_bases[it_mod_base].second) {
               really_modified = true;
            }
         }

         if (really_modified) {
            candidate_path_vector_tmp.push_back(candidate_path_vector_tmp_tmp[it_candidate]);
         }
      }

      candidate_path_vector_tmp_tmp = candidate_path_vector_tmp;

      //--------------------------------------------------
      // STEP 2-3: choose a final one in candidate_path_vector_tmp_tmp if possible
      //--------------------------------------------------
      // compare quality scores of candidate paths
      // if the number of paths in candidate_path_vector_tmp_tmp is larger than 1
      if (candidate_path_vector_tmp_tmp.size() > 1) {
         // each path
         std::vector<C_candidate_path>::iterator it_path;
         std::vector<C_candidate_path>::iterator it_path_1st;
         std::vector<C_candidate_path>::iterator it_path_2nd;

         std::size_t qs_1st(INIT_MIN_QS);
         std::size_t qs_2nd(INIT_MIN_QS);

         // find the first and second paths
         // each candidate path
         for (it_path = candidate_path_vector_tmp_tmp.begin(); it_path != candidate_path_vector_tmp_tmp.end(); it_path++) {
            // each modification
            for (std::size_t it_mod = 0; it_mod < (*it_path).modified_bases.size(); it_mod++) {
               // add quality scores of modified bases
               if (sequence_modified[(*it_path).modified_bases[it_mod].first] != (*it_path).modified_bases[it_mod].second) {
                  (*it_path).sum_qs += ((unsigned short int)quality_score[(*it_path).modified_bases[it_mod].first] - quality_score_offset);
               }
            }

            // compare quality scores of each path
            if ((*it_path).sum_qs <= qs_1st) {
               qs_2nd = qs_1st;
               qs_1st = (*it_path).sum_qs;

               it_path_2nd = it_path_1st;
               it_path_1st = it_path;
            }
            else if ((*it_path).sum_qs <= qs_2nd) {
               qs_2nd = (*it_path).sum_qs;

               it_path_2nd = it_path;
            }
         }

         // check whether too many bases are modified in the first path, which implies indel may exist in the original read
         // if too many modifications exist, the read is just trimmed
         //bool keep_going(true);
         bool too_many_corrections(false);

         // at least over MAX_MODIFICATION times of modifications from the first modified position
         if (((read_length - (*it_path_1st).modified_bases[0].first) >= min_check_length) &&
            ((*it_path_1st).modified_bases.size() > MAX_MODIFICATION)) {
            // each modified base in a right range
            //std::size_t partial_qs((*it_path_1st).sum_qs);

            for (unsigned int it_mod_base = 0; it_mod_base < ((*it_path_1st).modified_bases.size() - MAX_MODIFICATION); it_mod_base++) {
               // at least MAX_MODIFICATION times of modifications within min_check_length bases
               if (((*it_path_1st).modified_bases[it_mod_base + MAX_MODIFICATION].first - (*it_path_1st).modified_bases[it_mod_base].first) < min_check_length) {
                  //if (it_mod_base > 0) {
                  //   partial_qs -= ((unsigned short int)quality_score[(*it_path_1st).modified_bases[it_mod_base - 1].first] - quality_score_offset);
                  //}

                  // average quality score of the modified bases is too high
                  //if (1.0 * partial_qs / ((*it_path_1st).modified_bases.size() - it_mod_base) > quality_score_cutoff) {
                     // trim_5_end or trim_3_end
                     if ((read_length - (*it_path_1st).modified_bases[it_mod_base].first) <= max_trimmed_bases) {
                        trim_3_end = read_length - (*it_path_1st).modified_bases[it_mod_base].first;
                     }
                     else if (((*it_path_1st).modified_bases[it_mod_base].first + 1) <= max_trimmed_bases) {
                        trim_5_end = (*it_path_1st).modified_bases[it_mod_base].first + 1;
                     }

                     too_many_corrections = true;
                     break;
                  //}
               }
            }
         }

         // use the 1st path if not too many corrections are made AND if the 1st path has a sufficiently low score
         // if an indel exists using quality scores is not a good way to choose the best path
         if ((too_many_corrections == false) && (qs_1st + MIN_QS_DIFF <= qs_2nd)) {
            for (unsigned int it_base = 0; it_base < (*it_path_1st).modified_bases.size(); it_base++) {
               // filter out the bases that are equal to the original ones
               if (sequence_modified[(*it_path_1st).modified_bases[it_base].first] != (*it_path_1st).modified_bases[it_base].second) {
                  sequence_modification[(*it_path_1st).modified_bases[it_base].first] = (*it_path_1st).modified_bases[it_base].second;
                  sequence_modified[(*it_path_1st).modified_bases[it_base].first] = (*it_path_1st).modified_bases[it_base].second;

                  num_corrected_errors_local++;
               }
            }
         }
         // hard to choose one path
         else {
            // A AND B
            std::vector< std::pair<unsigned int, char> > base_vector_intersection;
            // A OR B
            std::vector< std::pair<unsigned int, char> > base_vector_union;

            // temporary vectors
            std::vector< std::pair<unsigned int, char> > base_vector_intersection_prev(candidate_path_vector_tmp_tmp[0].modified_bases);
            std::vector< std::pair<unsigned int, char> > base_vector_union_prev(candidate_path_vector_tmp_tmp[0].modified_bases);

            // each candidate path
            for (std::size_t it_p = 1; it_p < candidate_path_vector_tmp_tmp.size(); it_p++) {
               base_vector_intersection.clear();
               base_vector_union.clear();

               base_intersection(base_vector_intersection_prev.begin(), base_vector_intersection_prev.end(), candidate_path_vector_tmp_tmp[it_p].modified_bases.begin(), candidate_path_vector_tmp_tmp[it_p].modified_bases.end(), base_vector_intersection);
               base_union       (base_vector_union_prev.begin(),        base_vector_union_prev.end(),        candidate_path_vector_tmp_tmp[it_p].modified_bases.begin(), candidate_path_vector_tmp_tmp[it_p].modified_bases.end(), base_vector_union);

               base_vector_intersection_prev = base_vector_intersection;
               base_vector_union_prev        = base_vector_union;
            }

            // A - B
            std::vector< std::pair<unsigned int, char> > base_vector_difference;
            base_difference(base_vector_union.begin(), base_vector_union.end(), base_vector_intersection.begin(),base_vector_intersection.end(), base_vector_difference);

            // find trimmed region
            // correcting the 5'-end and 3'-end was not done
            // therefore the total number of trimmed bases is 0 yet
            if (base_vector_difference.size() > 0) {
               std::size_t vector_index_leftmost(0);
               std::size_t vector_index_rightmost(base_vector_difference.size() - 1);

               bool keep_going(true);

               while (keep_going) {
                  // # of trimmed bases at the 5' end: base_vector_difference[vector_index_leftmost].first + 1
                  // # of trimmed bases at the 3' end: read_length - base_vector_difference[vector_index_rightmost].first
                  // the 5'-end is smaller
                  if ((base_vector_difference[vector_index_leftmost].first + 1) < (read_length - base_vector_difference[vector_index_rightmost].first)) {
                     // check the total number of trimmed bases
                     if ((base_vector_difference[vector_index_leftmost].first + 1 + trim_3_end) <= max_trimmed_bases) {
                        trim_5_end = base_vector_difference[vector_index_leftmost].first + 1;

                        // two points are met
                        if (vector_index_leftmost == vector_index_rightmost) {
                           keep_going = false;
                        }
                        else {
                           vector_index_leftmost++;
                        }
                     }
                     // no need for more check
                     else {
                        keep_going = false;
                     }
                  }
                  // the 3'-end is smaller
                  else {
                     // check the total number of trimmed bases
                     if ((read_length - base_vector_difference[vector_index_rightmost].first) <= max_trimmed_bases) {
                        trim_3_end = read_length - base_vector_difference[vector_index_rightmost].first;

                        // two points are met
                        if (vector_index_leftmost == vector_index_rightmost) {
                           keep_going = false;
                        }
                        else {
                           vector_index_rightmost--;
                        }
                     }
                     // no need for more check
                     else {
                        keep_going = false;
                     }
                  }
               }
            }

            // find consensus modifications
            for (std::size_t it_inter = 0; it_inter < base_vector_intersection.size(); it_inter++) {
               // check whether the base is not in the trimmed regions
               if ((base_vector_intersection[it_inter].first <  (read_length - trim_3_end)) &&
                   (base_vector_intersection[it_inter].first >= trim_5_end)) {
                  // filter out the bases that are equal to the original ones
                  if (sequence_modified[base_vector_intersection[it_inter].first] != base_vector_intersection[it_inter].second) {
                     sequence_modification[base_vector_intersection[it_inter].first] = base_vector_intersection[it_inter].second;
                     sequence_modified[base_vector_intersection[it_inter].first] = base_vector_intersection[it_inter].second;

                     num_corrected_errors_local++;
                  }
               }
            }
         }
      }
      // only one path
      // correction succeeds
      else if (candidate_path_vector_tmp_tmp.size() == 1) {
         // check whether too many bases are modified in the first path, which implies indel may exist in the original read
         // if too many modifications exist, the read is just trimmed
         //bool keep_going(true);
         bool too_many_corrections(false);

         // at least over MAX_MODIFICATION times of modifications from the first modified position
         if (((read_length - candidate_path_vector_tmp_tmp[0].modified_bases[0].first) >= min_check_length) &&
            (candidate_path_vector_tmp_tmp[0].modified_bases.size() > MAX_MODIFICATION)) {
            // calculate sum_qs
            //for (unsigned int it_mod_base = 0; it_mod_base < candidate_path_vector_tmp_tmp[0].modified_bases.size(); it_mod_base++) {
            //   candidate_path_vector_tmp_tmp[0].sum_qs += ((unsigned short int)quality_score[candidate_path_vector_tmp_tmp[0].modified_bases[it_mod_base].first] - quality_score_offset);
            //}

            // each modified base in a right range
            //std::size_t partial_qs(candidate_path_vector_tmp_tmp[0].sum_qs);

            for (unsigned int it_mod_base = 0; it_mod_base < (candidate_path_vector_tmp_tmp[0].modified_bases.size() - MAX_MODIFICATION); it_mod_base++) {
               // at least MAX_MODIFICATION times of modifications within min_check_length bases
               if ((candidate_path_vector_tmp_tmp[0].modified_bases[it_mod_base + MAX_MODIFICATION].first - candidate_path_vector_tmp_tmp[0].modified_bases[it_mod_base].first) < min_check_length) {
                  //if (it_mod_base > 0) {
                  //   partial_qs -= ((unsigned short int)quality_score[candidate_path_vector_tmp_tmp[0].modified_bases[it_mod_base - 1].first] - quality_score_offset);
                  //}

                  // average quality score of the modified bases is too high
                  //if (1.0 * partial_qs / (candidate_path_vector_tmp_tmp[0].modified_bases.size() - it_mod_base) > quality_score_cutoff) {
                     // trim_5_end or trim_3_end
                     if ((read_length - candidate_path_vector_tmp_tmp[0].modified_bases[it_mod_base].first) <= max_trimmed_bases) {
                        trim_3_end = read_length - candidate_path_vector_tmp_tmp[0].modified_bases[it_mod_base].first;

                        // update sequence_modification for the non-trimmed corrections
                        if (it_mod_base > 0) {
                           for (unsigned int it_base = 0; it_base < (it_mod_base - 1); it_base++) {
                              sequence_modification[(candidate_path_vector_tmp_tmp[0]).modified_bases[it_base].first] = (candidate_path_vector_tmp_tmp[0]).modified_bases[it_base].second;
                              sequence_modified[(candidate_path_vector_tmp_tmp[0]).modified_bases[it_base].first] = (candidate_path_vector_tmp_tmp[0]).modified_bases[it_base].second;
                              num_corrected_errors_local++;
                           }
                        }
                     }
                     else if ((candidate_path_vector_tmp_tmp[0].modified_bases[it_mod_base].first + 1) <= max_trimmed_bases) {
                        trim_5_end = candidate_path_vector_tmp_tmp[0].modified_bases[it_mod_base].first + 1;

                        // update sequence_modification for the non-trimmed corrections
                        if (it_mod_base < (candidate_path_vector_tmp_tmp[0].modified_bases.size() - 1)) {
                           for (unsigned int it_base = (it_mod_base + 1); it_base < candidate_path_vector_tmp_tmp[0].modified_bases.size(); it_base++) {
                              sequence_modification[(candidate_path_vector_tmp_tmp[0]).modified_bases[it_base].first] = (candidate_path_vector_tmp_tmp[0]).modified_bases[it_base].second;
                              sequence_modified[(candidate_path_vector_tmp_tmp[0]).modified_bases[it_base].first] = (candidate_path_vector_tmp_tmp[0]).modified_bases[it_base].second;
                              num_corrected_errors_local++;
                           }
                        }
                     }

                     too_many_corrections = true;
                     break;
                  //}
               }
            }
         }

         // not too many modified bases
         if (too_many_corrections == false) {
            for (unsigned int it_base = 0; it_base < candidate_path_vector_tmp_tmp[0].modified_bases.size(); it_base++) {
               // filter out the bases that are equal to the original ones
               if (sequence_modified[candidate_path_vector_tmp_tmp[0].modified_bases[it_base].first] != candidate_path_vector_tmp_tmp[0].modified_bases[it_base].second) {
                  sequence_modification[candidate_path_vector_tmp_tmp[0].modified_bases[it_base].first] = candidate_path_vector_tmp_tmp[0].modified_bases[it_base].second;
                  sequence_modified[candidate_path_vector_tmp_tmp[0].modified_bases[it_base].first] = candidate_path_vector_tmp_tmp[0].modified_bases[it_base].second;

                  num_corrected_errors_local++;
               }
            }
         }
      }
   }
}



//----------------------------------------------------------------------
// correct_errors_between_solid_regions
//----------------------------------------------------------------------
inline void C_correct_errors::correct_errors_between_solid_regions(const std::string& org_sequence, std::string& sequence, const std::string& quality_score, const std::size_t& left_first, const std::size_t& index_start, const std::size_t& index_end, const std::size_t& right_second, const std::size_t& org_boundary_left, const std::size_t& org_boundary_right, std::string& sequence_modification, std::size_t& trim_5_end, std::size_t& trim_3_end, const std::size_t& num_solid_islands, const std::size_t& read_length, const std::size_t& max_trimmed_bases, const std::size_t& min_check_length, std::size_t& num_corrected_errors_local, const unsigned char* bit_vector, const std::vector<unsigned int>& hash_seed) {
   //--------------------------------------------------
   // from i-th region to (i + 1)-th region
   //--------------------------------------------------
   // i-th solid region | non-solid | (i+1)-th solid region
   // --------------------------------------- read
   //              |-----|                    (index_start)-th k-mer
   //                              |-----|    (index_end)-th k-mer: last non-solid k-mer
   //                         |-----|         (index_last_mod)-th k-mer = (index_end - kmer_length + 1)-th k-mer: last k-mer that can be modified
   //                               |----|    This regions should not be modified
   //--------------------------------------------------
   // list of candidate paths
   std::vector<C_candidate_path> candidate_path_vector_tmp;

   // index of the k-mer that can be modified
   // k-mers that are overlapped with a solid regioin cannot be modified
   std::size_t index_last_mod(index_end - kmer_length + 1);

   // make an initial k-mer
   std::string kmer_initial(sequence.substr(index_start, kmer_length));

   // each alternative neocletide
   run_exploration = true;

   for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
      // make a change
      kmer_initial[kmer_length - 1] = NEOCLEOTIDE[it_alter];

      // kmer_initial is solid
      if (query_text(kmer_initial, bit_vector, hash_seed) == true) {
         // generate a new path
         C_candidate_path candidate_path;

         if (sequence[index_start + kmer_length - 1] != NEOCLEOTIDE[it_alter]) {
            std::pair<unsigned int, char> pair_tmp;
            pair_tmp.first  = index_start + kmer_length - 1;
            pair_tmp.second = NEOCLEOTIDE[it_alter];

            candidate_path.modified_bases.push_back(pair_tmp);
         }

         // if this k-mer is the last k-mer that can be modified
         // running extend_a_kmer_right is not needed any more
         if (index_start == index_last_mod) {
            candidate_path_vector_tmp.push_back(candidate_path);
         }
         else {
            // trace  this kmer recursively and update candidate_path_vector_tmp
            extend_a_kmer(
                          kmer_initial,
                          sequence,
                          index_start,
                          index_last_mod,
                          candidate_path,
                          candidate_path_vector_tmp,
                          org_boundary_left,
                          org_boundary_right,
                          quality_score,
                          bit_vector,
                          hash_seed
                         );
         }
      }
   }

   // complete exploration was not done because there are too many candidata paths
   // remove all the paths in candidate_path_vector_tmp
   if (run_exploration == false) {
      candidate_path_vector_tmp.clear();
   }

   std::vector<C_candidate_path> candidate_path_vector;

   // check the solidness of k-mers between index_last_mod and index_end
   bool all_solid_wo_modification(false);

   // each candidate path
   for (std::vector<C_candidate_path>::iterator it_path = candidate_path_vector_tmp.begin(); it_path != candidate_path_vector_tmp.end(); it_path++) {
      if ((*it_path).modified_bases.size() == 0) {
         all_solid_wo_modification = true;
         break;
      }
      else {
         // checking is needed
         std::size_t index_last_modified_base((*it_path).modified_bases[(*it_path).modified_bases.size() - 1].first);

         if (index_last_modified_base > index_last_mod) {
            // generate a temporary sequence
            std::string sequence_tmp(sequence);
            for (unsigned int it_base = 0; it_base < (*it_path).modified_bases.size(); it_base++) {
               sequence_tmp[(*it_path).modified_bases[it_base].first] = (*it_path).modified_bases[it_base].second;
            }

            // check k-mers
            std::size_t num_success(0);
            std::size_t num_fails(0);
            for (std::size_t it_check = index_last_mod; it_check <= index_last_modified_base; it_check++) {
               std::string kmer_current(sequence_tmp.substr(it_check, kmer_length));

               if (query_text(kmer_current, bit_vector, hash_seed) == true) {
                  num_success++;
               }
               else {
                  num_fails++;

                  if (num_fails > NUM_ALLOWABLE_FAILS) {
                     break;
                  }
               }
            }

            if (num_success >= (index_last_modified_base - index_last_mod + 1 - NUM_ALLOWABLE_FAILS)) {
               candidate_path_vector.push_back(*it_path);
            }
         }
         // checking is not needed
         else {
            candidate_path_vector.push_back(*it_path);
         }
      }
   }

   // remain only really modified paths
   candidate_path_vector_tmp.clear();
   // each path
   for (std::size_t it_candidate = 0; it_candidate < candidate_path_vector.size(); it_candidate++) {
      // each modification
      bool really_modified(false);
      for (unsigned int it_mod_base = 0; it_mod_base < candidate_path_vector[it_candidate].modified_bases.size(); it_mod_base++) {
         if (org_sequence[candidate_path_vector[it_candidate].modified_bases[it_mod_base].first] != candidate_path_vector[it_candidate].modified_bases[it_mod_base].second) {
            really_modified = true;
         }
      }

      if (really_modified) {
         candidate_path_vector_tmp.push_back(candidate_path_vector[it_candidate]);
      }
   }

   candidate_path_vector = candidate_path_vector_tmp;

   // all k-mers are solid without any modification
   if (all_solid_wo_modification == true) {
   // do nothing
   }
   // compare quality scores of candidate paths
   // if the number of paths in candidate_path_vector is larger than 1
   else if (candidate_path_vector.size() > 1) {
      // each path
      std::vector<C_candidate_path>::iterator it_path;
      std::vector<C_candidate_path>::iterator it_path_1st;
      std::vector<C_candidate_path>::iterator it_path_2nd;

      std::size_t qs_1st(INIT_MIN_QS);
      std::size_t qs_2nd(INIT_MIN_QS);

      for (it_path = candidate_path_vector.begin(); it_path != candidate_path_vector.end(); it_path++) {
         // each modification
         for (std::vector< std::pair<unsigned int, char> >::iterator it_mod = (*it_path).modified_bases.begin(); it_mod != (*it_path).modified_bases.end(); it_mod++) {
            // add quality scores of modified bases
            (*it_path).sum_qs += ((unsigned short int)quality_score[(*it_mod).first] - quality_score_offset);
         }

         // compare quality scores of each path
         if ((*it_path).sum_qs <= qs_1st) {
            qs_2nd = qs_1st;
            qs_1st = (*it_path).sum_qs;

            it_path_2nd = it_path_1st;
            it_path_1st = it_path;
         }
         else if ((*it_path).sum_qs <= qs_2nd) {
            qs_2nd = (*it_path).sum_qs;

            it_path_2nd = it_path;
         }
      }

      // check whether too many bases are modified in the first path, which implies indel may exist in the original read
      // if too many modifications exist, the read is just trimmed
      bool too_many_corrections(false);

      // this is done only when the number of solid islands is two
      if (num_solid_islands == 2) {
         // the 1st island is small && the second island is big
         if (((index_start - left_first) > MAX_MODIFICATION) && ((right_second - index_end + 2) <= MAX_MODIFICATION)) {
            // at least over MAX_MODIFICATION times of modifications from the first modified position
            if ((((*it_path_1st).modified_bases[(*it_path_1st).modified_bases.size() - 1].first -kmer_length - index_start + 2) >= min_check_length) &&
                ((*it_path_1st).modified_bases.size() > MAX_MODIFICATION)) {
               // each modified base in a right range
               //std::size_t partial_qs((*it_path_1st).sum_qs);

               for (std::size_t it_mod_base = ((*it_path_1st).modified_bases.size() - 1); it_mod_base >= MAX_MODIFICATION; it_mod_base--) {
                  // at least MAX_MODIFICATION times of modifications within min_check_length bases
                  if (((*it_path_1st).modified_bases[it_mod_base].first - (*it_path_1st).modified_bases[it_mod_base - MAX_MODIFICATION].first) < min_check_length) {
                     //if (it_mod_base < ((*it_path_1st).modified_bases.size() - 1)) {
                     //   partial_qs -= ((unsigned short int)quality_score[(*it_path_1st).modified_bases[it_mod_base + 1].first] - quality_score_offset);
                     //}

                     // average quality score of the modified bases is too high
                     //if (1.0 * partial_qs / (it_mod_base + 1) > quality_score_cutoff) {
                        // trim_5_end or trim_3_end
                        if ((read_length - (*it_path_1st).modified_bases[it_mod_base].first) <= max_trimmed_bases) {
                           trim_3_end = read_length - (*it_path_1st).modified_bases[it_mod_base].first;
                        }
                        else if (((*it_path_1st).modified_bases[it_mod_base].first + 1) <= max_trimmed_bases) {
                           trim_5_end = (*it_path_1st).modified_bases[it_mod_base].first + 1;
                        }

                        too_many_corrections = true;
                        break;
                     //}
                  }
               }
            }
         }
         // the 1st island is big && the second island is small
         else if (((index_start - left_first) <= MAX_MODIFICATION) && ((right_second - index_end + 2) > MAX_MODIFICATION)) {
            // at least over MAX_MODIFICATION times of modifications from the first modified position
            if (((index_end - (*it_path_1st).modified_bases[0].first + 1) >= min_check_length) &&
               ((*it_path_1st).modified_bases.size() > MAX_MODIFICATION)) {
               // each modified base in a right range
               //std::size_t partial_qs((*it_path_1st).sum_qs);

               for (unsigned int it_mod_base = 0; it_mod_base < ((*it_path_1st).modified_bases.size() - MAX_MODIFICATION); it_mod_base++) {
                  // at least MAX_MODIFICATION times of modifications within min_check_length bases
                  if (((*it_path_1st).modified_bases[it_mod_base + MAX_MODIFICATION].first - (*it_path_1st).modified_bases[it_mod_base].first) < min_check_length) {
                     //if (it_mod_base > 0) {
                     //   partial_qs -= ((unsigned short int)quality_score[(*it_path_1st).modified_bases[it_mod_base - 1].first] - quality_score_offset);
                     //}

                     // average quality score of the modified bases is too high
                     //if (1.0 * partial_qs / ((*it_path_1st).modified_bases.size() - it_mod_base) > quality_score_cutoff) {
                        // trim_5_end or trim_3_end
                        if ((read_length - (*it_path_1st).modified_bases[it_mod_base].first) <= max_trimmed_bases) {
                           trim_3_end = read_length - (*it_path_1st).modified_bases[it_mod_base].first;
                        }
                        else if (((*it_path_1st).modified_bases[it_mod_base].first + 1) <= max_trimmed_bases) {
                           trim_5_end = (*it_path_1st).modified_bases[it_mod_base].first + 1;
                        }

                        too_many_corrections = true;
                        break;
                     //}
                  }
               }
            }
         }
      }

      // use the 1st path if not too many corrections are made AND if the 1st path has a sufficiently low score
      // if an indel exists using quality scores is not a good way to choose the best path
      if ((too_many_corrections == false) && (qs_1st + MIN_QS_DIFF <= qs_2nd)) {
         // each modification
         for (std::vector< std::pair<unsigned int, char> >::iterator it_base = (*it_path_1st).modified_bases.begin(); it_base != (*it_path_1st).modified_bases.end(); it_base++) {
            // update sequence_modification
            sequence_modification[(*it_base).first] = (*it_base).second;
            sequence[(*it_base).first] = (*it_base).second;
            num_corrected_errors_local++;
         }
      }
      // hard to choose one path
      else {
         // A AND B
         std::vector< std::pair<unsigned int, char> > base_vector_intersection;
         // A OR B
         std::vector< std::pair<unsigned int, char> > base_vector_union;

         // temporary vectors
         std::vector< std::pair<unsigned int, char> > base_vector_intersection_prev(candidate_path_vector[0].modified_bases);
         std::vector< std::pair<unsigned int, char> > base_vector_union_prev(candidate_path_vector[0].modified_bases);

         // each candidate path
         for (std::size_t it_p = 1; it_p < candidate_path_vector.size(); it_p++) {
            base_vector_intersection.clear();
            base_vector_union.clear();

            base_intersection(base_vector_intersection_prev.begin(), base_vector_intersection_prev.end(), candidate_path_vector[it_p].modified_bases.begin(), candidate_path_vector[it_p].modified_bases.end(), base_vector_intersection);
            base_union       (base_vector_union_prev.begin(),        base_vector_union_prev.end(),        candidate_path_vector[it_p].modified_bases.begin(), candidate_path_vector[it_p].modified_bases.end(), base_vector_union);

            base_vector_intersection_prev = base_vector_intersection;
            base_vector_union_prev        = base_vector_union;
         }

         // A - B
         std::vector< std::pair<unsigned int, char> > base_vector_difference;
         base_difference(base_vector_union.begin(), base_vector_union.end(), base_vector_intersection.begin(),base_vector_intersection.end(), base_vector_difference);

         // find trimmed region
         // correcting the 5'-end and 3'-end was not done
         // therefore the total number of trimmed bases is 0 yet
         if (base_vector_difference.size() > 0) {
            std::size_t vector_index_leftmost(0);
            std::size_t vector_index_rightmost(base_vector_difference.size() - 1);

            bool keep_going(true);

            while (keep_going) {
               // # of trimmed bases at the 5' end: base_vector_difference[vector_index_leftmost].first + 1
               // # of trimmed bases at the 3' end: read_length - base_vector_difference[vector_index_rightmost].first
               // the 5'-end is smaller
               if ((base_vector_difference[vector_index_leftmost].first + 1) < (read_length - base_vector_difference[vector_index_rightmost].first)) {
                  // check the total number of trimmed bases
                  if ((base_vector_difference[vector_index_leftmost].first + 1 + trim_3_end) <= max_trimmed_bases) {
                     if ((base_vector_difference[vector_index_leftmost].first + 1) > trim_5_end) {
                        trim_5_end = base_vector_difference[vector_index_leftmost].first + 1;
                     }

                     // two points are met
                     if (vector_index_leftmost == vector_index_rightmost) {
                        keep_going = false;
                     }
                     else {
                        vector_index_leftmost++;
                     }
                  }
                  // no need for more check
                  else {
                     keep_going = false;
                  }
               }
               // the 3'-end is smaller
               else {
                  // check the total number of trimmed bases
                  if ((read_length - base_vector_difference[vector_index_rightmost].first) <= max_trimmed_bases) {
                     if ((read_length - base_vector_difference[vector_index_rightmost].first) > trim_3_end) {
                        trim_3_end = read_length - base_vector_difference[vector_index_rightmost].first;
                     }

                     // two points are met
                     if (vector_index_leftmost == vector_index_rightmost) {
                        keep_going = false;
                     }
                     else {
                        vector_index_rightmost--;
                     }
                  }
                  // no need for more check
                  else {
                     keep_going = false;
                  }
               }
            }
         }

         // find consensus modifications
         for (std::size_t it_inter = 0; it_inter < base_vector_intersection.size(); it_inter++) {
            // check whether the base is not in the trimmed regions
            if ((base_vector_intersection[it_inter].first <  (read_length - trim_3_end)) &&
                (base_vector_intersection[it_inter].first >= trim_5_end)) {
               // filter out the bases that are equal to the original ones
               if (sequence[base_vector_intersection[it_inter].first] != base_vector_intersection[it_inter].second) {
                  sequence_modification[base_vector_intersection[it_inter].first] = base_vector_intersection[it_inter].second;
                  sequence[base_vector_intersection[it_inter].first] = base_vector_intersection[it_inter].second;

                  num_corrected_errors_local++;
               }
            }
         }
      }
   }
   // only one path
   else if (candidate_path_vector.size() == 1) {
      // check whether too many bases are modified in the first path, which implies indel may exist in the original read
      // if too many modifications exist, the read is just trimmed
      bool too_many_corrections(false);

      // this is done only when the number of solid islands is two
      if (num_solid_islands == 2) {
         // the 1st island is small && the second island is big
         if (((index_start - left_first) > MAX_MODIFICATION) && ((right_second - index_end + 2) <= MAX_MODIFICATION)) {
            // at least over MAX_MODIFICATION times of modifications from the first modified position
            if (((candidate_path_vector[0].modified_bases[candidate_path_vector[0].modified_bases.size() - 1].first -kmer_length - index_start + 2) >= min_check_length) &&
                (candidate_path_vector[0].modified_bases.size() > MAX_MODIFICATION)) {
               // calculate sum_qs
               //for (unsigned int it_mod_base = 0; it_mod_base < candidate_path_vector[0].modified_bases.size(); it_mod_base++) {
               //   candidate_path_vector[0].sum_qs += ((unsigned short int)quality_score[candidate_path_vector[0].modified_bases[it_mod_base].first] - quality_score_offset);
               //}

               // each modified base in a right range
               //std::size_t partial_qs(candidate_path_vector[0].sum_qs);

               for (std::size_t it_mod_base = (candidate_path_vector[0].modified_bases.size() - 1); it_mod_base >= MAX_MODIFICATION; it_mod_base--) {
                  // at least MAX_MODIFICATION times of modifications within min_check_length bases
                  if ((candidate_path_vector[0].modified_bases[it_mod_base].first - candidate_path_vector[0].modified_bases[it_mod_base - MAX_MODIFICATION].first) < min_check_length) {
                     //if (it_mod_base < (candidate_path_vector[0].modified_bases.size() - 1)) {
                     //   partial_qs -= ((unsigned short int)quality_score[candidate_path_vector[0].modified_bases[it_mod_base + 1].first] - quality_score_offset);
                     //}

                     // average quality score of the modified bases is too high
                     //if (1.0 * partial_qs / (it_mod_base + 1) > quality_score_cutoff) {
                        // trim_5_end or trim_3_end
                        if ((read_length - candidate_path_vector[0].modified_bases[it_mod_base].first) <= max_trimmed_bases) {
                           trim_3_end = read_length - candidate_path_vector[0].modified_bases[it_mod_base].first;

                           // update sequence_modification for the non-trimmed corrections
                           if (it_mod_base > 0) {
                              for (unsigned int it_base = 0; it_base < (it_mod_base - 1); it_base++) {
                                 sequence_modification[(candidate_path_vector[0]).modified_bases[it_base].first] = (candidate_path_vector[0]).modified_bases[it_base].second;
                                 sequence[(candidate_path_vector[0]).modified_bases[it_base].first] = (candidate_path_vector[0]).modified_bases[it_base].second;
                                 num_corrected_errors_local++;
                              }
                           }
                        }
                        else if ((candidate_path_vector[0].modified_bases[it_mod_base].first + 1) <= max_trimmed_bases) {
                           trim_5_end = candidate_path_vector[0].modified_bases[it_mod_base].first + 1;

                           // update sequence_modification for the non-trimmed corrections
                           if (it_mod_base < (candidate_path_vector[0].modified_bases.size() - 1)) {
                              for (unsigned int it_base = (it_mod_base + 1); it_base < candidate_path_vector[0].modified_bases.size(); it_base++) {
                                 sequence_modification[(candidate_path_vector[0]).modified_bases[it_base].first] = (candidate_path_vector[0]).modified_bases[it_base].second;
                                 sequence[(candidate_path_vector[0]).modified_bases[it_base].first] = (candidate_path_vector[0]).modified_bases[it_base].second;
                                 num_corrected_errors_local++;
                              }
                           }
                        }

                        too_many_corrections = true;
                        break;
                     //}
                  }
               }
            }
         }
         // the 1st island is big && the second island is small
         else if (((index_start - left_first) <= MAX_MODIFICATION) && ((right_second - index_end + 2) > MAX_MODIFICATION)) {
            // at least over MAX_MODIFICATION times of modifications from the first modified position
            if (((index_end - candidate_path_vector[0].modified_bases[0].first + 1) >= min_check_length) &&
               (candidate_path_vector[0].modified_bases.size() > MAX_MODIFICATION)) {
               // each modified base in a right range
               //std::size_t partial_qs(candidate_path_vector[0].sum_qs);

               for (unsigned int it_mod_base = 0; it_mod_base < (candidate_path_vector[0].modified_bases.size() - MAX_MODIFICATION); it_mod_base++) {
                  // at least MAX_MODIFICATION times of modifications within min_check_length bases
                  if ((candidate_path_vector[0].modified_bases[it_mod_base + MAX_MODIFICATION].first - candidate_path_vector[0].modified_bases[it_mod_base].first) < min_check_length) {
                     //if (it_mod_base > 0) {
                     //   partial_qs -= ((unsigned short int)quality_score[candidate_path_vector[0].modified_bases[it_mod_base - 1].first] - quality_score_offset);
                     //}

                     // average quality score of the modified bases is too high
                     //if (1.0 * partial_qs / (candidate_path_vector[0].modified_bases.size() - it_mod_base) > quality_score_cutoff) {
                        // trim_5_end or trim_3_end
                        if ((read_length - candidate_path_vector[0].modified_bases[it_mod_base].first) <= max_trimmed_bases) {
                           trim_3_end = read_length - candidate_path_vector[0].modified_bases[it_mod_base].first;
                        }
                        else if ((candidate_path_vector[0].modified_bases[it_mod_base].first + 1) <= max_trimmed_bases) {
                           trim_5_end = candidate_path_vector[0].modified_bases[it_mod_base].first + 1;
                        }

                        too_many_corrections = true;
                        break;
                     //}
                  }
               }
            }
         }
      }

      // not too many modified bases
      if (too_many_corrections == false) {
         // each modification
         for (unsigned int it_base = 0; it_base < candidate_path_vector[0].modified_bases.size(); it_base++) {
            // update sequence_modification
            sequence_modification[candidate_path_vector[0].modified_bases[it_base].first] = candidate_path_vector[0].modified_bases[it_base].second;
            sequence[candidate_path_vector[0].modified_bases[it_base].first] = candidate_path_vector[0].modified_bases[it_base].second;

            num_corrected_errors_local++;
         }
      }
   }
   // no path
   // if indels exist between two solid k-mer islands, checking the solidness of k-mers between index_last_mod and index_end will always fails
   // this kind of errors cannot be corrected without modifying existing solid k-mer islands
   // one of the both the sides should be trimmed
   else if (candidate_path_vector.size() == 0) {
      // trim all the bases to the right of the 1st solid k-mer island
      if ((read_length - (index_start + kmer_length - 1)) <= max_trimmed_bases) {
         trim_3_end = read_length - (index_start + kmer_length - 1);
      }
      // trim all the bases to the left of the 2nd solid k-mer island
      else if ((index_end + 1) <= max_trimmed_bases) {
         trim_5_end = index_end + 1;
      }
      // trim all the bases to the right of the 1st solid k-mer island
      else if ((read_length - (org_boundary_left + kmer_length - 1)) <= max_trimmed_bases) {
         trim_3_end = read_length - (org_boundary_left + kmer_length - 1);
      }
      // trim all the bases to the left of the 2nd solid k-mer island
      else if ((org_boundary_right + 1) <= max_trimmed_bases) {
         trim_5_end = org_boundary_right + 1;
      }
   }
}



//----------------------------------------------------------------------
// correct_errors_5_prime_end
//----------------------------------------------------------------------
inline void C_correct_errors::correct_errors_5_prime_end(const std::string& org_sequence, std::string& sequence, const std::string& quality_score, const std::size_t& index_start, std::string& sequence_modification, std::size_t& trim_5_end, std::size_t& trim_3_end, const std::size_t& org_boundary, const std::size_t& read_length, const std::size_t& max_trimmed_bases, const std::size_t& min_check_length, std::size_t& num_corrected_errors_local, const unsigned char* bit_vector, const std::vector<unsigned int>& hash_seed) {
   // |  non-solid  | 1st solid region
   // |--------------------------------------| read
   //         |-----|                          (index_start)-th k-mer
   //--------------------------------------------------
   // list of candidate paths
   std::vector<C_candidate_path> candidate_path_vector_tmp;

   // make an initial k-mer
   std::string kmer_initial(sequence.substr(index_start, kmer_length));

   // each alternative neocletide
   run_exploration = true;

   for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
      // make a change
      kmer_initial[0] = NEOCLEOTIDE[it_alter];

      // kmer_initial is solid
      if (query_text(kmer_initial, bit_vector, hash_seed) == true) {
         // generate a new path
         C_candidate_path candidate_path;

         if (sequence[index_start] != NEOCLEOTIDE[it_alter]) {
            std::pair<unsigned int, char> pair_tmp;
            pair_tmp.first  = index_start;
            pair_tmp.second = NEOCLEOTIDE[it_alter];

            candidate_path.modified_bases.push_back(pair_tmp);
         }

         // if this k-mer is the first k-mer in a read
         // running extend_a_kmer_5_prime_end is not needed any more
         if (index_start == 0) {
            candidate_path_vector_tmp.push_back(candidate_path);
         }
         else if (index_start > 0) {
            // trace  this kmer recursively and update candidate_path_vector_tmp
            extend_a_kmer_5_prime_end(
                                      kmer_initial,
                                      sequence,
                                      index_start,
                                      candidate_path,
                                      candidate_path_vector_tmp,
                                      org_boundary,
                                      quality_score,
                                      bit_vector,
                                      hash_seed
                                     );
         }
      }
   }

   // complete exploration was not done because there are too many candidata paths
   // remove all the paths in candidate_path_vector_tmp
   if (run_exploration == false) {
      candidate_path_vector_tmp.clear();
   }

   std::vector<C_candidate_path> candidate_path_vector_tmp_tmp;

   // check the solidness of the leftmost k-mers of each modified base
   bool all_solid_wo_modification(false);

   // each candidate path
   for (std::vector<C_candidate_path>::iterator it_path = candidate_path_vector_tmp.begin(); it_path != candidate_path_vector_tmp.end(); it_path++) {
      if ((*it_path).modified_bases.size() == 0) {
         all_solid_wo_modification = true;
         break;
      }
      else {
         std::string sequence_tmp(sequence);

         // index_smallest_modified
         std::size_t index_smallest_modified((*it_path).modified_bases[(*it_path).modified_bases.size() - 1].first);

         // number of bases that should be extended
         std::size_t extend_amount;

         // calculate extend_amount
         // no extension is needed
         // kmer_length = 11, max_extension = 5
         // |0|0|0|0|0|0|0|0|0|0|1|1|1|-
         // |0|1|2|3|4|5|6|7|8|9|0|1|2|-
         // |<------------------->|      k = 11
         // |--------------------------- read
         //                     |<------ index_smallest_modified >= 10
         if (index_smallest_modified >= kmer_length - 1) {
            candidate_path_vector_tmp_tmp.push_back(*it_path);
         }
         // extension is needed
         else {
            // applied the modified bases to sequence_tmp
            for (unsigned int it_base = 0; it_base < (*it_path).modified_bases.size(); it_base++) {
               // modify sequence_tmp
               sequence_tmp[(*it_path).modified_bases[it_base].first] = (*it_path).modified_bases[it_base].second;
            }

            // determine the number of extensions
            // extension amount = kmer_length - index_smallest_modified - 1
            // kmer_length = 11, max_extension = 5
            // |0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|1|1|1|-
            // |5|4|3|2|1|0|1|2|3|4|5|6|7|8|9|0|1|2|-
            //           |<------------------->|      k = 11
            //           |--------------------------- read
            //                     |<------->|        (index_smallest_modified < 10) AND (index_smallest_modified >= 5)
            //     |<------------------->|            index_smallest_modified = 7 -> extend_amount = 3
            if (index_smallest_modified >= kmer_length - max_extension - 1) {
               extend_amount = kmer_length - index_smallest_modified - 1;
            }
            // kmer_length = 11, max_extension = 5
            // |0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|1|1|1|-
            // |5|4|3|2|1|0|1|2|3|4|5|6|7|8|9|0|1|2|-
            //           |<------------------->|      k = 11
            //           |--------------------------- read
            //           |<------->|                  index_smallest_modified < 5
            else {
               extend_amount = max_extension;
            }

            bool extension_success(false);

            // generate an initial k-mer
            std::string kmer_initial(sequence_tmp.substr(0, kmer_length - 1));
            kmer_initial = '0' + kmer_initial;

            // each alternative neocletide
            for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
               // make a change
               kmer_initial[0] = NEOCLEOTIDE[it_alter];

               // kmer_initial is solid
               if (query_text(kmer_initial, bit_vector, hash_seed) == true) {
                  // if extend_amount == 1
                  // running extend_out_left is not needed any more
                  if (extend_amount == 1) {
                     extension_success = true;
                     break;
                  }
                  else {
                     // trace  this kmer recursively and update candidate_path_vector_tmp
                     extend_out_left(
                                     kmer_initial,
                                     1,
                                     extend_amount,
                                     extension_success,
                                     bit_vector,
                                     hash_seed
                                    );
                  }
               }
            }

            if (extension_success == true) {
               candidate_path_vector_tmp_tmp.push_back(*it_path);
            }
         }
      }
   }

   // remain only really modified paths
   candidate_path_vector_tmp.clear();
   // each path
   for (std::size_t it_candidate = 0; it_candidate < candidate_path_vector_tmp_tmp.size(); it_candidate++) {
      // each modification
      bool really_modified(false);
      for (unsigned int it_mod_base = 0; it_mod_base < candidate_path_vector_tmp_tmp[it_candidate].modified_bases.size(); it_mod_base++) {
         if (org_sequence[candidate_path_vector_tmp_tmp[it_candidate].modified_bases[it_mod_base].first] != candidate_path_vector_tmp_tmp[it_candidate].modified_bases[it_mod_base].second) {
            really_modified = true;
         }
      }

      if (really_modified) {
         candidate_path_vector_tmp.push_back(candidate_path_vector_tmp_tmp[it_candidate]);
      }
   }

   candidate_path_vector_tmp_tmp = candidate_path_vector_tmp;

   // all k-mers are solid without any modification
   // do nothing
   if (all_solid_wo_modification == true) {
   }
   // compare quality scores of candidate paths
   // if the number of paths in candidate_path_vector_tmp_tmp is larger than 1
   else if (candidate_path_vector_tmp_tmp.size() > 1) {
      // each path
      std::vector<C_candidate_path>::iterator it_path;
      std::vector<C_candidate_path>::iterator it_path_1st;
      std::vector<C_candidate_path>::iterator it_path_2nd;

      std::size_t qs_1st(INIT_MIN_QS);
      std::size_t qs_2nd(INIT_MIN_QS);

      // each candidate path
      for (it_path = candidate_path_vector_tmp_tmp.begin(); it_path != candidate_path_vector_tmp_tmp.end(); it_path++) {
         // each modification
         for (std::size_t it_mod = 0; it_mod < (*it_path).modified_bases.size(); it_mod++) {
            // add quality scores of modified bases
            if (sequence[(*it_path).modified_bases[it_mod].first] != (*it_path).modified_bases[it_mod].second) {
               (*it_path).sum_qs += ((unsigned short int)quality_score[(*it_path).modified_bases[it_mod].first] - quality_score_offset);
            }
         }

         // compare quality scores of each path
         if ((*it_path).sum_qs <= qs_1st) {
            qs_2nd = qs_1st;
            qs_1st = (*it_path).sum_qs;

            it_path_2nd = it_path_1st;
            it_path_1st = it_path;
         }
         else if ((*it_path).sum_qs <= qs_2nd) {
            qs_2nd = (*it_path).sum_qs;

            it_path_2nd = it_path;
         }
      }

      // check whether too many bases are modified in the first path, which implies indel may exist in the original read
      // if too many modifications exist, the read is just trimmed
      //bool keep_going(true);
      bool too_many_corrections(false);

      // at least over MAX_MODIFICATION times of modifications from the first modified position
      if ((((*it_path_1st).modified_bases[(*it_path_1st).modified_bases.size() - 1].first + 1) >= min_check_length) &&
         ((*it_path_1st).modified_bases.size() > MAX_MODIFICATION)) {
         // each modified base in a right range
         //std::size_t partial_qs((*it_path_1st).sum_qs);

         for (std::size_t it_mod_base = ((*it_path_1st).modified_bases.size() - 1); it_mod_base >= MAX_MODIFICATION; it_mod_base--) {
            // at least MAX_MODIFICATION times of modifications within min_check_length bases
            if (((*it_path_1st).modified_bases[it_mod_base].first - (*it_path_1st).modified_bases[it_mod_base - MAX_MODIFICATION].first) < min_check_length) {
               //if (it_mod_base < ((*it_path_1st).modified_bases.size() - 1)) {
               //   partial_qs -= ((unsigned short int)quality_score[(*it_path_1st).modified_bases[it_mod_base + 1].first] - quality_score_offset);
               //}

               // average quality score of the modified bases is too high
               //if (1.0 * partial_qs / (it_mod_base + 1) > quality_score_cutoff) {
                  // trim_5_end or trim_3_end
                  if ((read_length - (*it_path_1st).modified_bases[it_mod_base].first) <= max_trimmed_bases) {
                     trim_3_end = read_length - (*it_path_1st).modified_bases[it_mod_base].first;
                  }
                  else if (((*it_path_1st).modified_bases[it_mod_base].first + 1) <= max_trimmed_bases) {
                     trim_5_end = (*it_path_1st).modified_bases[it_mod_base].first + 1;
                  }

                  too_many_corrections = true;
                  break;
               //}
            }
         }
      }

      // use the 1st path if not too many corrections are made AND if the 1st path has a sufficiently low score
      // if an indel exists using quality scores is not a good way to choose the best path
      if ((too_many_corrections == false) && (qs_1st + MIN_QS_DIFF <= qs_2nd)) {
         // update sequence_modification
         for (unsigned int it_base = 0; it_base < (*it_path_1st).modified_bases.size(); it_base++) {
            sequence_modification[(*it_path_1st).modified_bases[it_base].first] = (*it_path_1st).modified_bases[it_base].second;

            sequence[(*it_path_1st).modified_bases[it_base].first] = (*it_path_1st).modified_bases[it_base].second;
            num_corrected_errors_local++;
         }
      }
      // hard to choose one path
      else {
         // A AND B
         std::vector< std::pair<unsigned int, char> > base_vector_intersection;
         // A OR B
         std::vector< std::pair<unsigned int, char> > base_vector_union;

         // temporary vectors
         std::vector< std::pair<unsigned int, char> > base_vector_intersection_prev(candidate_path_vector_tmp_tmp[0].modified_bases);
         std::vector< std::pair<unsigned int, char> > base_vector_union_prev(candidate_path_vector_tmp_tmp[0].modified_bases);

         // each candidate path
         for (std::size_t it_p = 1; it_p < candidate_path_vector_tmp_tmp.size(); it_p++) {
            base_vector_intersection.clear();
            base_vector_union.clear();

            base_intersection(base_vector_intersection_prev.begin(), base_vector_intersection_prev.end(), candidate_path_vector_tmp_tmp[it_p].modified_bases.begin(), candidate_path_vector_tmp_tmp[it_p].modified_bases.end(), base_vector_intersection);
            base_union       (base_vector_union_prev.begin(),        base_vector_union_prev.end(),        candidate_path_vector_tmp_tmp[it_p].modified_bases.begin(), candidate_path_vector_tmp_tmp[it_p].modified_bases.end(), base_vector_union);

            base_vector_intersection_prev = base_vector_intersection;
            base_vector_union_prev        = base_vector_union;
         }

         // A - B
         std::vector< std::pair<unsigned int, char> > base_vector_difference;
         base_difference(base_vector_union.begin(), base_vector_union.end(), base_vector_intersection.begin(),base_vector_intersection.end(), base_vector_difference);

         // find trimmed region
         // correcting the 5'-end and 3'-end was not done
         // therefore the total number of trimmed bases is 0 yet
         if (base_vector_difference.size() > 0) {
            std::size_t vector_index_leftmost(0);
            std::size_t vector_index_rightmost(base_vector_difference.size() - 1);

            bool keep_going(true);

            while (keep_going) {
               // # of trimmed bases at the 5' end: base_vector_difference[vector_index_leftmost].first + 1
               // # of trimmed bases at the 3' end: read_length - base_vector_difference[vector_index_rightmost].first
               // the 5'-end is smaller
               if ((base_vector_difference[vector_index_leftmost].first + 1) < (read_length - base_vector_difference[vector_index_rightmost].first)) {
                  // check the total number of trimmed bases
                  if ((base_vector_difference[vector_index_leftmost].first + 1 + trim_3_end) <= max_trimmed_bases) {
                     if ((base_vector_difference[vector_index_leftmost].first + 1) > trim_5_end) {
                        trim_5_end = base_vector_difference[vector_index_leftmost].first + 1;
                     }

                     // two points are met
                     if (vector_index_leftmost == vector_index_rightmost) {
                        keep_going = false;
                     }
                     else {
                        vector_index_leftmost++;
                     }
                  }
                  // no need for more check
                  else {
                     keep_going = false;
                  }
               }
               // the 3'-end is smaller
               else {
                  // check the total number of trimmed bases
                  if ((read_length - base_vector_difference[vector_index_rightmost].first) <= max_trimmed_bases) {
                     if ((read_length - base_vector_difference[vector_index_rightmost].first) > trim_3_end) {
                        trim_3_end = read_length - base_vector_difference[vector_index_rightmost].first;
                     }

                     // two points are met
                     if (vector_index_leftmost == vector_index_rightmost) {
                        keep_going = false;
                     }
                     else {
                        vector_index_rightmost--;
                     }
                  }
                  // no need for more check
                  else {
                     keep_going = false;
                  }
               }
            }
         }

         // find consensus modifications
         for (std::size_t it_inter = 0; it_inter < base_vector_intersection.size(); it_inter++) {
            // check whether the base is not in the trimmed regions
            if ((base_vector_intersection[it_inter].first <  (read_length - trim_3_end)) &&
                (base_vector_intersection[it_inter].first >= trim_5_end)) {
               // filter out the bases that are equal to the original ones
               if (sequence[base_vector_intersection[it_inter].first] != base_vector_intersection[it_inter].second) {
                  sequence_modification[base_vector_intersection[it_inter].first] = base_vector_intersection[it_inter].second;
                  sequence[base_vector_intersection[it_inter].first] = base_vector_intersection[it_inter].second;

                  num_corrected_errors_local++;
               }
            }
         }
      }
   }
   // only one path
   else if (candidate_path_vector_tmp_tmp.size() == 1) {
      // check whether too many bases are modified in the first path, which implies indel may exist in the original read
      // if too many modifications exist, the read is just trimmed
      //bool keep_going(true);
      bool too_many_corrections(false);

      // at least over MAX_MODIFICATION times of modifications from the first modified position
      if (((candidate_path_vector_tmp_tmp[0].modified_bases[0].first + 1) >= min_check_length) &&
         (candidate_path_vector_tmp_tmp[0].modified_bases.size() > MAX_MODIFICATION)) {
         // calculate sum_qs
         //for (unsigned int it_mod_base = 0; it_mod_base < candidate_path_vector_tmp_tmp[0].modified_bases.size(); it_mod_base++) {
         //   candidate_path_vector_tmp_tmp[0].sum_qs += ((unsigned short int)quality_score[candidate_path_vector_tmp_tmp[0].modified_bases[it_mod_base].first] - quality_score_offset);
         //}

         // each modified base in a right range
         //std::size_t partial_qs(candidate_path_vector_tmp_tmp[0].sum_qs);

         for (std::size_t it_mod_base = (candidate_path_vector_tmp_tmp[0].modified_bases.size() - 1); it_mod_base >= MAX_MODIFICATION; it_mod_base--) {
            // at least MAX_MODIFICATION times of modifications within min_check_length bases
            if ((candidate_path_vector_tmp_tmp[0].modified_bases[it_mod_base].first - candidate_path_vector_tmp_tmp[0].modified_bases[it_mod_base - MAX_MODIFICATION].first) < min_check_length) {
               //if (it_mod_base < (candidate_path_vector_tmp_tmp[0].modified_bases.size() - 1)) {
               //   partial_qs -= ((unsigned short int)quality_score[candidate_path_vector_tmp_tmp[0].modified_bases[it_mod_base + 1].first] - quality_score_offset);
               //}

               // average quality score of the modified bases is too high
               //if (1.0 * partial_qs / (it_mod_base + 1) > quality_score_cutoff) {
                  // trim_5_end or trim_3_end
                  if ((read_length - candidate_path_vector_tmp_tmp[0].modified_bases[it_mod_base].first) <= max_trimmed_bases) {
                     trim_3_end = read_length - candidate_path_vector_tmp_tmp[0].modified_bases[it_mod_base].first;

                     // update sequence_modification for the non-trimmed corrections
                     if (it_mod_base > 0) {
                        for (unsigned int it_base = 0; it_base < (it_mod_base - 1); it_base++) {
                           sequence_modification[(candidate_path_vector_tmp_tmp[0]).modified_bases[it_base].first] = (candidate_path_vector_tmp_tmp[0]).modified_bases[it_base].second;
                           sequence[(candidate_path_vector_tmp_tmp[0]).modified_bases[it_base].first] = (candidate_path_vector_tmp_tmp[0]).modified_bases[it_base].second;
                           num_corrected_errors_local++;
                        }
                     }
                  }
                  else if ((candidate_path_vector_tmp_tmp[0].modified_bases[it_mod_base].first + 1) <= max_trimmed_bases) {
                     trim_5_end = candidate_path_vector_tmp_tmp[0].modified_bases[it_mod_base].first + 1;

                     // update sequence_modification for the non-trimmed corrections
                     if (it_mod_base < (candidate_path_vector_tmp_tmp[0].modified_bases.size() - 1)) {
                        for (unsigned int it_base = (it_mod_base + 1); it_base < candidate_path_vector_tmp_tmp[0].modified_bases.size(); it_base++) {
                           sequence_modification[(candidate_path_vector_tmp_tmp[0]).modified_bases[it_base].first] = (candidate_path_vector_tmp_tmp[0]).modified_bases[it_base].second;
                           sequence[(candidate_path_vector_tmp_tmp[0]).modified_bases[it_base].first] = (candidate_path_vector_tmp_tmp[0]).modified_bases[it_base].second;
                           num_corrected_errors_local++;
                        }
                     }
                  }

                  too_many_corrections = true;
                  break;
               //}
            }
         }
      }

      // not too many modified bases
      if (too_many_corrections == false) {
         // update sequence_modification
         for (unsigned int it_base = 0; it_base < candidate_path_vector_tmp_tmp[0].modified_bases.size(); it_base++) {
            sequence_modification[candidate_path_vector_tmp_tmp[0].modified_bases[it_base].first] = candidate_path_vector_tmp_tmp[0].modified_bases[it_base].second;

            sequence[candidate_path_vector_tmp_tmp[0].modified_bases[it_base].first] = candidate_path_vector_tmp_tmp[0].modified_bases[it_base].second;

            num_corrected_errors_local++;
         }
      }
   }
}



//----------------------------------------------------------------------
// correct_errors_3_prime_end
//----------------------------------------------------------------------
inline void C_correct_errors::correct_errors_3_prime_end(const std::string& org_sequence, std::string& sequence, const std::string& quality_score, const std::size_t& index_start, std::string& sequence_modification, std::size_t& trim_5_end, std::size_t& trim_3_end, const std::size_t& org_boundary, const std::size_t& read_length, const std::size_t& max_trimmed_bases, const std::size_t& min_check_length, std::size_t& num_corrected_errors_local, const unsigned char* bit_vector, const std::vector<unsigned int>& hash_seed) {
   //  last solid region | non-solid region |
   // --------------------------------------| read
   //               |-----|                   (index_start)-th k-mer
   //--------------------------------------------------
   // list of candidate paths
   std::vector<C_candidate_path> candidate_path_vector_tmp;

   // make an initial k-mer
   std::string kmer_initial(sequence.substr(index_start, kmer_length));

   // each alternative neocletide
   run_exploration = true;

   for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
      // make a change
      kmer_initial[kmer_length - 1] = NEOCLEOTIDE[it_alter];

      // kmer_initial is solid
      if (query_text(kmer_initial, bit_vector, hash_seed) == true) {
         // generate a new path
         C_candidate_path candidate_path;

         if (sequence[index_start + kmer_length - 1] != NEOCLEOTIDE[it_alter]) {
            std::pair<unsigned int, char> pair_tmp;
            pair_tmp.first  = index_start + kmer_length - 1;
            pair_tmp.second = NEOCLEOTIDE[it_alter];

            candidate_path.modified_bases.push_back(pair_tmp);
         }

         // if this k-mer is the last k-mer in a read
         // running extend_a_kmer_3_prime_end is not needed any more
         if (index_start == (sequence.length() - kmer_length)) {
            candidate_path_vector_tmp.push_back(candidate_path);
         }
         else if (index_start < (sequence.length() - kmer_length)) {
            // trace  this kmer recursively and update candidate_path_vector_tmp
            extend_a_kmer_3_prime_end(
                                      kmer_initial,
                                      sequence,
                                      index_start,
                                      candidate_path,
                                      candidate_path_vector_tmp,
                                      org_boundary,
                                      quality_score,
                                      bit_vector,
                                      hash_seed
                                     );
         }
      }
   }

   // complete exploration was not done because there are too many candidata paths
   // remove all the paths in candidate_path_vector_tmp
   if (run_exploration == false) {
      candidate_path_vector_tmp.clear();
   }

   std::vector<C_candidate_path> candidate_path_vector_tmp_tmp;

   // check the solidness of the rightmost k-mers of each modified base
   bool all_solid_wo_modification(false);

   // each candidate path
   for (std::vector<C_candidate_path>::iterator it_path = candidate_path_vector_tmp.begin(); it_path != candidate_path_vector_tmp.end(); it_path++) {
      if ((*it_path).modified_bases.size() == 0) {
         all_solid_wo_modification = true;
         break;
      }
      else {
         std::string sequence_tmp(sequence);

         // index_largest_modified
         std::size_t index_largest_modified((*it_path).modified_bases[(*it_path).modified_bases.size() - 1].first);

         // number of bases that should be extended
         std::size_t extend_amount;

         // calculate extend_amount
         // no extension is needed
         // sequence.length() = 20, kmer_length = 11, max_extension = 5
         // |0|0|0|1|1|1|1|1|1|1|1|1|1|
         // |7|8|9|0|1|2|3|4|5|6|7|8|9|
         //     |<------------------->| k = 11
         // --------------------------| read
         // ----->|                     index_largest_modified <= 9
         if (index_largest_modified <= sequence.length() - kmer_length) {
            candidate_path_vector_tmp_tmp.push_back(*it_path);
         }
         // extension is needed
         else {
            // applied the modified bases to sequence_tmp
            for (unsigned int it_base = 0; it_base < (*it_path).modified_bases.size(); it_base++) {
               // modify sequence_tmp
               sequence_tmp[(*it_path).modified_bases[it_base].first] = (*it_path).modified_bases[it_base].second;
            }

            // determine the number of extensions
            // sequence.length() = 20, kmer_length = 11, max_extension = 5
            // |0|0|0|1|1|1|1|1|1|1|1|1|1|0|0|0|0|0|
            // |7|8|9|0|1|2|3|4|5|6|7|8|9|1|2|3|4|5|
            //     |<------------------->|           k = 11
            // --------------------------|           read
            //       |<------->|                     (index_largest_modified > 10) AND (index_largest_modified <= 14)
            //           |<------------------->|     index_largest_modified = 12 -> extend_amout = 3
            if (index_largest_modified <= sequence.length() + max_extension - kmer_length) {
               extend_amount = kmer_length - (sequence.length() - index_largest_modified);
            }
            // sequence.length() = 20, kmer_length = 11, max_extension = 5
            // |0|0|0|1|1|1|1|1|1|1|1|1|1|0|0|0|0|0|
            // |7|8|9|0|1|2|3|4|5|6|7|8|9|1|2|3|4|5|
            //     |<------------------->|           k = 11
            // --------------------------|           read
            //                 |<------->|           index_largest_modified > 15
            else {
               extend_amount = max_extension;
            }

            bool extension_success(false);

            // generate an initial k-mer
            // sequence.length() = 20, kmer_length = 11
            // |0|0|0|1|1|1|1|1|1|1|1|1|1|
            // |7|8|9|0|1|2|3|4|5|6|7|8|9|
            //       |<----------------->| kmer_length - 1 = 10
            // --------------------------| read
            //       |-|                   20 - 11 + 1 = 10
            std::string kmer_initial(sequence_tmp.substr(sequence.length() - kmer_length + 1, kmer_length - 1));
            kmer_initial = kmer_initial + '0';

            // each alternative neocletide
            for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
               // make a change
               kmer_initial[kmer_length - 1] = NEOCLEOTIDE[it_alter];

               // kmer_initial is solid
               if (query_text(kmer_initial, bit_vector, hash_seed) == true) {
                  // if extend_amount == 1
                  // running extend_out_right is not needed any more
                  if (extend_amount == 1) {
                     extension_success = true;
                     break;
                  }
                  else {
                     // trace  this kmer recursively and update candidate_path_vector_tmp
                     extend_out_right(
                                      kmer_initial,
                                      1,
                                      extend_amount,
                                      extension_success,
                                      bit_vector,
                                      hash_seed
                                     );
                  }
               }
            }

            if (extension_success == true) {
               candidate_path_vector_tmp_tmp.push_back(*it_path);
            }
         }
      }
   }

   // remain only really modified paths
   candidate_path_vector_tmp.clear();
   // each path
   for (std::size_t it_candidate = 0; it_candidate < candidate_path_vector_tmp_tmp.size(); it_candidate++) {
      // each modification
      bool really_modified(false);
      for (unsigned int it_mod_base = 0; it_mod_base < candidate_path_vector_tmp_tmp[it_candidate].modified_bases.size(); it_mod_base++) {
         if (org_sequence[candidate_path_vector_tmp_tmp[it_candidate].modified_bases[it_mod_base].first] != candidate_path_vector_tmp_tmp[it_candidate].modified_bases[it_mod_base].second) {
            really_modified = true;
         }
      }

      if (really_modified) {
         candidate_path_vector_tmp.push_back(candidate_path_vector_tmp_tmp[it_candidate]);
      }
   }

   candidate_path_vector_tmp_tmp = candidate_path_vector_tmp;

   // all k-mers are solid without any modification
   // do nothing
   if (all_solid_wo_modification == true) {
   }
   // compare quality scores of candidate paths
   // if the number of paths in candidate_path_vector_tmp_tmp is larger than 1
   else if (candidate_path_vector_tmp_tmp.size() > 1) {
      // each path
      std::vector<C_candidate_path>::iterator it_path;
      std::vector<C_candidate_path>::iterator it_path_1st;
      std::vector<C_candidate_path>::iterator it_path_2nd;

      std::size_t qs_1st(INIT_MIN_QS);
      std::size_t qs_2nd(INIT_MIN_QS);

      // each candidate path
      for (it_path = candidate_path_vector_tmp_tmp.begin(); it_path != candidate_path_vector_tmp_tmp.end(); it_path++) {
         // each modification
         for (std::size_t it_mod = 0; it_mod < (*it_path).modified_bases.size(); it_mod++) {
            // add quality scores of modified bases
            (*it_path).sum_qs += ((unsigned short int)quality_score[(*it_path).modified_bases[it_mod].first] - quality_score_offset);
         }

         // compare quality scores of each path
         if ((*it_path).sum_qs <= qs_1st) {
            qs_2nd = qs_1st;
            qs_1st = (*it_path).sum_qs;

            it_path_2nd = it_path_1st;
            it_path_1st = it_path;
         }
         else if ((*it_path).sum_qs <= qs_1st) {
            qs_2nd = (*it_path).sum_qs;

            it_path_2nd = it_path;
         }
      }

      // check whether too many bases are modified in the first path, which implies indel may exist in the original read
      // if too many modifications exist, the read is just trimmed
      //bool keep_going(true);
      bool too_many_corrections(false);

      // at least over MAX_MODIFICATION times of modifications from the first modified position
      if (((read_length - (*it_path_1st).modified_bases[0].first) >= min_check_length) &&
         ((*it_path_1st).modified_bases.size() > MAX_MODIFICATION)) {
         // each modified base in a right range
         //std::size_t partial_qs((*it_path_1st).sum_qs);

         for (unsigned int it_mod_base = 0; it_mod_base < ((*it_path_1st).modified_bases.size() - MAX_MODIFICATION); it_mod_base++) {
            // at least MAX_MODIFICATION times of modifications within min_check_length bases
            if (((*it_path_1st).modified_bases[it_mod_base + MAX_MODIFICATION].first - (*it_path_1st).modified_bases[it_mod_base].first) < min_check_length) {
               //if (it_mod_base > 0) {
               //   partial_qs -= ((unsigned short int)quality_score[(*it_path_1st).modified_bases[it_mod_base - 1].first] - quality_score_offset);
               //}

               // average quality score of the modified bases is too high
               //if (1.0 * partial_qs / ((*it_path_1st).modified_bases.size() - it_mod_base) > quality_score_cutoff) {
                  // trim_5_end or trim_3_end
                  if ((read_length - (*it_path_1st).modified_bases[it_mod_base].first) <= max_trimmed_bases) {
                     trim_3_end = read_length - (*it_path_1st).modified_bases[it_mod_base].first;
                  }
                  else if (((*it_path_1st).modified_bases[it_mod_base].first + 1) <= max_trimmed_bases) {
                     trim_5_end = (*it_path_1st).modified_bases[it_mod_base].first + 1;
                  }

                  too_many_corrections = true;
                  break;
               //}
            }
         }
      }

      // use the 1st path if not too many corrections are made AND if the 1st path has a sufficiently low score
      // if an indel exists using quality scores is not a good way to choose the best path
      if ((too_many_corrections == false) && (qs_1st + MIN_QS_DIFF <= qs_2nd)) {
         // update sequence_modification
         for (unsigned int it_base = 0; it_base < (*it_path_1st).modified_bases.size(); it_base++) {
            sequence_modification[(*it_path_1st).modified_bases[it_base].first] = (*it_path_1st).modified_bases[it_base].second;
            sequence[(*it_path_1st).modified_bases[it_base].first] = (*it_path_1st).modified_bases[it_base].second;

            num_corrected_errors_local++;
         }
      }
      // hard to choose one path
      else {
         // A AND B
         std::vector< std::pair<unsigned int, char> > base_vector_intersection;
         // A OR B
         std::vector< std::pair<unsigned int, char> > base_vector_union;

         // temporary vectors
         std::vector< std::pair<unsigned int, char> > base_vector_intersection_prev(candidate_path_vector_tmp_tmp[0].modified_bases);
         std::vector< std::pair<unsigned int, char> > base_vector_union_prev(candidate_path_vector_tmp_tmp[0].modified_bases);

         // each candidate path
         for (std::size_t it_p = 1; it_p < candidate_path_vector_tmp_tmp.size(); it_p++) {
            base_vector_intersection.clear();
            base_vector_union.clear();

            base_intersection(base_vector_intersection_prev.begin(), base_vector_intersection_prev.end(), candidate_path_vector_tmp_tmp[it_p].modified_bases.begin(), candidate_path_vector_tmp_tmp[it_p].modified_bases.end(), base_vector_intersection);
            base_union       (base_vector_union_prev.begin(),        base_vector_union_prev.end(),        candidate_path_vector_tmp_tmp[it_p].modified_bases.begin(), candidate_path_vector_tmp_tmp[it_p].modified_bases.end(), base_vector_union);

            base_vector_intersection_prev = base_vector_intersection;
            base_vector_union_prev        = base_vector_union;
         }

         // A - B
         std::vector< std::pair<unsigned int, char> > base_vector_difference;
         base_difference(base_vector_union.begin(), base_vector_union.end(), base_vector_intersection.begin(),base_vector_intersection.end(), base_vector_difference);

         // find trimmed region
         // correcting the 5'-end and 3'-end was not done
         // therefore the total number of trimmed bases is 0 yet
         if (base_vector_difference.size() > 0) {
            std::size_t vector_index_leftmost(0);
            std::size_t vector_index_rightmost(base_vector_difference.size() - 1);

            bool keep_going(true);

            while (keep_going) {
               // # of trimmed bases at the 5' end: base_vector_difference[vector_index_leftmost].first + 1
               // # of trimmed bases at the 3' end: read_length - base_vector_difference[vector_index_rightmost].first
               // the 5'-end is smaller
               if ((base_vector_difference[vector_index_leftmost].first + 1) < (read_length - base_vector_difference[vector_index_rightmost].first)) {
                  // check the total number of trimmed bases
                  if ((base_vector_difference[vector_index_leftmost].first + 1 + trim_3_end) <= max_trimmed_bases) {
                     if ((base_vector_difference[vector_index_leftmost].first + 1) > trim_5_end) {
                        trim_5_end = base_vector_difference[vector_index_leftmost].first + 1;
                     }

                     // two points are met
                     if (vector_index_leftmost == vector_index_rightmost) {
                        keep_going = false;
                     }
                     else {
                        vector_index_leftmost++;
                     }
                  }
                  // no need for more check
                  else {
                     keep_going = false;
                  }
               }
               // the 3'-end is smaller
               else {
                  // check the total number of trimmed bases
                  if ((read_length - base_vector_difference[vector_index_rightmost].first) <= max_trimmed_bases) {
                     if ((read_length - base_vector_difference[vector_index_rightmost].first) > trim_3_end) {
                        trim_3_end = read_length - base_vector_difference[vector_index_rightmost].first;
                     }

                     // two points are met
                     if (vector_index_leftmost == vector_index_rightmost) {
                        keep_going = false;
                     }
                     else {
                        vector_index_rightmost--;
                     }
                  }
                  // no need for more check
                  else {
                     keep_going = false;
                  }
               }
            }
         }

         // find consensus modifications
         for (std::size_t it_inter = 0; it_inter < base_vector_intersection.size(); it_inter++) {
            // check whether the base is not in the trimmed regions
            if ((base_vector_intersection[it_inter].first <  (read_length - trim_3_end)) &&
                (base_vector_intersection[it_inter].first >= trim_5_end)) {
               // filter out the bases that are equal to the original ones
               if (sequence[base_vector_intersection[it_inter].first] != base_vector_intersection[it_inter].second) {
                  sequence_modification[base_vector_intersection[it_inter].first] = base_vector_intersection[it_inter].second;
                  sequence[base_vector_intersection[it_inter].first] = base_vector_intersection[it_inter].second;

                  num_corrected_errors_local++;
               }
            }
         }
      }
   }
   // only one path
   else if (candidate_path_vector_tmp_tmp.size() == 1) {
      // check whether too many bases are modified in the first path, which implies indel may exist in the original read
      // if too many modifications exist, the read is just trimmed
      //bool keep_going(true);
      bool too_many_corrections(false);

      // at least over MAX_MODIFICATION times of modifications from the first modified position
      if (((read_length - (candidate_path_vector_tmp_tmp[0]).modified_bases[0].first) >= min_check_length) &&
         ((candidate_path_vector_tmp_tmp[0]).modified_bases.size() > MAX_MODIFICATION)) {
         // calculate sum_qs
         //for (unsigned int it_mod_base = 0; it_mod_base < candidate_path_vector_tmp_tmp[0].modified_bases.size(); it_mod_base++) {
         //   candidate_path_vector_tmp_tmp[0].sum_qs += ((unsigned short int)quality_score[candidate_path_vector_tmp_tmp[0].modified_bases[it_mod_base].first] - quality_score_offset);
         //}

         // each modified base in a right range
         //std::size_t partial_qs((candidate_path_vector_tmp_tmp[0]).sum_qs);

         for (unsigned int it_mod_base = 0; it_mod_base < ((candidate_path_vector_tmp_tmp[0]).modified_bases.size() - MAX_MODIFICATION); it_mod_base++) {
            // at least MAX_MODIFICATION times of modifications within min_check_length bases
            if (((candidate_path_vector_tmp_tmp[0]).modified_bases[it_mod_base + MAX_MODIFICATION].first - (candidate_path_vector_tmp_tmp[0]).modified_bases[it_mod_base].first) < min_check_length) {
               //if (it_mod_base > 0) {
               //   partial_qs -= ((unsigned short int)quality_score[(candidate_path_vector_tmp_tmp[0]).modified_bases[it_mod_base - 1].first] - quality_score_offset);
               //}

               // average quality score of the modified bases is too high
               //if (1.0 * partial_qs / ((candidate_path_vector_tmp_tmp[0]).modified_bases.size() - it_mod_base) > quality_score_cutoff) {
                  // trim_5_end or trim_3_end
                  if ((read_length - (candidate_path_vector_tmp_tmp[0]).modified_bases[it_mod_base].first) <= max_trimmed_bases) {
                     trim_3_end = read_length - (candidate_path_vector_tmp_tmp[0]).modified_bases[it_mod_base].first;

                     // update sequence_modification for the non-trimmed corrections
                     if (it_mod_base > 0) {
                        for (unsigned int it_base = 0; it_base < (it_mod_base - 1); it_base++) {
                           sequence_modification[(candidate_path_vector_tmp_tmp[0]).modified_bases[it_base].first] = (candidate_path_vector_tmp_tmp[0]).modified_bases[it_base].second;
                           sequence[(candidate_path_vector_tmp_tmp[0]).modified_bases[it_base].first] = (candidate_path_vector_tmp_tmp[0]).modified_bases[it_base].second;

                           num_corrected_errors_local++;
                        }
                     }
                  }
                  else if (((candidate_path_vector_tmp_tmp[0]).modified_bases[it_mod_base].first + 1) <= max_trimmed_bases) {
                     trim_5_end = (candidate_path_vector_tmp_tmp[0]).modified_bases[it_mod_base].first + 1;

                     // update sequence_modification for the non-trimmed corrections
                     if (it_mod_base < (candidate_path_vector_tmp_tmp[0].modified_bases.size() - 1)) {
                        for (unsigned int it_base = (it_mod_base + 1); it_base < candidate_path_vector_tmp_tmp[0].modified_bases.size(); it_base++) {
                           sequence_modification[(candidate_path_vector_tmp_tmp[0]).modified_bases[it_base].first] = (candidate_path_vector_tmp_tmp[0]).modified_bases[it_base].second;
                           sequence[(candidate_path_vector_tmp_tmp[0]).modified_bases[it_base].first] = (candidate_path_vector_tmp_tmp[0]).modified_bases[it_base].second;

                           num_corrected_errors_local++;
                        }
                     }
                  }

                  too_many_corrections = true;
                  break;
               //}
            }
         }
      }

      // not too many modified bases
      if (too_many_corrections == false) {
         // update sequence_modification
         for (unsigned int it_base = 0; it_base < candidate_path_vector_tmp_tmp[0].modified_bases.size(); it_base++) {
            sequence_modification[(candidate_path_vector_tmp_tmp[0]).modified_bases[it_base].first] = (candidate_path_vector_tmp_tmp[0]).modified_bases[it_base].second;
            sequence[(candidate_path_vector_tmp_tmp[0]).modified_bases[it_base].first] = (candidate_path_vector_tmp_tmp[0]).modified_bases[it_base].second;

            num_corrected_errors_local++;
         }
      }
   }
}



//----------------------------------------------------------------------
// correct_errors_first_kmer
//----------------------------------------------------------------------
inline void C_correct_errors::correct_errors_first_kmer(const std::string& sequence, const std::string& quality_score, std::string& sequence_modification, std::vector<C_candidate_path>& candidate_path_vector, const unsigned char* bit_vector, const std::vector<unsigned int>& hash_seed) {
   std::string first_kmer(sequence.substr(0, kmer_length));

   std::vector<unsigned int> low_qs_indexes;

   for (unsigned int it_bases = 0; it_bases < kmer_length; it_bases++) {
      if (((unsigned short int)quality_score[it_bases] - quality_score_offset) < quality_score_cutoff) {
         low_qs_indexes.push_back(it_bases);
      }
   }

   // low quality bases exist
   if (low_qs_indexes.size() == 0) {
      // the first k-mer is solid
      if (query_text(first_kmer, bit_vector, hash_seed) == true) {
         C_candidate_path candidate_path;

         candidate_path_vector.push_back(candidate_path);
      }
      else {
         // change each base in the first k-mer and check whether it is solid
         for (unsigned int it_bases = 0; it_bases < kmer_length; it_bases++) {
            std::string kmer_tmp(first_kmer);

            // each alternative neocletide
            for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
               // not equal to the original character
               if (first_kmer[it_bases] != NEOCLEOTIDE[it_alter]) {
                  // generate a new k-mer
                  kmer_tmp[it_bases] = NEOCLEOTIDE[it_alter];

                  // add kmer_tmp to candidate_path_tmp if it is solid
                  if (query_text(kmer_tmp, bit_vector, hash_seed) == true) {
                     // generate a new candidate path
                     C_candidate_path candidate_path;

                     std::pair<unsigned int, char> pair_tmp;
                     pair_tmp.first = it_bases;
                     pair_tmp.second = NEOCLEOTIDE[it_alter];
                     candidate_path.modified_bases.push_back(pair_tmp);

                     candidate_path_vector.push_back(candidate_path);
                  }
               }
            }
         }
      }
   }
   // low quality bases exist
   else if (low_qs_indexes.size() > 0) {
      // the number of low-quality bases is smaller than the threshold
      if (low_qs_indexes.size() <= MAX_LOW_QS_BASES) {
         C_candidate_path candidate_path;

         check_first_kmer(
                          first_kmer,
                          candidate_path,
                          low_qs_indexes,
                          candidate_path_vector,
                          0,
                          bit_vector,
                          hash_seed
                         );

         // no candidate path is found
         if (candidate_path_vector.size() == 0) {
            // the first k-mer is solid
            if (query_text(first_kmer, bit_vector, hash_seed) == true) {
               candidate_path.clear_path();

               candidate_path_vector.push_back(candidate_path);
            }
            else {
               // change each base in the first k-mer and check whether it is solid
               for (unsigned int it_bases = 0; it_bases < kmer_length; it_bases++) {
                  std::string kmer_tmp(first_kmer);

                  // each alternative neocletide
                  for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
                     // not equal to the original character
                     if (first_kmer[it_bases] != NEOCLEOTIDE[it_alter]) {
                        // generate a new k-mer
                        kmer_tmp[it_bases] = NEOCLEOTIDE[it_alter];

                        // add kmer_tmp to candidate_path_tmp if it is solid
                        if (query_text(kmer_tmp, bit_vector, hash_seed) == true) {
                           // generate a new candidate path
                           candidate_path.clear_path();

                           std::pair<unsigned int, char> pair_tmp;
                           pair_tmp.first = it_bases;
                           pair_tmp.second = NEOCLEOTIDE[it_alter];
                           candidate_path.modified_bases.push_back(pair_tmp);

                           candidate_path_vector.push_back(candidate_path);
                        }
                     }
                  }
               }
            }
         }
      }
      // too many low-quality bases
      else {
         // change each base in the first k-mer and check whether it is solid
         C_candidate_path candidate_path;
         for (unsigned int it_bases = 0; it_bases < kmer_length; it_bases++) {
            std::string kmer_tmp(first_kmer);

            // each alternative neocletide
            for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
               // not equal to the original character
               if (first_kmer[it_bases] != NEOCLEOTIDE[it_alter]) {
                  // generate a new k-mer
                  kmer_tmp[it_bases] = NEOCLEOTIDE[it_alter];

                  // add kmer_tmp to candidate_path_tmp if it is solid
                  if (query_text(kmer_tmp, bit_vector, hash_seed) == true) {
                     // generate a new candidate path
                     candidate_path.clear_path();

                     std::pair<unsigned int, char> pair_tmp;
                     pair_tmp.first = it_bases;
                     pair_tmp.second = NEOCLEOTIDE[it_alter];
                     candidate_path.modified_bases.push_back(pair_tmp);

                     candidate_path_vector.push_back(candidate_path);
                  }
               }
            }
         }

         // the first k-mer is solid
         if (query_text(first_kmer, bit_vector, hash_seed) == true) {
            candidate_path.clear_path();

            candidate_path_vector.push_back(candidate_path);
         }
      }
   }
}



//----------------------------------------------------------------------
// check_first_kmer
//----------------------------------------------------------------------
inline void C_correct_errors::check_first_kmer(const std::string& kmer, const C_candidate_path& candidate_path_in, const std::vector<unsigned int>& low_qs_indexes, std::vector<C_candidate_path>& candidate_path_vector, const std::size_t& index, const unsigned char* bit_vector, const std::vector<unsigned int>& hash_seed) {
   std::string new_kmer(kmer);

   for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
      // make a new k-mer
      new_kmer[low_qs_indexes[index]] = NEOCLEOTIDE[it_alter];

      C_candidate_path candidate_path_next(candidate_path_in);

      if (kmer[low_qs_indexes[index]] != NEOCLEOTIDE[it_alter]) {
         std::pair<unsigned int, char> pair_tmp;
         pair_tmp.first = low_qs_indexes[index];
         pair_tmp.second = NEOCLEOTIDE[it_alter];
         candidate_path_next.modified_bases.push_back(pair_tmp);
      }

      // the rightmost low quality base is reached
      // add the path to the vector
      // the original k-mer is also included
      if (index == low_qs_indexes.size() - 1) {
         if (query_text(new_kmer, bit_vector, hash_seed) == true) {
            candidate_path_vector.push_back(candidate_path_next);
         }
      }
      // the rightmost low quality base is not reached
      // do recursively
      else {
         check_first_kmer(
                          new_kmer,
                          candidate_path_next,
                          low_qs_indexes,
                          candidate_path_vector,
                          index + 1,
                          bit_vector,
                          hash_seed
                         );
      }
   }
}



//----------------------------------------------------------------------
// solid_first_kmer
//----------------------------------------------------------------------
inline void C_correct_errors::solid_first_kmer(const C_candidate_path& candidate_path, const std::string& sequence, bool& extension_success, const unsigned char* bit_vector, const std::vector<unsigned int>& hash_seed) {
   // index_smallest_modified
   std::size_t index_smallest_modified(candidate_path.modified_bases[0].first);

   // number of bases that should be extended
   std::size_t extend_amount;

   // applied the modified bases to first_kmer
   std::string first_kmer(sequence.substr(0, kmer_length));
   for (unsigned int it_base = 0; it_base < candidate_path.modified_bases.size(); it_base++) {
      first_kmer[candidate_path.modified_bases[it_base].first] = candidate_path.modified_bases[it_base].second;
   }

   // determine the number of extensions
   // extension amount = kmer_length - index_smallest_modified - 1
   // kmer_length = 11, max_extension = 5
   // |0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|1|1|1|-
   // |5|4|3|2|1|0|1|2|3|4|5|6|7|8|9|0|1|2|-
   //           |<------------------->|      k = 11
   //           |--------------------------- read
   //                     |<------->|        (index_smallest_modified < 10) AND (index_smallest_modified >= 5)
   //     |<------------------->|            index_smallest_modified = 7 -> extend_amount = 3
   if (index_smallest_modified >= kmer_length - max_extension - 1) {
      extend_amount = kmer_length - index_smallest_modified - 1;
   }
   // kmer_length = 11, max_extension = 5
   // |0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|1|1|1|-
   // |5|4|3|2|1|0|1|2|3|4|5|6|7|8|9|0|1|2|-
   //           |<------------------->|      k = 11
   //           |--------------------------- read
   //           |<------->|                  index_smallest_modified < 5
   else {
      extend_amount = max_extension;
   }

   // generate an initial k-mer
   std::string kmer_initial(first_kmer.substr(0, kmer_length - 1));
   kmer_initial = '0' + kmer_initial;

   // each alternative neocletide
   for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
      // make a change
      kmer_initial[0] = NEOCLEOTIDE[it_alter];

      // kmer_initial is solid
      if (query_text(kmer_initial, bit_vector, hash_seed) == true) {
         // if extend_amount == 1
         // running extend_out_left is not needed any more
         if (extend_amount == 1) {
            extension_success = true;
            break;
         }
         else {
            // trace  this kmer recursively and update candidate_path_vector_tmp
            extend_out_left(
                            kmer_initial,
                            1,
                            extend_amount,
                            extension_success,
                            bit_vector,
                            hash_seed
                           );
         }
      }
   }
}



//----------------------------------------------------------------------
// extend_first_kmer_to_right
//----------------------------------------------------------------------
inline void C_correct_errors::extend_first_kmer_to_right(const std::string& sequence, const std::string& quality_score, C_candidate_path& candidate_path_in, std::vector<C_candidate_path>& candidate_path_vector_all, const std::size_t& read_length, const unsigned char* bit_vector, const std::vector<unsigned int>& hash_seed) {
   // generate the first k-mer
   std::string first_kmer(sequence.substr(0, kmer_length + 1));
   for (unsigned int it_base = 0; it_base < candidate_path_in.modified_bases.size(); it_base++) {
      first_kmer[candidate_path_in.modified_bases[it_base].first] = candidate_path_in.modified_bases[it_base].second;
   }

   // generate the second k-mer
   std::string second_kmer(first_kmer.substr(1, kmer_length));

   // list of candidate paths
   std::vector<C_candidate_path> candidate_path_vector_tmp;

   // second_kmer is solid
   if (query_text(second_kmer, bit_vector, hash_seed) == true) {
      // if this k-mer is the last k-mer in a read
      // running extend_a_kmer_3_prime_end is not needed any more
      if ((sequence.length() - kmer_length) == 1) {
         candidate_path_vector_tmp.push_back(candidate_path_in);
      }
      else if ((sequence.length() - kmer_length) > 1) {
         // trace  this kmer recursively and update candidate_path_vector_tmp
         // org_boundary is not needed in this case: read_length is being used as a dummy number
         extend_a_kmer_3_prime_end(
                                   second_kmer,
                                   sequence,
                                   1,
                                   candidate_path_in,
                                   candidate_path_vector_tmp,
                                   read_length,
                                   quality_score,
                                   bit_vector,
                                   hash_seed
                                  );
      }
   }
   // second_kmer is not solid
   else {
      // each alternative neocletide
      for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
         // not equal to the original character
         if (sequence[kmer_length] != NEOCLEOTIDE[it_alter]) {
            // make a change
            second_kmer[kmer_length - 1] = NEOCLEOTIDE[it_alter];

            // new second_kmer is solid
            if (query_text(second_kmer, bit_vector, hash_seed) == true) {
               // generate a new path
               C_candidate_path new_candidate_path(candidate_path_in);

               std::pair<unsigned int, char> pair_tmp;
               pair_tmp.first  = kmer_length;
               pair_tmp.second = NEOCLEOTIDE[it_alter];

               new_candidate_path.modified_bases.push_back(pair_tmp);

               // if this k-mer is the last k-mer in a read
               // running extend_a_kmer_3_prime_end is not needed any more
               if ((sequence.length() - kmer_length) == 1) {
                  candidate_path_vector_tmp.push_back(new_candidate_path);
               }
               else if ((sequence.length() - kmer_length) > 1) {
                  // trace  this kmer recursively and update candidate_path_vector_tmp
                  // org_boundary is not needed in this case: read_length is being used as a dummy number
                  extend_a_kmer_3_prime_end(
                                            second_kmer,
                                            sequence,
                                            1,
                                            new_candidate_path,
                                            candidate_path_vector_tmp,
                                            read_length,
                                            quality_score,
                                            bit_vector,
                                            hash_seed
                                           );
               }
            }
         }
      }
   }

   // complete exploration was not done because there are too many candidata paths
   // remove all the paths in candidate_path_vector_tmp
   if (run_exploration == false) {
      candidate_path_vector_tmp.clear();
   }

   // each candidate path
   for (std::vector<C_candidate_path>::iterator it_path = candidate_path_vector_tmp.begin(); it_path != candidate_path_vector_tmp.end(); it_path++) {
      std::string sequence_tmp(sequence);

      // index_largest_modified
      std::size_t index_largest_modified((*it_path).modified_bases[(*it_path).modified_bases.size() - 1].first);

      // number of bases that should be extended
      std::size_t extend_amount;

      // calculate extend_amount
      // no extension is needed
      // sequence.length() = 20, kmer_length = 11, max_extension = 5
      // |0|0|0|1|1|1|1|1|1|1|1|1|1|
      // |7|8|9|0|1|2|3|4|5|6|7|8|9|
      //     |<------------------->| k = 11
      // --------------------------| read
      // ----->|                     index_largest_modified <= 9
      if (index_largest_modified <= sequence.length() - kmer_length) {
         candidate_path_vector_all.push_back(*it_path);
      }
      // extension is needed
      else {
         // applied the modified bases to sequence_tmp
         for (unsigned int it_base = 0; it_base < (*it_path).modified_bases.size(); it_base++) {
            // modify sequence_tmp
            sequence_tmp[(*it_path).modified_bases[it_base].first] = (*it_path).modified_bases[it_base].second;
         }

         // determine the number of extensions
         // sequence.length() = 20, kmer_length = 11, max_extension = 5
         // |0|0|0|1|1|1|1|1|1|1|1|1|1|0|0|0|0|0|
         // |7|8|9|0|1|2|3|4|5|6|7|8|9|1|2|3|4|5|
         //     |<------------------->|           k = 11
         // --------------------------|           read
         //       |<------->|                     (index_largest_modified > 10) AND (index_largest_modified <= 14)
         //           |<------------------->|     index_largest_modified = 12 -> extend_amout = 3
         if (index_largest_modified <= sequence.length() + max_extension - kmer_length) {
            extend_amount = kmer_length - (sequence.length() - index_largest_modified);
         }
         // sequence.length() = 20, kmer_length = 11, max_extension = 5
         // |0|0|0|1|1|1|1|1|1|1|1|1|1|0|0|0|0|0|
         // |7|8|9|0|1|2|3|4|5|6|7|8|9|1|2|3|4|5|
         //     |<------------------->|           k = 11
         // --------------------------|           read
         //                 |<------->|           index_largest_modified > 15
         else {
            extend_amount = max_extension;
         }

         bool extension_success(false);

         // generate an initial k-mer
         // sequence.length() = 20, kmer_length = 11
         // |0|0|0|1|1|1|1|1|1|1|1|1|1|
         // |7|8|9|0|1|2|3|4|5|6|7|8|9|
         //       |<----------------->| kmer_length - 1 = 10
         // --------------------------| read
         //       |-|                   20 - 11 + 1 = 10
         std::string kmer_initial(sequence_tmp.substr(sequence.length() - kmer_length + 1, kmer_length - 1));
         kmer_initial = kmer_initial + '0';

         // each alternative neocletide
         for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
            // make a change
            kmer_initial[kmer_length - 1] = NEOCLEOTIDE[it_alter];

            // kmer_initial is solid
            if (query_text(kmer_initial, bit_vector, hash_seed) == true) {
               // if extend_amount == 1
               // running extend_out_right is not needed any more
               if (extend_amount == 1) {
                  extension_success = true;
                  break;
               }
               else {
                  // trace  this kmer recursively and update candidate_path_vector_tmp
                  extend_out_right(
                                   kmer_initial,
                                   1,
                                   extend_amount,
                                   extension_success,
                                   bit_vector,
                                   hash_seed
                                  );
               }
            }
         }

         if (extension_success == true) {
            candidate_path_vector_all.push_back(*it_path);
         }
      }
   }
}



//----------------------------------------------------------------------
// extend_a_kmer
//----------------------------------------------------------------------
inline void C_correct_errors::extend_a_kmer(const std::string& kmer, const std::string& sequence, const std::size_t& index_kmer, const std::size_t& index_last_mod, C_candidate_path& current_path, std::vector<C_candidate_path>& candidate_path_vector, const std::size_t& org_boundary_left, const std::size_t& org_boundary_right, const std::string& quality_score, const unsigned char* bit_vector, const std::vector<unsigned int>& hash_seed) {
   if (run_exploration == true) {
      // generate a new k-mer
      std::string kmer_new(kmer.substr(1, kmer_length - 1));
      kmer_new.push_back(sequence[index_kmer + kmer_length]);

      // this was the real boundary between solid k-mers and weak k-mers
      // check all possible cases
      if ((index_kmer == (org_boundary_left - 1)) || (index_kmer == (org_boundary_right - 1)) || (((unsigned short int)quality_score[index_kmer + kmer_length] - quality_score_offset) <= extremely_low_quality_score)) {
         // each alternative neocletide
         for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
            // make a change
            kmer_new[kmer_length - 1] = NEOCLEOTIDE[it_alter];

            // kmer_new is solid
            if (query_text(kmer_new, bit_vector, hash_seed) == true) {
               // generate a new path
               C_candidate_path temporary_path(current_path);

               if (sequence[index_kmer + kmer_length] != NEOCLEOTIDE[it_alter]) {
                  std::pair<unsigned int, char> pair_tmp;
                  pair_tmp.first  = index_kmer + kmer_length;
                  pair_tmp.second = NEOCLEOTIDE[it_alter];

                  temporary_path.modified_bases.push_back(pair_tmp);
               }

               // if this k-mer is the last k-mer that can be modified
               // running extend_a_kmer_right is not needed any more
               if ((index_kmer + 1) == index_last_mod) {
                  candidate_path_vector.push_back(temporary_path);

                  // too many candidate paths
                  if (candidate_path_vector.size() > MAX_CANDIDATE_PATHS) {
                     run_exploration = false;
                  }
               }
               else {
                  // trace  this kmer recursively and update candidate_path_vector
                  extend_a_kmer(
                                kmer_new,
                                sequence,
                                index_kmer + 1,
                                index_last_mod,
                                temporary_path,
                                candidate_path_vector,
                                org_boundary_left,
                                org_boundary_right,
                                quality_score,
                                bit_vector,
                                hash_seed
                               );
               }
            }
         }
      }
      else {
         // kmer_new is a solid k-mer
         if (query_text(kmer_new, bit_vector, hash_seed) == true) {
            // if this k-mer is the last k-mer that can be modified
            // running extend_a_kmer_right is not needed any more
            if ((index_kmer + 1) == index_last_mod) {
               candidate_path_vector.push_back(current_path);

               // too many candidate paths
               if (candidate_path_vector.size() > MAX_CANDIDATE_PATHS) {
                  run_exploration = false;
               }
            }
            else {
               extend_a_kmer(
                             kmer_new,
                             sequence,
                             index_kmer + 1,
                             index_last_mod,
                             current_path,
                             candidate_path_vector,
                             org_boundary_left,
                             org_boundary_right,
                             quality_score,
                             bit_vector,
                             hash_seed
                            );
            }
         }
         else {
            // each alternative neocletide
            for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
               // not equal to the original character
               if (sequence[index_kmer + kmer_length] != NEOCLEOTIDE[it_alter]) {
                  // make a change
                  kmer_new[kmer_length - 1] = NEOCLEOTIDE[it_alter];

                  // kmer_new is solid
                  if (query_text(kmer_new, bit_vector, hash_seed) == true) {
                     // generate a new path
                     C_candidate_path temporary_path(current_path);

                     std::pair<unsigned int, char> pair_tmp;
                     pair_tmp.first  = index_kmer + kmer_length;
                     pair_tmp.second = NEOCLEOTIDE[it_alter];

                     temporary_path.modified_bases.push_back(pair_tmp);

                     // if this k-mer is the last k-mer that can be modified
                     // running extend_a_kmer_right is not needed any more
                     if ((index_kmer + 1) == index_last_mod) {
                        candidate_path_vector.push_back(temporary_path);

                        // too many candidate paths
                        if (candidate_path_vector.size() > MAX_CANDIDATE_PATHS) {
                           run_exploration = false;
                        }
                     }
                     else {
                        // trace  this kmer recursively and update candidate_path_vector
                        extend_a_kmer(
                                      kmer_new,
                                      sequence,
                                      index_kmer + 1,
                                      index_last_mod,
                                      temporary_path,
                                      candidate_path_vector,
                                      org_boundary_left,
                                      org_boundary_right,
                                      quality_score,
                                      bit_vector,
                                      hash_seed
                                     );
                     }
                  }
               }
            }
         }
      }
   }
}



//----------------------------------------------------------------------
// extend_a_kmer_5_prime_end
//----------------------------------------------------------------------
inline void C_correct_errors::extend_a_kmer_5_prime_end(const std::string& kmer, const std::string& sequence, const std::size_t& index_kmer, C_candidate_path& current_path, std::vector<C_candidate_path>& candidate_path_vector, const std::size_t& org_boundary, const std::string& quality_score, const unsigned char* bit_vector, const std::vector<unsigned int>& hash_seed) {
   if (run_exploration == true) {
      // generate a new k-mer
      std::string kmer_new(kmer.substr(0, kmer_length - 1));
      kmer_new = sequence[index_kmer - 1] + kmer_new;

      // this was the real boundary between solid k-mers and weak k-mers
      // check all possible cases
      if ((index_kmer == (org_boundary + 1)) || (((unsigned short int)quality_score[index_kmer - 1] - quality_score_offset) <= extremely_low_quality_score)) {
         // each alternative neocletide
         for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
            // make a change
            kmer_new[0] = NEOCLEOTIDE[it_alter];

            // kmer_new is solid
            if (query_text(kmer_new, bit_vector, hash_seed) == true) {
               // generate a new path
               C_candidate_path temporary_path(current_path);

               // not equal to the original character
               if (sequence[index_kmer - 1] != NEOCLEOTIDE[it_alter]) {
                  std::pair<unsigned int, char> pair_tmp;
                  pair_tmp.first  = index_kmer - 1;
                  pair_tmp.second = NEOCLEOTIDE[it_alter];

                  temporary_path.modified_bases.push_back(pair_tmp);
               }

               // if this k-mer is the first k-mer in a read
               // running extend_a_kmer_5_prime_end is not needed any more
               if ((index_kmer - 1) == 0) {
                  candidate_path_vector.push_back(temporary_path);

                  // too many candidate paths
                  if (candidate_path_vector.size() > MAX_CANDIDATE_PATHS) {
                     run_exploration = false;
                  }
               }
               else if ((index_kmer - 1) > 0) {
                  // trace  this kmer recursively and update candidate_path_vector
                  extend_a_kmer_5_prime_end(
                                            kmer_new,
                                            sequence,
                                            index_kmer - 1,
                                            temporary_path,
                                            candidate_path_vector,
                                            org_boundary,
                                            quality_score,
                                            bit_vector,
                                            hash_seed
                                           );
               }
            }
         }
      }
      else {
         // kmer_new is a solid k-mer
         if (query_text(kmer_new, bit_vector, hash_seed) == true) {
            // if this k-mer is the first k-mer in a read
            // running extend_a_kmer_5_prime_end is not needed any more
            if ((index_kmer - 1) == 0) {
               candidate_path_vector.push_back(current_path);

               // too many candidate paths
               if (candidate_path_vector.size() > MAX_CANDIDATE_PATHS) {
                  run_exploration = false;
               }
            }
            else if ((index_kmer - 1) > 0) {
               extend_a_kmer_5_prime_end(
                                         kmer_new,
                                         sequence,
                                         index_kmer - 1,
                                         current_path,
                                         candidate_path_vector,
                                         org_boundary,
                                         quality_score,
                                         bit_vector,
                                         hash_seed
                                        );
            }
         }
         else {
            // each alternative neocletide
            for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
               // not equal to the original character
               if (sequence[index_kmer - 1] != NEOCLEOTIDE[it_alter]) {
                  // make a change
                  kmer_new[0] = NEOCLEOTIDE[it_alter];

                  // kmer_new is solid
                  if (query_text(kmer_new, bit_vector, hash_seed) == true) {
                     // generate a new path
                     C_candidate_path temporary_path(current_path);

                     std::pair<unsigned int, char> pair_tmp;
                     pair_tmp.first  = index_kmer - 1;
                     pair_tmp.second = NEOCLEOTIDE[it_alter];

                     temporary_path.modified_bases.push_back(pair_tmp);

                     // if this k-mer is the first k-mer in a read
                     // running extend_a_kmer_5_prime_end is not needed any more
                     if ((index_kmer - 1) == 0) {
                        candidate_path_vector.push_back(temporary_path);

                        // too many candidate paths
                        if (candidate_path_vector.size() > MAX_CANDIDATE_PATHS) {
                           run_exploration = false;
                        }
                     }
                     else if ((index_kmer - 1) > 0) {
                        // trace  this kmer recursively and update candidate_path_vector
                        extend_a_kmer_5_prime_end(
                                                  kmer_new,
                                                  sequence,
                                                  index_kmer - 1,
                                                  temporary_path,
                                                  candidate_path_vector,
                                                  org_boundary,
                                                  quality_score,
                                                  bit_vector,
                                                  hash_seed
                                                 );
                     }
                  }
               }
            }
         }
      }
   }
}



//----------------------------------------------------------------------
// extend_a_kmer_3_prime_end
//----------------------------------------------------------------------
inline void C_correct_errors::extend_a_kmer_3_prime_end(const std::string& kmer, const std::string& sequence, const std::size_t& index_kmer, C_candidate_path& current_path, std::vector<C_candidate_path>& candidate_path_vector, const std::size_t& org_boundary, const std::string& quality_score, const unsigned char* bit_vector, const std::vector<unsigned int>& hash_seed) {
   if (run_exploration == true) {
      // generate a new k-mer
      std::string kmer_new(kmer.substr(1, kmer_length - 1));
      kmer_new = kmer_new + sequence[index_kmer + kmer_length];

      // this was the real boundary between solid k-mers and weak k-mers
      // check all possible cases
      if ((index_kmer == (org_boundary - 1)) || (((unsigned short int)quality_score[index_kmer + kmer_length] - quality_score_offset) <= extremely_low_quality_score)) {
         // each alternative neocletide
         for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
            // make a change
            kmer_new[kmer_length - 1] = NEOCLEOTIDE[it_alter];

            // kmer_new is solid
            if (query_text(kmer_new, bit_vector, hash_seed) == true) {
               // generate a new path
               C_candidate_path temporary_path(current_path);

               // not equal to the original character
               if (sequence[index_kmer + kmer_length] != NEOCLEOTIDE[it_alter]) {
                  std::pair<unsigned int, char> pair_tmp;
                  pair_tmp.first  = index_kmer + kmer_length;
                  pair_tmp.second = NEOCLEOTIDE[it_alter];

                  temporary_path.modified_bases.push_back(pair_tmp);
               }

               // if this k-mer is the last k-mer in a read
               // running extend_a_kmer_3_prime_end is not needed any more
               if ((index_kmer + 1) == (sequence.length() - kmer_length)) {
                  candidate_path_vector.push_back(temporary_path);

                  // too many candidate paths
                  if (candidate_path_vector.size() > MAX_CANDIDATE_PATHS) {
                     run_exploration = false;
                  }
               }
               else if ((index_kmer + 1) < (sequence.length() - kmer_length)) {
                  // trace  this kmer recursively and update candidate_path_vector
                  extend_a_kmer_3_prime_end(
                                            kmer_new,
                                            sequence,
                                            index_kmer + 1,
                                            temporary_path,
                                            candidate_path_vector,
                                            org_boundary,
                                            quality_score,
                                            bit_vector,
                                            hash_seed
                                           );
               }
            }
         }
      }
      else {
         // kmer_new is a solid k-mer
         if (query_text(kmer_new, bit_vector, hash_seed) == true) {
            // if this k-mer is the last k-mer in a read
            // running extend_a_kmer_3_prime_end is not needed any more
            if ((index_kmer + 1) == (sequence.length() - kmer_length)) {
               candidate_path_vector.push_back(current_path);

               // too many candidate paths
               if (candidate_path_vector.size() > MAX_CANDIDATE_PATHS) {
                  run_exploration = false;
               }
            }
            else if ((index_kmer + 1) < (sequence.length() - kmer_length)) {
               extend_a_kmer_3_prime_end(
                                         kmer_new,
                                         sequence,
                                         index_kmer + 1,
                                         current_path,
                                         candidate_path_vector,
                                         org_boundary,
                                         quality_score,
                                         bit_vector,
                                         hash_seed
                                        );
            }
         }
         else {
            // each alternative neocletide
            for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
               // not equal to the original character
               if (sequence[index_kmer + kmer_length] != NEOCLEOTIDE[it_alter]) {
                  // make a change
                  kmer_new[kmer_length - 1] = NEOCLEOTIDE[it_alter];

                  // kmer_new is solid
                  if (query_text(kmer_new, bit_vector, hash_seed) == true) {
                     // generate a new path
                     C_candidate_path temporary_path(current_path);

                     std::pair<unsigned int, char> pair_tmp;
                     pair_tmp.first  = index_kmer + kmer_length;
                     pair_tmp.second = NEOCLEOTIDE[it_alter];

                     temporary_path.modified_bases.push_back(pair_tmp);

                     // if this k-mer is the last k-mer in a read
                     // running extend_a_kmer_3_prime_end is not needed any more
                     if ((index_kmer + 1) == (sequence.length() - kmer_length)) {
                        candidate_path_vector.push_back(temporary_path);

                        // too many candidate paths
                        if (candidate_path_vector.size() > MAX_CANDIDATE_PATHS) {
                           run_exploration = false;
                        }
                     }
                     else if ((index_kmer + 1) < (sequence.length() - kmer_length)) {
                        // trace  this kmer recursively and update candidate_path_vector
                        extend_a_kmer_3_prime_end(
                                                  kmer_new,
                                                  sequence,
                                                  index_kmer + 1,
                                                  temporary_path,
                                                  candidate_path_vector,
                                                  org_boundary,
                                                  quality_score,
                                                  bit_vector,
                                                  hash_seed
                                                 );
                     }
                  }
               }
            }
         }
      }
   }
}



//----------------------------------------------------------------------
// extend_out_left
//----------------------------------------------------------------------
inline void C_correct_errors::extend_out_left(const std::string& kmer, const std::size_t& num_extend, const std::size_t& extend_amount, bool& extension_success, const unsigned char* bit_vector, const std::vector<unsigned int>& hash_seed) {
   // generate a new k-mer
   std::string kmer_new(kmer.substr(0, kmer_length - 1));
   kmer_new = '0' + kmer_new;

   // each alternative neocletide
   for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
      // generate kmer_new
      kmer_new[0] = NEOCLEOTIDE[it_alter];

      // kmer_new is solid
      if (query_text(kmer_new, bit_vector, hash_seed) == true) {
         // if current num_extend = extend_amount
         // running extend_out_left is not needed any more
         if ((num_extend + 1) == extend_amount) {
            extension_success = true;
            break;
         }
         else {
            // trace  this kmer recursively
            extend_out_left(
                            kmer_new,
                            num_extend + 1,
                            extend_amount,
                            extension_success,
                            bit_vector,
                            hash_seed
                           );
         }
      }
   }
}



//----------------------------------------------------------------------
// extend_out_right
//----------------------------------------------------------------------
inline void C_correct_errors::extend_out_right(const std::string& kmer, const std::size_t& num_extend, const std::size_t& extend_amount, bool& extension_success, const unsigned char* bit_vector, const std::vector<unsigned int>& hash_seed) {
   // generate a new k-mer
   std::string kmer_new(kmer.substr(1, kmer_length - 1));
   kmer_new = kmer_new + '0';

   // each alternative neocletide
   for (unsigned short int it_alter = A; it_alter <= T; it_alter++) {
      // generate kmer_new
      kmer_new[kmer_length - 1] = NEOCLEOTIDE[it_alter];

      // kmer_new is solid
      if (query_text(kmer_new, bit_vector, hash_seed) == true) {
         // if current num_extend = extend_amount
         // running extend_out_right is not needed any more
         if ((num_extend + 1) == extend_amount) {
            extension_success = true;
            break;
         }
         else {
            // trace  this kmer recursively
            extend_out_right(
                             kmer_new,
                             num_extend + 1,
                             extend_amount,
                             extension_success,
                             bit_vector,
                             hash_seed
                            );
         }
      }
   }
}



//----------------------------------------------------------------------
// query_text
//----------------------------------------------------------------------
inline bool C_correct_errors::query_text(const std::string& kmer, const unsigned char*& bit_vector, const std::vector<unsigned int>& hash_seed) {
   //--------------------------------------------------
   // reverse complement or not
   //--------------------------------------------------
   std::string kmer_internal;
   std::string kmer_rc;

   // lexicographical comparison
   //  1: kmer is faster than its reverse complement in the lexicographical order
   //  0: same
   // -1: reverse complement of kmer is faster than kmer in the lexicographical order
   int comparison_result(0);
   for (unsigned int it = 0, it_rc = (kmer_length - 1); it < kmer_length; it++, it_rc--) {
      if (kmer[it] == 'A') {
         if (kmer[it_rc] == 'T') {
         }
         else {
            comparison_result = 1;
            break;
         }
      }
      else if (kmer[it] == 'C') {
         if (kmer[it_rc] == 'G') {
         }
         else if (kmer[it_rc] == 'T') {
            comparison_result = -1;
            break;
         }
         else {
            comparison_result = 1;
            break;
         }
      }
      else if (kmer[it] == 'G') {
         if (kmer[it_rc] == 'A') {
            comparison_result = 1;
            break;
         }
         else if (kmer[it_rc] == 'C') {
         }
         else {
            comparison_result = -1;
            break;
         }
      }
      else if (kmer[it] == 'T') {
         if (kmer[it_rc] == 'A') {
         }
         else {
            comparison_result = -1;
            break;
         }
      }
      else {
         std::cout << std::endl << "ERROR: Illegal character " << kmer[it] << " (query_text)" << std::endl << std::endl;
         MPI_Abort(MPI_COMM_WORLD, 226);
      }
   }

   // kmer is faster
   if (comparison_result >= 0) {
      kmer_internal = kmer;
   }
   // reverse complement is faster
   else {
      // reverse complement
      kmer_internal.resize(kmer_length);

      for (unsigned int it = 0, rc_index = kmer_length - 1; it < kmer_length; it++, rc_index--) {
         switch (kmer[it]) {
            case 'A' :
               kmer_internal[rc_index] = 'T';
               break;
            case 'C' :
               kmer_internal[rc_index] = 'G';
               break;
            case 'G' :
               kmer_internal[rc_index] = 'C';
               break;
            case 'T' :
               kmer_internal[rc_index] = 'A';
               break;
            default :
               std::cout << std::endl << "ERROR: Illegal character " << kmer[it] << " (query_text)" << std::endl << std::endl;
               MPI_Abort(MPI_COMM_WORLD, 227);
               break;
         }
      }
   }

   bloom_type original_index;

   bloom_type hash [2];

   unsigned short int bit_index;

   for (unsigned short int it_hash_func = 0; it_hash_func < num_hash_func; it_hash_func++) {
      MurmurHash3_x64_128(kmer_internal.c_str(), kmer_length, hash_seed[it_hash_func], hash);

      original_index = hash[0] % bit_vector_width;
      bit_index      = original_index % BITS_PER_CHAR;

      if ((bit_vector[original_index / BITS_PER_CHAR] & BIT_MASK[bit_index]) != BIT_MASK[bit_index]) {
         return false;
      }
   }

   return true;
}



//----------------------------------------------------------------------
// base_intersection
//----------------------------------------------------------------------
inline void C_correct_errors::base_intersection(std::vector< std::pair<unsigned int, char> >::iterator in1_begin, const std::vector< std::pair<unsigned int, char> >::iterator& in1_end, std::vector< std::pair<unsigned int, char> >::iterator in2_begin, const std::vector< std::pair<unsigned int, char> >::iterator& in2_end, std::vector< std::pair<unsigned int, char> >& out) {
   while ((in1_begin != in1_end) && (in2_begin != in2_end)) {
      if ((*in1_begin).first < (*in2_begin).first) {
         in1_begin++;
      }
      else if ((*in2_begin).first < (*in1_begin).first) {
         in2_begin++;
      }
      else {
         if ((*in1_begin).second == (*in2_begin).second) {
            out.push_back(*in1_begin);
         }

         in1_begin++;
         in2_begin++;
      }
   }
}



//----------------------------------------------------------------------
// base_union
//----------------------------------------------------------------------
inline void C_correct_errors::base_union(std::vector< std::pair<unsigned int, char> >::iterator in1_begin, const std::vector< std::pair<unsigned int, char> >::iterator& in1_end, std::vector< std::pair<unsigned int, char> >::iterator in2_begin, const std::vector< std::pair<unsigned int, char> >::iterator& in2_end, std::vector< std::pair<unsigned int, char> >& out) {
   while (true) {
      if (in1_begin == in1_end) {
         out.insert(out.end(), in2_begin, in2_end);
         return;
      }
      else if (in2_begin == in2_end) {
         out.insert(out.end(), in1_begin, in1_end);
         return;
      }
      else {
         if (((*in1_begin).first) < ((*in2_begin).first)) {
            out.push_back(*in1_begin);
            in1_begin++;
         }
         else if (((*in2_begin).first) < ((*in1_begin).first)) {
            out.push_back(*in2_begin);
            in2_begin++;
         }
         else {
            out.push_back(*in1_begin);
            in1_begin++;
            in2_begin++;
         }
      }
   }
}



//----------------------------------------------------------------------
// base_difference
//----------------------------------------------------------------------
inline void C_correct_errors::base_difference(std::vector< std::pair<unsigned int, char> >::iterator in1_begin, const std::vector< std::pair<unsigned int, char> >::iterator& in1_end, std::vector< std::pair<unsigned int, char> >::iterator in2_begin, const std::vector< std::pair<unsigned int, char> >::iterator& in2_end, std::vector< std::pair<unsigned int, char> >& out) {
   while ((in1_begin != in1_end) && (in2_begin != in2_end)) {
      if ((*in1_begin).first < (*in2_begin).first) {
         out.push_back(*in1_begin);
         in1_begin++;
      }
      else if ((*in2_begin).first < (*in1_begin).first) {
         in2_begin++;
      }
      else {
         in1_begin++;
         in2_begin++;
      }
   }

   out.insert(out.end(), in1_begin, in1_end);
}



//----------------------------------------------------------------------
// generate_hash_seed
//----------------------------------------------------------------------
void C_correct_errors::generate_hash_seed(const bloom_type random_seed, std::vector<unsigned int>& hash_seed) {
   hash_seed.resize(num_hash_func);

   srand(static_cast<unsigned int>(random_seed));

   for (unsigned short int it_seed = 0; it_seed < num_hash_func; it_seed++) {
      hash_seed[it_seed] = static_cast<unsigned int>(rand());
   }
}



//----------------------------------------------------------------------
// clear_path
//----------------------------------------------------------------------
void C_candidate_path::clear_path() {
   modified_bases.clear();
   sum_qs = 0;
}



//----------------------------------------------------------------------
// correct_errors_in_reads_single_fastq_gzipped
//----------------------------------------------------------------------
void C_correct_errors::correct_errors_in_reads_single_fastq_gzipped(const C_arg& c_inst_args) {
   //----------------------------------------------------------------------
   // load bloom filter data
   //----------------------------------------------------------------------
   // generate unique hash seeds
   std::vector<unsigned int> hash_seed;
   // no bloom filter data load
   // use a command line random seed
   if (c_inst_args.load_bf == false) {
      generate_hash_seed(c_inst_args.random_seed, hash_seed);
   }
   // load bloom filter data
   // use a loaded random seed
   else {
      generate_hash_seed(random_seed, hash_seed);
   }

   // generate a bloom filter
   unsigned char* bit_vector(new unsigned char[static_cast<std::size_t>(bit_vector_width_byte)]);
   std::fill_n(bit_vector, bit_vector_width_byte, 0x00);

   // open the bloom filter data file
   std::string bf_data_file_name;
   std::string bf_size_file_name;
   if (c_inst_args.load_bf == false) {
      bf_data_file_name = c_inst_args.bf_data_file_name;
      bf_size_file_name = c_inst_args.bf_size_file_name;
   }
   else {
      bf_data_file_name = c_inst_args.loaded_bf_data_file_name;
      bf_size_file_name = c_inst_args.loaded_bf_size_file_name;
   }

   std::ifstream f_bf_dump_data;
   f_bf_dump_data.open(bf_data_file_name.c_str(), std::ios::binary);
   if (f_bf_dump_data.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << bf_data_file_name << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 211);
   }

   // load bit vector data
   f_bf_dump_data.read(reinterpret_cast<char*>(bit_vector), bit_vector_width_byte);

   f_bf_dump_data.close();

   if (rank_node == 0) {
      // calculate the size of the bloom filter in megabyte
      bloom_type bit_vector_size_mb(bit_vector_width_byte / (1024 * 1024));

      if (c_inst_args.load_bf == true) {
         std::cout << "     Bloom filter data file    : " << c_inst_args.bf_data_file_name << std::endl;
         std::cout << "     Bloom filter size file    : " << c_inst_args.bf_size_file_name << std::endl;
         std::cout << "     Number of keys            : " << num_unique_solid_kmers << std::endl;
         std::cout << "     Bit-vector size           : " << bit_vector_size_mb << " MB" << std::endl;
         std::cout << "     Number of hash functions  : " << num_hash_func << std::endl;
         std::cout << "     k-mer threshold           : " << kmer_occurrence_threshold << std::endl;
      }
   }

   //----------------------------------------------------------------------
   // correct errors
   //----------------------------------------------------------------------
   //
   // variables
   //
   std::vector<std::string> read_vector;

   std::string sequence_modification;
   std::string trimmed_seq;
   std::string trimmed_qs;
   std::string buffer_not_needed;

   std::size_t read_length;
   std::size_t min_check_length;
   std::size_t max_allowed_ns;
   std::size_t max_trimmed_bases;
   std::size_t num_corrected_errors_local;
   std::size_t trim_5_end;
   std::size_t trim_3_end;
   std::size_t num_corrected_errors_tmp;
   std::size_t num_corrected_reads_tmp;
   std::size_t num_trimmed_bases_tmp;
   std::size_t read_vector_index;
   std::size_t current_read_index;
   std::size_t current_read_index_write;
   std::size_t end_read_prev_rank;
   std::size_t num_reads_local;

   bool too_many_errors;

   std::regex regex_non_acgtn("[^ACGTN]");

   read_vector.resize(read_block_size * 3);

   num_corrected_errors_tmp = 0;
   num_corrected_reads_tmp  = 0;
   num_trimmed_bases_tmp    = 0;

   // initialize variables
   num_reads_local          = 0;
   read_vector_index        = 0;
   current_read_index_write = 0;

   // calculate end_read_prev_rank
   // the current node rank == 0
   if (rank_node == 0) {
      end_read_prev_rank = 0;
   }
   // the current node rank > 0
   else {
      // if no read is assigned to this node
      // finding a start point is not needed
      if (num_reads_vector[rank_node] == 0) {
         end_read_prev_rank = 0;
      }
      else {
         for (int it_rank = 0; it_rank < rank_node; it_rank++) {
            end_read_prev_rank += num_reads_vector[it_rank];
         }
      }
   }

   // open a temporary output file
   std::string corrected_read_file_name(c_inst_args.corrected_read_file_name + '.' + rank_node_text);

   FILE* f_corrected_read(fopen(corrected_read_file_name.c_str(), "w"));
   if (f_corrected_read == NULL) {
      std::cout << std::endl << "ERROR: Cannot open " << corrected_read_file_name << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 212);
   }

   // open an input read file
   std::ifstream f_read;
   f_read.open(c_inst_args.read_file_name.c_str(), std::ios_base::binary);

   if (!f_read.is_open()) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.read_file_name << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 214);
   }

   // set a stream filter for gzipped files
   boost::iostreams::filtering_istream f_read_filter;

   f_read_filter.push(boost::iostreams::gzip_decompressor());
   f_read_filter.push(f_read);

   // skip unnecessary reads
   for (std::size_t it_read_order = 0; it_read_order < end_read_prev_rank; it_read_order++) {
      std::getline(f_read_filter, buffer_not_needed);
      std::getline(f_read_filter, buffer_not_needed);
      std::getline(f_read_filter, buffer_not_needed);
      std::getline(f_read_filter, buffer_not_needed);
   }

   // process reads
   for (std::size_t it_read_order = 0; it_read_order < num_reads_vector[rank_node]; it_read_order++) {
      // calculate the index for read_vector
      current_read_index_write = read_vector_index * 3;

      //
      // header
      //
      std::getline(f_read_filter, read_vector[current_read_index_write]);

      //
      // sequence
      //
      std::getline(f_read_filter, read_vector[current_read_index_write + 1]);

      //
      // connector
      //
      std::getline(f_read_filter, buffer_not_needed);

      //
      // quality score
      //
      std::getline(f_read_filter, read_vector[current_read_index_write + 2]);

      read_vector_index++;
      num_reads_local++;

      // read_vector is full
      if (read_vector_index == read_block_size) {
         //--------------------------------------------------
         // correct reads in a block
         //--------------------------------------------------
         #pragma omp parallel shared(bit_vector, hash_seed, read_vector) private(current_read_index, read_length, min_check_length, max_allowed_ns, max_trimmed_bases, num_corrected_errors_local, trim_5_end, trim_3_end, sequence_modification, too_many_errors) reduction(+: num_corrected_errors_tmp, num_corrected_reads_tmp, num_trimmed_bases_tmp)
         {
            // iterate reads
            #pragma omp for schedule(dynamic, OPENMP_CHUNK_SIZE)
            for (std::size_t it_read = 0; it_read < read_block_size; it_read++) {
               // calculate the current index
               current_read_index = it_read * 3;

               // change sequences to upper case
               std::transform(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), read_vector[current_read_index + 1].begin(), static_cast<int(*)(int)>(std::toupper));

               // substitute non-standard characters with Ns
               read_vector[current_read_index + 1] = std::regex_replace(read_vector[current_read_index + 1], regex_non_acgtn, "N");

               // set various thresholds
               read_length      = read_vector[current_read_index + 1].length();
               min_check_length = read_length * CHECK_RANGE_RATIO;
               max_allowed_ns   = read_length * MAX_N_RATIO;

               // forward
               // too short read: no trimming
               if (read_length <= MIN_BASES_AFTER_TRIMMING) {
                  max_trimmed_bases = 0;
               }
               else {
                  max_trimmed_bases = std::min((std::size_t)(read_length * MAX_TRIMMING_RATE), read_length - MIN_BASES_AFTER_TRIMMING);
               }

               // check the number of Ns in the read
               if (((std::size_t)std::count(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), 'N') <= max_allowed_ns) &&
                   (read_length >= kmer_length)) {
                  // substitute Ns other characters
                  std::replace(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), 'N', SUBST_CHAR);

                  // initialize variables
                  num_corrected_errors_local = 0;
                  trim_5_end                 = 0;
                  trim_3_end                 = 0;

                  // storage for the modification of the reads
                  // # of entries = read length
                  sequence_modification.assign(read_length, '0');

                  //----------------------------------------------------------------------
                  // correct errors in a read
                  //----------------------------------------------------------------------
                  // forward read
                  // errors cannot be corrected if k is equal to read length
                  if (read_length > kmer_length) {
                     correct_errors_in_a_read_fastq(
                                                    read_vector[current_read_index + 1],
                                                    sequence_modification,
                                                    read_vector[current_read_index + 2],
                                                    trim_5_end,
                                                    trim_3_end,
                                                    read_length,
                                                    max_trimmed_bases,
                                                    min_check_length,
                                                    num_corrected_errors_local,
                                                    bit_vector,
                                                    hash_seed
                                                   );
                  }

                  // no trim
                  if (c_inst_args.notrim == true) {
                     trim_5_end = 0;
                     trim_3_end = 0;
                  }
                  // adjust the number of trimmed bases
                  else {
                     if ((trim_5_end + trim_3_end) > max_trimmed_bases) {
                        if (trim_3_end <= max_trimmed_bases) {
                           trim_5_end = 0;
                        }
                        else if (trim_5_end <= max_trimmed_bases) {
                           trim_3_end = 0;
                        }
                        else {
                           trim_5_end = 0;
                           trim_3_end = max_trimmed_bases;
                        }
                     }
                  }

                  num_trimmed_bases_tmp += (trim_5_end + trim_3_end);

                  // update num_corrected_reads
                  too_many_errors = false;
                  if (num_corrected_errors_local > (read_vector[current_read_index + 1].length() * MAX_ERROR_RATE)) {
                     too_many_errors = true;
                  }
                  else if (num_corrected_errors_local > 0) {
                     num_corrected_errors_tmp += num_corrected_errors_local;

                     num_corrected_reads_tmp++;
                  }
                  else if (c_inst_args.notrim == false) {
                     if ((trim_5_end > 0) || (trim_3_end > 0)) {
                        num_corrected_reads_tmp++;
                     }
                  }

                  // make a corrected read
                  if (too_many_errors == false) {
                     // apply modifications to the read
                     for (unsigned int it_base = trim_5_end; it_base < (read_vector[current_read_index + 1].length() - trim_3_end); it_base++) {
                        if (sequence_modification[it_base] != '0') {
                           read_vector[current_read_index + 1][it_base] = sequence_modification[it_base];
                        }
                     }
                  }

                  read_vector[current_read_index] = read_vector[current_read_index] + '\n';

                  // make a trimmed read
                  if ((trim_5_end + trim_3_end) > 0) {
                     read_vector[current_read_index + 1] = read_vector[current_read_index + 1].substr(trim_5_end, read_vector[current_read_index + 1].length() - trim_5_end - trim_3_end) + '\n';
                     read_vector[current_read_index + 2] = read_vector[current_read_index + 2].substr(trim_5_end, read_vector[current_read_index + 2].length() - trim_5_end - trim_3_end) + '\n';
                  }
                  else {
                     read_vector[current_read_index + 1] = read_vector[current_read_index + 1] + '\n';
                     read_vector[current_read_index + 2] = read_vector[current_read_index + 2] + '\n';
                  }
               }
               // too many Ns
               // write reads without modification
               else {
                  // sequence
                  read_vector[current_read_index + 1] += '\n';

                  // quality score
                  read_vector[current_read_index + 2] += '\n';
               }
            // it_read
            }
         // omp parallel
         }

         // write corrected_reads
         for (std::size_t it_write = 0; it_write < read_block_size; it_write++) {
            current_read_index = it_write * 3;

            // header
            fwrite(read_vector[current_read_index].c_str(), 1, read_vector[current_read_index].length(), f_corrected_read);
            // sequence
            fwrite(read_vector[current_read_index + 1].c_str(), 1, read_vector[current_read_index + 1].length(), f_corrected_read);
            // connector
            fwrite("+\n", 1, 2, f_corrected_read);
            // quality score
            fwrite(read_vector[current_read_index + 2].c_str(), 1, read_vector[current_read_index + 2].length(), f_corrected_read);
         }

         read_vector_index = 0;
      }
   }

   f_read.close();

   // correct errors in remaining reads
   if (read_vector_index > 0) {
      #pragma omp parallel shared(bit_vector, hash_seed, read_vector) private(current_read_index, read_length, min_check_length, max_allowed_ns, max_trimmed_bases, num_corrected_errors_local, trim_5_end, trim_3_end, sequence_modification, too_many_errors) reduction(+: num_corrected_errors_tmp, num_corrected_reads_tmp, num_trimmed_bases_tmp)
      {
         // iterate reads
         #pragma omp for schedule(dynamic, OPENMP_CHUNK_SIZE)
         for (std::size_t it_read = 0; it_read < read_vector_index; it_read++) {
            // calculate the current index
            current_read_index = it_read * 3;

            // change sequences to upper case
            std::transform(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), read_vector[current_read_index + 1].begin(), static_cast<int(*)(int)>(std::toupper));

            // substitute non-standard characters with Ns
            read_vector[current_read_index + 1] = std::regex_replace(read_vector[current_read_index + 1], regex_non_acgtn, "N");

            // set various thresholds
            read_length      = read_vector[current_read_index + 1].length();
            min_check_length = read_length * CHECK_RANGE_RATIO;
            max_allowed_ns   = read_length * MAX_N_RATIO;

            // forward
            // too short read: no trimming
            if (read_length <= MIN_BASES_AFTER_TRIMMING) {
               max_trimmed_bases = 0;
            }
            else {
               max_trimmed_bases = std::min((std::size_t)(read_length * MAX_TRIMMING_RATE), read_length - MIN_BASES_AFTER_TRIMMING);
            }

            // check the number of Ns in the read
            if (((std::size_t)std::count(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), 'N') <= max_allowed_ns) &&
                (read_length >= kmer_length)) {
               // substitute Ns other characters
               std::replace(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), 'N', SUBST_CHAR);

               // initialize variables
               num_corrected_errors_local = 0;
               trim_5_end                 = 0;
               trim_3_end                 = 0;

               // storage for the modification of the reads
               // # of entries = read length
               sequence_modification.assign(read_length, '0');

               //----------------------------------------------------------------------
               // correct errors in a read
               //----------------------------------------------------------------------
               // forward read
               // errors cannot be corrected if k is equal to read length
               if (read_length > kmer_length) {
                  correct_errors_in_a_read_fastq(
                                                 read_vector[current_read_index + 1],
                                                 sequence_modification,
                                                 read_vector[current_read_index + 2],
                                                 trim_5_end,
                                                 trim_3_end,
                                                 read_length,
                                                 max_trimmed_bases,
                                                 min_check_length,
                                                 num_corrected_errors_local,
                                                 bit_vector,
                                                 hash_seed
                                                );
               }

               // no trim
               if (c_inst_args.notrim == true) {
                  trim_5_end = 0;
                  trim_3_end = 0;
               }
               // adjust the number of trimmed bases
               else {
                  if ((trim_5_end + trim_3_end) > max_trimmed_bases) {
                     if (trim_3_end <= max_trimmed_bases) {
                        trim_5_end = 0;
                     }
                     else if (trim_5_end <= max_trimmed_bases) {
                        trim_3_end = 0;
                     }
                     else {
                        trim_5_end = 0;
                        trim_3_end = max_trimmed_bases;
                     }
                  }
               }

               num_trimmed_bases_tmp += (trim_5_end + trim_3_end);

               // update num_corrected_reads
               too_many_errors = false;
               if (num_corrected_errors_local > (read_vector[current_read_index + 1].length() * MAX_ERROR_RATE)) {
                  too_many_errors = true;
               }
               else if (num_corrected_errors_local > 0) {
                  num_corrected_errors_tmp += num_corrected_errors_local;

                  num_corrected_reads_tmp++;
               }
               else if (c_inst_args.notrim == false) {
                  if ((trim_5_end > 0) || (trim_3_end > 0)) {
                     num_corrected_reads_tmp++;
                  }
               }

               // make a corrected read
               if (too_many_errors == false) {
                  // apply modifications to the read
                  for (unsigned int it_base = trim_5_end; it_base < (read_vector[current_read_index + 1].length() - trim_3_end); it_base++) {
                     if (sequence_modification[it_base] != '0') {
                        read_vector[current_read_index + 1][it_base] = sequence_modification[it_base];
                     }
                  }
               }

               read_vector[current_read_index] = read_vector[current_read_index] + '\n';

               // make a trimmed quality score
               if ((trim_5_end + trim_3_end) > 0) {
                  read_vector[current_read_index + 1] = read_vector[current_read_index + 1].substr(trim_5_end, read_vector[current_read_index + 1].length() - trim_5_end - trim_3_end) + '\n';
                  read_vector[current_read_index + 2] = read_vector[current_read_index + 2].substr(trim_5_end, read_vector[current_read_index + 2].length() - trim_5_end - trim_3_end) + '\n';
               }
               else {
                  read_vector[current_read_index + 1] = read_vector[current_read_index + 1] + '\n';
                  read_vector[current_read_index + 2] = read_vector[current_read_index + 2] + '\n';
               }
            }
            // too many Ns
            // write reads without modification
            else {
               // sequence
               read_vector[current_read_index + 1] += '\n';

               // quality score
               read_vector[current_read_index + 2] += '\n';
            }
         // it_read
         }

      // omp parallel
      }

      // write corrected_reads
      for (std::size_t it_read = 0; it_read < read_vector_index; it_read++) {
         current_read_index = it_read * 3;

         // header
         fwrite(read_vector[current_read_index].c_str(), 1, read_vector[current_read_index].length(), f_corrected_read);
         // sequence
         fwrite(read_vector[current_read_index + 1].c_str(), 1, read_vector[current_read_index + 1].length(), f_corrected_read);
         // connector
         fwrite("+\n", 1, 2, f_corrected_read);
         // quality score
         fwrite(read_vector[current_read_index + 2].c_str(), 1, read_vector[current_read_index + 2].length(), f_corrected_read);
      }
   }

   fclose(f_corrected_read);

   read_vector.clear();

   num_corrected_errors = num_corrected_errors_tmp;
   num_corrected_reads  = num_corrected_reads_tmp;
   num_trimmed_bases    = num_trimmed_bases_tmp;
}



//----------------------------------------------------------------------
// correct_errors_in_reads_paired_fastq_gzipped
//----------------------------------------------------------------------
void C_correct_errors::correct_errors_in_reads_paired_fastq_gzipped(const C_arg& c_inst_args) {
   //----------------------------------------------------------------------
   // load bloom filter data
   //----------------------------------------------------------------------
   // generate unique hash seeds
   std::vector<unsigned int> hash_seed;
   // no bloom filter data load
   // use a command line random seed
   if (c_inst_args.load_bf == false) {
      generate_hash_seed(c_inst_args.random_seed, hash_seed);
   }
   // load bloom filter data
   // use a loaded random seed
   else {
      generate_hash_seed(random_seed, hash_seed);
   }

   // generate a bloom filter
   unsigned char* bit_vector(new unsigned char[static_cast<std::size_t>(bit_vector_width_byte)]);
   std::fill_n(bit_vector, bit_vector_width_byte, 0x00);

   // open the bloom filter data file
   std::string bf_data_file_name;
   std::string bf_size_file_name;
   if (c_inst_args.load_bf == false) {
      bf_data_file_name = c_inst_args.bf_data_file_name;
      bf_size_file_name = c_inst_args.bf_size_file_name;
   }
   else {
      bf_data_file_name = c_inst_args.loaded_bf_data_file_name;
      bf_size_file_name = c_inst_args.loaded_bf_size_file_name;
   }

   std::ifstream f_bf_dump_data;
   f_bf_dump_data.open(bf_data_file_name.c_str(), std::ios::binary);
   if (f_bf_dump_data.is_open() == false) {
      std::cout << std::endl << "ERROR: Cannot open " << bf_data_file_name << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 216);
   }

   // load bit vector data
   f_bf_dump_data.read(reinterpret_cast<char*>(bit_vector), bit_vector_width_byte);

   f_bf_dump_data.close();

   if (rank_node == 0) {
      // calculate the size of the bloom filter in megabyte
      bloom_type bit_vector_size_mb(bit_vector_width_byte / (1024 * 1024));

      if (c_inst_args.load_bf == true) {
         std::cout << "     Bloom filter data file    : " << c_inst_args.bf_data_file_name << std::endl;
         std::cout << "     Bloom filter size file    : " << c_inst_args.bf_size_file_name << std::endl;
         std::cout << "     Number of keys            : " << num_unique_solid_kmers << std::endl;
         std::cout << "     Bit-vector size           : " << bit_vector_size_mb << " MB" << std::endl;
         std::cout << "     Number of hash functions  : " << num_hash_func << std::endl;
         std::cout << "     k-mer threshold           : " << kmer_occurrence_threshold << std::endl;
      }
   }

   //--------------------------------------------------
   // correct errors
   //--------------------------------------------------
   //
   // variables
   //
   std::ifstream f_read;

   std::vector<std::string> read_vector;

   std::string sequence_modification;
   std::string trimmed_seq;
   std::string trimmed_qs;
   std::string buffer_not_needed;

   std::size_t read_length;
   std::size_t min_check_length;
   std::size_t max_allowed_ns;
   std::size_t max_trimmed_bases;
   std::size_t num_corrected_errors_local;
   std::size_t trim_5_end;
   std::size_t trim_3_end;
   std::size_t num_corrected_errors_tmp;
   std::size_t num_corrected_reads_tmp;
   std::size_t num_trimmed_bases_tmp;
   std::size_t read_vector_index;
   std::size_t current_read_index;
   std::size_t current_read_index_write;
   std::size_t end_read_prev_rank;
   std::size_t num_reads_local;

   bool too_many_errors;

   std::regex regex_non_acgtn("[^ACGTN]");

   boost::iostreams::filtering_istream f_read_filter;

   // initialize variables
   read_vector.resize(read_block_size * 3);

   num_corrected_errors_tmp = 0;
   num_corrected_reads_tmp  = 0;
   num_trimmed_bases_tmp    = 0;

   //
   // forward
   //
   // open a temporary output file
   std::string corrected_read_file_name1(c_inst_args.corrected_read_file_name1 + '.' + rank_node_text);

   FILE* f_corrected_read1(fopen(corrected_read_file_name1.c_str(), "w"));
   if (f_corrected_read1 == NULL) {
      std::cout << std::endl << "ERROR: Cannot open " << corrected_read_file_name1 << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 217);
   }

   // open an input read file
   // open the file
   f_read.open(c_inst_args.read_file_name1.c_str(), std::ios_base::binary);

   if (!f_read.is_open()) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.read_file_name1 << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 214);
   }

   // set a stream filter for gzipped files
   f_read_filter.reset();
   f_read_filter.push(boost::iostreams::gzip_decompressor());
   f_read_filter.push(f_read);

   // initialize variables
   num_reads_local          = 0;
   read_vector_index        = 0;
   current_read_index_write = 0;

   // calculate end_read_prev_rank
   // the current node rank == 0
   if (rank_node == 0) {
      end_read_prev_rank = 0;
   }
   // the current node rank > 0
   else {
      // if no read is assigned to this node
      // finding a start point is not needed
      if (num_reads_vector1[rank_node] == 0) {
         end_read_prev_rank = 0;
      }
      else {
         for (int it_rank = 0; it_rank < rank_node; it_rank++) {
            end_read_prev_rank += num_reads_vector1[it_rank];
         }
      }
   }

   // skip unnecessary reads
   for (std::size_t it_read_order = 0; it_read_order < end_read_prev_rank; it_read_order++) {
      std::getline(f_read_filter, buffer_not_needed);
      std::getline(f_read_filter, buffer_not_needed);
      std::getline(f_read_filter, buffer_not_needed);
      std::getline(f_read_filter, buffer_not_needed);
   }

   // process reads
   for (std::size_t it_read_order = 0; it_read_order < num_reads_vector1[rank_node]; it_read_order++) {
      // calculate the index for read_vector
      current_read_index_write = read_vector_index * 3;

      //
      // header
      //
      std::getline(f_read_filter, read_vector[current_read_index_write]);

      //
      // sequence
      //
      std::getline(f_read_filter, read_vector[current_read_index_write + 1]);

      //
      // connector
      //
      std::getline(f_read_filter, buffer_not_needed);

      //
      // quality score
      //
      std::getline(f_read_filter, read_vector[current_read_index_write + 2]);

      read_vector_index++;
      num_reads_local++;

      // read_vector is full
      if (read_vector_index == read_block_size) {
         //--------------------------------------------------
         // correct reads in a block
         //--------------------------------------------------
         #pragma omp parallel shared(bit_vector, hash_seed, read_vector) private(current_read_index, read_length, min_check_length, max_allowed_ns, max_trimmed_bases, num_corrected_errors_local, trim_5_end, trim_3_end, sequence_modification, too_many_errors) reduction(+: num_corrected_errors_tmp, num_corrected_reads_tmp, num_trimmed_bases_tmp)
         {
            // iterate reads
            #pragma omp for schedule(dynamic, OPENMP_CHUNK_SIZE)
            for (std::size_t it_read = 0; it_read < read_block_size; it_read++) {
               // calculate the current index
               current_read_index = it_read * 3;

               // change sequences to upper case
               std::transform(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), read_vector[current_read_index + 1].begin(), static_cast<int(*)(int)>(std::toupper));

               // substitute non-standard characters with Ns
               read_vector[current_read_index + 1] = std::regex_replace(read_vector[current_read_index + 1], regex_non_acgtn, "N");

               // set various thresholds
               read_length      = read_vector[current_read_index + 1].length();
               min_check_length = read_length * CHECK_RANGE_RATIO;
               max_allowed_ns   = read_length * MAX_N_RATIO;

               // too short read: no trimming
               if (read_length <= MIN_BASES_AFTER_TRIMMING) {
                  max_trimmed_bases = 0;
               }
               else {
                  max_trimmed_bases = std::min((std::size_t)(read_length * MAX_TRIMMING_RATE), read_length - MIN_BASES_AFTER_TRIMMING);
               }

               // check the number of Ns in the read
               if (((std::size_t)std::count(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), 'N') <= max_allowed_ns) &&
                   (read_length >= kmer_length)) {
                  // substitute Ns other characters
                  std::replace(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), 'N', SUBST_CHAR);

                  // initialize variables
                  num_corrected_errors_local = 0;
                  trim_5_end                 = 0;
                  trim_3_end                 = 0;

                  // storage for the modification of the reads
                  // # of entries = read length
                  sequence_modification.assign(read_length, '0');

                  //----------------------------------------------------------------------
                  // correct errors in a read
                  //----------------------------------------------------------------------
                  // forward read
                  // errors cannot be corrected if k is equal to read length
                  if (read_length > kmer_length) {
                     correct_errors_in_a_read_fastq(
                                                    read_vector[current_read_index + 1],
                                                    sequence_modification,
                                                    read_vector[current_read_index + 2],
                                                    trim_5_end,
                                                    trim_3_end,
                                                    read_length,
                                                    max_trimmed_bases,
                                                    min_check_length,
                                                    num_corrected_errors_local,
                                                    bit_vector,
                                                    hash_seed
                                                   );
                  }

                  // no trim
                  if (c_inst_args.notrim == true) {
                     trim_5_end = 0;
                     trim_3_end = 0;
                  }
                  // adjust the number of trimmed bases
                  else {
                     if ((trim_5_end + trim_3_end) > max_trimmed_bases) {
                        if (trim_3_end <= max_trimmed_bases) {
                           trim_5_end = 0;
                        }
                        else if (trim_5_end <= max_trimmed_bases) {
                           trim_3_end = 0;
                        }
                        else {
                           trim_5_end = 0;
                           trim_3_end = max_trimmed_bases;
                        }
                     }
                  }

                  num_trimmed_bases_tmp += (trim_5_end + trim_3_end);

                  // update num_corrected_reads
                  too_many_errors = false;
                  if (num_corrected_errors_local > (read_vector[current_read_index + 1].length() * MAX_ERROR_RATE)) {
                     too_many_errors = true;
                  }
                  else if (num_corrected_errors_local > 0) {
                     num_corrected_errors_tmp += num_corrected_errors_local;

                     num_corrected_reads_tmp++;
                  }
                  else if (c_inst_args.notrim == false) {
                     if ((trim_5_end > 0) || (trim_3_end > 0)) {
                        num_corrected_reads_tmp++;
                     }
                  }

                  // make a corrected read
                  if (too_many_errors == false) {
                     // apply modifications to the read
                     for (unsigned int it_base = trim_5_end; it_base < (read_vector[current_read_index + 1].length() - trim_3_end); it_base++) {
                        if (sequence_modification[it_base] != '0') {
                           read_vector[current_read_index + 1][it_base] = sequence_modification[it_base];
                        }
                     }
                  }

                  read_vector[current_read_index] = read_vector[current_read_index] + '\n';

                  // make a trimmed read
                  if ((trim_5_end + trim_3_end) > 0) {
                     read_vector[current_read_index + 1] = read_vector[current_read_index + 1].substr(trim_5_end, read_vector[current_read_index + 1].length() - trim_5_end - trim_3_end) + '\n';
                     read_vector[current_read_index + 2] = read_vector[current_read_index + 2].substr(trim_5_end, read_vector[current_read_index + 2].length() - trim_5_end - trim_3_end) + '\n';
                  }
                  else {
                     read_vector[current_read_index + 1] = read_vector[current_read_index + 1] + '\n';
                     read_vector[current_read_index + 2] = read_vector[current_read_index + 2] + '\n';
                  }
               }
               // too many Ns
               // write reads without modification
               else {
                  // sequence
                  read_vector[current_read_index + 1] += '\n';

                  // quality score
                  read_vector[current_read_index + 2] += '\n';
               }
            // it_read
            }
         // omp parallel
         }

         // write corrected_reads
         for (std::size_t it_write = 0; it_write < read_block_size; it_write++) {
            current_read_index = it_write * 3;

            // header
            fwrite(read_vector[current_read_index].c_str(), 1, read_vector[current_read_index].length(), f_corrected_read1);
            // sequence
            fwrite(read_vector[current_read_index + 1].c_str(), 1, read_vector[current_read_index + 1].length(), f_corrected_read1);
            // connector
            fwrite("+\n", 1, 2, f_corrected_read1);
            // quality score
            fwrite(read_vector[current_read_index + 2].c_str(), 1, read_vector[current_read_index + 2].length(), f_corrected_read1);
         }

         read_vector_index = 0;
      }
   }

   f_read.close();

   // correct errors in remaining reads
   if (read_vector_index > 0) {
      #pragma omp parallel shared(bit_vector, hash_seed, read_vector) private(current_read_index, read_length, min_check_length, max_allowed_ns, max_trimmed_bases, num_corrected_errors_local, trim_5_end, trim_3_end, sequence_modification, too_many_errors) reduction(+: num_corrected_errors_tmp, num_corrected_reads_tmp, num_trimmed_bases_tmp)
      {
         // iterate reads
         #pragma omp for schedule(dynamic, OPENMP_CHUNK_SIZE)
         for (std::size_t it_read = 0; it_read < read_vector_index; it_read++) {
            // calculate the current index
            current_read_index = it_read * 3;

            // change sequences to upper case
            std::transform(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), read_vector[current_read_index + 1].begin(), static_cast<int(*)(int)>(std::toupper));

            // substitute non-standard characters with Ns
            read_vector[current_read_index + 1] = std::regex_replace(read_vector[current_read_index + 1], regex_non_acgtn, "N");

            // set various thresholds
            read_length      = read_vector[current_read_index + 1].length();
            min_check_length = read_length * CHECK_RANGE_RATIO;
            max_allowed_ns   = read_length * MAX_N_RATIO;

            // too short read: no trimming
            if (read_length <= MIN_BASES_AFTER_TRIMMING) {
               max_trimmed_bases = 0;
            }
            else {
               max_trimmed_bases = std::min((std::size_t)(read_length * MAX_TRIMMING_RATE), read_length - MIN_BASES_AFTER_TRIMMING);
            }

            // check the number of Ns in the read
            if (((std::size_t)std::count(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), 'N') <= max_allowed_ns) &&
                (read_length >= kmer_length)) {
               // substitute Ns other characters
               std::replace(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), 'N', SUBST_CHAR);

               // initialize variables
               num_corrected_errors_local = 0;
               trim_5_end                 = 0;
               trim_3_end                 = 0;

               // storage for the modification of the reads
               // # of entries = read length
               sequence_modification.assign(read_length, '0');

               //----------------------------------------------------------------------
               // correct errors in a read
               //----------------------------------------------------------------------
               // forward read
               // errors cannot be corrected if k is equal to read length
               if (read_length > kmer_length) {
                  correct_errors_in_a_read_fastq(
                                                 read_vector[current_read_index + 1],
                                                 sequence_modification,
                                                 read_vector[current_read_index + 2],
                                                 trim_5_end,
                                                 trim_3_end,
                                                 read_length,
                                                 max_trimmed_bases,
                                                 min_check_length,
                                                 num_corrected_errors_local,
                                                 bit_vector,
                                                 hash_seed
                                                );
               }

               // no trim
               if (c_inst_args.notrim == true) {
                  trim_5_end = 0;
                  trim_3_end = 0;
               }
               // adjust the number of trimmed bases
               else {
                  if ((trim_5_end + trim_3_end) > max_trimmed_bases) {
                     if (trim_3_end <= max_trimmed_bases) {
                        trim_5_end = 0;
                     }
                     else if (trim_5_end <= max_trimmed_bases) {
                        trim_3_end = 0;
                     }
                     else {
                        trim_5_end = 0;
                        trim_3_end = max_trimmed_bases;
                     }
                  }
               }

               num_trimmed_bases_tmp += (trim_5_end + trim_3_end);

               // update num_corrected_reads
               too_many_errors = false;
               if (num_corrected_errors_local > (read_vector[current_read_index + 1].length() * MAX_ERROR_RATE)) {
                  too_many_errors = true;
               }
               else if (num_corrected_errors_local > 0) {
                  num_corrected_errors_tmp += num_corrected_errors_local;

                  num_corrected_reads_tmp++;
               }
               else if (c_inst_args.notrim == false) {
                  if ((trim_5_end > 0) || (trim_3_end > 0)) {
                     num_corrected_reads_tmp++;
                  }
               }

               // make a corrected read
               if (too_many_errors == false) {
                  // apply modifications to the read
                  for (unsigned int it_base = trim_5_end; it_base < (read_vector[current_read_index + 1].length() - trim_3_end); it_base++) {
                     if (sequence_modification[it_base] != '0') {
                        read_vector[current_read_index + 1][it_base] = sequence_modification[it_base];
                     }
                  }
               }

               read_vector[current_read_index] = read_vector[current_read_index] + '\n';

               // make a trimmed quality score
               if ((trim_5_end + trim_3_end) > 0) {
                  read_vector[current_read_index + 1] = read_vector[current_read_index + 1].substr(trim_5_end, read_vector[current_read_index + 1].length() - trim_5_end - trim_3_end) + '\n';
                  read_vector[current_read_index + 2] = read_vector[current_read_index + 2].substr(trim_5_end, read_vector[current_read_index + 2].length() - trim_5_end - trim_3_end) + '\n';
               }
               else {
                  read_vector[current_read_index + 1] = read_vector[current_read_index + 1] + '\n';
                  read_vector[current_read_index + 2] = read_vector[current_read_index + 2] + '\n';
               }
            }
            // too many Ns
            // write reads without modification
            else {
               // sequence
               read_vector[current_read_index + 1] += '\n';

               // quality score
               read_vector[current_read_index + 2] += '\n';
            }
         // it_read
         }

      // omp parallel
      }

      // write corrected_reads
      for (std::size_t it_read = 0; it_read < read_vector_index; it_read++) {
         current_read_index = it_read * 3;

         // header
         fwrite(read_vector[current_read_index].c_str(), 1, read_vector[current_read_index].length(), f_corrected_read1);
         // sequence
         fwrite(read_vector[current_read_index + 1].c_str(), 1, read_vector[current_read_index + 1].length(), f_corrected_read1);
         // connector
         fwrite("+\n", 1, 2, f_corrected_read1);
         // quality score
         fwrite(read_vector[current_read_index + 2].c_str(), 1, read_vector[current_read_index + 2].length(), f_corrected_read1);
      }
   }

   fclose(f_corrected_read1);

   //
   // reverse
   //
   // open a temporary output file
   std::string corrected_read_file_name2(c_inst_args.corrected_read_file_name2 + '.' + rank_node_text);

   FILE* f_corrected_read2(fopen(corrected_read_file_name2.c_str(), "w"));
   if (f_corrected_read2 == NULL) {
      std::cout << std::endl << "ERROR: Cannot open " << corrected_read_file_name2 << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 221);
   }

   // open an input read file
   // open the file
   f_read.open(c_inst_args.read_file_name2.c_str(), std::ios_base::binary);

   if (!f_read.is_open()) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.read_file_name2 << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 214);
   }

   // set a stream filter for gzipped files
   f_read_filter.reset();
   f_read_filter.push(boost::iostreams::gzip_decompressor());
   f_read_filter.push(f_read);

   // initialize variables
   num_reads_local          = 0;
   read_vector_index        = 0;
   current_read_index_write = 0;

   // calculate end_read_prev_rank
   // the current node rank == 0
   if (rank_node == 0) {
      end_read_prev_rank = 0;
   }
   // the current node rank > 0
   else {
      // if no read is assigned to this node
      // finding a start point is not needed
      if (num_reads_vector2[rank_node] == 0) {
         end_read_prev_rank = 0;
      }
      else {
         for (int it_rank = 0; it_rank < rank_node; it_rank++) {
            end_read_prev_rank += num_reads_vector2[it_rank];
         }
      }
   }

   // skip unnecessary reads
   for (std::size_t it_read_order = 0; it_read_order < end_read_prev_rank; it_read_order++) {
      std::getline(f_read_filter, buffer_not_needed);
      std::getline(f_read_filter, buffer_not_needed);
      std::getline(f_read_filter, buffer_not_needed);
      std::getline(f_read_filter, buffer_not_needed);
   }

   // process reads
   for (std::size_t it_read_order = 0; it_read_order < num_reads_vector2[rank_node]; it_read_order++) {
      // calculate the index for read_vector
      current_read_index_write = read_vector_index * 3;

      //
      // header
      //
      std::getline(f_read_filter, read_vector[current_read_index_write]);

      //
      // sequence
      //
      std::getline(f_read_filter, read_vector[current_read_index_write + 1]);

      //
      // connector
      //
      std::getline(f_read_filter, buffer_not_needed);

      //
      // quality score
      //
      std::getline(f_read_filter, read_vector[current_read_index_write + 2]);

      read_vector_index++;
      num_reads_local++;

      // read_vector is full
      if (read_vector_index == read_block_size) {
         //--------------------------------------------------
         // correct reads in a block
         //--------------------------------------------------
         #pragma omp parallel shared(bit_vector, hash_seed, read_vector) private(current_read_index, read_length, min_check_length, max_allowed_ns, max_trimmed_bases, num_corrected_errors_local, trim_5_end, trim_3_end, sequence_modification, too_many_errors) reduction(+: num_corrected_errors_tmp, num_corrected_reads_tmp, num_trimmed_bases_tmp)
         {
            // iterate reads
            #pragma omp for schedule(dynamic, OPENMP_CHUNK_SIZE)
            for (std::size_t it_read = 0; it_read < read_block_size; it_read++) {
               // calculate the current index
               current_read_index = it_read * 3;

               // change sequences to upper case
               std::transform(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), read_vector[current_read_index + 1].begin(), static_cast<int(*)(int)>(std::toupper));

               // substitute non-standard characters with Ns
               read_vector[current_read_index + 1] = std::regex_replace(read_vector[current_read_index + 1], regex_non_acgtn, "N");

               // set various thresholds
               read_length      = read_vector[current_read_index + 1].length();
               min_check_length = read_length * CHECK_RANGE_RATIO;
               max_allowed_ns   = read_length * MAX_N_RATIO;

               // too short read: no trimming
               if (read_length <= MIN_BASES_AFTER_TRIMMING) {
                  max_trimmed_bases = 0;
               }
               else {
                  max_trimmed_bases = std::min((std::size_t)(read_length * MAX_TRIMMING_RATE), read_length - MIN_BASES_AFTER_TRIMMING);
               }

               // check the number of Ns in the read
               if (((std::size_t)std::count(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), 'N') <= max_allowed_ns) &&
                   (read_length >= kmer_length)) {
                  // substitute Ns other characters
                  std::replace(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), 'N', SUBST_CHAR);

                  // initialize variables
                  num_corrected_errors_local = 0;
                  trim_5_end                 = 0;
                  trim_3_end                 = 0;

                  // storage for the modification of the reads
                  // # of entries = read length
                  sequence_modification.assign(read_length, '0');

                  //----------------------------------------------------------------------
                  // correct errors in a read
                  //----------------------------------------------------------------------
                  // forward read
                  // errors cannot be corrected if k is equal to read length
                  if (read_length > kmer_length) {
                     correct_errors_in_a_read_fastq(
                                                    read_vector[current_read_index + 1],
                                                    sequence_modification,
                                                    read_vector[current_read_index + 2],
                                                    trim_5_end,
                                                    trim_3_end,
                                                    read_length,
                                                    max_trimmed_bases,
                                                    min_check_length,
                                                    num_corrected_errors_local,
                                                    bit_vector,
                                                    hash_seed
                                                   );
                  }

                  // no trim
                  if (c_inst_args.notrim == true) {
                     trim_5_end = 0;
                     trim_3_end = 0;
                  }
                  // adjust the number of trimmed bases
                  else {
                     if ((trim_5_end + trim_3_end) > max_trimmed_bases) {
                        if (trim_3_end <= max_trimmed_bases) {
                           trim_5_end = 0;
                        }
                        else if (trim_5_end <= max_trimmed_bases) {
                           trim_3_end = 0;
                        }
                        else {
                           trim_5_end = 0;
                           trim_3_end = max_trimmed_bases;
                        }
                     }
                  }

                  num_trimmed_bases_tmp += (trim_5_end + trim_3_end);

                  // update num_corrected_reads
                  too_many_errors = false;
                  if (num_corrected_errors_local > (read_vector[current_read_index + 1].length() * MAX_ERROR_RATE)) {
                     too_many_errors = true;
                  }
                  else if (num_corrected_errors_local > 0) {
                     num_corrected_errors_tmp += num_corrected_errors_local;

                     num_corrected_reads_tmp++;
                  }
                  else if (c_inst_args.notrim == false) {
                     if ((trim_5_end > 0) || (trim_3_end > 0)) {
                        num_corrected_reads_tmp++;
                     }
                  }

                  // make a corrected read
                  if (too_many_errors == false) {
                     // apply modifications to the read
                     for (unsigned int it_base = trim_5_end; it_base < (read_vector[current_read_index + 1].length() - trim_3_end); it_base++) {
                        if (sequence_modification[it_base] != '0') {
                           read_vector[current_read_index + 1][it_base] = sequence_modification[it_base];
                        }
                     }
                  }

                  read_vector[current_read_index] = read_vector[current_read_index] + '\n';

                  // make a trimmed read
                  if ((trim_5_end + trim_3_end) > 0) {
                     read_vector[current_read_index + 1] = read_vector[current_read_index + 1].substr(trim_5_end, read_vector[current_read_index + 1].length() - trim_5_end - trim_3_end) + '\n';
                     read_vector[current_read_index + 2] = read_vector[current_read_index + 2].substr(trim_5_end, read_vector[current_read_index + 2].length() - trim_5_end - trim_3_end) + '\n';
                  }
                  else {
                     read_vector[current_read_index + 1] = read_vector[current_read_index + 1] + '\n';
                     read_vector[current_read_index + 2] = read_vector[current_read_index + 2] + '\n';
                  }
               }
               // too many Ns
               // write reads without modification
               else {
                  // sequence
                  read_vector[current_read_index + 1] += '\n';

                  // quality score
                  read_vector[current_read_index + 2] += '\n';
               }
            // it_read
            }
         // omp parallel
         }

         // write corrected_reads
         for (std::size_t it_write = 0; it_write < read_block_size; it_write++) {
            current_read_index = it_write * 3;

            // header
            fwrite(read_vector[current_read_index].c_str(), 1, read_vector[current_read_index].length(), f_corrected_read2);
            // sequence
            fwrite(read_vector[current_read_index + 1].c_str(), 1, read_vector[current_read_index + 1].length(), f_corrected_read2);
            // connector
            fwrite("+\n", 1, 2, f_corrected_read2);
            // quality score
            fwrite(read_vector[current_read_index + 2].c_str(), 1, read_vector[current_read_index + 2].length(), f_corrected_read2);
         }

         read_vector_index = 0;
      }
   }

   f_read.close();

   // correct errors in remaining reads
   if (read_vector_index > 0) {
      #pragma omp parallel shared(bit_vector, hash_seed, read_vector) private(current_read_index, read_length, min_check_length, max_allowed_ns, max_trimmed_bases, num_corrected_errors_local, trim_5_end, trim_3_end, sequence_modification, too_many_errors) reduction(+: num_corrected_errors_tmp, num_corrected_reads_tmp, num_trimmed_bases_tmp)
      {
         // iterate reads
         #pragma omp for schedule(dynamic, OPENMP_CHUNK_SIZE)
         for (std::size_t it_read = 0; it_read < read_vector_index; it_read++) {
            // calculate the current index
            current_read_index = it_read * 3;

            // change sequences to upper case
            std::transform(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), read_vector[current_read_index + 1].begin(), static_cast<int(*)(int)>(std::toupper));

            // substitute non-standard characters with Ns
            read_vector[current_read_index + 1] = std::regex_replace(read_vector[current_read_index + 1], regex_non_acgtn, "N");

            // set various thresholds
            read_length      = read_vector[current_read_index + 1].length();
            min_check_length = read_length * CHECK_RANGE_RATIO;
            max_allowed_ns   = read_length * MAX_N_RATIO;

            // too short read: no trimming
            if (read_length <= MIN_BASES_AFTER_TRIMMING) {
               max_trimmed_bases = 0;
            }
            else {
               max_trimmed_bases = std::min((std::size_t)(read_length * MAX_TRIMMING_RATE), read_length - MIN_BASES_AFTER_TRIMMING);
            }

            // check the number of Ns in the read
            if (((std::size_t)std::count(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), 'N') <= max_allowed_ns) &&
                (read_length >= kmer_length)) {
               // substitute Ns other characters
               std::replace(read_vector[current_read_index + 1].begin(), read_vector[current_read_index + 1].end(), 'N', SUBST_CHAR);

               // initialize variables
               num_corrected_errors_local = 0;
               trim_5_end                 = 0;
               trim_3_end                 = 0;

               // storage for the modification of the reads
               // # of entries = read length
               sequence_modification.assign(read_length, '0');

               //----------------------------------------------------------------------
               // correct errors in a read
               //----------------------------------------------------------------------
               // forward read
               // errors cannot be corrected if k is equal to read length
               if (read_length > kmer_length) {
                  correct_errors_in_a_read_fastq(
                                                 read_vector[current_read_index + 1],
                                                 sequence_modification,
                                                 read_vector[current_read_index + 2],
                                                 trim_5_end,
                                                 trim_3_end,
                                                 read_length,
                                                 max_trimmed_bases,
                                                 min_check_length,
                                                 num_corrected_errors_local,
                                                 bit_vector,
                                                 hash_seed
                                                );
               }

               // no trim
               if (c_inst_args.notrim == true) {
                  trim_5_end = 0;
                  trim_3_end = 0;
               }
               // adjust the number of trimmed bases
               else {
                  if ((trim_5_end + trim_3_end) > max_trimmed_bases) {
                     if (trim_3_end <= max_trimmed_bases) {
                        trim_5_end = 0;
                     }
                     else if (trim_5_end <= max_trimmed_bases) {
                        trim_3_end = 0;
                     }
                     else {
                        trim_5_end = 0;
                        trim_3_end = max_trimmed_bases;
                     }
                  }
               }

               num_trimmed_bases_tmp += (trim_5_end + trim_3_end);

               // update num_corrected_reads
               too_many_errors = false;
               if (num_corrected_errors_local > (read_vector[current_read_index + 1].length() * MAX_ERROR_RATE)) {
                  too_many_errors = true;
               }
               else if (num_corrected_errors_local > 0) {
                  num_corrected_errors_tmp += num_corrected_errors_local;

                  num_corrected_reads_tmp++;
               }
               else if (c_inst_args.notrim == false) {
                  if ((trim_5_end > 0) || (trim_3_end > 0)) {
                     num_corrected_reads_tmp++;
                  }
               }

               // make a corrected read
               if (too_many_errors == false) {
                  // apply modifications to the read
                  for (unsigned int it_base = trim_5_end; it_base < (read_vector[current_read_index + 1].length() - trim_3_end); it_base++) {
                     if (sequence_modification[it_base] != '0') {
                        read_vector[current_read_index + 1][it_base] = sequence_modification[it_base];
                     }
                  }
               }

               read_vector[current_read_index] = read_vector[current_read_index] + '\n';

               // make a trimmed quality score
               if ((trim_5_end + trim_3_end) > 0) {
                  read_vector[current_read_index + 1] = read_vector[current_read_index + 1].substr(trim_5_end, read_vector[current_read_index + 1].length() - trim_5_end - trim_3_end) + '\n';
                  read_vector[current_read_index + 2] = read_vector[current_read_index + 2].substr(trim_5_end, read_vector[current_read_index + 2].length() - trim_5_end - trim_3_end) + '\n';
               }
               else {
                  read_vector[current_read_index + 1] = read_vector[current_read_index + 1] + '\n';
                  read_vector[current_read_index + 2] = read_vector[current_read_index + 2] + '\n';
               }
            }
            // too many Ns
            // write reads without modification
            else {
               // sequence
               read_vector[current_read_index + 1] += '\n';

               // quality score
               read_vector[current_read_index + 2] += '\n';
            }
         // it_read
         }

      // omp parallel
      }

      // write corrected_reads
      for (std::size_t it_read = 0; it_read < read_vector_index; it_read++) {
         current_read_index = it_read * 3;

         // header
         fwrite(read_vector[current_read_index].c_str(), 1, read_vector[current_read_index].length(), f_corrected_read2);
         // sequence
         fwrite(read_vector[current_read_index + 1].c_str(), 1, read_vector[current_read_index + 1].length(), f_corrected_read2);
         // connector
         fwrite("+\n", 1, 2, f_corrected_read2);
         // quality score
         fwrite(read_vector[current_read_index + 2].c_str(), 1, read_vector[current_read_index + 2].length(), f_corrected_read2);
      }
   }

   fclose(f_corrected_read2);

   read_vector.clear();

   num_corrected_errors = num_corrected_errors_tmp;
   num_corrected_reads  = num_corrected_reads_tmp;
   num_trimmed_bases    = num_trimmed_bases_tmp;
}
