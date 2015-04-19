#include "check_inputs.hpp"



KSEQ_INIT(gzFile, gzread)



//----------------------------------------------------------------------
// check_read_file
//----------------------------------------------------------------------
void C_check_read::check_read_file(const C_arg& c_inst_args, C_time& c_inst_time) {
   // start measuring run-time
   time_t rawtime;
   time(&rawtime);
   c_inst_time.start_check_read_file = asctime(localtime(&rawtime));

   std::cout << "Checking input read files" << std::endl;

   // paired
   if (c_inst_args.paired_read == true) {
      check_read_file_fastq_paired(c_inst_args);
   }
   // single
   else {
      check_read_file_fastq_single(c_inst_args);
   }

   time(&rawtime);
   c_inst_time.end_check_read_file = asctime(localtime(&rawtime));
}



//----------------------------------------------------------------------
// check_read_file_fastq_single
//----------------------------------------------------------------------
void C_check_read::check_read_file_fastq_single(const C_arg& c_inst_args) {
   //----------------------------------------------------------------------
   // construct a quality score histogram
   //----------------------------------------------------------------------
   qs_histo_vec.resize(QS_HISTOGRAM_MAX + 1, 0);

   if (c_inst_args.gzipped_input_read) {
      construct_quality_score_histogram_single_gzipped(c_inst_args);
   }
   else {
      construct_quality_score_histogram_single_unzipped(c_inst_args);
   }

   //----------------------------------------------------------------
   // find the max/min quality scores
   //----------------------------------------------------------------
   max_quality_score = -1000;
   min_quality_score =  1000;

   // max
   for (int it_histo = 0; it_histo <= QS_HISTOGRAM_MAX; it_histo++) {
      if (qs_histo_vec[it_histo] > 0) {
         max_quality_score = it_histo;
      }
   }

   // min
   for (int it_histo = QS_HISTOGRAM_MAX; it_histo >= 0; it_histo--) {
      if (qs_histo_vec[it_histo] > 0) {
         min_quality_score = it_histo;
      }
   }

   //----------------------------------------------------------------
   // determine the quality score offset
   //----------------------------------------------------------------
   // http://en.wikipedia.org/wiki/FASTQ_format
   if (min_quality_score < 33) {
      std::cout << std::endl << "ERROR: There is a quality score < 33" << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 100);
   }
   // 33 <= min_quality_score <= 58
   else if (min_quality_score <= 58) {
      if (max_quality_score <= 74) {
         quality_score_offset = PHRED33;
      }
      else {
         std::cout << std::endl << "ERROR: Irregular quality score range " << min_quality_score << "-" << max_quality_score << std::endl << std::endl;
         MPI_Abort(MPI_COMM_WORLD, 101);
      }
   }
   // 59 <= min_quality_score <= 74
   else if (min_quality_score <= 74) {
      if (max_quality_score <= 74) {
         std::cout << std::endl << "WARNING: Hard to determine the quality score offset (min: " << min_quality_score << ", max: " << max_quality_score << ")" << std::endl;
         std::cout <<              "         Phred+33 will be used" << std::endl << std::endl;
         quality_score_offset = PHRED33;
      }
      // max_quality_score >= 75
      else {
         quality_score_offset = PHRED64;
      }
   }
   // min_quality_score >= 75
   else {
      quality_score_offset = PHRED64;
   }

   //----------------------------------------------------------------------
   // write a histogram file
   //----------------------------------------------------------------------
   std::size_t total_bases(0);

   std::ofstream f_histo;

   f_histo.open(c_inst_args.qs_histo_file_name.c_str());

   if (f_histo.is_open()) {
      for (std::size_t it_histo = quality_score_offset; it_histo <= quality_score_offset + MAX_PHRED; it_histo++) {
         f_histo << std::setw(4) << it_histo - quality_score_offset << ": " << std::setw(14) << qs_histo_vec[it_histo] << std::endl;

         total_bases += qs_histo_vec[it_histo];
      }
   }

   f_histo.close();

   // find the threshold
   std::size_t partial_sum(0);

   bool set_quality_score_cutoff(false);
   bool set_extremely_low_quality_score(false);

   for (std::size_t it_histo = 0; it_histo <= QS_HISTOGRAM_MAX; it_histo++) {
      partial_sum += qs_histo_vec[it_histo];

      if (((1.0 * partial_sum / total_bases) >= QS_CUTOFF_RATIO) && (set_quality_score_cutoff == false)) {
         quality_score_cutoff = it_histo - quality_score_offset;
         set_quality_score_cutoff = true;
      }

      if (((1.0 * partial_sum / total_bases) >= QS_EXTREMELY_LOW_RATIO) && (set_extremely_low_quality_score == false)) {
         extremely_low_quality_score = it_histo - quality_score_offset;
         set_extremely_low_quality_score = true;

         if ((1.0 * partial_sum / total_bases) >= (QS_EXTREMELY_LOW_RATIO + 0.05)) {
            std::cout << std::endl;
            std::cout << "WARNING: Some quality score thresholds are set to a high value" << std::endl;
            std::cout << "         because overall quality scores are too bad." << std::endl;
            std::cout << "         It may cause long runtime and large memory usage." << std::endl << std::endl;
         }
      }
   }

   std::cout << "     Number of reads           : " << num_reads << std::endl;
   std::cout << "     Quality score offset      : " << quality_score_offset << std::endl;
   std::cout << "     Quality score threshold   : " << quality_score_cutoff << std::endl;
   std::cout << "     Low quality score threhold: " << extremely_low_quality_score << std::endl;
   std::cout << "     Checking input read files : done" << std::endl;
   std::cout << std::endl;
}



//----------------------------------------------------------------------
// check_read_file_fastq_paired
//----------------------------------------------------------------------
void C_check_read::check_read_file_fastq_paired(const C_arg& c_inst_args) {
   //----------------------------------------------------------------------
   // construct a quality score histogram
   //----------------------------------------------------------------------
   qs_histo_vec.resize(QS_HISTOGRAM_MAX + 1, 0);

   if (c_inst_args.gzipped_input_read) {
      construct_quality_score_histogram_paired_gzipped(c_inst_args);
   }
   else {
      construct_quality_score_histogram_paired_unzipped(c_inst_args);
   }

   //----------------------------------------------------------------
   // find the max/min quality scores
   //----------------------------------------------------------------
   max_quality_score = -1000;
   min_quality_score =  1000;

   // max
   for (int it_histo = 0; it_histo <= QS_HISTOGRAM_MAX; it_histo++) {
      if (qs_histo_vec[it_histo] > 0) {
         max_quality_score = it_histo;
      }
   }

   // min
   for (int it_histo = QS_HISTOGRAM_MAX; it_histo >= 0; it_histo--) {
      if (qs_histo_vec[it_histo] > 0) {
         min_quality_score = it_histo;
      }
   }

   //----------------------------------------------------------------
   // determine the quality score offset
   //----------------------------------------------------------------
   // http://en.wikipedia.org/wiki/FASTQ_format
   if (min_quality_score < 33) {
      std::cout << std::endl << "ERROR: There is a quality score < 33" << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 102);
   }
   // 33 <= min_quality_score <= 58
   else if (min_quality_score <= 58) {
      if (max_quality_score <= 74) {
         quality_score_offset = PHRED33;
      }
      else {
         std::cout << std::endl << "ERROR: Irregular quality score range " << min_quality_score << "-" << max_quality_score << std::endl << std::endl;
         MPI_Abort(MPI_COMM_WORLD, 103);
      }
   }
   // 59 <= min_quality_score <= 74
   else if (min_quality_score <= 74) {
      if (max_quality_score <= 74) {
         std::cout << std::endl << "WARNING: Hard to determine the quality score offset (min: " << min_quality_score << ", max: " << max_quality_score << ")" << std::endl;
         std::cout <<              "         Phred+33 will be used" << std::endl << std::endl;
         quality_score_offset = PHRED33;
      }
      // max_quality_score >= 75
      else {
         quality_score_offset = PHRED64;
      }
   }
   // min_quality_score >= 75
   else {
      quality_score_offset = PHRED64;
   }

   //----------------------------------------------------------------------
   // write a histogram file
   //----------------------------------------------------------------------
   std::size_t total_bases(0);

   std::ofstream f_histo;

   f_histo.open(c_inst_args.qs_histo_file_name.c_str());

   if (f_histo.is_open()) {
      for (std::size_t it_histo = quality_score_offset; it_histo <= quality_score_offset + MAX_PHRED; it_histo++) {
         f_histo << std::setw(4) << it_histo - quality_score_offset << ": " << std::setw(14) << qs_histo_vec[it_histo] << std::endl;

         total_bases += qs_histo_vec[it_histo];
      }
   }

   f_histo.close();

   // find the threshold
   std::size_t partial_sum(0);

   bool set_quality_score_cutoff(false);
   bool set_extremely_low_quality_score(false);

   for (std::size_t it_histo = 0; it_histo <= QS_HISTOGRAM_MAX; it_histo++) {
      partial_sum += qs_histo_vec[it_histo];

      if (((1.0 * partial_sum / total_bases) >= QS_CUTOFF_RATIO) && (set_quality_score_cutoff == false)) {
         quality_score_cutoff = it_histo - quality_score_offset;
         set_quality_score_cutoff = true;
      }

      if (((1.0 * partial_sum / total_bases) >= QS_EXTREMELY_LOW_RATIO) && (set_extremely_low_quality_score == false)) {
         extremely_low_quality_score = it_histo - quality_score_offset;
         set_extremely_low_quality_score = true;

         if ((1.0 * partial_sum / total_bases) >= (QS_EXTREMELY_LOW_RATIO + 0.05)) {
            std::cout << std::endl;
            std::cout << "WARNING: Some quality score thresholds are set to a high value" << std::endl;
            std::cout << "         because overall quality scores are too bad." << std::endl;
            std::cout << "         It may cause long runtime and large memory usage." << std::endl << std::endl;
         }
      }
   }

   std::cout << "     Number of pairs           : " << num_reads << std::endl;
   std::cout << "     Quality score offset      : " << quality_score_offset << std::endl;
   std::cout << "     Quality score threshold   : " << quality_score_cutoff << std::endl;
   std::cout << "     Low quality score threhold: " << extremely_low_quality_score << std::endl;
   std::cout << "     Checking input read files : done" << std::endl << std::endl;
}



//----------------------------------------------------------------------
// construct_quality_score_histogram_single_unzipped
//----------------------------------------------------------------------
void C_check_read::construct_quality_score_histogram_single_unzipped(const C_arg& c_inst_args) {

   // open an input read file
   boost::iostreams::mapped_file_source f_read;

   // find out the size of read files
   f_read.open(c_inst_args.read_file_name.c_str());

   if (!f_read.is_open()) {
      std::cout << std::endl << "ERROR: Cannot open " << f_read << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 104);
   }

   read_file_size_byte      = f_read.size();
   read_file_unit_size_byte = f_read.alignment();

   std::size_t num_units(ceil(1.0 * read_file_size_byte / read_file_unit_size_byte));

   // check whether the size of the input file is 0
   if (read_file_size_byte == 0) {
      std::cout << std::endl << "ERROR: The size of " << c_inst_args.read_file_name << " is 0" << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 105);
   }

   // check whether the alignment unit offset is larger than MMAP_FILE_SIZE
   if (MMAP_FILE_SIZE <= read_file_unit_size_byte) {
      std::cout << std::endl << "ERROR: The block size for " << f_read << " should be larger than " << f_read.alignment() << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 106);
   }

   // check whether the size of the input file < read_file_unit_size_byte X <number of nodes>
   // if it is true, some cores may not have one alignment unit in the input file
   if (num_units < (std::size_t)size_node) {
      std::cout << std::endl << "ERROR: The size of " << c_inst_args.read_file_name << " is too small. Run BLESS again without MPI." << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 107);
   }

   f_read.close();

   //----------------------------------------------------------------
   // define variables
   //----------------------------------------------------------------
   bool flag_keep_going;

   std::size_t line_length;

   short int line_order;

   long long int processed_bytes;
   long long int alignment_offset;
   long long int remaining_bytes;
   long long int prev_remaining_bytes;
   long long int this_iteration_bytes;
   long long int num_units_per_node;
   long long int num_bytes_per_node;
   long long int prev_num_reads;

   int current_node_index;

   const char* pt_current;
   const char* pt_last;
   const char* pt_current_read_start;
   const char* pt_current_line_start;

   // initialize variables
   max_quality_score  = -1000;
   min_quality_score  = 1000;
   num_units_per_node = ceil(1.0 *  num_units / size_node);
   num_bytes_per_node = num_units_per_node * read_file_unit_size_byte;

   //----------------------------------------------------------------
   // construct the quality score histogram
   //----------------------------------------------------------------
   // initialize variables (for each iteration)
   flag_keep_going      = true;
   processed_bytes      = 0;
   alignment_offset     = 0;
   remaining_bytes      = read_file_size_byte;
   prev_remaining_bytes = read_file_size_byte;
   current_node_index   = 1;
   prev_num_reads       = 0;

   starting_point_vector.push_back(0);

   // read loop
   while (flag_keep_going) {
      // calculate the number of bytes that will be processed in this iteration
      if (remaining_bytes <= MMAP_FILE_SIZE) {
         this_iteration_bytes = remaining_bytes;
      }
      else {
         this_iteration_bytes = MMAP_FILE_SIZE;
      }

      // open the file
      //open(<file name>, <length>, <offset>)
      f_read.open(c_inst_args.read_file_name.c_str(), this_iteration_bytes, processed_bytes);

      if (!f_read.is_open()) {
         std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.read_file_name << std::endl << std::endl;
         MPI_Abort(MPI_COMM_WORLD, 108);
      }

      // initialize variables
      pt_current            = f_read.data() + alignment_offset;
      pt_last               = f_read.data() + f_read.size();
      pt_current_read_start = pt_current;
      pt_current_line_start = pt_current;

      processed_bytes += alignment_offset;

      line_order = 0;

      while (pt_current && (pt_current != pt_last)) {
         if ((pt_current = static_cast<const char*>(memchr(pt_current, '\n', pt_last - pt_current)))) {
            // quality score
            if (line_order == 3) {
               // find the max and min values
               line_length = pt_current - pt_current_line_start;

               // check read length
               if (line_length < c_inst_args.kmer_length) {
                  std::cout << std::endl << "ERROR: There are reads that are shorter than given k" << std::endl << std::endl;
                  MPI_Abort(MPI_COMM_WORLD, 109);
               }

               for (unsigned int it = 0; it < line_length; it++) {
                  qs_histo_vec[(int)pt_current_line_start[it]]++;
               }

               line_order            = 0;
               processed_bytes      += (pt_current - pt_current_read_start + 1);
               pt_current_read_start = pt_current + 1;
               num_reads++;

               if (processed_bytes >= num_bytes_per_node * current_node_index) {
                  num_reads_vector.push_back(num_reads - prev_num_reads);
                  prev_num_reads = num_reads;

                  starting_point_vector.push_back(processed_bytes);

                  current_node_index++;
               }
            }
            // header or sequence or connector
            else {
               line_order++;
            }

            pt_current++;

            pt_current_line_start = pt_current;
         }
      }

      remaining_bytes  = read_file_size_byte - processed_bytes;

      if (remaining_bytes <= 0) {
         flag_keep_going = false;
      }
      else if (remaining_bytes == prev_remaining_bytes) {
         std::cout << std::endl << "ERROR: The number of lines in " << c_inst_args.read_file_name << " is wrong" << std::endl << std::endl;
         MPI_Abort(MPI_COMM_WORLD, 110);
      }
      else {
         // align processedbytes to the multiple of unit_alignment_offset
         alignment_offset     = processed_bytes % read_file_unit_size_byte;
         processed_bytes     -= alignment_offset;
         prev_remaining_bytes = remaining_bytes;
         remaining_bytes     += alignment_offset;
      }

      f_read.close();
   }

   // if the number of bytes of each node is exactly same
   // the last element would have already been added to the vector
   if (num_reads_vector.size() < (unsigned int)size_node) {
      num_reads_vector.push_back(num_reads - prev_num_reads);
   }
}



//----------------------------------------------------------------------
// construct_quality_score_histogram_paired_unzipped
//----------------------------------------------------------------------
void C_check_read::construct_quality_score_histogram_paired_unzipped(const C_arg& c_inst_args) {
   // open input read files
   boost::iostreams::mapped_file_source f_read_1;
   boost::iostreams::mapped_file_source f_read_2;

   // find out the size of read files
   f_read_1.open(c_inst_args.read_file_name1.c_str());
   f_read_2.open(c_inst_args.read_file_name2.c_str());

   if (!f_read_1.is_open()) {
      std::cout << std::endl << "ERROR: Cannot open " << f_read_1 << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 111);
   }
   if (!f_read_2.is_open()) {
      std::cout << std::endl << "ERROR: Cannot open " << f_read_2 << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 112);
   }

   read_file_size_byte1      = f_read_1.size();
   read_file_size_byte2      = f_read_2.size();
   read_file_unit_size_byte1 = f_read_1.alignment();
   read_file_unit_size_byte2 = f_read_2.alignment();

   std::size_t num_units1(ceil(1.0 * read_file_size_byte1 / read_file_unit_size_byte1));
   std::size_t num_units2(ceil(1.0 * read_file_size_byte2 / read_file_unit_size_byte2));

   // check whether the size of the input file is 0
   if (read_file_size_byte1 == 0) {
      std::cout << std::endl << "ERROR: The size of " << c_inst_args.read_file_name1 << " is 0" << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 113);
   }
   if (read_file_size_byte2 == 0) {
      std::cout << std::endl << "ERROR: The size of " << c_inst_args.read_file_name2 << " is 0" << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 114);
   }

   // check whether the alignment unit offset is larger than MMAP_FILE_SIZE
   if (MMAP_FILE_SIZE <= read_file_unit_size_byte1) {
      std::cout << std::endl << "ERROR: The block size for " << f_read_1 << " should be larger than " << f_read_1.alignment() << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 115);
   }
   if (MMAP_FILE_SIZE <= read_file_unit_size_byte2) {
      std::cout << std::endl << "ERROR: The block size for " << f_read_2 << " should be larger than " << f_read_2.alignment() << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 116);
   }

   // check whether the size of the input file < read_file_unit_size_byte X <number of nodes>
   // if it is true, some cores may not have one alignment unit in the input file
   if (num_units1 < (std::size_t)size_node) {
      std::cout << std::endl << "ERROR: The size of " << c_inst_args.read_file_name1 << " is too small. Run BLESS again without MPI." << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 117);
   }
   if (num_units2 < (std::size_t)size_node) {
      std::cout << std::endl << "ERROR: The size of " << c_inst_args.read_file_name2 << " is too small. Run BLESS again without MPI." << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 118);
   }

   f_read_1.close();
   f_read_2.close();

   //----------------------------------------------------------------
   // define variables
   //----------------------------------------------------------------
   std::size_t num_reads_1(0);
   std::size_t num_reads_2(0);

   bool flag_keep_going;

   std::size_t line_length;

   short int line_order;

   long long int processed_bytes;
   long long int alignment_offset;
   long long int remaining_bytes;
   long long int prev_remaining_bytes;
   long long int this_iteration_bytes;
   long long int num_units_per_node1;
   long long int num_units_per_node2;
   long long int num_bytes_per_node1;
   long long int num_bytes_per_node2;
   long long int prev_num_reads;

   int current_node_index;

   const char* pt_current;
   const char* pt_last;
   const char* pt_current_read_start;
   const char* pt_current_line_start;

   // initialize variables
   max_quality_score   = -1000;
   min_quality_score   = 1000;
   num_units_per_node1 = ceil(1.0 * num_units1 / size_node);
   num_units_per_node2 = ceil(1.0 * num_units2 / size_node);
   num_bytes_per_node1 = num_units_per_node1 * read_file_unit_size_byte1;
   num_bytes_per_node2 = num_units_per_node2 * read_file_unit_size_byte2;

   //----------------------------------------------------------------
   // find out the max/min values of quality scores
   //----------------------------------------------------------------
   //
   // forward read file
   //
   // initialize variables
   flag_keep_going      = true;
   processed_bytes      = 0;
   alignment_offset     = 0;
   remaining_bytes      = read_file_size_byte1;
   prev_remaining_bytes = read_file_size_byte1;
   current_node_index   = 1;
   prev_num_reads       = 0;

   starting_point_vector1.push_back(0);

   // read loop
   while (flag_keep_going) {
      // calculate the number of bytes that will be processed in this iteration
      if (remaining_bytes <= MMAP_FILE_SIZE) {
         this_iteration_bytes = remaining_bytes;
      }
      else {
         this_iteration_bytes = MMAP_FILE_SIZE;
      }

      // open the file
      //open(<file name>, <length>, <offset>)
      f_read_1.open(c_inst_args.read_file_name1.c_str(), this_iteration_bytes, processed_bytes);

      if (!f_read_1.is_open()) {
         std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.read_file_name1 << std::endl << std::endl;
         MPI_Abort(MPI_COMM_WORLD, 119);
      }

      // initialize variables
      pt_current            = f_read_1.data() + alignment_offset;
      pt_last               = f_read_1.data() + f_read_1.size();
      pt_current_read_start = pt_current;
      pt_current_line_start = pt_current;

      processed_bytes += alignment_offset;

      line_order = 0;

      while (pt_current && (pt_current != pt_last)) {
         if ((pt_current = static_cast<const char*>(memchr(pt_current, '\n', pt_last - pt_current)))) {
            // quality score
            if (line_order == 3) {
               // find the max and min values
               line_length = pt_current - pt_current_line_start;

               // check read length
               if (line_length < c_inst_args.kmer_length) {
                  std::cout << std::endl << "ERROR: There are reads that are shorter than given k" << std::endl << std::endl;
                  MPI_Abort(MPI_COMM_WORLD, 120);
               }

               for (unsigned int it = 0; it < line_length; it++) {
                  qs_histo_vec[(int)pt_current_line_start[it]]++;
               }

               line_order            = 0;
               processed_bytes      += (pt_current - pt_current_read_start + 1);
               pt_current_read_start = pt_current + 1;
               num_reads_1++;

               if (processed_bytes >= num_bytes_per_node1 * current_node_index) {
                  num_reads_vector1.push_back(num_reads_1 - prev_num_reads);
                  prev_num_reads = num_reads_1;

                  starting_point_vector1.push_back(processed_bytes);

                  current_node_index++;
               }
            }
            // header or sequence or connector
            else {
               line_order++;
            }

            pt_current++;

            pt_current_line_start = pt_current;
         }
      }

      remaining_bytes  = read_file_size_byte1 - processed_bytes;

      if (remaining_bytes <= 0) {
         flag_keep_going = false;
      }
      else if (remaining_bytes == prev_remaining_bytes) {
         std::cout << std::endl << "ERROR: The number of lines in " << c_inst_args.read_file_name1 << " is wrong" << std::endl << std::endl;
         MPI_Abort(MPI_COMM_WORLD, 121);
      }
      else {
         // align processedbytes to the multiple of unit_alignment_offset
         alignment_offset     = processed_bytes % read_file_unit_size_byte1;
         processed_bytes     -= alignment_offset;
         prev_remaining_bytes = remaining_bytes;
         remaining_bytes     += alignment_offset;
      }

      f_read_1.close();
   }

   // if the number of bytes of each node is exactly same
   // the last element would have already been added to the vector
   if (num_reads_vector1.size() < (unsigned int)size_node) {
      num_reads_vector1.push_back(num_reads_1 - prev_num_reads);
   }

   //
   // reverse read file
   //
   // initialize variables
   flag_keep_going           = true;
   processed_bytes      = 0;
   alignment_offset     = 0;
   remaining_bytes      = read_file_size_byte2;
   prev_remaining_bytes = read_file_size_byte2;
   current_node_index   = 1;
   prev_num_reads       = 0;

   starting_point_vector2.push_back(0);

   // read loop
   while (flag_keep_going) {
      // calculate the number of bytes that will be processed in this iteration
      if (remaining_bytes <= MMAP_FILE_SIZE) {
         this_iteration_bytes = remaining_bytes;
      }
      else {
         this_iteration_bytes = MMAP_FILE_SIZE;
      }

      // open the file
      //open(<file name>, <length>, <offset>)
      f_read_2.open(c_inst_args.read_file_name2.c_str(), this_iteration_bytes, processed_bytes);

      if (!f_read_2.is_open()) {
         std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.read_file_name2 << std::endl << std::endl;
         MPI_Abort(MPI_COMM_WORLD, 122);
      }

      // initialize variables
      pt_current            = f_read_2.data() + alignment_offset;
      pt_last               = f_read_2.data() + f_read_2.size();
      pt_current_read_start = pt_current;
      pt_current_line_start = pt_current;

      processed_bytes += alignment_offset;

      line_order = 0;

      while (pt_current && (pt_current != pt_last)) {
         if ((pt_current = static_cast<const char*>(memchr(pt_current, '\n', pt_last - pt_current)))) {
            // quality score
            if (line_order == 3) {
               // find the max and min values
               line_length = pt_current - pt_current_line_start;

               // check read length
               if (line_length < c_inst_args.kmer_length) {
                  std::cout << std::endl << "ERROR: There are reads that are shorter than given k" << std::endl << std::endl;
                  MPI_Abort(MPI_COMM_WORLD, 123);
               }

               for (unsigned int it = 0; it < line_length; it++) {
                  qs_histo_vec[(int)pt_current_line_start[it]]++;
               }

               line_order            = 0;
               processed_bytes      += (pt_current - pt_current_read_start + 1);
               pt_current_read_start = pt_current + 1;
               num_reads_2++;

               if (processed_bytes >= num_bytes_per_node2 * current_node_index) {
                  num_reads_vector2.push_back(num_reads_2 - prev_num_reads);
                  prev_num_reads = num_reads_2;

                  starting_point_vector2.push_back(processed_bytes);

                  current_node_index++;
               }
            }
            // header or sequence or connector
            else {
               line_order++;
            }

            pt_current++;

            pt_current_line_start = pt_current;
         }
      }

      remaining_bytes  = read_file_size_byte2 - processed_bytes;

      if (remaining_bytes <= 0) {
         flag_keep_going = false;
      }
      else if (remaining_bytes == prev_remaining_bytes) {
         std::cout << std::endl << "ERROR: The number of lines in " << c_inst_args.read_file_name2 << " is wrong" << std::endl << std::endl;
         MPI_Abort(MPI_COMM_WORLD, 124);
      }
      else {
         // align processedbytes to the multiple of unit_alignment_offset
         alignment_offset     = processed_bytes % read_file_unit_size_byte2;
         processed_bytes     -= alignment_offset;
         prev_remaining_bytes = remaining_bytes;
         remaining_bytes     += alignment_offset;
      }

      f_read_2.close();
   }

   // if the number of bytes of each node is exactly same
   // the last element would have already been added to the vector
   if (num_reads_vector2.size() < (unsigned int)size_node) {
      num_reads_vector2.push_back(num_reads_2 - prev_num_reads);
   }

   // compair the number of reads in two input files
   if (num_reads_1 == num_reads_2) {
      num_reads = num_reads_1;
   }
   else {
      std::cout << std::endl << "ERROR: Number of lines in two input files are not same" << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 125);
   }
}



//----------------------------------------------------------------------
// construct_quality_score_histogram_single_gzipped
//----------------------------------------------------------------------
void C_check_read::construct_quality_score_histogram_single_gzipped(const C_arg& c_inst_args) {
   gzFile f_read;

   kseq_t* each_read;

   // open an input read
   f_read = gzopen(c_inst_args.read_file_name.c_str(), "r");

   if (f_read == Z_NULL) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.read_file_name << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // initialize each_read
   each_read = kseq_init(f_read);

   // build the quality score histogram
   int length;

   num_reads = 0;

   while ((length = kseq_read(each_read)) >= 0) {
      for (unsigned int it_qual = 0; it_qual < each_read->qual.l; it_qual++) {
         qs_histo_vec[(int)each_read->qual.s[it_qual]]++;
      }

      num_reads++;
   }

   // close the file
   kseq_destroy(each_read);
   gzclose(f_read);

   // set num_reads_vector
   std::size_t num_reads_per_node(ceil(1.0 * num_reads / size_node));
   std::size_t num_reads_tmp(num_reads);

   for (int it = 0; it < size_node; it++) {
      if (num_reads_tmp >= num_reads_per_node) {
         num_reads_vector.push_back(num_reads_per_node);
         num_reads_tmp = num_reads_tmp - num_reads_per_node;
      }
      else {
         num_reads_vector.push_back(num_reads_tmp);
         num_reads_tmp = 0;
      }
   }
}



//----------------------------------------------------------------------
// construct_quality_score_histogram_paired_gzipped
//----------------------------------------------------------------------
void C_check_read::construct_quality_score_histogram_paired_gzipped(const C_arg& c_inst_args) {
   int length;

   //
   // forward
   //
   kseq_t* each_read_1;

   // open an input read
   gzFile f_read_1;
   f_read_1 = gzopen(c_inst_args.read_file_name1.c_str(), "r");

   if (f_read_1 == Z_NULL) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.read_file_name1 << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // initialize each_read
   each_read_1 = kseq_init(f_read_1);

   // build the quality score histogram
   std::size_t num_reads_1(0);

   while ((length = kseq_read(each_read_1)) >= 0) {
      for (unsigned int it_qual = 0; it_qual < each_read_1->qual.l; it_qual++) {
         qs_histo_vec[(int)each_read_1->qual.s[it_qual]]++;
      }

      num_reads_1++;
   }

   // close the file
   kseq_destroy(each_read_1);
   gzclose(f_read_1);

   //
   // reverse
   //
   kseq_t* each_read_2;

   // open an input read
   gzFile f_read_2;
   f_read_2 = gzopen(c_inst_args.read_file_name2.c_str(), "r");

   if (f_read_2 == Z_NULL) {
      std::cout << std::endl << "ERROR: Cannot open " << c_inst_args.read_file_name2 << std::endl << std::endl;
      exit(EXIT_FAILURE);
   }

   // initialize each_read
   each_read_2 = kseq_init(f_read_2);

   // build the quality score histogram
   std::size_t num_reads_2(0);

   while ((length = kseq_read(each_read_2)) >= 0) {
      for (unsigned int it_qual = 0; it_qual < each_read_2->qual.l; it_qual++) {
         qs_histo_vec[(int)each_read_2->qual.s[it_qual]]++;
      }

      num_reads_2++;
   }

   // close the file
   kseq_destroy(each_read_2);
   gzclose(f_read_2);

   //
   // set num_reads_vector_*
   //
   std::size_t num_reads_per_node_1(ceil(1.0 * num_reads_1 / size_node));
   std::size_t num_reads_per_node_2(ceil(1.0 * num_reads_2 / size_node));
   std::size_t num_reads_tmp_1(num_reads_1);
   std::size_t num_reads_tmp_2(num_reads_2);

   for (int it = 0; it < size_node; it++) {
      // forward
      if (num_reads_tmp_1 >= num_reads_per_node_1) {
         num_reads_vector1.push_back(num_reads_per_node_1);
         num_reads_tmp_1 = num_reads_tmp_1 - num_reads_per_node_1;
      }
      else {
         num_reads_vector1.push_back(num_reads_tmp_1);
         num_reads_tmp_1 = 0;
      }

      // reverse
      if (num_reads_tmp_2 >= num_reads_per_node_2) {
         num_reads_vector2.push_back(num_reads_per_node_2);
         num_reads_tmp_2 = num_reads_tmp_2 - num_reads_per_node_2;
      }
      else {
         num_reads_vector2.push_back(num_reads_tmp_2);
         num_reads_tmp_2 = 0;
      }
   }

   // compare num_reads_1 and num_reads_2
   if (num_reads_1 == num_reads_2) {
      num_reads = num_reads_1;
   }
   else {
      std::cout << std::endl << "ERROR: The number of lines in the two input read files is not matched" << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 134);
   }
}
