#include "parse_args.hpp"



//----------------------------------------------------------------------
// read_args
//----------------------------------------------------------------------
void C_arg::read_args() {
   if (num_args == 1) {
      if (rank_global == 0) {
         print_usage();
      }

      exit(EXIT_SUCCESS);
   }

   //--------------------------------------------------
   // parse arguments
   //--------------------------------------------------
   std::string bf_prefix("");

   for (int it_arg = 1; it_arg < num_args; it_arg++) {
      if (strcmp(args[it_arg], "-help") == 0) {
         print_usage();
         exit(EXIT_SUCCESS);
      }
      else if (strcmp(args[it_arg], "-read") == 0) {
         if (it_arg <= num_args - 2) {
            read_file_name = args[it_arg + 1];
            it_arg++;
         }
         else {
            std::cout << std::endl << "ERROR: The read file name is not specified after the -read option" << std::endl << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 400);
         }
      }
      else if (strcmp(args[it_arg], "-read1") == 0) {
         if (it_arg <= num_args - 2) {
            read_file_name1 = args[it_arg + 1];
            it_arg++;
         }
         else {
            std::cout << std::endl << "ERROR: The first read file name is not specified after the -read1 option" << std::endl << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 401);
         }
      }
      else if (strcmp(args[it_arg], "-read2") == 0) {
         if (it_arg <= num_args - 2) {
            read_file_name2 = args[it_arg + 1];
            it_arg++;
         }
         else {
            std::cout << std::endl << "ERROR: The second read file name is not specified after the -read2 option" << std::endl << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 402);
         }
      }
      else if (strcmp(args[it_arg], "-prefix") == 0) {
         if (it_arg <= num_args - 2) {
            prefix = args[it_arg + 1];
            it_arg++;
         }
         else {
            std::cout << std::endl << "ERROR: The prefix is not specified after the -prefix option" << std::endl << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 403);
         }
      }
      else if (strcmp(args[it_arg], "-load") == 0) {
         if (it_arg <= num_args - 2) {
            bf_prefix = args[it_arg + 1];
            load_bf = true;
            load_bf_text = "On";
            it_arg++;
         }
         else {
            std::cout << std::endl << "ERROR: The Bloom filter file name prefix is not specified after the -load option" << std::endl << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 405);
         }
      }
      else if (strcmp(args[it_arg], "-kmerlength") == 0) {
         if (it_arg <= num_args - 2) {
            kmer_length = atoi(args[it_arg + 1]);
            it_arg++;
         }
         else {
            std::cout << std::endl << "ERROR: The k-mer length is not specified after the -kmerlength option" << std::endl << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 406);
         }
      }
      else if (strcmp(args[it_arg], "-fpr") == 0) {
         if (it_arg <= num_args - 2) {
            target_false_positive_prob = atof(args[it_arg + 1]);
            it_arg++;
         }
         else {
            std::cout << std::endl << "ERROR: The false positvie probability is not specified after the -fpr option" << std::endl << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 407);
         }
      }
      else if (strcmp(args[it_arg], "-seed") == 0) {
         if (it_arg <= num_args - 2) {
            random_seed = atol(args[it_arg + 1]);
            it_arg++;
         }
         else {
            std::cout << std::endl << "ERROR: The seed is not specified after the -seed option" << std::endl << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 408);
         }
      }
      else if (strcmp(args[it_arg], "-count") == 0) {
         if (it_arg <= num_args - 2) {
            kmer_occurrence_threshold = atoi(args[it_arg + 1]);
            set_kmer_occurrence_threshold = true;
            it_arg++;
         }
         else {
            std::cout << std::endl << "ERROR: The k-mer occurrence threshold is not specified after the -count option" << std::endl << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 409);
         }
      }
      else if (strcmp(args[it_arg], "-extend") == 0) {
         if (it_arg <= num_args - 2) {
            extend = atoi(args[it_arg + 1]);
            it_arg++;
         }
         else {
            std::cout << std::endl << "ERROR: The max extension is not specified after the -extend option" << std::endl << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 410);
         }
      }
      else if (strcmp(args[it_arg], "-max_mem") == 0) {
         if (it_arg <= num_args - 2) {
            max_mem = atoi(args[it_arg + 1]);
            it_arg++;
         }
         else {
            std::cout << std::endl << "ERROR: The max memory usage is not specified after the -max_mem option" << std::endl << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 411);
         }
      }
      else if (strcmp(args[it_arg], "-smpthread") == 0) {
         if (it_arg <= num_args - 2) {
            smpthread = atoi(args[it_arg + 1]);
            it_arg++;
         }
         else {
            std::cout << std::endl << "ERROR: The number of threads in a SMP node is not specified after the -smpthread option" << std::endl << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 412);
         }
      }
      else if (strcmp(args[it_arg], "-debug") == 0) {
         debug = true;
      }
      else if (strcmp(args[it_arg], "-notrim") == 0) {
         notrim = true;
         notrim_text = "On";
      }
      else {
         std::cout << std::endl << "ERROR: Illegal option " << args[it_arg] << std::endl << std::endl;
         MPI_Abort(MPI_COMM_WORLD, 413);
      }
   }


   //--------------------------------------------------
   // check options
   //--------------------------------------------------
   std::smatch smatch1;

   std::regex rx_gz_ext("\\.gz$");

   // paried-end input
   if (read_file_name.empty()) {
      paired_read = true;

      if (read_file_name1.empty()) {
         std::cout << std::endl << "ERROR: The first read file name is not specified" << std::endl << std::endl;
         MPI_Abort(MPI_COMM_WORLD, 414);
      }
      else {
         std::ifstream f_tmp;
         f_tmp.open(read_file_name1.c_str());
         if (f_tmp.is_open() == false) {
            std::cout << std::endl << "ERROR: Cannot open " << read_file_name1 << std::endl << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 415);
         }
         f_tmp.close();
      }

      if (read_file_name2.empty()) {
         std::cout << std::endl << "ERROR: The second read file name is not specified" << std::endl << std::endl;
         MPI_Abort(MPI_COMM_WORLD, 416);
      }
      else {
         std::ifstream f_tmp;
         f_tmp.open(read_file_name2.c_str());
         if (f_tmp.is_open() == false) {
            std::cout << std::endl << "ERROR: Cannot open " << read_file_name2 << std::endl << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 417);
         }
         f_tmp.close();
      }

      // check whether both the files are gzipped
      // the first file is gzipped
      if (std::regex_search(read_file_name1, smatch1, rx_gz_ext)) {
         // the second file is gzipped
         if (std::regex_search(read_file_name2, smatch1, rx_gz_ext)) {
            gzipped_input_read = true;
         }
         // the second file is not gzipped
         else {
            std::cout <<  std::endl << "ERROR: Only one out of two input read files is gzipped" << read_file_name2 << std::endl << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 417);
         }
      }
      // the first file is not gzipped
      else {
         // the second file is gzipped
         if (std::regex_search(read_file_name1, smatch1, rx_gz_ext)) {
            std::cout <<  std::endl << "ERROR: Only one out of two input read files is gzipped" << read_file_name2 << std::endl << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 417);
         }
      }
   }
   // single-end input
   else {
      paired_read = false;
      paired_read_text = "Off";

      std::ifstream f_tmp;
      f_tmp.open(read_file_name.c_str());
      if (f_tmp.is_open() == false) {
         std::cout << std::endl << "ERROR: Cannot open " << read_file_name << std::endl << std::endl;
         MPI_Abort(MPI_COMM_WORLD, 418);
      }
      f_tmp.close();

      if (read_file_name1.empty() == false) {
         std::cout << std::endl << "ERROR: -read1 cannot be used with -read" << std::endl;
         MPI_Abort(MPI_COMM_WORLD, 419);
      }
      if (read_file_name2.empty() == false) {
         std::cout << std::endl << "ERROR: -read2 cannot be used with -read" << std::endl;
         MPI_Abort(MPI_COMM_WORLD, 420);
      }

      // check whether both the files are gzipped
      if (std::regex_search(read_file_name, smatch1, rx_gz_ext)) {
         gzipped_input_read = true;
      }
   }

   // prefix
   if (prefix.empty()) {
      std::cout << std::endl << "ERROR: The prefix is not specified" << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 421);
   }

   boost::filesystem::path prefix_tmp(prefix);
   boost::filesystem::path prefix_tmp_absolute(boost::filesystem::complete(prefix_tmp));
   boost::filesystem::path prefix_dir(prefix_tmp_absolute.parent_path());

   // kmc
   boost::filesystem::path bless(args[0]);
   boost::filesystem::path bless_absolute(boost::filesystem::complete(bless));
   boost::filesystem::path bless_dir(bless_absolute.parent_path());

   kmc_binary = bless_dir.string() + "/" + KMC_DIR + "/" + KMC_BINARY;

   if (kmer_length == 0) {
      std::cout << std::endl << "ERROR: k-mer length not specified" << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 423);
   }
   else if (kmer_length < MIN_KMER_LENGTH) {
      std::cout << std::endl << "ERROR: k-mer length should be >= " << MIN_KMER_LENGTH << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 424);
   }
   else if (kmer_length > MAX_KMER_LENGTH) {
      std::cout << std::endl << "ERROR: k-mer length should be < " << MAX_KMER_LENGTH << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 425);
   }

   if (target_false_positive_prob == 0.0) {
      std::cout << std::endl << "ERROR: Target false positive probability not specified" << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 426);
   }
   else if (target_false_positive_prob < 0 || target_false_positive_prob > 1) {
      std::cout << std::endl << "ERROR: False positive rate should be between 0 and 1" << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 427);
   }

   if (set_kmer_occurrence_threshold == true) {
      if (kmer_occurrence_threshold < 2) {
         std::cout << std::endl << "ERROR: The k-mer occurrence threshold should be >= 2" << std::endl << std::endl;
         MPI_Abort(MPI_COMM_WORLD, 428);
      }
      else if (kmer_occurrence_threshold > MAX_KMER_THRESHOLD) {
         std::cout << std::endl << "ERROR: The k-mer occurrence threshold cannot exceed " << MAX_KMER_THRESHOLD << std::endl << std::endl;
         MPI_Abort(MPI_COMM_WORLD, 429);
      }

      if (load_bf == true) {
         std::cout << std::endl << "ERROR: The k-mer occurrence threshold cannot be used with the -count option" << std::endl << std::endl;
         MPI_Abort(MPI_COMM_WORLD, 430);
      }
   }

   // extend amount
   if (extend < 1) {
      std::cout << std::endl << "ERROR: The max extension should be >= 1" << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 431);
   }

   // max memory usage for kmc
   if (max_mem < MIN_MAX_MEM) {
      std::cout << std::endl << "ERROR: The max memeory usage should be >= " << MIN_MAX_MEM << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 432);
   }
   else if (max_mem > MAX_MAX_MEM) {
      std::cout << std::endl << "ERROR: The max memeory usage should be <= " << MAX_MAX_MEM << std::endl << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 433);
   }

   // number of cores in a smp node
   // the number is not given
   if (smpthread == 0) {
      smpthread = omp_get_max_threads();
   }
   else if (smpthread > (unsigned int)omp_get_max_threads()) {
      if (rank_global == 0) {
         std::cout << std::endl << "WARNING: The option for the number of threads is larger than the number of cores in a SMP (" << omp_get_max_threads() << ")" << std::endl << std::endl;
      }
   }

   omp_set_num_threads(smpthread);

   //--------------------------------------------------
   // output file names
   //--------------------------------------------------
   // set trim_file_name
   trim_file_name  = prefix + ".trim";
   trim_file_name1 = prefix + ".1.trim";
   trim_file_name2 = prefix + ".2.trim";

   // set a hisgram file name
   kmer_histo_file_name = prefix + ".histo.k-mer";
   qs_histo_file_name   = prefix + ".histo.qs";

   // set error_corrected read file names
   corrected_read_file_name  = prefix + ".corrected.fastq";
   corrected_read_file_name1 = prefix + ".1.corrected.fastq";
   corrected_read_file_name2 = prefix + ".2.corrected.fastq";

   // set bloom filter dump file names
   bf_data_file_name        = prefix    + ".bf.data";
   bf_size_file_name        = prefix    + ".bf.size";
   loaded_bf_data_file_name = bf_prefix + ".bf.data";
   loaded_bf_size_file_name = bf_prefix + ".bf.size";

   // set an input list file
   kmc_prefix      = prefix + ".kmc";
   input_list_file = kmc_prefix + ".input-list";

   //--------------------------------------------------
   // print options
   //--------------------------------------------------
   if (rank_global == 0) {
      // stdout
      std::cout << std::endl;
      std::cout << "----------------------------------------------------------------------" << std::endl;
      std::cout << "           BLESS: Bloom-filter based error correction tool" << std::endl;
      std::cout << "----------------------------------------------------------------------" << std::endl;
      std::cout << "CONTACT : Yun Heo (yunheo1@illinois.edu)" << std::endl;
      std::cout << "VERSION: " VERSION << std::endl;
      std::cout << "DATE   : " DATE << std::endl;
      std::cout << "----------------------------------------------------------------------" << std::endl;
      std::cout << std::endl;
      std::cout << "Parsing arguments" << std::endl;

      std::cout << "     Number of nodes              : " << size_node << std::endl;
      std::cout << "     Number of threads per node   : " << smpthread << std::endl;
      std::cout << "     Max memory usage for KMC     : " << max_mem << " GB" << std::endl;

      if (paired_read == true) {
         std::cout << "     1st read file name           : " << read_file_name1 << std::endl;
         std::cout << "     2nd read file name           : " << read_file_name2 << std::endl;
      }
      else {
         std::cout << "     Read file name               : " << read_file_name << std::endl;
      }

      std::cout << "     k-mer length                 : " << kmer_length << std::endl;
      std::cout << "     Target false positive prob.  : " << target_false_positive_prob << std::endl;
      std::cout << "     Random seed                  : " << random_seed << std::endl;

      if (set_kmer_occurrence_threshold == true) {
         std::cout << "     k-mer occurence threshold    : " << kmer_occurrence_threshold << std::endl;
      }
      else {
         std::cout << "     k-mer occurence threshold    : Not specified" << std::endl;
      }

      std::cout << "     Read extension amount        : " << extend << std::endl;
      std::cout << "     Load given Bloom filter data : " << load_bf_text << std::endl;
      std::cout << "     No trim reads                : " << notrim_text << std::endl;
      std::cout << "     Parsing arguments            : done" << std::endl;
      std::cout << std::endl;
   }
}

//----------------------------------------------------------------------
// print_usage
//----------------------------------------------------------------------
void C_arg::print_usage() {
   std::cout << std::endl;
   std::cout << "----------------------------------------------------------------------" << std::endl;
   std::cout << "           BLESS: Bloom-filter based error correction tool" << std::endl;
   std::cout << "----------------------------------------------------------------------" << std::endl;
   std::cout << "VERSION: " VERSION << std::endl;
   std::cout << "DATE   : " DATE << std::endl;
   std::cout << "----------------------------------------------------------------------" << std::endl;
   std::cout << std::endl;
   std::cout << "1. USAGE" << std::endl;
   std::cout << "     " << args[0] << " <OPTIONS>" << std::endl;
   std::cout << std::endl;
   std::cout << "2. OPTIONS" << std::endl;
   std::cout << "     1) REQUIRED" << std::endl;
   std::cout << "     -read <file>: Unpaired input fastq file. One of \"-read\" or" << std::endl;
   std::cout << "          \"-read1/-read2\" should be used." << std::endl;
   std::cout << "     -read1 <file>: Forward fastq file of paired-end reads. It should" << std::endl;
   std::cout << "          be used with \"-read2\"." << std::endl;
   std::cout << "     -read2 <file>: Reverse fastq file of paired-end reads. It should" << std::endl;
   std::cout << "          be used with \"-read1\"." << std::endl;
   std::cout << "     -kmerlength <number>: Length of k-mers." << std::endl;
   std::cout << "     -prefix <string>: Output file prefix." << std::endl;
   std::cout << std::endl;
   std::cout << "     2) OPTIONAL" << std::endl;
   std::cout << "     -count <integer>: Minimum occurrences for solid k-mers." << std::endl;
   std::cout << "          This number is automatically determined when it is not" << std::endl;
   std::cout << "          given." << std::endl;
   std::cout << "     -extend <integer>: Read extension amount. Default: " << extend << "." << std::endl;
   std::cout << "     -fpr <float>: Target false positive probability for the Bloom" << std::endl;
   std::cout << "          filter. Default: " << DEFAULT_FPR << "." << std::endl;
   std::cout << "     -smpthread <integer>: Number of threads used in a SMP node." << std::endl;
   std::cout << "          Default: number of cores in each SMP node." << std::endl;
   std::cout << "     -load <prefix>: Skip the solid k-mer finding step and load" << std::endl;
   std::cout << "          pre-built Bloom filter data. When BLESS is executed with" << std::endl;
   std::cout << "          \"-prefix <prefix>\" option <prefix>.bf.data and" << std::endl;
   std::cout << "          <prefix>.bf.size are generated. If \"-load <prefix>\" option" << std::endl;
   std::cout << "          is used, BLESS does not try to fiile solid k-mers in inputs" << std::endl;
   std::cout << "          and load the existing Bloom filter data files." << std::endl;
   std::cout << "     -max_mem <4-1024>: Set maximum memory usage for KMC. Default: " << MIN_MAX_MEM << "." << std::endl;
   std::cout << "     -notrim: Turn trimming off." << std::endl;
   std::cout << "     -seed <integer>: Set a seed for random number generation." << std::endl;
   std::cout << "          Default: " << DEFAULT_SEED << "." << std::endl;
   std::cout << std::endl;
   std::cout << "3. EXAMPLES" << std::endl;
   std::cout << "     1) PAIRED INPUT READS" << std::endl;
   std::cout << "     bless -read1 in1.fastq -read2 in2.fastq -prefix \\" << std::endl;
   std::cout << "          directory/prefix -kmerlength 31" << std::endl;
   std::cout << std::endl;
   std::cout << "     2) UNPAIRED INPUT READ" << std::endl;
   std::cout << "     bless -read in.fastq -prefix directory/prefix -kmerlength 31" << std::endl;
   std::cout << std::endl;
   std::cout << "4. CONTACT" << std::endl;
   std::cout << "     Yun Heo (yunheo1@illinois.edu)" << std::endl;
   std::cout << std::endl;
}
