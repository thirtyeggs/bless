#ifndef _PARSE_ARGS_H
#define _PARSE_ARGS_H



#include "define.hpp"
#include "time.hpp"



//----------------------------------------------------------------------
// C_arg
//----------------------------------------------------------------------
class C_arg {
public:
   // variables
   std::string bf_data_file_name;
   std::string bf_size_file_name;
   std::string corrected_read_file_name;
   std::string corrected_read_file_name1;
   std::string corrected_read_file_name2;
   std::string input_list_file;
   std::string kmc_binary;
   std::string kmc_prefix;
   std::string kmer_histo_file_name;
   std::string loaded_bf_data_file_name;
   std::string loaded_bf_size_file_name;
   std::string load_bf_text;
   std::string notrim_text;
   std::string paired_read_text;
   std::string prefix;
   std::string qs_histo_file_name;
   std::string read_file_name;
   std::string read_file_name1;
   std::string read_file_name2;
   std::string trim_file_name;
   std::string trim_file_name1;
   std::string trim_file_name2;

   bloom_type random_seed;

   double target_false_positive_prob;

   int rank_global;
   int size_node;

   std::size_t extend;
   std::size_t kmer_occurrence_threshold;
   std::size_t max_mem;
   std::size_t smpthread;

   unsigned int kmer_length;

   bool debug;
   bool load_bf;
   bool notrim;
   bool paired_read;
   bool set_kmer_occurrence_threshold;
   bool gzipped_input_read;

   // constructors
   explicit C_arg(
                  int argc,
                  char** argv,
                  C_time& c_inst_time,
                  int in_rank_global,
                  int in_size_node
                 ) :
                     load_bf_text("Off"),
                     notrim_text("Off"),
                     paired_read_text("On"),
                     random_seed(DEFAULT_SEED),
                     target_false_positive_prob(DEFAULT_FPR),
                     rank_global(in_rank_global),
                     size_node(in_size_node),
                     extend(MAX_EXTENSION),
                     kmer_occurrence_threshold(0),
                     max_mem(MIN_MAX_MEM),
                     smpthread(0),
                     kmer_length(0),
                     debug(false),
                     load_bf(false),
                     notrim(false),
                     paired_read(true),
                     set_kmer_occurrence_threshold(false),
                     gzipped_input_read(false),
                     num_args(argc),
                     args(argv)
                    {
      time_t rawtime;
      time(&rawtime);
      c_inst_time.start_parse_args = asctime(localtime(&rawtime));

      read_args();

      time(&rawtime);
      c_inst_time.end_parse_args = asctime(localtime(&rawtime));
   };

private:
   // variables
   int num_args;

   char** args;

   // functions
   void print_usage();
   void read_args();
};



#endif
