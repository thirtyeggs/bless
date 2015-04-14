#ifndef _COUNT_SOLID_KMERS_H
#define _COUNT_SOLID_KMERS_H



#include "parse_args.hpp"
#include "check_inputs.hpp"
#include "kmc_file.h"



//----------------------------------------------------------------------
// C_count_solid_kmers
//----------------------------------------------------------------------
class C_count_solid_kmers {
public:
   // variables
   bloom_type num_unique_solid_kmers;

   bloom_type* num_occurrences_histogram;
   bloom_type* reduced_num_occurrences_histogram;

   std::size_t kmer_occurrence_threshold;
   std::size_t valey_point;

   unsigned int kmer_length;

   int rank_node;
   int size_node;

   bool get_valey_point;

   MPI_Comm comm_node;

   // constructors
   C_count_solid_kmers() :
                          num_unique_solid_kmers(0),
                          kmer_occurrence_threshold(0),
                          valey_point(0),
                          kmer_length(0),
                          rank_node(-1),
                          size_node(0),
                          get_valey_point(false)
                         {};

   // functions
   void count_kmers(const C_arg& c_inst_args, C_time& c_inst_time);

private:
   CKMCFile kmer_database;

   std::string rank_node_text;
   std::string kmc_prefix;

   //functions
   void run_kmc(const C_arg& c_inst_args);
   void determine_kmer_occurrence_threshold(const C_arg& c_inst_args);
   void generate_kmer_occurrence_histogram(const C_arg& c_inst_args);
   void count_unique_solid_kmers(const C_arg& c_inst_args);
   void write_kmer_histogram(const C_arg& c_inst_args);
};



#endif
