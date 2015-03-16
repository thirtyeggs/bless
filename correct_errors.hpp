#ifndef _CORRECT_H
#define _CORRECT_H



#include "parse_args.hpp"
#include "check_inputs.hpp"
#include "kmc_file.h"
#include "murmurhash3/MurmurHash3.h"



//----------------------------------------------------------------------
// C_candidate_path
//----------------------------------------------------------------------
class C_candidate_path {
public:
   // modified positions and modification results
   std::vector< std::pair<unsigned int, char> > modified_bases;

   // sum of quality scores of modified bases
   std::size_t sum_qs;

   // constructors
   C_candidate_path() : sum_qs(0) {};

   void clear_path();
};



//----------------------------------------------------------------------
// C_correct_errors
//----------------------------------------------------------------------
class C_correct_errors {
public:
   // variables
   bloom_type num_unique_solid_kmers;
   bloom_type random_seed;

   std::size_t bit_vector_width;
   std::size_t bit_vector_width_byte;
   std::size_t extremely_low_quality_score;
   std::size_t kmer_occurrence_threshold;
   std::size_t max_extension;
   std::size_t num_reads;
   std::size_t num_trimmed_bases;
   std::size_t quality_score_cutoff;
   std::size_t quality_score_offset;

   long long int read_file_size_byte;
   long long int read_file_size_byte1;
   long long int read_file_size_byte2;

   unsigned int kmer_length;

   int rank_node;
   int rank_smp;
   int read_file_unit_size_byte;
   int read_file_unit_size_byte1;
   int read_file_unit_size_byte2;
   int size_node;

   unsigned short int num_hash_func;

   std::vector<std::size_t> num_reads_vector;
   std::vector<std::size_t> num_reads_vector1;
   std::vector<std::size_t> num_reads_vector2;
   std::vector<std::size_t> starting_point_vector;
   std::vector<std::size_t> starting_point_vector1;
   std::vector<std::size_t> starting_point_vector2;

   MPI_Comm comm_smp;
   MPI_Comm comm_node;

   // constructors
   C_correct_errors() :
                       num_unique_solid_kmers(0),
                       random_seed(0),
                       bit_vector_width(0),
                       bit_vector_width_byte(0),
                       extremely_low_quality_score(0),
                       kmer_occurrence_threshold(0),
                       max_extension(0),
                       num_reads(0),
                       num_trimmed_bases(0),
                       quality_score_cutoff(0),
                       quality_score_offset(0),
                       read_file_size_byte(0),
                       read_file_size_byte1(0),
                       read_file_size_byte2(0),
                       kmer_length(0),
                       rank_node(-1),
                       rank_smp(-1),
                       read_file_unit_size_byte(0),
                       read_file_unit_size_byte1(0),
                       read_file_unit_size_byte2(0),
                       size_node(-1),
                       num_hash_func(0),
                       num_corrected_errors(0),
                       num_corrected_reads(0),
                       num_wrongly_corrected_errors(0),
                       num_wrongly_corrected_errors_check(0),
                       read_block_size(0),
                       run_exploration(true)
                      {};

   // functions
   void correct_errors_in_reads(const C_arg& c_inst_args, C_time& c_inst_time);
   void determine_bloom_filter_parameters(const C_arg& c_inst_args, C_time& c_inst_time);
   void program_kmers(const C_arg& c_inst_args, C_time& c_inst_time);
   void read_bloom_filter_parameters(const C_arg& c_inst_args, C_time& c_inst_time);
   void summarize_outputs(const C_arg& c_inst_args, C_time& c_inst_time);

private:

   std::size_t num_corrected_errors;
   std::size_t num_corrected_reads;
   std::size_t num_wrongly_corrected_errors;
   std::size_t num_wrongly_corrected_errors_check;
   std::size_t read_block_size;

   std::string rank_node_text;
   std::string kmc_prefix;

   bool run_exploration;

   // functions
   void base_difference(std::vector< std::pair<unsigned int, char> >::iterator in1_begin, const std::vector< std::pair<unsigned int, char> >::iterator& in1_end, std::vector< std::pair<unsigned int, char> >::iterator in2_begin, const std::vector< std::pair<unsigned int, char> >::iterator& in2_end, std::vector< std::pair<unsigned int, char> >& out);
   void base_intersection(std::vector< std::pair<unsigned int, char> >::iterator in1_begin, const std::vector< std::pair<unsigned int, char> >::iterator& in1_end, std::vector< std::pair<unsigned int, char> >::iterator in2_begin, const std::vector< std::pair<unsigned int, char> >::iterator& in2_end, std::vector< std::pair<unsigned int, char> >& out);
   void base_union(std::vector< std::pair<unsigned int, char> >::iterator in1_begin, const std::vector< std::pair<unsigned int, char> >::iterator& in1_end, std::vector< std::pair<unsigned int, char> >::iterator in2_begin, const std::vector< std::pair<unsigned int, char> >::iterator& in2_end, std::vector< std::pair<unsigned int, char> >& out);
   void check_first_kmer(const std::string& kmer, const C_candidate_path& candidate_path_in, const std::vector<unsigned int>& low_qs_indexes, std::vector<C_candidate_path>& candidate_path_vector, const std::size_t& index, const unsigned char* bit_vector, const std::vector<unsigned int>& hash_seed);
   void correct_errors_in_reads_single_fastq(const C_arg& c_inst_args);
   void correct_errors_in_reads_paired_fastq(const C_arg& c_inst_args);
   void correct_errors_in_a_read_fastq(const std::string& sequence, std::string& sequence_modification, const std::string& quality_score, std::size_t& trim_5_end, std::size_t& trim_3_end, const std::size_t& read_length, const std::size_t& max_trimmed_bases, const std::size_t& min_check_length, std::size_t& num_corrected_errors_local, const unsigned char* bit_vector, const std::vector<unsigned int>& hash_seed);
   void correct_errors_between_solid_regions(const std::string& org_sequence, std::string& sequence, const std::string& quality_score, const std::size_t& left_first, const std::size_t& index_start, const std::size_t& index_end, const std::size_t& right_second, const std::size_t& org_boundary_left, const std::size_t& org_boundary_right, std::string& sequence_modification, std::size_t& trim_5_end, std::size_t& trim_3_end, const std::size_t& num_solid_islands, const std::size_t& read_length, const std::size_t& max_trimmed_bases, const std::size_t& min_check_length, std::size_t& num_corrected_errors_local, const unsigned char* bit_vector, const std::vector<unsigned int>& hash_seed);
   void correct_errors_5_prime_end(const std::string& org_sequence, std::string& sequence, const std::string& quality_score, const std::size_t& index_start, std::string& sequence_modification, std::size_t& trim_5_end, std::size_t& trim_3_end, const std::size_t& org_boundary, const std::size_t& read_length, const std::size_t& max_trimmed_bases, const std::size_t& min_check_length, std::size_t& num_corrected_errors_local, const unsigned char* bit_vector, const std::vector<unsigned int>& hash_seed);
   void correct_errors_3_prime_end(const std::string& org_sequence, std::string& sequence, const std::string& quality_score, const std::size_t& index_start, std::string& sequence_modification, std::size_t& trim_5_end, std::size_t& trim_3_end, const std::size_t& org_boundary, const std::size_t& read_length, const std::size_t& max_trimmed_bases, const std::size_t& min_check_length, std::size_t& num_corrected_errors_local, const unsigned char* bit_vector, const std::vector<unsigned int>& hash_seed);
   void correct_errors_first_kmer(const std::string& sequence, const std::string& quality_score, std::string& sequence_modification, std::vector<C_candidate_path>& candidate_path_vector, const unsigned char* bit_vector, const std::vector<unsigned int>& hash_seed);
   void extend_a_kmer(const std::string& kmer, const std::string& sequence, const std::size_t& index_kmer, const std::size_t& index_last_mod, C_candidate_path& current_path, std::vector<C_candidate_path>& candidate_path_vector, const std::size_t& org_boundary_left, const std::size_t& org_boundary_right, const std::string& quality_score, const unsigned char* bit_vector, const std::vector<unsigned int>& hash_seed);
   void extend_a_kmer_5_prime_end(const std::string& kmer, const std::string& sequence, const std::size_t& index_kmer, C_candidate_path& current_path, std::vector<C_candidate_path>& candidate_path_vector, const std::size_t& org_boundary, const std::string& quality_score, const unsigned char* bit_vector, const std::vector<unsigned int>& hash_seed);
   void extend_a_kmer_3_prime_end(const std::string& kmer, const std::string& sequence, const std::size_t& index_kmer, C_candidate_path& current_path, std::vector<C_candidate_path>& candidate_path_vector, const std::size_t& org_boundary, const std::string& quality_score, const unsigned char* bit_vector, const std::vector<unsigned int>& hash_seed);
   void extend_first_kmer_to_right(const std::string& sequence, const std::string& quality_score, C_candidate_path& candidate_path_in, std::vector<C_candidate_path>& candidate_path_vector_all, const std::size_t& read_length, const unsigned char* bit_vector, const std::vector<unsigned int>& hash_seed);
   void extend_out_left(const std::string& kmer, const std::size_t& num_extend, const std::size_t& extend_amount, bool& extension_success, const unsigned char* bit_vector, const std::vector<unsigned int>& hash_seed);
   void extend_out_right(const std::string& kmer, const std::size_t& num_extend, const std::size_t& extend_amount, bool& extension_success, const unsigned char* bit_vector, const std::vector<unsigned int>& hash_seed);
   void generate_hash_seed(const bloom_type random_seed, std::vector<unsigned int>& hash_seed);
   void solid_first_kmer(const C_candidate_path& candidate_path, const std::string& sequence, bool& extension_success, const unsigned char* bit_vector, const std::vector<unsigned int>& hash_seed);

   bool query_text(const std::string& kmer, const unsigned char*& bit_vector, const std::vector<unsigned int>& hash_seed);

   std::string remove_new_line(std::string in_string);
};



#endif
