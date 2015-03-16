#include "time.hpp"
#include "parse_args.hpp"
#include "check_inputs.hpp"
#include "count_solid_kmers.hpp"
#include "correct_errors.hpp"



//----------------------------------------------------------------------
// main
//----------------------------------------------------------------------
int main(int argc, char** argv) {
   // variables
   int size_global;
   int size_smp;
   int size_node;

   int rank_global;
   int rank_smp;
   int rank_node;

   //----------------------------------------------------------------------
   // MPI initialization
   //----------------------------------------------------------------------
   MPI_Init(&argc, &argv);

   // the real relationship between rank_global and the other ranks may be different from this
   //    0                    1                    2                : rank_node
   // -------------------  -------------------  -------------------
   // | [ ] [ ] [ ] [ ] |  | [ ] [ ] [ ] [ ] |  | [ ] [ ] [ ] [ ] |
   // |  0   1   2   3  |  |  0   1   2   3  |  |  0   1   2   3  | : rank_smp
   // -------------------  -------------------  -------------------
   //    0   1   2   3        4   5   6   7        8   9  10  11    : rank_global

   // global
   MPI_Comm_size(MPI_COMM_WORLD, &size_global);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank_global);

   // each smp
   MPI_Comm comm_smp;
   MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm_smp);
   MPI_Comm_rank(comm_smp, &rank_smp);
   MPI_Comm_size(comm_smp, &size_smp);

   // each node
   MPI_Comm comm_node;
   MPI_Comm_split(MPI_COMM_WORLD, rank_smp, 0, &comm_node);
   MPI_Comm_size(comm_node, &size_node);

   if (rank_smp == 0) {
      MPI_Comm_rank(comm_node, &rank_node);
   }

   //----------------------------------------------------------------------
   // parse arguments
   //----------------------------------------------------------------------
   C_time c_inst_time;
   C_arg  c_inst_args(argc, argv, c_inst_time, rank_global, size_node);

   //----------------------------------------------------------------------
   // check input read files
   //----------------------------------------------------------------------
   C_check_read c_inst_check_reads;

   c_inst_check_reads.size_node = size_node;

   // check input read files
   if (rank_global == 0) {
      c_inst_check_reads.check_read_file(c_inst_args, c_inst_time);
   }
   else {
      c_inst_check_reads.num_reads_vector.resize(size_node);
      c_inst_check_reads.num_reads_vector1.resize(size_node);
      c_inst_check_reads.num_reads_vector2.resize(size_node);

      c_inst_check_reads.starting_point_vector.resize(size_node);
      c_inst_check_reads.starting_point_vector1.resize(size_node);
      c_inst_check_reads.starting_point_vector2.resize(size_node);
   }

   // broadcast variables from global rank 0 to the others
   MPI_Bcast(&c_inst_check_reads.quality_score_offset,        1,         MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
   MPI_Bcast(&c_inst_check_reads.extremely_low_quality_score, 1,         MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
   MPI_Bcast(&c_inst_check_reads.quality_score_cutoff,        1,         MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
   MPI_Bcast(&c_inst_check_reads.num_reads,                   1,         MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);

   if (c_inst_args.paired_read) {
      MPI_Bcast(&c_inst_check_reads.read_file_size_byte1,      1,         MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
      MPI_Bcast(&c_inst_check_reads.read_file_size_byte2,      1,         MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
      MPI_Bcast(&c_inst_check_reads.read_file_unit_size_byte1, 1,         MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
      MPI_Bcast(&c_inst_check_reads.read_file_unit_size_byte2, 1,         MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
      MPI_Bcast(&c_inst_check_reads.num_reads_vector1[0],      size_node, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
      MPI_Bcast(&c_inst_check_reads.num_reads_vector2[0],      size_node, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
      MPI_Bcast(&c_inst_check_reads.starting_point_vector1[0], size_node, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
      MPI_Bcast(&c_inst_check_reads.starting_point_vector2[0], size_node, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
   }
   else {
      MPI_Bcast(&c_inst_check_reads.read_file_size_byte,      1,         MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
      MPI_Bcast(&c_inst_check_reads.read_file_unit_size_byte, 1,         MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
      MPI_Bcast(&c_inst_check_reads.num_reads_vector[0],      size_node, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
      MPI_Bcast(&c_inst_check_reads.starting_point_vector[0], size_node, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
   }

   //----------------------------------------------------------------------
   // count k-mers
   //----------------------------------------------------------------------
   C_count_solid_kmers c_inst_count_solid_kmers;

   c_inst_count_solid_kmers.kmer_length = c_inst_args.kmer_length;
   c_inst_count_solid_kmers.size_node   = size_node;
   c_inst_count_solid_kmers.rank_node   = rank_node;
   c_inst_count_solid_kmers.comm_node   = comm_node;

   if (rank_smp == 0) {
      if (c_inst_args.load_bf == false) {
         c_inst_count_solid_kmers.kmer_occurrence_threshold = c_inst_args.kmer_occurrence_threshold;

         // count k-mers
         c_inst_count_solid_kmers.count_kmers(c_inst_args, c_inst_time);
      }
   }

   // broadcast variables
   MPI_Bcast(&c_inst_count_solid_kmers.num_unique_solid_kmers,    1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
   MPI_Bcast(&c_inst_count_solid_kmers.kmer_occurrence_threshold, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);

   //----------------------------------------------------------------------
   // correct errors
   //----------------------------------------------------------------------
   // construct a class to correct errors
   C_correct_errors c_inst_correct_errors;

   c_inst_correct_errors.rank_node                   = rank_node;
   c_inst_correct_errors.rank_smp                    = rank_smp;
   c_inst_correct_errors.size_node                   = size_node;
   c_inst_correct_errors.comm_smp                    = comm_smp;
   c_inst_correct_errors.comm_node                   = comm_node;
   c_inst_correct_errors.max_extension               = c_inst_args.extend;
   c_inst_correct_errors.quality_score_offset        = c_inst_check_reads.quality_score_offset;
   c_inst_correct_errors.quality_score_cutoff        = c_inst_check_reads.quality_score_cutoff;
   c_inst_correct_errors.extremely_low_quality_score = c_inst_check_reads.extremely_low_quality_score;
   c_inst_correct_errors.num_reads                   = c_inst_check_reads.num_reads;

   if (c_inst_args.paired_read) {
      c_inst_correct_errors.read_file_size_byte1      = c_inst_check_reads.read_file_size_byte1;
      c_inst_correct_errors.read_file_size_byte2      = c_inst_check_reads.read_file_size_byte2;
      c_inst_correct_errors.read_file_unit_size_byte1 = c_inst_check_reads.read_file_unit_size_byte1;
      c_inst_correct_errors.read_file_unit_size_byte2 = c_inst_check_reads.read_file_unit_size_byte2;
      c_inst_correct_errors.num_reads_vector1         = c_inst_check_reads.num_reads_vector1;
      c_inst_correct_errors.num_reads_vector2         = c_inst_check_reads.num_reads_vector2;
      c_inst_correct_errors.starting_point_vector1    = c_inst_check_reads.starting_point_vector1;
      c_inst_correct_errors.starting_point_vector2    = c_inst_check_reads.starting_point_vector2;
   }
   else {
      c_inst_correct_errors.read_file_size_byte      = c_inst_check_reads.read_file_size_byte;
      c_inst_correct_errors.read_file_unit_size_byte = c_inst_check_reads.read_file_unit_size_byte;
      c_inst_correct_errors.num_reads_vector         = c_inst_check_reads.num_reads_vector;
      c_inst_correct_errors.starting_point_vector    = c_inst_check_reads.starting_point_vector;
   }

   c_inst_correct_errors.kmer_length = c_inst_args.kmer_length;

   // construct a bloom filter
   if (rank_smp == 0) {
      // newly generate bloom filter data
      if (c_inst_args.load_bf == false) {
         c_inst_correct_errors.kmer_occurrence_threshold = c_inst_count_solid_kmers.kmer_occurrence_threshold;
         c_inst_correct_errors.num_unique_solid_kmers    = (std::size_t)(c_inst_count_solid_kmers.num_unique_solid_kmers);

         // program k-mers into a Bloom filter
         c_inst_correct_errors.determine_bloom_filter_parameters(c_inst_args, c_inst_time);

         // program solid k-mers into the bloom filter
         c_inst_correct_errors.program_kmers(c_inst_args, c_inst_time);
      }
      // reuse existing bloom filter data
      else {
         // read the bloom filter parameters
         // data will be loaded in load_split_bloom_filter later
         // but this is still needed to load some parameters
         c_inst_correct_errors.read_bloom_filter_parameters(c_inst_args, c_inst_time);
      }
   }

   // broadcast variables
   MPI_Bcast(&c_inst_correct_errors.bit_vector_width,          1, MPI_UNSIGNED_LONG_LONG, 0, comm_node);
   MPI_Bcast(&c_inst_correct_errors.bit_vector_width_byte,     1, MPI_UNSIGNED_LONG_LONG, 0, comm_node);
   MPI_Bcast(&c_inst_correct_errors.num_hash_func,             1, MPI_UNSIGNED_SHORT,     0, comm_node);
   MPI_Bcast(&c_inst_correct_errors.kmer_occurrence_threshold, 1, MPI_UNSIGNED_LONG_LONG, 0, comm_node);
   MPI_Bcast(&c_inst_correct_errors.random_seed,               1, MPI_UNSIGNED_LONG_LONG, 0, comm_node);
   MPI_Bcast(&c_inst_correct_errors.num_unique_solid_kmers,    1, MPI_UNSIGNED_LONG_LONG, 0, comm_node);

   // correct errors
   if (rank_smp == 0) {
      c_inst_correct_errors.correct_errors_in_reads(c_inst_args, c_inst_time);
   }

   // summarize output results
   if ((rank_node == 0) && (rank_smp == 0)) {
      c_inst_correct_errors.summarize_outputs(c_inst_args, c_inst_time);
   }

   MPI_Finalize();

   return (EXIT_SUCCESS);
}
