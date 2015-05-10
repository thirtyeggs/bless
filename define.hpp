#ifndef _DEFINE_H
#define _DEFINE_H



// c++ libraries
#include <algorithm>
#include <bitset>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
#include <queue>
#include <regex>
#include <sstream>
#include <stack>
#include <string>
#include <vector>

// c posix libraries
#include <sys/stat.h>
#include <sys/types.h>

// boost libraries
#include <boost/array.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/functional/hash/extensions.hpp>
//#include <boost/regex.hpp>

// google sparse hash
#include <google/sparse_hash_map>

// openmp
#include <omp.h>

// mpi
#include <mpi.h>

// klib
#include <zlib.h>
#include <kseq.h>



// definitions
#define HEX_UC(x)                 std::setw(2) << std::setfill('0') << std::hex << (int)(x)
#define VERSION                   "1.01"
#define DATE                      "May. 10, 2015"
#define BITS_PER_CHAR             8
#define MAX_KMER_THRESHOLD        65535
#define NUM_NEOCLEOTIDE           4
#define A                         0
#define C                         1
#define G                         2
#define T                         3
#define MIN_KMER_LENGTH           10
#define MAX_KMER_LENGTH           255 // max value of KMC
#define ZERO                      0
#define BIT1                      1
#define BIT8                      128
#define BPS_PER_BYTE              4
#define KMER_HISTOGRAM_SIZE       256 // KMC default
#define DEFAULT_FPR               0.001
#define PHRED33                   33
#define PHRED64                   64
#define MIN_NON_SOLID_LENGTH      2
#define MIN_SOLID_LENGTH          2
#define INIT_MIN_QS               1000000
#define MIN_QS_DIFF               10
#define MAX_EXTENSION             5
#define MAX_LOW_QS_BASES          3
#define FP_SUSPECT_LENGTH         1
#define SOLID_REGION_ADJUST_RANGE 4
#define SUBST_CHAR                'A'
#define MAX_ERROR_RATE            0.25
#define MAX_TRIMMING_RATE         0.6
#define MIN_BASES_AFTER_TRIMMING  30
#define DEFAULT_SEED              0
#define CHECK_RANGE_RATIO         0.07
#define NUM_ALLOWABLE_FAILS       2
#define MAX_AVG_QUALITY_SCORE     29
#define MAX_MODIFICATION          4
#define QS_HISTOGRAM_MAX          126
#define QS_CUTOFF_RATIO           0.05
#define QS_EXTREMELY_LOW_RATIO    0.01
#define MAX_N_RATIO               0.1
#define MAX_CANDIDATE_PATHS       74
#define KMC_DIR                   "kmc/bin"
#define KMC_BINARY                "kmc"
#define MIN_MAX_MEM               4
#define MAX_MAX_MEM               1024
#define KMC_DEFAULT_MIN_COUNT     2
#define PIGZ_DIR                  "pigz/pigz-2.3.3"
#define PIGZ_BINARY               "pigz"
#define DEFAULT_QS_CUTOFF         -10
#define BIT_VEC_INC_RATIO         1.2
#define OPENMP_CHUNK_SIZE         10
#define READ_BLOCK_SIZE_RATIO     1000
#define MMAP_FILE_SIZE            209715200 // 200 * 1024 * 1024 = 200 MB
#define BLOOM_RCV_BUFFER_SIZE     104857600 // 100 * 1024 * 1024 = 100 MB
#define MAX_PHRED                 41
#define KMER_BLOCK_SIZE           100000

static const char NEOCLEOTIDE[NUM_NEOCLEOTIDE] = {'A', 'C', 'G', 'T'};

static const unsigned char BIT_MASK[BITS_PER_CHAR] = {
                                                      0x01, //00000001
                                                      0x02, //00000010
                                                      0x04, //00000100
                                                      0x08, //00001000
                                                      0x10, //00010000
                                                      0x20, //00100000
                                                      0x40, //01000000
                                                      0x80  //10000000
                                                     };



typedef unsigned long long int bloom_type;
const unsigned int PREDEF_NUM_HASH_FUNC = 128;
static const bloom_type PREDEF_HASH_FUNC[PREDEF_NUM_HASH_FUNC] = {
                                                                    0xAAAAAAAAAAAAAAAAULL, 0x5555555555555555ULL, 0x3333333333333333ULL, 0xCCCCCCCCCCCCCCCCULL,
                                                                    0x6666666666666666ULL, 0x9999999999999999ULL, 0xB5B5B5B5B5B5B5B5ULL, 0x4B4B4B4B4B4B4B4BULL,
                                                                    0xAA55AA55AA55AA55ULL, 0x5533553355335533ULL, 0x33CC33CC33CC33CCULL, 0xCC66CC66CC66CC66ULL,
                                                                    0x6699669966996699ULL, 0x99B599B599B599B5ULL, 0xB54BB54BB54BB54BULL, 0x4BAA4BAA4BAA4BAAULL,
                                                                    0xAA33AA33AA33AA33ULL, 0x55CC55CC55CC55CCULL, 0x3366336633663366ULL, 0xCC99CC99CC99CC99ULL,
                                                                    0x66B566B566B566B5ULL, 0x994B994B994B994BULL, 0xB5AAB5AAB5AAB5AAULL, 0xAAAAAA33AAAAAA33ULL,
                                                                    0x555555CC555555CCULL, 0x3333336633333366ULL, 0xCCCCCC99CCCCCC99ULL, 0x666666B5666666B5ULL,
                                                                    0x9999994B9999994BULL, 0xB5B5B5AAB5B5B5AAULL, 0xFFFFFFFFFFFFFFFFULL, 0xFFFF0000FFFF0000ULL,
                                                                    0xB823D5EBB823D5EBULL, 0xC1191CDFC1191CDFULL, 0xF623AEB3F623AEB3ULL, 0xDB58499FDB58499FULL,
                                                                    0xC8D42E70C8D42E70ULL, 0xB173F616B173F616ULL, 0xA91A5967A91A5967ULL, 0xDA427D63DA427D63ULL,
                                                                    0xB1E8A2EAB1E8A2EAULL, 0xF6C0D155F6C0D155ULL, 0x4909FEA34909FEA3ULL, 0xA68CC6A7A68CC6A7ULL,
                                                                    0xC395E782C395E782ULL, 0xA26057EBA26057EBULL, 0x0CD5DA280CD5DA28ULL, 0x467C5492467C5492ULL,
                                                                    0xF15E6982F15E6982ULL, 0x61C6FAD361C6FAD3ULL, 0x9615E3529615E352ULL, 0x6E9E355A6E9E355AULL,
                                                                    0x689B563E689B563EULL, 0x0C9831A80C9831A8ULL, 0x6753C18B6753C18BULL, 0xA622689BA622689BULL,
                                                                    0x8CA63C478CA63C47ULL, 0x42CC288442CC2884ULL, 0x8E89919B8E89919BULL, 0x6EDBD7D36EDBD7D3ULL,
                                                                    0x15B6796C15B6796CULL, 0x1D6FDFE41D6FDFE4ULL, 0x63FF909263FF9092ULL, 0xE7401432E7401432ULL,
                                                                    0xEFFE9412EFFE9412ULL, 0xAEAEDF79AEAEDF79ULL, 0x9F245A319F245A31ULL, 0x83C136FC83C136FCULL,
                                                                    0xC3DA4A8CC3DA4A8CULL, 0xA5112C8CA5112C8CULL, 0x5271F4915271F491ULL, 0x9A948DAB9A948DABULL,
                                                                    0xCEE59A8DCEE59A8DULL, 0xB5F525ABB5F525ABULL, 0x59D1321759D13217ULL, 0x24E7C33124E7C331ULL,
                                                                    0x697C2103697C2103ULL, 0x84B0A46084B0A460ULL, 0x86156DA986156DA9ULL, 0xAEF2AC68AEF2AC68ULL,
                                                                    0x23243DA523243DA5ULL, 0x3F6496433F649643ULL, 0x5FA495A85FA495A8ULL, 0x67710DF867710DF8ULL,
                                                                    0x9A6C499E9A6C499EULL, 0xDCFB0227DCFB0227ULL, 0x46A4343346A43433ULL, 0x1832B07A1832B07AULL,
                                                                    0xC46AFF3CC46AFF3CULL, 0xB9C8FFF0B9C8FFF0ULL, 0xC9500467C9500467ULL, 0x34431BDF34431BDFULL,
                                                                    0xB652432BB652432BULL, 0xE367F12BE367F12BULL, 0x427F4C1B427F4C1BULL, 0x224C006E224C006EULL,
                                                                    0x2E7E5A892E7E5A89ULL, 0x96F99AA596F99AA5ULL, 0x0BEB452A0BEB452AULL, 0x2FD87C392FD87C39ULL,
                                                                    0x74B2E1FB74B2E1FBULL, 0x222EFD24222EFD24ULL, 0xF357F60CF357F60CULL, 0x440FCB1E440FCB1EULL,
                                                                    0x8BBE030F8BBE030FULL, 0x6704DC296704DC29ULL, 0x1144D12F1144D12FULL, 0x948B1355948B1355ULL,
                                                                    0x6D8FD7E96D8FD7E9ULL, 0x1C11A0141C11A014ULL, 0xADD1592FADD1592FULL, 0xFB3C712EFB3C712EULL,
                                                                    0xFC77642FFC77642FULL, 0xF9C4CE8CF9C4CE8CULL, 0x31312FB931312FB9ULL, 0x08B0DD7908B0DD79ULL,
                                                                    0x318FA6E7318FA6E7ULL, 0xC040D23DC040D23DULL, 0xC0589AA7C0589AA7ULL, 0x0CA5C0750CA5C075ULL,
                                                                    0xF874B172F874B172ULL, 0x0CF914D50CF914D5ULL, 0x784D3280784D3280ULL, 0x4E8CFEBC4E8CFEBCULL,
                                                                    0xC569F575C569F575ULL, 0xCDB2A091CDB2A091ULL, 0x2CC016B42CC016B4ULL, 0x5C5F44215C5F4421ULL
                                                                 };



#endif
