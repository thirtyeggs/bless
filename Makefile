CC=mpicxx
PROG=bless
CFLAGS=-Wall -O3 -I ./kmc/kmc_api -fopenmp -std=c++11
LDFLAGS=-lboost_filesystem -lboost_system -lboost_iostreams -lboost_regex -fopenmp -std=c++11

SRCS=kmc/kmc_api/kmc_file.cpp kmc/kmc_api/kmer_api.cpp kmc/kmc_api/mmer.cpp murmurhash3/MurmurHash3.cpp check_inputs.cpp correct_errors.cpp count_solid_kmers.cpp main.cpp parse_args.cpp
OBJS=$(SRCS:.cpp=.o)

$(PROG): $(OBJS)
	$(CC) $(OBJS) -o $@ $(LDFLAGS)
	cd kmc; make CC=$(CC)

kmc/kmc_api/%.o: kmc/kmc_api/%.cpp
	$(CC) $(CFLAGS) -c $(DEF) $? -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) -c $(DEF) $? -o $@

clean:
	rm -rf $(PROG) $(OBJS)
	cd kmc; make clean
