
KMER_DIR := kmer
FUNC_DIR := function
UTIL_DIR := utils
HASH_DIR := hash
IO_DIR := io
KMC_API_DIR := $(KMER_DIR)/kmc_api
ALGLIB_DIR := $(UTIL_DIR)/alglib-3.15.0/src
CUR_DIR := $(notdir $(shell pwd))

INC := -I ./
LIBS := -L ./
CXX := g++
CXXFLAGS := -std=c++11 -O3 $(INC)
CXXFLAGS2 := -std=c++11 -O3
CXXLINK := -O3 -std=c++11 -lm -lz -lpthread -lboost_system -lboost_thread

KMER_OBJS := $(KMER_DIR)/Kmer.o $(KMER_DIR)/KmerIntPair.o $(KMER_DIR)/kmc_reader.o
UTIL_OBJS := $(UTIL_DIR)/fqstq.o
HASH_OBJS := $(HASH_DIR)/hash.o $(HASH_DIR)/lshash.o
IO_OBJS := $(IO_DIR)/ioHT.o $(IO_DIR)/ioMatrix.o 
FUNC_OBJS := $(FUNC_DIR)/cluster.o $(FUNC_DIR)/distance.o $(FUNC_DIR)/funcAB.o

KMC_API_OBJS := \
$(KMC_API_DIR)/mmer.o \
$(KMC_API_DIR)/kmc_file.o \
$(KMC_API_DIR)/kmer_api.o

ALGLIB_OBJS := \
$(ALGLIB_DIR)/ap.o \
$(ALGLIB_DIR)/alglibinternal.o \
$(ALGLIB_DIR)/alglibmisc.o \
$(ALGLIB_DIR)/linalg.o \
$(ALGLIB_DIR)/specialfunctions.o \
$(ALGLIB_DIR)/statistics.o

$(KMER_DIR)/%.o: $(KMER_DIR)/%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@
$(FUNC_DIR)/%.o: $(FUNC_DIR)/%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@
$(UTIL_DIR)/%.o: $(UTIL_DIR)/%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@
$(HASH_DIR)/%.o: $(HASH_DIR)/%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@
$(IO_DIR)/%.o: $(IO_DIR)/%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@
$(KMC_API_DIR)/%.o: $(KMC_API_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ 
$(ALGLIB_DIR)/%.o: $(ALGLIB_DIR)/%.cpp
	$(CXX) $(CXXFLAGS2) -c $< -o $@

kmerLSH: app/kmerLSH.cc $(KMER_OBJS) $(HASH_OBJS) $(IO_OBJS) $(UTIL_OBJS) $(FUNC_OBJS) $(KMC_API_OBJS) $(ALGLIB_OBJS)
	$(CXX) -o $@ $^ $(INC) $(LIBS) $(CXXLINK)

clean: 
	-rm -f $(KMER_OBJS)/*.o
	-rm -f $(HASH_OBJS)/*.o
	-rm -f $(IO_OBJS)/*.o
	-rm -f $(UTIL_OBJS)/*.o
	-rm -f $(FUNC_OBJS)/*.o
	-rm -f $(KMC_API_DIR)/*.o
	-rm -f $(ALGLIB_DIR)/*.o
	-rm -f kmerLSH

all: kmerLSH




