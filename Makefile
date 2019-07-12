# Backward-cpp (https://github.com/bombela/backward-cpp) is useful for stack tracing. Optional.
DBACKWARD = #-DBACKWARD_HAS_BFD

CXX      = g++
CC       = gcc
CPPFLAGS = -g -O2
CFLAGS   = -g -O2
CXXFLAGS = -c -std=c++11
LIBS     = -pthread -lstdc++ -lcurl -lz -lbz2
DFLAGS   = $(DBACKWARD)
INCLUDE  = -I third_party/htslib/htslib-1.9 -I third_party/eigen/eigen-3.3.7 -I third_party/boost/boost_1_69_0 -I src -I src/stats -I src/kallisto -I third_party/gzip-hpp/gzip-hpp-master/include -I third_party/zlib/zlib-1.2.11

EXEC         = anaquin
SOURCES      = $(wildcard src/*.cpp src/RNA/*.cpp src/Meta/*.cpp src/kallisto/*.cpp src/tools/*.cpp src/analyzers/*.cpp src/Genomics/*.cpp src/data/*.cpp src/parsers/*.cpp src/writers/*.cpp src/stats/*.cpp)
OBJECTS      = $(SOURCES:.cpp=.o)
SOURCES_LIB  = $(wildcard src/tools/*.c third_party/htslib/htslib-1.9/*.c third_party/htslib/htslib-1.9/cram/*.c third_party/zlib/zlib-1.2.11/*.c)
OBJECTS_LIB  = $(SOURCES_LIB:.c=.o)

$(EXEC): $(OBJECTS) $(OBJECTS_LIB)
	$(CXX) $(OBJECTS) $(OBJECTS_LIB) $(CFLAGS) $(DFLAGS) -o $(EXEC) $(LIBS)

%.o: %.c
	$(CC) $(DFLAGS) $(CFLAGS) -c $(INCLUDE) $< -o $@

%.o: %.cpp
	$(CXX) $(DFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDE) $< -o $@

clean:
	rm -f $(EXEC) $(OBJECTS) $(OBJECTS_LIB)