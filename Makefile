CXX := g++
CXXFLAGS := -O3 -march=native -g -ggdb
all: mcsp

mcsp: solve_mcs.cc mcsp.cc graph.cc solve_mcs.hh graph.hh
	$(CXX) $(CXXFLAGS) -Wall -std=c++11 -o mcsp solve_mcs.cc graph.cc mcsp.cc -pthread
