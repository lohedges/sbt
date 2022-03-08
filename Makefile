.DEFAULT_GOAL := all

CXX := g++
CXXFLAGS := -O3 -Isrc

ABIFLAG := 1
POPLARFLAGS :=-std=c++17 -L/opt/poplar/lib -lpoplar -lpopops -lpoputil

OMPFLAGS := -fopenmp

all:
	$(CXX) $(CXXFLAGS) -D_GLIBCXX_USE_CXX11_ABI=$(ABIFLAG) $(POPLARFLAGS) $(OMPFLAGS) $(OPTFLAGS) src/main.cpp src/EventReader.cpp src/SearchByTripletIPU.cpp -o search_by_triplet

clean:
	rm -f search_by_triplet
