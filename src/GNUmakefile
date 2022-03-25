HOME = ./..
PROTO = $(HOME)/proto/include

DIM = 2
USE_BINARY = 0
CXX = clang++
CXXFLAGS = -std=c++11
CXXFLAGS += -g -O0 # debug flags
CXXFLAGS += -D DIM=$(DIM)  -D VISIT_USE_BINARY=$(USE_BINARY) -I$(PROTO)
LIBS = matrix.cpp velocities.cpp interpolate.cpp w.cpp writers.cpp

tests: GNUmakefile main.cpp
	$(CXX) main.cpp $(LIBS) $(CXXFLAGS) -o tests-$(DIM)D.exe

clean:
	rm -r *.exe *.vtk
