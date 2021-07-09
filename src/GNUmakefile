HOME = ./..
PROTO = $(HOME)/proto/include

DIM = 2
USE_BINARY = 0
CXX = g++
CXXFLAGS = -std=c++11
CXXFLAGS += -g -O0
CXXFLAGS += -D DIM=$(DIM)  -D VISIT_USE_BINARY=$(USE_BINARY) -I$(PROTO)
LIBS = matrix.cpp velocities.cpp interpolate.cpp w.cpp

new_project: GNUmakefile main.cpp
	$(CXX) main.cpp $(LIBS) $(CXXFLAGS) -g -o new.exe