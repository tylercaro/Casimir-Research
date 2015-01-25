#################################
# Simple Makefile with Variables
#################################

# Compiler
CXX = g++

# Compiler flags
CXXFLAGS += -Wall
#CXXFLAGS += -O3
CXXFLAGS += -O0
CXXFLAGS += -ggdb

# Include Directories
INC += -I/opt/local/include

# Libraries
LIBS += -L/opt/local/lib/
LIBS += -lgsl
LIBS += -lgslcblas
LIBS += -lm

# Build Rules
solved: solved.o langevin.o
	$(CXX) -o solved $(CXXFLAGS) $(LIBS) solved.o langevin.o

solved.o: solved.cpp odesolver.h
	$(CXX) -c -o solved.o $(CXXFLAGS) $(INC) solved.cpp

langevin.o: langevin.cpp odesolver.h
	$(CXX) -c -o langevin.o $(CXXFLAGS) $(INC) langevin.cpp

clean:
	rm solved solved.o langevin.o
