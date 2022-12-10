
# -----------------------------------------------------
# On local Ubuntu
# -----------------------------------------------------
#
# Be sure to use MPI compiler wrappers.
# Whenever compiling an MPI program, you should use the MPI wrappers:
#
#    C - mpicc
#    C++ - mpiCC, mpicxx, mpic++
#    FORTRAN - mpifort, mpif77, mpif90
#
# These wrappers do all of the dirty work for you of making sure that
# all of the appropriate compiler flags, libraries, include directories,
# library directories, etc. are included when you compile your program.


# -----------------------------------------------------
# On Summit:
# -----------------------------------------------------
#
# module load intel
# module load impi

# -----------------------------------------------------
# Set compiler based on platform
# -----------------------------------------------------

# Default compiler is Intel MPI compiler

CXX      := mpicxx
EXE		 := bin/raytracing
CXXFLAGS := -g -O0 -std=c++2a
SOURCES  := $(addprefix src/,ovenWalls.cpp Cylinder.hpp geom.hpp arrayUtils.hpp)
CPPOBJECTS  := $(addprefix bin/,$(basename $(SOURCES:.cpp=.o)))
HPPOBJECTS  := $(addprefix bin/,$(basename $(SOURCES:.hpp=.o)))
# -----------------------------------------------------
# Make esPIC
# -----------------------------------------------------


$(EXE): src/ovenWalls.cpp $(SOURCES)
	@mkdir -p bin
	@mkdir -p animation
	@mkdir -p data
	$(CXX) $(CXXFLAGS) src/ovenWalls.cpp -o $(EXE)

.PHONY: clean
clean: clean-artifacts
	rm -f $(EXE)

.PHONY: clean-artifacts
clean-artifacts:
	rm -f data/* animation/*