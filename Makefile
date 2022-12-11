
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

CXX      := pgc++
EXE		 := bin/raytracing
# CXXFLAGS := -g -O3 -std=c++17 -Minit-msg -Minfo=accel -Minline -finline-functions -Mipa=acc
CXXFLAGS := -g -O3 -acc -std=c++23 -Minit-msg -Minfo=accel -Minline -finline-functions
SOURCES  := $(addprefix src/,ovenWalls.cpp Cylinder.hpp geom.hpp arrayUtils.hpp)
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
	rm -f $(EXE) *.o

.PHONY: clean-artifacts
clean-artifacts:
	rm -f data/* animation/*