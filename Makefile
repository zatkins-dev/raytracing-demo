

CXX      := pgc++
EXE		 := bin/raytracing
EXE_NOMP := bin/raytracing_nomp
CXXFLAGS := -g -fast -mp -tp=zen3 -acc=multicore -acc=gpu -std=c++23 -Minit-msg -Minfo=accel -Minline -finline-functions -Midiom -Munroll -gpu=fastmath
SOURCES  := $(addprefix src/,ovenWalls.cpp Cylinder.hpp geom.hpp arrayUtils.hpp ovenWalls.h)
# -----------------------------------------------------
# Make esPIC
# -----------------------------------------------------

.PHONY: all
all: $(EXE) $(EXE_NOMP)


$(EXE): src/ovenWalls.cpp $(SOURCES)
	@mkdir -p bin
	@mkdir -p animation
	@mkdir -p data
	$(CXX) $(CXXFLAGS)  src/ovenWalls.cpp -o $(EXE)

$(EXE_NOMP): src/ovenWalls.cpp $(SOURCES)
	$(CXX) $(CXXFLAGS) -Mnoopenmp src/ovenWalls.cpp -o $(EXE_NOMP)

.PHONY: clean
clean: clean-artifacts
	rm -f $(EXE) $(EXE_NOMP) *.o

.PHONY: clean-artifacts
clean-artifacts:
	rm -f data/* animation/* plot/*