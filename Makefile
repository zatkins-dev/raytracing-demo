

CXX      := pgc++
EXE		 := bin/raytracing
CXXFLAGS := -g -fast -mp -tp=zen3 -acc=multicore -acc=gpu -std=c++23 -Minit-msg -Minfo=accel -Minline -finline-functions -Midiom -Munroll -gpu=fastmath
SOURCES  := $(addprefix src/,ovenWalls.cpp Cylinder.hpp geom.hpp arrayUtils.hpp ovenWalls.h)
# -----------------------------------------------------
# Make esPIC
# -----------------------------------------------------

.PHONY: all
all: $(EXE)


$(EXE): src/ovenWalls.cpp $(SOURCES)
	@mkdir -p bin
	@mkdir -p animation
	@mkdir -p data
	$(CXX) $(CXXFLAGS) src/ovenWalls.cpp -o $(EXE)

.PHONY: clean
clean: clean-artifacts
	rm -f $(EXE)  *.o

.PHONY: clean-artifacts
clean-artifacts:
	rm -f data/* animation/* plot/*