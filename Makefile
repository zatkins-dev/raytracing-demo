

CXX      := pgc++
EXE		 := bin/raytracing
CXXFLAGS := -g -fast -fastsse -O3 -mp=multicore -acc=gpu -std=c++23 -Minit-msg -Minfo=accel -Minline -finline-functions -Midiom -Munroll -mavx2 -mcmodel=medium -gpu=fastmath -gpu=nvlamath -gpu=fma -nostdpar
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