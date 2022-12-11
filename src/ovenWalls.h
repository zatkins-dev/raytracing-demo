//  ================================================================================
//  || ||
//  ||              ovenWalls ||
//  ||              ------------------------------------------------------ ||
//  ||              T H E R M A L   R A D I A T I O N ||
//  || ||
//  ||              D E M O N S T R A T I O N   C O D E ||
//  ||              ------------------------------------------------------ ||
//  || ||
//  ||       Developed by: Scott R. Runnels, Ph.D. ||
//  ||                     University of Colorado Boulder ||
//  || ||
//  ||                For: CU Boulder CSCI 4576/5576 and associated labs ||
//  || ||
//  ||           Copyright 2020 Scott Runnels ||
//  || ||
//  ||                     Not for distribution or use outside of the ||
//  ||                     this course. ||
//  || ||
//  ================================================================================

#include <chrono>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <openacc.h>
#include <sstream>
#include <string>
#include <vector>

#include "Cylinder.hpp"
#include "geom.hpp"

using std::cout;
using std::endl;
using std::string;
using std::stringstream;
using std::vector;

typedef vector<double> VD;
typedef vector<vector<double>> VDD;
typedef vector<int> VI;
typedef vector<vector<int>> VII;

void FatalError(string msg) {

  cout << " " << endl;
  cout << " " << endl;
  cout << "Fatal Error: " << msg << endl;
  cout << " " << endl;
  exit(0);
}

class timingInfo {
public:
  // Name of the intenty for which timing is desired
  struct time_result_t {
    double wall_s = 0;
    double cpu_s = 0;
  };

  string name = "default";

  // Two methods of recording time:

  std::chrono::steady_clock::time_point startSeconds, endSeconds; // (1)
  struct timespec startHighResTime, endHighResTime;               // (2)

  timingInfo() = default;

  timingInfo(string _name) : name(_name) {}

  void start() {
    // Record start time using two methods

    startSeconds = std::chrono::steady_clock::now();            // (1)
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &startHighResTime); // (2)
  }

  time_result_t finish() {

    // Record finish time using two methods

    endSeconds = std::chrono::steady_clock::now();            // (1)
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &endHighResTime); // (2)

    // Compute elapsed time
    std::chrono::duration<double> wall = endSeconds - startSeconds;
    double elapsedWall = wall.count();
    double elapsedCPU =
        timespecsToDuration(startHighResTime, endHighResTime).count(); // (2)

    return time_result_t{.wall_s = elapsedWall, .cpu_s = elapsedCPU};
  }

private:
  std::chrono::duration<double> timespecsToDuration(struct timespec ts0,
                                                    struct timespec ts1) {
    using std::chrono::nanoseconds;
    using std::chrono::seconds;
    return seconds{ts1.tv_sec - ts0.tv_sec} +
           nanoseconds{ts1.tv_nsec - ts0.tv_nsec};
  }
};

inline void dump(std::string dest, int **intersections, int imax, int jmax) {
  std::ofstream out(dest);
  for (int i = 0; i < imax; i++) {
    for (int j = 0; j < jmax; j++) {
      out << intersections[i][j] << " ";
    }
    out << "\n";
  }
}
