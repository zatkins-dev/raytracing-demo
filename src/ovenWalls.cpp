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

#include "ovenWalls.h"
#include <any>
#include <functional>
#include <unordered_map>

//  ==
//  ||
//  ||
//  ||  Main Program
//  ||
//  ||
//  ==

int main(int argc, char *argv[]) {
#pragma acc init
  int outer_nCellh = 10, outer_nCellc = 10;
  int inner_nCellh = 10, inner_nCellc = 10;
  double outer_r = 10., inner_r = outer_r / 2;
  double outer_h = 10., inner_h = outer_h;
  bool save = false, quiet = false;
  Mode mode = Mode::gpu_acc;

  std::unordered_map<std::string, std::function<void()>> flags{
      {"-save", [&]() { save = true; }},
      {"-q", [&]() { quiet = true; }},
      {"-quiet", [&]() { quiet = true; }},
  };

  std::unordered_map<std::string, std::function<void(std::string)>> arguments{
      {"-nCellh_o", [&](auto a) { outer_nCellh = std::stoi(a); }},
      {"-nCellc_o", [&](auto a) { outer_nCellc = std::stoi(a); }},
      {"-nCellh_i", [&](auto a) { inner_nCellh = std::stoi(a); }},
      {"-nCellc_i", [&](auto a) { inner_nCellc = std::stoi(a); }},
      {"-r_o", [&](auto a) { outer_r = std::stod(a); }},
      {"-r_i", [&](auto a) { inner_r = std::stod(a); }},
      {"-h_o", [&](auto a) { outer_h = std::stod(a); }},
      {"-h_i", [&](auto a) { inner_h = std::stod(a); }},
      {"-mode", [&](auto a) { mode = str2mode(a); }},
  };

  for (int count = 0; count < argc; ++count) {
    if (flags.find(argv[count]) != flags.end())
      flags[argv[count]]();
    else if (arguments.find(argv[count]) != arguments.end()) {
      arguments[argv[count]](argv[count + 1]);
      ++count;
    }
  }

  if (!quiet) {
    cout << "Ray Tracing Demo Code" << endl;
    cout << "Input Summary: " << endl;
    cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" << endl;
    cout << "Outer Cylinder" << endl;
    cout << " │Cells in circular dir: " << outer_nCellc << endl;
    cout << " │         height   dir: " << outer_nCellh << endl;
    cout << " │Radius               : " << outer_r << endl;
    cout << " │Height               : " << outer_h << endl << endl;
    cout << "Inner Cylinder" << endl;
    cout << " │Cells in circular dir: " << inner_nCellc << endl;
    cout << " │         height   dir: " << inner_nCellh << endl;
    cout << " │Radius               : " << inner_r << endl;
    cout << " │Height               : " << inner_h << endl << endl;
    cout << "Saving all plot data?  : " << std::boolalpha << save << endl;
    cout << "Mode: " << mode2str(mode) << endl;
    cout << endl << endl;
  }

  timingInfo timer;
  timer.start();

  const Cylinder CylA(outer_nCellh, outer_nCellc, outer_r, outer_h, false);
  const Cylinder CylB(inner_nCellh, inner_nCellc, inner_r, inner_h, true);
  auto intersections = array2d<int>(CylA.nFaces, CylA.nFaces);
  auto points = array2d<Point3>(CylA.nFaces, CylA.nFaces);
  const int nFacesA = CylA.nFaces;

  switch (mode) {
  case Mode::cpu_serial: { // Test CPU Performance, no parallelism
    for (int src_face = 0; src_face < nFacesA; src_face++) {
      for (int dst_face = 0; dst_face < nFacesA; ++dst_face) {
        // Ray to target
        Ray ray{CylA.faces[src_face].c, CylA.faces[dst_face].c};
        // Find blockers of ray
        auto [face, point] = CylB & ray;
        intersections[src_face][dst_face] = face;
        points[src_face][dst_face] = face >= 0 ? point : CylA.faces[dst_face].c;
      }
    }
    break;
  }

  case Mode::cpu_omp: { // Test CPU Performance, multithreading with OpenMP
#pragma omp parallel for collapse(2) shared(intersections) shared(points)
    for (int src_face = 0; src_face < nFacesA; src_face++) {
      for (int dst_face = 0; dst_face < nFacesA; ++dst_face) {
        // Ray to target
        Ray ray{CylA.faces[src_face].c, CylA.faces[dst_face].c};
        // Find blockers of ray
        auto [face, point] = CylB & ray;
        intersections[src_face][dst_face] = face;
        points[src_face][dst_face] = face >= 0 ? point : CylA.faces[dst_face].c;
      }
    }
    break;
  }

  case Mode::gpu_acc: // Test GPU Performance with OpenACC
  default: {
    CylA.todev();
    CylB.todev();
#pragma acc parallel loop independent collapse(2) gang vector,                 \
    present(CylA, CylB),                                                       \
    copyout(intersections [0:CylA.nFaces] [0:CylA.nFaces]),                    \
    copyout(points [0:CylA.nFaces] [0:CylA.nFaces]),
    for (int src_face = 0; src_face < nFacesA; src_face++) {
      for (int dst_face = 0; dst_face < nFacesA; dst_face++) {
        // Ray to target
        Ray ray{CylA.faces[src_face].c, CylA.faces[dst_face].c};
        // Find blockers of ray
        auto [face, point] = CylB & ray;
        intersections[src_face][dst_face] = face;
        points[src_face][dst_face] = face >= 0 ? point : CylA.faces[dst_face].c;
      }
    }
    break;
  }
  }

  auto result = timer.stop();

  cout << mode2str(mode) << " Time: " << result.cpu_s << " s (CPU Time), "
       << result.wall_s << " s (Wall Time)\n";

  if (save) {
    plot(CylA, "data/CylA");
    plot(CylB, "data/CylB");
    for (int i = 0; i < CylA.nFaces; i++) {
      VI facesHit, facesBlocking;
      for (int j = 0; j < CylA.nFaces; j++) {
        if (intersections[i][j] < 0)
          facesHit.push_back(j);
        else
          facesBlocking.push_back(intersections[i][j]);
      }

      char src_face_str[7];
      sprintf(src_face_str, "%06d", i);
      const string src_face{src_face_str};
      plotFacesInList(CylA, "data/hit_" + src_face, facesHit);
      plotFacesInList(CylB, "data/blockers_" + src_face, facesBlocking);
      plotFace(CylA, "data/source_" + src_face, i);

      std::ofstream hit("data/rays_hit_" + src_face + ".plt");
      std::ofstream blocked("data/rays_blocked_" + src_face + ".plt");
      for (int j = 0; j < CylA.nFaces; j++) {
        if (i == j)
          continue;
        if (intersections[i][j] < 0)
          hit << CylA.faces[i].c << '\n' << points[i][j] << "\n\n\n\n";
        else
          blocked << CylA.faces[i].c << '\n' << points[i][j] << "\n\n\n\n";
      }
    }
  }
  freeArray(intersections);
  freeArray(points);

  return 0;
}
