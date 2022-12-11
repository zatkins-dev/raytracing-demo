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
  int outer_nCellh = 10, outer_nCellc = 10;
  int inner_nCellh = 10, inner_nCellc = 10;
  double outer_r = 10., inner_r = outer_r / 2;
  double outer_h = 10., inner_h = outer_h;
  int source_face_id = 0, target_face_id = -1;

  std::unordered_map<std::string, std::function<void(std::string)>> arguments{
      {"-nCellh_o", [&](auto a) { outer_nCellh = std::stoi(a); }},
      {"-nCellc_o", [&](auto a) { outer_nCellc = std::stoi(a); }},
      {"-nCellh_i", [&](auto a) { inner_nCellh = std::stoi(a); }},
      {"-nCellc_i", [&](auto a) { inner_nCellc = std::stoi(a); }},
      {"-r_o", [&](auto a) { outer_r = std::stod(a); }},
      {"-r_i", [&](auto a) { inner_r = std::stod(a); }},
      {"-h_o", [&](auto a) { outer_h = std::stod(a); }},
      {"-h_i", [&](auto a) { inner_h = std::stod(a); }},
      {"-source", [&](auto a) { source_face_id = std::stoi(a); }},
      {"-target", [&](auto a) { target_face_id = std::stoi(a); }},
  };

  for (int count = 0; count < argc; ++count) {
    if (arguments.find(argv[count]) != arguments.end()) {
      arguments[argv[count]](argv[count + 1]);
      ++count;
    }
  }

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
  cout << "Source face id         : " << source_face_id << endl;
  cout << "Target face id         : " << target_face_id << endl;
  cout << endl << endl;

  const Cylinder CylA(outer_nCellh, outer_nCellc, outer_r, outer_h, 2);
  const Cylinder CylB(inner_nCellh, inner_nCellc, inner_r, inner_h, 1);

  plot(CylA, "data/CylA");
  plot(CylB, "data/CylB");

  // CPU Test
  {
    timingInfo cpu;
    cpu.start();
    // Set up parameters for ray tracing
    auto intersections = array2d<int>(CylA.nFaces, CylA.nFaces);

    plotFace(CylA, "data/source", source_face_id);

    for (int sourceFaceID = 0; sourceFaceID < CylA.nFaces; sourceFaceID++) {
      // (1) Point on source face
      const Point3 sourceCenter = CylA.centers[source_face_id];

      for (int targetFaceID = 0; targetFaceID < CylA.nFaces; ++targetFaceID) {
        //  (2) Ray to target
        Ray ray{sourceCenter, CylA.centers[targetFaceID]};

        // (4) Find blockers of ray
        intersections[sourceFaceID][targetFaceID] = CylB & ray;
      }
    }

    auto result = cpu.finish();
    cout << "CPU Test: " << result.cpu_s << " s (CPU Time), " << result.wall_s
         << " s (Wall Time)\n";

    VI facesHit;
    VI facesBlocking;
    for (int i = 0; i < CylA.nFaces; i++) {
      if (intersections[source_face_id][i] < 0)
        facesHit.push_back(i);
      else
        facesBlocking.push_back(intersections[source_face_id][i]);
    }

    plotFacesInList(CylA, "data/hit", facesHit);
    plotFacesInList(CylB, "data/blockers", facesBlocking);
  }

  // GPU Test
  {
    timingInfo gpu;
    gpu.start();

    // Set up parameters for ray tracing
    auto intersections = array2d<int>(CylA.nFaces, CylA.nFaces);

    plotFace(CylA, "data/source", source_face_id);

    // (1) Point on source face
    const Point3 sourceCenter = CylA.centers[source_face_id];
#pragma acc data pcopyin(CylA), pcopyin(CylB),                                 \
    pcopyout(intersections [0:CylA.nFaces] [0:CylA.nFaces]),                   \
    pcopyin(sourceCenter)
    {
#pragma acc parallel loop collapse(2)
      for (int sourceFaceID = 0; sourceFaceID < CylA.nFaces; sourceFaceID++) {
        for (int targetFaceID = 0; targetFaceID < CylA.nFaces; ++targetFaceID) {
          //  (2) Ray to target
          Ray ray{CylA.centers[sourceFaceID], CylA.centers[targetFaceID]};

          // (4) Find blockers of ray
          intersections[sourceFaceID][targetFaceID] = CylB & ray;
        }
      }
    }

    auto result = gpu.finish();
    cout << "GPU Test: " << result.cpu_s << " s (CPU Time), " << result.wall_s
         << " s (Wall Time)\n";

    VI facesHit;
    VI facesBlocking;
    for (int i = 0; i < CylA.nFaces; i++) {
      if (intersections[source_face_id][i] < 0)
        facesHit.push_back(i);
      else
        facesBlocking.push_back(intersections[source_face_id][i]);
    }

    plotFacesInList(CylA, "data/hit", facesHit);
    plotFacesInList(CylB, "data/blockers", facesBlocking);
  }

  // {
  //   std::ofstream file("data/blockers_points.plt");
  //   for (auto &p : pointsBlocking) {
  //     file << sourceCenter << endl;
  //     file << p << endl;
  //     file << "\n\n\n";
  //   }
  // }

  // {
  //   std::ofstream file("data/hit_points.plt");
  //   for (auto &p : pointsHit) {
  //     file << sourceCenter << endl;
  //     file << p << endl;
  //     file << "\n\n\n";
  //   }
  // }

  return 0;
}
