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
    if (arguments.contains(argv[count])) {
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

  Cylinder CylA(outer_nCellh, outer_nCellc, outer_r, outer_h, 2);
  Cylinder CylB(inner_nCellh, inner_nCellc, inner_r, inner_h, 1);

  CylA.plot("data/CylA");
  CylB.plot("data/CylB");

  // Set up parameters for ray tracing

  VI facesBlocking;
  VI facesHit;
  std::vector<Point3> pointsBlocking;
  std::vector<Point3> pointsHit;

  CylA.plotFace("data/source", source_face_id);

  // (1) Point on source face
  const Point3 sourceCenter = CylA.centers[source_face_id];
  for (int targetFaceID = 0; targetFaceID < CylA.nFaces; ++targetFaceID) {
    if (target_face_id >= 0)
      targetFaceID = target_face_id;
    //  (2) Ray to target
    Ray ray(sourceCenter, CylA.centers[targetFaceID]);
    // if (fabs(dot(ray.n, CylA.normals[sourceFaceID])) < 1e-10)
    //   continue;

    // (3) Plot ray
    CylA.plotRay("data/ray", ray);

    // (4) Find blockers of ray
    auto maybe_intersection = CylB & ray;
    if (maybe_intersection) {
      facesBlocking.push_back(maybe_intersection->face_id);
      pointsBlocking.push_back(maybe_intersection->p);
    } else {
      auto intersect = CylA & ray;
      facesHit.push_back(targetFaceID);
      if (intersect) {
        pointsHit.push_back(intersect->p);
      }
    }
    if (target_face_id >= 0)
      break;
  }

  // Plot hits and blockers

  CylA.plotFacesInList("data/hit", facesHit);
  CylB.plotFacesInList("data/blockers", facesBlocking);

  {
    std::ofstream file("data/blockers_points.plt");
    for (auto &p : pointsBlocking) {
      file << sourceCenter << endl;
      file << p << endl;
      file << "\n\n\n";
    }
  }

  {
    std::ofstream file("data/hit_points.plt");
    for (auto &p : pointsHit) {
      file << sourceCenter << endl;
      file << p << endl;
      file << "\n\n\n";
    }
  }

  return 0;
}
