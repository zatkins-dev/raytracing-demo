#pragma once

#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using std::cout, std::endl;
using std::string;

#include "arrayUtils.hpp"
#include "geom.hpp"

class Cylinder {
public:
  const int nCellh;    // Number of cells in the height direction
  const int nCellc;    // Number of cells in the circular direction
  const int nRealh;    // Number of nodes in the height direction
  const int nRealc;    // Number of nodes in the circular direction
  const int nFaces;    // Number of faces
  const double radius; // Cylinder radius
  const double height; // Cylinder height
  const Mode mode;
  Face *faces; // Faces: faces[f] is a Face containing the vertices and normal

  Cylinder(int _nCellh, int _nCellc, double _radius, double _length, bool out,
           Mode _mode)
      : nCellh(_nCellh), nCellc(_nCellc), radius(_radius), height(_length),
        nRealh(_nCellh + 1), nRealc(_nCellc), nFaces(_nCellc * _nCellh),
        mode(_mode) {
    // Form mesh
    faces = new Face[nFaces];

    double dtheta = 360. / nCellc;
    double dz = height / nCellh;
    auto compute_coord = [=](int i, int j) -> Point3 {
      return {
          _radius * std::cos(dtheta * j * M_PI / 180.),
          _radius * std::sin(dtheta * j * M_PI / 180.),
          dz * i,
      };
    };
    for (int j = 0; j < _nCellc; ++j) {
      for (int i = 0; i < _nCellh; ++i) {
        // Form face, i.e., for this face number (which happens to be p)
        // collect the pids of the 4 nodes that comprise this face.
        const Point3 q1 = compute_coord(i, j), q2 = compute_coord(i + 1, j),
                     q3 = compute_coord(i + 1, (j + 1) % nRealc),
                     q4 = compute_coord(i, (j + 1) % nRealc);

        const Point3 avg = 0.25 * (q1 + q2 + q3 + q4);
        const double mag = sqrt(avg.x * avg.x + avg.y * avg.y);
        const Vec3 normal = (2 * out - 1) * Vec3{avg.x / mag, avg.y / mag, 0.};
        faces[i + j * nCellh] = Face(q1, q2, q3, q4, normal);
      }
    }
    if (mode == Mode::gpu_acc) {
#pragma acc enter data copyin(this [0:1])
      const int n = nFaces;
#pragma acc enter data copyin(faces [0:n])
    }
  }

  ~Cylinder() { delete[] faces; }

#pragma acc routine seq
  constexpr auto operator&(const Ray &ray) const noexcept {
    for (int f = 0; f < nFaces; f++) { // Check for front-face intersections
      auto [intersects, p] = ray & faces[f];
      if (intersects)
        return std::make_pair(f, p);
    }
    for (int f = 0; f < nFaces; f++) {
      // Check for back-face intersections ONLY IF no front-face intersections
      auto [intersects, p] = ray & -faces[f];
      if (intersects)
        return std::make_pair(f, p);
    }
    return std::make_pair(-1, Point3{});
  }

#pragma acc routine seq
  int pid(int i, int j) { return (i + j * nRealh); }
};

void plot(const Cylinder &c, string descriptor) {
  // Open plot file
  std::ofstream file(descriptor + ".plt");

  // Write to plot file
  for (int f = 0; f < c.nFaces; ++f) {
    const auto &fc = c.faces[f];
    file << fc.q1 << '\n';
    file << fc.q2 << '\n';
    file << fc.q3 << '\n';
    file << fc.q4 << '\n';
    file << fc.q1 << '\n';
    file << "\n\n\n";
  }
}

void plotFacesInList(const Cylinder &c, string descriptor,
                     std::vector<int> &facesToPlot) {
  // Open plot file
  std::ofstream file(descriptor + ".plt");

  // Write to plot file
  for (int i = 0; i < facesToPlot.size(); ++i) {
    const auto &fc = c.faces[facesToPlot[i]];
    file << fc.q1 << '\n';
    file << fc.q2 << '\n';
    file << fc.q3 << '\n';
    file << fc.q4 << '\n';
    file << fc.q1 << '\n';
    file << "\n\n\n";
  }
}

void plotFace(const Cylinder &c, string descriptor, int faceID) {
  // Open plot file
  std::ofstream file(descriptor + ".plt");
  const auto &fc = c.faces[faceID];
  // Write to plot file
  file << fc.q1 << '\n';
  file << fc.q2 << '\n';
  file << fc.q3 << '\n';
  file << fc.q4 << '\n';
  file << fc.q1 << '\n';
  file << "\n\n\n";
}

void plotRay(string descriptor, const Ray &ray) {
  // Open plot file
  std::ofstream file(descriptor + ".plt", std::ios::app);

  file << ray.p << endl;
  file << ray.p + ray.n << endl;
  file << ray.p << endl;
  file << endl;
}

void plotPoint(string descriptor, const Point3 &p) {
  // Open plot file
  std::ofstream file(descriptor + ".plt");
  file << p << endl;
}
