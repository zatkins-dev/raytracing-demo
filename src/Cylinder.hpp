#pragma once

#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using std::cout, std::endl;
using std::sin, std::cos;
using std::string;

#include "arrayUtils.hpp"
#include "geom.hpp"

class Cylinder {
public:
  const int nFaces; // Number of faces
  const double r;   // Cylinder radius
  const double h;   // Cylinder height
  Face *faces; // Faces: faces[f] is a Face containing the vertices and normal

  Cylinder(int nCellh, int nCellc, double radius, double height, bool outward_n)
      : r(radius), h(height), nFaces(nCellc * nCellh) {
    faces = new Face[nFaces];
    const double dt = 2. * M_PI / nCellc;
    const double dz = h / nCellh;
    for (int j = 0; j < nCellc; ++j) {
      for (int i = 0; i < nCellh; ++i) {
        const int jp1 = (j + 1) % nCellc;
        const Point3 q1 = {r * cos(dt * j), r * sin(dt * j), dz * i},
                     q2 = {r * cos(dt * j), r * sin(dt * j), dz * (i + 1)},
                     q3 = {r * cos(dt * jp1), r * sin(dt * jp1), dz * (i + 1)},
                     q4 = {r * cos(dt * jp1), r * sin(dt * jp1), dz * i};

        const Point3 avg = 0.25 * (q1 + q2 + q3 + q4);
        const double mag = sqrt(avg.x * avg.x + avg.y * avg.y);
        const int scale_normal = (2 * outward_n - 1);
        const Vec3 normal = scale_normal * Vec3{avg.x / mag, avg.y / mag, 0.};
        faces[i + j * nCellh] = Face(q1, q2, q3, q4, normal);
      }
    }
  }

  ~Cylinder() { delete[] faces; }

  inline void todev() const {
#pragma acc enter data copyin(this [0:1])
#pragma acc enter data copyin(faces [0:nFaces])
  }

  inline void fromdev() {
#pragma acc exit data copyout(faces [0:nFaces])
#pragma acc exit data copyout(this [0:1])
  }

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
