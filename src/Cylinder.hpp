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
  int nCellh;      // Number of cells in the height direction
  int nCellc;      // Number of cells in the circular direction
  int nRealh;      // Number of nodes in the height direction
  int nRealc;      // Number of nodes in the circular direction
  int nField;      // Number of nodes
  double dtheta;   // Angular spacing of nodes, in degrees
  double dz;       // Axial spacing of nodes
  double radius;   // Cylinder radius
  double length;   // Cylinder length
  int **face;      // Faces:  face[f][i] gives the nodes of face f
  int nFaces;      // Number of faces
  Point3 *coord;   // Vertices
  Vec3 *normals;   // Face normals
  Point3 *centers; // Face centers

  Cylinder(int _nCellh, int _nCellc, double _radius, double _length,
           int Out1In2)
      : nCellh(_nCellh), nCellc(_nCellc), radius(_radius), length(_length) {
    // Store user inputs
    nRealh = nCellh + 1;
    nRealc = nCellc;

    nField = nRealh * nRealc;
    nFaces = nCellh * nCellc;

    dtheta = 360. / nCellc;
    dz = length / nCellh;

    // Form mesh

    coord = array<Point3>(nField);
    face = array2d<int>(nFaces, 4);
    normals = array<Vec3>(nFaces);
    centers = array<Point3>(nFaces);

    int faceCount = 0;

    for (int j = 1; j <= nRealc; ++j) {
      for (int i = 1; i <= nRealh; ++i) {
        int p = pid(i, j);

        double theta = dtheta * j;
        double zval = dz * (i - 1);

        coord[p] = {
            radius * std::cos(theta * M_PI / 180.),
            radius * std::sin(theta * M_PI / 180.),
            zval,
        };

        // Form face, i.e., for this face number (which happens to be p)
        // collect the pids of the 4 nodes that comprise this face.

        const int point[4]{
            pid(i, j),
            pid(i + 1, j),
            j == nRealc ? pid(i + 1, 1) : pid(i + 1, j % nRealc + 1),
            j == nRealc ? pid(i, 1) : pid(i, j % nRealc + 1),
        };

        // Correct point for when we have wrapped completely around the
        // circle

        // Store the point values in this face, p (but not for the last point in
        // an i-row

        if (i < nRealh) {
          for (int k = 0; k < 4; ++k)
            face[faceCount][k] = point[k];
          faceCount++;
        }
      }
    }

    assert(faceCount == nFaces);

    // Compute normals for each face

    for (int f = 0; f < nFaces; ++f) {
      Point3 avg;
      for (int k = 0; k < 4; ++k) {
        avg += coord[face[f][k]];
      }
      avg *= 0.25;

      double mag = sqrt(avg.x * avg.x + avg.y * avg.y);

      normals[f] = {
          avg.x / mag,
          avg.y / mag,
          0.,
      };

      centers[f] = avg;

      if (Out1In2 == 2)
        normals[f] *= -1.;
    }
#pragma acc enter data copyin(this[:1])
#pragma acc enter data copyin(                                                 \
    coord[:nField], face[:nFaces][:4], normals[:nFaces], centers[:nFaces])
  }

  void updatehost() { // update host copy of data
#pragma acc update self(coord [0:nField], face [0:nFaces] [0:4],               \
                        normals [0:nFaces], centers [0:nFaces])
  }
  void updatedev() { // update device copy of data
#pragma acc update device(coord [0:nField], face [0:nFaces] [0:4],             \
                          normals [0:nFaces], centers [0:nFaces])
  }

#pragma acc routine seq
  constexpr int operator&(const Ray &ray) const noexcept {
    int intersection = -1;

#pragma acc loop reduction(max : intersection)
    for (int f = 0; f < nFaces; f++) {
      // Vertices of potential blocker
      const Point3 c[]{
          coord[face[f][0]],
          coord[face[f][1]],
          coord[face[f][2]],
          coord[face[f][3]],
      };
      Face t_face{c, normals[f]};

      // Look for intersection
      // ignore the source face
      if (normsq(t_face.c - ray.p) < 1e-12)
        continue;
      intersection = ray & t_face ? std::max(f, intersection) : intersection;
    }
    return intersection;
  }

  int pid(int i, int j) { return (i - 1 + (j - 1) * nRealh); }
};

void plot(const Cylinder &c, string descriptor) {
  // Open plot file
  std::ofstream file(descriptor + ".plt");

  // Write to plot file
  for (int f = 0; f < c.nFaces; ++f) {
    const std::array<int, 5> pts{
        c.face[f][0], c.face[f][1], c.face[f][2], c.face[f][3], c.face[f][0],
    };

    for (auto &p : pts)
      file << c.coord[p] << endl;

    file << "\n\n\n";
  }

  std::ofstream file1(descriptor + "_centers.plt");

  // Write to plot file
  for (int f = 0; f < c.nFaces; ++f) {
    file1 << c.centers[f] << endl;
  }
}

void plotFacesInList(const Cylinder &c, string descriptor,
                     std::vector<int> &facesToPlot) {
  // Open plot file
  std::ofstream file(descriptor + ".plt");

  // Write to plot file
  for (int i = 0; i < facesToPlot.size(); ++i) {
    const std::array<int, 5> pts{
        c.face[facesToPlot[i]][0], c.face[facesToPlot[i]][1],
        c.face[facesToPlot[i]][2], c.face[facesToPlot[i]][3],
        c.face[facesToPlot[i]][0],
    };

    for (auto &p : pts)
      file << c.coord[p] << endl;

    file << "\n\n\n";
  }
}

void plotFace(const Cylinder &c, string descriptor, int faceID) {
  // Open plot file
  std::ofstream file(descriptor + ".plt");

  // Write to plot file
  std::array<int, 5> pts{
      c.face[faceID][0], c.face[faceID][1], c.face[faceID][2],
      c.face[faceID][3], c.face[faceID][0],
  };
  for (auto &p : pts)
    file << c.coord[p] << endl;

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
