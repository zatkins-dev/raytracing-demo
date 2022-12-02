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
  struct PointOnFace {
    Point3 p;
    int face_id;
  };

public:
  using intersect_result_t = std::optional<PointOnFace>;
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

  Cylinder() {}

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
            .x = radius * std::cos(theta * M_PI / 180.),
            .y = radius * std::sin(theta * M_PI / 180.),
            .z = zval,
        };

        // Form face, i.e., for this face number (which happens to be p)
        // collect the pids of the 4 nodes that comprise this face.

        int point[4];
        point[0] = pid(i, j);
        point[1] = pid(i + 1, j);
        point[2] = pid(i + 1, j % nRealc + 1);
        point[3] = pid(i, j % nRealc + 1);

        // Correct point for when we have wrapped completely around the circle

        if (j == nRealc) {
          point[2] = pid(i + 1, 1);
          point[3] = pid(i, 1);
        }

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
          .x = avg.x / mag,
          .y = avg.y / mag,
          .z = 0.,
      };

      centers[f] = avg;

      if (Out1In2 == 2)
        normals[f] *= -1.;
    }
  }

  constexpr intersect_result_t operator&(const Ray &ray) const {
    int blockerFaceID = 0;
    bool blocked = false;
    double min_dist = INFINITY;
    Point3 closest_point{};
    int closest_face = -1;

    for (int f = 0; f < nFaces; f++) {
      // Vertices of potential blocker
      Face t_face{
          {
              coord[face[f][0]],
              coord[face[f][1]],
              coord[face[f][2]],
              coord[face[f][3]],
          },
          normals[f],
      };

      // Look for intersection
      // ignore the source face
      if (normsq(t_face.c - ray.p) < 1e-12)
        continue;
      auto intersection = ray & t_face;
      if (intersection) {
        double dist = normsq(*intersection - ray.p);
        if (dist < min_dist) {
          min_dist = dist;
          closest_face = f;
          closest_point = *intersection;
        }
      }
    }
    if (closest_face > 0)
      return PointOnFace{closest_point, closest_face};
    else
      return std::nullopt;
  }

  int pid(int i, int j) { return (i - 1 + (j - 1) * nRealh); }

  void plot(string descriptor) {
    // Open plot file
    std::ofstream file(descriptor + ".plt");

    // Write to plot file
    for (int f = 0; f < nFaces; ++f) {
      const std::array<int, 5> pts{
          face[f][0], face[f][1], face[f][2], face[f][3], face[f][0],
      };

      for (auto &p : pts)
        file << coord[p] << endl;

      file << "\n\n\n";
    }

    std::ofstream file1(descriptor + "_centers.plt");

    // Write to plot file
    for (int f = 0; f < nFaces; ++f) {
      file1 << centers[f] << endl;
    }
  }

  void plotFacesInList(string descriptor, std::vector<int> &facesToPlot) {
    // Open plot file
    std::ofstream file(descriptor + ".plt");

    // Write to plot file
    for (int i = 0; i < facesToPlot.size(); ++i) {
      const std::array<int, 5> pts{
          face[facesToPlot[i]][0], face[facesToPlot[i]][1],
          face[facesToPlot[i]][2], face[facesToPlot[i]][3],
          face[facesToPlot[i]][0],
      };

      for (auto &p : pts)
        file << coord[p] << endl;

      file << "\n\n\n";
    }
  }

  void plotFace(string descriptor, int faceID) {
    // Open plot file
    std::ofstream file(descriptor + ".plt");

    // Write to plot file
    std::array<int, 5> pts{
        face[faceID][0], face[faceID][1], face[faceID][2],
        face[faceID][3], face[faceID][0],
    };
    for (auto &p : pts)
      file << coord[p] << endl;

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
};