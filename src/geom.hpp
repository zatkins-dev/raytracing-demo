/**
 * @file geom.hpp
 * @author Zach Atkins (zach.atkins@colorado.edu)
 * @brief Geometric utility functions and types
 * @date 2022-11-18
 *
 * @copyright Copyright (c) 2022
 *
 **/
#pragma once

#include <array>
#include <cmath>
#include <iostream>
#include <openacc.h>
#include <optional>

struct Vec3 {
  double x = 0., y = 0., z = 0.;

  template <class T> constexpr const Vec3 operator+(const T &w) const {
    return {x + w.x, y + w.y, z + w.z};
  }

  constexpr const Vec3 operator*(double a) const {
    return {a * x, a * y, a * z};
  }

  constexpr const Vec3 operator-() const { return *this * -1; }

  constexpr Vec3 &operator*=(double a) {
    *this = *this * a;
    return *this;
  }
};

struct Point3 {
  double x = 0., y = 0., z = 0.;

  template <class T> constexpr const Point3 operator+(const T &w) const {
    return {x + w.x, y + w.y, z + w.z};
  }

  constexpr const Point3 operator*(double a) const {
    return {a * x, a * y, a * z};
  }

  constexpr const Point3 operator-() const { return *this * -1; }

  constexpr const Vec3 operator-(const Point3 &v) const {
    return {x - v.x, y - v.y, z - v.z};
  }

  template <class T> constexpr Point3 &operator+=(T v) {
    *this = *this + v;
    return *this;
  }
  constexpr Point3 &operator*=(double a) {
    *this = *this * a;
    return *this;
  }
};

template <class T> constexpr T operator*(double a, const T &v) { return v * a; }

std::ostream &operator<<(std::ostream &os, const Vec3 &v) {
  return os << v.x << " " << v.y << " " << v.z;
}

std::ostream &operator<<(std::ostream &os, const Point3 &v) {
  return os << v.x << " " << v.y << " " << v.z;
}

constexpr Vec3 cross(const Vec3 &v, const Vec3 &w) {
  return {
      v.y * w.z - v.z * w.y,
      v.z * w.x - v.x * w.z,
      v.x * w.y - v.y * w.x,
  };
}

constexpr double dot(const Vec3 &v, const Vec3 &w) {
  return v.x * w.x + v.y * w.y + v.z * w.z;
}

constexpr double normsq(const Vec3 &v) {
  return v.x * v.x + v.y * v.y + v.z * v.z;
}

constexpr double norm(const Vec3 &v) { return std::sqrt(normsq(v)); }

constexpr Vec3 normalized(const Vec3 &v) { return (1. / norm(v)) * v; }

class Face {
public:
  Point3 verts[4];
  Point3 c;
  Vec3 n;
  const Point3 &q1 = verts[0], &q2 = verts[1], &q3 = verts[2], &q4 = verts[3];

  constexpr Face(const Point3 verts_[4], Vec3 n_) : n(n_) {
    for (int i = 0; i < 4; i++) {
      verts[i] = verts_[i];
      c += verts[i];
    }
    c *= 0.25;
  }
};

struct Plane {
  Point3 p;
  Vec3 n;
};

class Ray {
public:
  Point3 p;
  Vec3 n;

  constexpr Ray(Point3 p_, Vec3 n_) : p(p_), n(normalized(n_)) {}

  constexpr Ray(Point3 p1, Point3 p2) {
    p = p1;
    n = p2 - p1;
  }
};

/**
 * @brief Determines whether P is within the corner CAB
 *
 * (1) Form vector pointing from A to B: vAB
 * (2) Form vector pointing from C to A: vCA
 * (3) Form vector pointing from A to P: vAP
 * (4) P is inside the C-A-B corner if the cross products
 *     of these two vectors are aligned, i.e.,
 *     dot(cA, cB) > 0, cA = vAB x vAP  and cB = vCA x vAP
 *
 * @param P
 * @param A
 * @param B
 * @param C
 * @return true if dot(vAB x vAP, vCA x vAP) > 0
 **/
#pragma acc routine seq
constexpr bool insideCorner(const Point3 &P, const Point3 &C, const Point3 &A,
                            const Point3 &B) noexcept {
  return dot(cross(B - A, P - A), cross(A - C, P - A)) >= 0;
}

#pragma acc routine seq
constexpr bool onEdge(const Point3 &P, const Point3 &A, const Point3 &B) {
  return std::abs(norm(P - A) + norm(B - P) - norm(A - B)) < 1e-12;
}

#pragma acc routine seq
constexpr bool onVertex(const Point3 &P, const Face &face) {
  return normsq(P - face.q1) < 1e-12 || normsq(P - face.q2) < 1e-12 ||
         normsq(P - face.q3) < 1e-12 || normsq(P - face.q4) < 1e-12;
}

#pragma acc routine seq
constexpr bool onEdge(const Point3 &P, const Face &face) {
  return onEdge(P, face.q1, face.q2) || onEdge(P, face.q2, face.q3) ||
         onEdge(P, face.q3, face.q4) || onEdge(P, face.q4, face.q1);
}

#pragma acc routine seq
constexpr bool onFace(const Point3 &P, const Face &face) {
  return insideCorner(P, face.q1, face.q4, face.q3) &&
         insideCorner(P, face.q2, face.q1, face.q4) &&
         insideCorner(P, face.q3, face.q2, face.q1) &&
         insideCorner(P, face.q4, face.q3, face.q2);
}
//  ==
//  ||
//  ||  LineHitsFace: Returns 1 if the line specified by a point L0 and vector L
//  ||                intersect a face specified by four points Q0 - Q3.
//  ||
//  ==
#pragma acc routine seq
constexpr auto operator&(const Ray &ray, const Face &face) noexcept {
  const auto l_dot_n = dot(ray.n, face.n);
  if (l_dot_n > 1.e-14)
    return std::make_pair(false, Point3{});
  const auto d = dot(face.c - ray.p, face.n) / l_dot_n;
  const Point3 p = ray.p + d * ray.n;
  return std::make_pair(onFace(p, face), p);
  // return onVertex(p, face) || onEdge(p, face) || onFace(p, face);
}
