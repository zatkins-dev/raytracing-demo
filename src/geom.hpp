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
#include <optional>

struct Vec3 {
  double x = 0., y = 0., z = 0.;

  template <class T> constexpr const Vec3 operator+(const T &w) const {
    return {.x = x + w.x, .y = y + w.y, .z = z + w.z};
  }

  constexpr const Vec3 operator*(double a) const {
    return {.x = a * x, .y = a * y, .z = a * z};
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
    return {.x = x + w.x, .y = y + w.y, .z = z + w.z};
  }

  constexpr const Point3 operator*(double a) const {
    return {.x = a * x, .y = a * y, .z = a * z};
  }

  constexpr const Point3 operator-() const { return *this * -1; }

  constexpr const Vec3 operator-(const Point3 &v) const {
    return {.x = x - v.x, .y = y - v.y, .z = z - v.z};
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
      .x = v.y * w.z - v.z * w.y,
      .y = v.z * w.x - v.x * w.z,
      .z = v.x * w.y - v.y * w.x,
  };
}

constexpr double dot(const Vec3 &v, const Vec3 &w) {
  return v.x * w.x + v.y * w.y + v.z * w.z;
}

constexpr double normsq(const Vec3 &v) { return dot(v, v); }

constexpr Vec3 normalized(const Vec3 &v) { return (1. / sqrt(normsq(v))) * v; }

struct Face {
  std::array<Point3, 4> verts;
  Point3 c;
  Vec3 n;
  const Point3 &q1 = verts[0], &q2 = verts[1], &q3 = verts[2], &q4 = verts[3];

  constexpr Face(std::array<Point3, 4> verts_, Vec3 n_) : verts(verts_), n(n_) {
    for (auto &q : verts) {
      c += q;
    }
    c *= 0.25;
  }
};

struct Plane {
  Point3 p;
  Vec3 n;
};

struct Ray {
  Point3 p;
  Vec3 n;

  constexpr explicit Ray(const Point3 &p_, const Vec3 &n_)
      : p(p_), n(normalized(n_)) {}

  constexpr Ray(const Point3 &p1, const Point3 &p2)
      : p(p1), n(normalized(p2 - p1)) {}
};

constexpr std::optional<Point3> operator&(const Ray &ray, const Plane &plane) {
  auto l_dot_n = dot(ray.n, plane.n);
  if (fabs(l_dot_n) < 1.e-10)
    return {};
  auto d = dot(plane.p - ray.p, plane.n) / l_dot_n;
  return ray.p + d * ray.n;
}

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
constexpr bool insideCorner(const Point3 &P, const Point3 &C, const Point3 &A,
                            const Point3 &B) {
  const Vec3 vAP = P - A;
  return dot(cross(B - A, vAP), cross(A - C, vAP)) >= 0;
}

//  ==
//  ||
//  ||  insideQuad: Given a quadrilateral as points in x,y,z arrays
//  ||              determine if point P is inside the quadrilateral.
//  ||
//  ||              It is inside the quadrilateral if it is inside
//  ||              each of its four corners.
//  ||
//  ==
constexpr std::optional<Point3> operator&(const Point3 &p, const Face &face) {
  const bool inside = insideCorner(p, face.q1, face.q4, face.q3) &&
                      insideCorner(p, face.q2, face.q1, face.q4) &&
                      insideCorner(p, face.q3, face.q2, face.q1) &&
                      insideCorner(p, face.q4, face.q3, face.q2);
  return inside ? std::optional(p) : std::nullopt;
}

//  ==
//  ||
//  ||  LineHitsFace: Returns 1 if the line specified by a point L0 and vector L
//  ||                intersect a face specified by four points Q0 - Q3.
//  ||
//  ==
constexpr std::optional<Point3> operator&(const Ray &line, const Face &face) {
  auto intersects = line & Plane{face.c, face.n};
  if (!intersects)
    return std::nullopt;
  return intersects.value() & face;
}
