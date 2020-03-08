#ifndef __SCENE_OBJECTS__
#define __SCENE_OBJECTS__

#include "glm/vec2.hpp"
#include "glm/vec3.hpp"
#include "Material.h"
#include <deque>
#include <vector>

using namespace glm;

struct Sphere {
    vec3 center;
    float radius;
    Material material;

    Sphere(const vec3 &center, const float &radius): center(center), radius(radius) {}
    Sphere(const vec3 &c, const float &r, const Material &m): center(c), radius(r), material(m) {}

    bool rayIntersect(const vec3 &orig, const vec3 &dir, float &dist) const;
};

struct Triangle {
    vec3 v0, v1, v2;
    Material material;

    Triangle(const vec3 &vert0, const vec3 &vert1, const vec3 &vert2) :
        v0(vert0), v1(vert1), v2(vert2) {} 
    Triangle(const vec3 &vert0, const vec3 &vert1, const vec3 &vert2, const Material &m) :
        v0(vert0), v1(vert1), v2(vert2), material(m) {}

    bool rayIntersect(const vec3 &orig, const vec3 &dir, float &dist) const;

};

struct Quadrangle {
    vec3 v0, v1, v2, v3;
    Material material;

    Quadrangle(const vec3 &vert0, const vec3 &vert1, const vec3 &vert2, const vec3 &vert3) :
        v0(vert0), v1(vert1), v2(vert2), v3(vert3) {}
    Quadrangle(const vec3 &vert0, const vec3 &vert1, const vec3 &vert2, const vec3 &vert3, const Material &m) :
        v0(vert0), v1(vert1), v2(vert2), v3(vert3), material(m) {}

    std::deque<Triangle> toTriangles() const;
    bool rayIntersect(const vec3 &orig, const vec3 &dir, float &dist) const;
};

struct Light {
    vec3 position;
    float intensity;
    Light(const vec3 &pos, const float &i) : position(pos), intensity(i) {}
};

struct Model {
    std::deque<Triangle> triangles;
    float scale;
    Material material;
    Model(const std::string &filename, const float &scale,
          const Material &m);
};

#endif