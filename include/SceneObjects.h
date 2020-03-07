#ifndef __SCENE_OBJECTS__
#define __SCENE_OBJECTS__

#include "glm/vec2.hpp"
#include "glm/vec3.hpp"
#include "glm/gtx/intersect.hpp"
#include "Material.h"

using namespace glm;

struct Sphere {
    vec3 center;
    float radius;
    Material material;

    Sphere(const vec3 &center, const float &radius): center(center), radius(radius) {}
    Sphere(const vec3 &c, const float &r, const Material &m): center(c), radius(r), material(m) {}

    bool rayIntersect(const vec3 &orig, const vec3 &dir, float &dist) const {
        return intersectRaySphere(orig, dir, center, radius, dist);
    }
};

struct Triangle {
    vec3 v0, v1, v2;
    Material material;
    Triangle(const vec3 &vert0, const vec3 &vert1, const vec3 &vert2) :
        v0(vert0), v1(vert1), v2(vert2) {} 
    Triangle(const vec3 &vert0, const vec3 &vert1, const vec3 &vert2, const Material &m) :
        v0(vert0), v1(vert1), v2(vert2), material(m) {}

    bool rayIntersect(const vec3 &orig, const vec3 &dir, float &dist) const {
        vec2 baryPosition;
        return intersectRayTriangle(orig, dir, v0, v1, v2, baryPosition, dist);
    }
};

struct Light {
    vec3 position;
    float intensity;
    Light(const vec3 &pos, const float &i) : position(pos), intensity(i) {}
};

#endif