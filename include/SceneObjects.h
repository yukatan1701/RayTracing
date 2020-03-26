#ifndef __SCENE_OBJECTS__
#define __SCENE_OBJECTS__

#include "glm/vec2.hpp"
#include "glm/vec3.hpp"
#include "Material.h"
#include <deque>
#include <array>
#include <unordered_set>
#include <iostream>
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
        v0(vert0), v1(vert1), v2(vert2) {  } 
    Triangle(const vec3 &vert0, const vec3 &vert1, const vec3 &vert2, const Material &m) :
        v0(vert0), v1(vert1), v2(vert2), material(m) {  }
    
    std::array<vec3, 3> getVerts() const {
        std::array<vec3, 3> v;
        v[0] = v0, v[1] = v1, v[2]= v2;
        v[0] = v0, v[1] = v1, v[2]= v2;
        return v;
    }

    bool rayIntersect(const vec3 &orig, const vec3 &dir, float &dist) const;

};

struct Quadrangle {
    vec3 v0, v1, v2, v3;
    Material material;

    Quadrangle(const vec3 &vert0, const vec3 &vert1, const vec3 &vert2, const vec3 &vert3) :
        v0(vert0), v1(vert1), v2(vert2), v3(vert3) {}
    Quadrangle(const vec3 &vert0, const vec3 &vert1, const vec3 &vert2, const vec3 &vert3, const Material &m) :
        v0(vert0), v1(vert1), v2(vert2), v3(vert3), material(m) {}
    Quadrangle(const std::vector<vec3> &v, const Material &m) :
        v0(v[0]), v1(v[1]), v2(v[2]), v3(v[3]), material(m) {}

    std::deque<Triangle> toTriangles() const;
    bool rayIntersect(const vec3 &orig, const vec3 &dir, float &dist, vec3 &n) const;
    void print() const {
        printf("(%.3f, %.3f, %.3f) ", v0.x, v0.y, v0.z);
        printf("(%.3f, %.3f, %.3f) ", v1.x, v1.y, v1.z);
        printf("(%.3f, %.3f, %.3f) ", v2.x, v2.y, v2.z);
        printf("(%.3f, %.3f, %.3f)\n", v3.x, v3.y, v3.z);
    }
};

struct Light {
    vec3 position;
    float intensity;
    Light(const vec3 &pos, const float &i) : position(pos), intensity(i) {}
};

struct Cube {
    vec3 minPoint, maxPoint;
    std::pair<float, float> dx, dy, dz;
    Material material;
    Cube() {}
    Cube(const vec3 &leftBottom, const vec3 &rightTop, const Material &m);
    void loadFaces(const vec3 &leftBottom, const vec3 &rightTop);
    bool rayIntersect(const vec3 &orig, const vec3 &dir, float &dist) const;
    bool rayIntersect(const vec3 &orig, const vec3 &dir, float &dist, vec3 &n) const;
};

struct BoundingBox : public Cube {
    static const int size = 5;
    Cube grid[size][size][size];
    std::unordered_set<const Triangle *> tgrid[size][size][size];
    BoundingBox() : Cube() {}
    BoundingBox(const vec3 &leftBottom, const vec3 &rightTop) :
        Cube(leftBottom, rightTop, Material()) {}
    void init();
    void initTriangles(const std::unordered_set<const Triangle *> &trs);
};

struct Model {
    BoundingBox box;
    std::unordered_set<const Triangle *> triangles;
    float scale;
    Material material;
    Model(const std::string &filename, const float &scale,
          const vec3 &offset, const Material &m);
    ~Model();
    bool boxIntersect(const vec3 &orig, const vec3 &dir, float &dist, vec3 &n) const;
    bool boxIntersect(const vec3 &orig, const vec3 &dir, float &dist) const;
};

struct SceneObjects {
    const std::deque<Sphere> &spheres;
    const std::deque<Triangle> &triangles;
    const std::deque<Light> &lights;
    const std::deque<Model> &models;
    const std::deque<Cube> &cubes;
    SceneObjects(const std::deque<Sphere> &s, const std::deque<Triangle> &t,
                 const std::deque<Light> &l, const std::deque<Model> &m,
                 const std::deque<Cube> &c) :
                 spheres(s), triangles(t), lights(l), models(m), cubes(c) {}
                 
};

#endif