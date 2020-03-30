#ifndef __SCENE_OBJECTS__
#define __SCENE_OBJECTS__

#include "Material.h"
#include "includes.h"
#include "glm/gtx/intersect.hpp"
#include <deque>
#include <array>
#include <iostream>
#include <vector>

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
    void load(const vec3 &leftBottom, const vec3 &rightTop);
    bool rayIntersect(const vec3 &orig, const vec3 &dir, float &dist) const;
    bool rayIntersect(const vec3 &orig, const vec3 &dir, float &dist, vec3 &n) const;
};

struct BoundingBox : public Cube {
    static const int size = 7;
    Cube grid[size][size][size];
    objset<const Triangle *> tgrid[size][size][size];
    BoundingBox() : Cube() {}
    BoundingBox(const vec3 &leftBottom, const vec3 &rightTop) :
        Cube(leftBottom, rightTop, Material()) {}
    void init();
    void initTriangles(const objset<const Triangle *> &trs);
};

struct Model {
    BoundingBox box;
    objset<const Triangle *> triangles;
    float scale;
    Material material;
    Model(const std::string &filename, const float &scale,
          const vec3 &offset, const Material &m);
    ~Model();
    bool boxIntersect(const vec3 &orig, const vec3 &dir, float &dist, vec3 &n) const;
    bool boxIntersect(const vec3 &orig, const vec3 &dir, float &dist) const;
};

struct Island {
    vec3 center;
    float majorAxis, minorAxis;
    vec3 normal = vec3(0.0f, 1.0f, 0.0f);
    const vec3 lightGreen = vec3(136, 204, 0) * 0.3f / 255.0f;
    const vec3 darkGreen = vec3(102, 153, 0) * 0.3f / 255.0f;
    Material material;
    Island(const vec3 &center, const float &mjAx, const float &mnAx,
           const Material &material) :
        center(center), majorAxis(mjAx), minorAxis(mnAx),
        material(material) {}
    bool rayIntersect(const vec3 &orig, const vec3 &dir, float &dist);
};

struct SceneObjects {
    const objset<Sphere *> spheres;
    const objset<Triangle *> triangles;
    const objset<Light *> lights;
    const objset<Model *> models;
    const objset<Cube *> cubes;
    const objset<Island *> islands;
    SceneObjects(const objset<Sphere *> &s,
                 const objset<Triangle *> &t,
                 const objset<Light *> &l,
                 const objset<Model *> &m,
                 const objset<Cube *> &c,
                 const objset<Island *> &i) :
                 spheres(s), triangles(t), lights(l), models(m), cubes(c), islands(i) {}
    void spheresIntersect(const vec3 &orig, const vec3 &dir, vec3 &hit,
                          vec3 &N, Material &material, float &minDist) const;
    void trianglesIntersect(const vec3 &orig, const vec3 &dir, vec3 &hit, vec3 &N,
                            Material &material, float &minDist) const;
    void modelsIntersect(const vec3 &orig, const vec3 &dir, vec3 &hit, vec3 &N,
                         Material &material, float &minDist) const;
    void cubesIntersect(const vec3 &orig, const vec3 &dir, vec3 &hit, vec3 &N,
                        Material &material, float &minDist) const;
    void islandsIntersect(const vec3 &orig, const vec3 &dir, vec3 &hit, vec3 &N,
                          Material &material, float &minDist) const;
                 
};

#endif