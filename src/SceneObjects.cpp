#include "SceneObjects.h"
#include "glm/detail/func_geometric.inl"
#include <vector>

using namespace glm;

bool Sphere::rayIntersect(const vec3 &orig, const vec3 &dir, float &dist) const {
    float Epsilon = std::numeric_limits<float>::epsilon();
    vec3 diff = center - orig;
    float t0 = dot(diff, dir);
    float dSquared = dot(diff, diff) - t0 * t0;
    float radius2 = radius * radius;
    if(dSquared > radius2)
        return false;
    float t1 = sqrt(radius2 - dSquared);
    dist = t0 > t1 + Epsilon ? t0 - t1 : t0 + t1;
    return dist > Epsilon;
}

bool Triangle::rayIntersect(const vec3 &orig, const vec3 &dir, float &dist) const {
    vec3 edge1 = v1 - v0;
    vec3 edge2 = v2 - v0;
    vec3 pvec = cross(dir, edge2);
    float det = dot(edge1, pvec);
    const float eps = 0.00001f;
    if (det < eps)
        return false;

    vec3 tvec = orig - v0;
    float u = dot(tvec, pvec);
    if (u < 0 || u > det)
        return false;

    vec3 qvec = cross(tvec, edge1);
    float v = dot(dir, qvec);
    if (v < 0 || u + v > det)
        return false;

    dist = dot(edge2, qvec) * (1.0f / det);
    return dist > eps;
}

std::vector<Triangle> Quadrangle::toTriangles() const {
    std::vector<Triangle> tr;
    tr.push_back(Triangle(v0, v1, v2, material));
    tr.push_back(Triangle(v0, v2, v3, material));
    return tr;
}

bool Quadrangle::rayIntersect(const vec3 &orig, const vec3 &dir, float &dist) const {
        Triangle t0(v0, v1, v2), t1(v0, v2, v3);
        float dist0, dist1;
        if (t0.rayIntersect(orig, dir, dist0)) {
            dist = dist0;
            return true;
        }
        if (t1.rayIntersect(orig, dir, dist1)) {
            dist = dist1;
            return true;
        }
        return false;
    }