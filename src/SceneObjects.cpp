#include "SceneObjects.h"
#include "glm/gtx/normal.hpp"
#include <iostream>
#include <fstream>
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


std::deque<Triangle> Quadrangle::toTriangles() const {
    std::deque<Triangle> tr;
    tr.push_front(Triangle(v0, v1, v2, material));
    tr.push_front(Triangle(v0, v2, v3, material));
    return tr;
}

bool Quadrangle::rayIntersect(const vec3 &orig, const vec3 &dir, float &dist, vec3 &n) const {
    Triangle t0(v0, v1, v2), t1(v0, v2, v3);
    float dist0, dist1;
    if (t0.rayIntersect(orig, dir, dist0)) {
        n = triangleNormal(t0.v0, t0.v1, t0.v2);
        dist = dist0;
        return true;
    }
    if (t1.rayIntersect(orig, dir, dist1)) {
        n = triangleNormal(t1.v0, t1.v1, t1.v2);
        dist = dist1;
        return true;
    }
    return false;
}

Cube::Cube(const vec3 &leftBottom, const vec3 &rightTop,
           const Material &m = Material()) : material(m){
    const vec3 &lb = leftBottom;
    const vec3 &rt = rightTop;
    float minX = std::min(lb.x, rt.x), minY = std::min(lb.y, rt.y),
          minZ = std::min(lb.z, rt.z);
    float maxX = std::max(lb.x, rt.x), maxY = std::max(lb.y, rt.y),
          maxZ = std::max(lb.z, rt.z);
    std::vector<vec3> &bv = bottomVerts;
    std::vector<vec3> &tv = topVerts;
    bv.push_back(vec3(minX, minY, minZ));
    bv.push_back(vec3(minX, maxY, minZ));
    bv.push_back(vec3(maxX, maxY, minZ));
    bv.push_back(vec3(maxX, minY, minZ));
    tv.push_back(vec3(minX, minY, maxZ));
    tv.push_back(vec3(maxX, minY, maxZ));
    tv.push_back(vec3(maxX, maxY, maxZ));
    tv.push_back(vec3(minX, maxY, maxZ));
    
    std::vector<vec3> f1 { bv[0], bv[3], tv[1], tv[0] };
    std::vector<vec3> f2 { bv[3], bv[2], tv[2], tv[1] };
    std::vector<vec3> f3 { bv[2], bv[1], tv[3], tv[2] };
    std::vector<vec3> f4 { bv[1], bv[0], tv[0], tv[3] };
    
    faces.push_front(Quadrangle(bottomVerts, material));
    faces.push_front(Quadrangle(topVerts, material));
    faces.push_front(Quadrangle(f1, material));
    faces.push_front(Quadrangle(f2, material));
    faces.push_front(Quadrangle(f3, material));
    faces.push_front(Quadrangle(f4, material));
}

bool Cube::rayIntersect(const vec3 &orig, const vec3 &dir, float &dist, vec3 &n) const {
    for (auto f: faces) {
        if (f.rayIntersect(orig, dir, dist, n)) {
            return true;
        }
    }
    return false;
}

Model::Model(const std::string &filename, const float &scale, const vec3 &offset,
             const Material &m = Material()) :
             scale(scale), material(m) {
    std::ifstream file(filename, std::ios::in);
    if (!file.is_open()) {
        std::cerr << "Failed to open '" << filename << "' file." << std::endl;
    }
    int vertNum = 0, faceNum = 0;
    file >> vertNum >> faceNum;
    std::vector<vec3> vertices;
    int n = 10;
    float min_x, min_y, min_z, max_x, max_y, max_z;
    min_x = min_y = min_z = 1000;
    max_x = max_y = max_z = -1000;
    while (vertNum--) {
        double x = 0.0, y = 0.0, z = 0.0;
        char ch = 0;
        file >> ch; // v
        file >> x >> y >> z;
        vec3 v = vec3((float) x, (float) y, (float) z) * scale + offset;
        min_x = std::min(min_x, v.x);
        min_y = std::min(min_y, v.y);
        min_z = std::min(min_z, v.z);
        max_x = std::max(max_x, v.x);
        max_y = std::max(max_y, v.y);
        max_z = std::max(max_z, v.z);
        vertices.push_back(v);
    }
    std::cout << min_x << " " << min_y << " " << min_z << std::endl;
    std::cout << max_x << " " << max_y << " " << max_z << std::endl;
    boundingCube = Cube(vec3(min_x, min_y, min_z), vec3(max_x, max_y, max_z));
    //boundingCube.printCube();
    while (faceNum--) {
        int ix, iy, iz;
        char ch = 0;
        file >> ch; // f
        file >> ix >> iy >> iz;
        Triangle tr(vertices[ix-1], vertices[iy-1], vertices[iz-1], material);
        triangles.push_front(tr);
    } 
    file.close();
}

bool Model::cubeIntersect(const vec3 &orig, const vec3 &dir, float &dist, vec3 &n) const {
    return boundingCube.rayIntersect(orig, dir, dist, n);
}