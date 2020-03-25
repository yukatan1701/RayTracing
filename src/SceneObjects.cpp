#include "SceneObjects.h"
#include "glm/gtx/normal.hpp"
#include <iostream>
#include <fstream>
#include <vector>

using namespace glm;

void pvec(const vec3 &vec) {
    printf("(%.3f, %.3f, %.3f)", vec.x, vec.y, vec.z);
}

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

void Cube::loadFaces(const vec3 &leftBottom, const vec3 &rightTop) {
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

    dx = std::pair<float, float>(minX, maxX);
    dy = std::pair<float, float>(minY, maxY);
    dz = std::pair<float, float>(minZ, maxZ);
    //printf("[%.3f][%.3f][%.3f]\n", dx.first, dy.first, dz.first);

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

Cube::Cube(const vec3 &leftBottom, const vec3 &rightTop,
           const Material &m = Material()) : material(m){
    loadFaces(leftBottom, rightTop);
}

bool Cube::rayIntersect(const vec3 &orig, const vec3 &dir, float &dist, vec3 &n) const {
    for (auto f: faces) {
        if (f.rayIntersect(orig, dir, dist, n)) {
            return true;
        }
    }
    return false;
}

bool Cube::isInCube(const vec3 &v) const {

}

void BoundingBox::init() {
    float lx = dx.second - dx.first, ly = dy.second - dy.first, lz = dz.second - dz.first;
    float px = lx / float(size), py = ly / float(size), pz = lz / float(size);
    float x = dx.first, y = dy.first, z = dz.first;
    vec3 leftBottom, rightTop;
    for (int i = 0; i < size; i++) {
        leftBottom.x = x, rightTop.x = x + px;
        for (int j = 0; j < size; j++) {
            leftBottom.y = y, rightTop.y = y + py;
            for (int k = 0; k < size; k++) {
                leftBottom.z = z, rightTop.z = z + pz;
                grid[i][j][k].loadFaces(leftBottom, rightTop);
                z += pz;
            }
            y += py;
        }
        x += px;
    }
}

void BoundingBox::initTriangles(const std::deque<Triangle> &trs) {
    float lx = (dx.second - dx.first) * size;
    float ly = (dy.second - dy.first) * size;
    float lz = (dz.second - dz.first) * size;
    for (auto &t: trs) {
        const vec3 &v0 = t.v0, &v1 = t.v1, &v2 = t.v2;
        float fx0 = (v0.x - dx.first) / lx;
        float fy0 = (v0.y - dy.first) / ly;
        float fz0 = (v0.z - dz.first) / lz;
        float fx1 = (v1.x - dx.first) / lx;
        float fy1 = (v1.y - dy.first) / ly;
        float fz1 = (v1.z - dz.first) / lz;
        float fx2 = (v2.x - dx.first) / lx;
        float fy2 = (v2.y - dy.first) / ly;
        float fz2 = (v2.z - dz.first) / lz;
        int i0 = int(fx0), j0 = int(fy0), k0 = int(fz0);
        int i1 = int(fx1), j1 = int(fy1), k1 = int(fz1);
        int i2 = int(fx2), j2 = int(fy2), k2 = int(fz2);
        //if (i > size-1 || j > size-1 || k > size-1)
        //    printf("%d, %d, %d\n", i, j, k);
        if (*tgrid[i0][j0][k0].begin() != &t)
            tgrid[i0][j0][k0].push_front(&t);
        if (*tgrid[i1][j1][k1].begin() != &t)
            tgrid[i1][j1][k1].push_front(&t);
        if (*tgrid[i2][j2][k2].begin() != &t)
            tgrid[i2][j2][k2].push_front(&t);
    }
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
    float eps = 0.0001f;
    vec3 veps(eps);
    box.loadFaces(vec3(min_x, min_y, min_z) - veps, vec3(max_x, max_y, max_z) + veps);
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
    box.init();
    box.initTriangles(triangles);
}

bool Model::boxIntersect(const vec3 &orig, const vec3 &dir, float &dist, vec3 &n) const {
    return box.rayIntersect(orig, dir, dist, n);
}