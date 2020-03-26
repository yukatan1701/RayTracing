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

void Cube::loadFaces(const vec3 &minPoint, const vec3 &maxPoint) {
    const vec3 &minP = minPoint;
    const vec3 &maxP = maxPoint;
    float minX = std::min(minP.x, maxP.x), minY = std::min(minP.y, maxP.y),
          minZ = std::min(minP.z, maxP.z);
    float maxX = std::max(minP.x, maxP.x), maxY = std::max(minP.y, maxP.y),
          maxZ = std::max(minP.z, maxP.z);
    Cube::minPoint = vec3(minX, minY, minZ);
    Cube::maxPoint = vec3(maxX, maxY, maxZ);

    dx = std::pair<float, float>(minX, maxX);
    dy = std::pair<float, float>(minY, maxY);
    dz = std::pair<float, float>(minZ, maxZ);
    //printf("[%.3f][%.3f][%.3f]\n", dx.first, dy.first, dz.first);
}

Cube::Cube(const vec3 &leftBottom, const vec3 &rightTop,
           const Material &m = Material()) : material(m){
    loadFaces(leftBottom, rightTop);
}

bool Cube::rayIntersect(const vec3 &orig, const vec3 &dir, float &dist) const {
    // orig point inside the cude
    if (orig.x >= minPoint.x && orig.x <= maxPoint.x &&
        orig.y >= minPoint.y && orig.y <= maxPoint.y &&
        orig.z >= minPoint.z && orig.z <= maxPoint.z) {
            return true;
        }
    float t_near = std::numeric_limits<float>::min();
    float t_far = std::numeric_limits<float>::max();
    float t1, t2;
    float eps = std::numeric_limits<float>::epsilon();

    for (int i = 0; i < 3; ++i) {
        if (abs(dir[i]) >= eps) {
            t1 = (minPoint[i] - orig[i]) / dir[i];
            t2 = (maxPoint[i] - orig[i]) / dir[i];

            if (t1 > t2)
                std::swap(t1, t2);
            if (t1 > t_near)
                t_near = t1;
            if (t2 < t_far)
                t_far = t2;

            if (t_near > t_far)
                return false;
            if (t_far < std::numeric_limits<float>::epsilon())
                return false;
        } else {
            if (orig[i] < minPoint[i] || orig[i] > maxPoint[i])
                return false;
        }
    }
    if (!(t_near <= t_far && t_far >= eps))
        return false;
    dist = t_near;
    return true;
}

bool Cube::rayIntersect(const vec3 &orig, const vec3 &dir, float &dist, vec3 &n) const {
    // orig point inside the cude
    if (orig.x >= minPoint.x && orig.x <= maxPoint.x &&
        orig.y >= minPoint.y && orig.y <= maxPoint.y &&
        orig.z >= minPoint.z && orig.z <= maxPoint.z) {
            return true;
        }
    float t_near = std::numeric_limits<float>::min();
    float t_far = std::numeric_limits<float>::max();
    float t1, t2;
    float eps = std::numeric_limits<float>::epsilon();

    for (int i = 0; i < 3; ++i) {
        if (abs(dir[i]) >= eps) {
            t1 = (minPoint[i] - orig[i]) / dir[i];
            t2 = (maxPoint[i] - orig[i]) / dir[i];
            if (t1 > t2)
                std::swap(t1, t2);
            if (t1 > t_near)
                t_near = t1;
            if (t2 < t_far)
                t_far = t2;

            if (t_near > t_far)
                return false;
            if (t_far < std::numeric_limits<float>::epsilon())
                return false;
        } else {
            if (orig[i] < minPoint[i] || orig[i] > maxPoint[i])
                return false;
        }
    }
    if (!(t_near <= t_far && t_far >= eps))
        return false;
    dist = t_near;
    vec3 hit = orig + dir * dist;
    vec3 c = (minPoint + maxPoint) * 0.5f;
    vec3 d = (minPoint - maxPoint) * 0.5f;
    vec3 p = hit - c;
    float bias = 1.000001f;
    n = normalize(vec3(trunc(p.x / abs(d.x) * bias),
                       trunc(p.y / abs(d.y) * bias),
                       trunc(p.z / abs(d.z) * bias)));
    return true;
}

void BoundingBox::init() {
    float lx = dx.second - dx.first, ly = dy.second - dy.first, lz = dz.second - dz.first;
    std::cout << lx << " " << ly << " " << lz << std::endl;
    std::cout << size << std::endl;
    float px = lx / float(size), py = ly / float(size), pz = lz / float(size);
    std::cout << px << " " << py << " " << pz << std::endl;
    float x = dx.first, y = dy.first, z = dz.first;
    vec3 leftBottom, rightTop;
    for (int i = 0; i < size; i++) {
        leftBottom.x = dx.first + px * i;
        rightTop.x = dx.first + px * (i + 1);
        for (int j = 0; j < size; j++) {
            leftBottom.y = dy.first + py * j;
            rightTop.y = dy.first + py * (j + 1);
            for (int k = 0; k < size; k++) {
                leftBottom.z = dz.first + pz * k;
                rightTop.z = dz.first + pz * (k + 1);
                grid[i][j][k].loadFaces(leftBottom, rightTop);
            }
        }
    }
}

void BoundingBox::initTriangles(const std::unordered_set<const Triangle *> &trs) {
    float lx = (dx.second - dx.first);
    float ly = (dy.second - dy.first);
    float lz = (dz.second - dz.first);
    for (auto &t: trs) {
        for (const vec3 &v : t->getVerts()) {
            float fx = (v.x - dx.first) * size / lx;
            float fy = (v.y - dy.first) * size / ly;
            float fz = (v.z - dz.first) * size / lz;
            int i = int(fx), j = int(fy), k = int(fz);
            tgrid[i][j][k].insert(t);
        }
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
    box.init();
    //boundingCube.printCube();
    while (faceNum--) {
        int ix, iy, iz;
        char ch = 0;
        file >> ch; // f
        file >> ix >> iy >> iz;
        triangles.insert(new Triangle(vertices[ix-1], vertices[iy-1], vertices[iz-1], material));
    } 
    file.close();
    box.initTriangles(triangles);
}

bool Model::boxIntersect(const vec3 &orig, const vec3 &dir, float &dist, vec3 &n) const {
    return box.rayIntersect(orig, dir, dist, n);
}

bool Model::boxIntersect(const vec3 &orig, const vec3 &dir, float &dist) const {
    return box.rayIntersect(orig, dir, dist);
}

Model::~Model() {
    for (auto tr: triangles) {
       // delete tr;
    }
}