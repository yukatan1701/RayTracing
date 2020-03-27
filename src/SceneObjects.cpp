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
    const float eps = std::numeric_limits<float>::epsilon();
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

void Cube::load(const vec3 &minPoint, const vec3 &maxPoint) {
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
    load(leftBottom, rightTop);
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
    float t1 = 1.0f, t2 = 1.0f;
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
    if (!rayIntersect(orig, dir, dist)) {
        n = vec3(0.0f);
        return false;
    }
    vec3 hit = orig + dir * dist;
    vec3 c = (minPoint + maxPoint) * 0.5f;
    vec3 d = (minPoint - maxPoint) * 0.5f;
    vec3 p = hit - c;
    float bias = 1.000001f;
    float v1 = 0.0f, v2 = 0.0f, v3 = 0.0f;
    float eps = std::numeric_limits<float>::epsilon();
    v1 = p.x / abs(d.x) * bias;
    v2 = p.y / abs(d.y) * bias;
    v3 = p.z / abs(d.z) * bias;
    int iv1(v1), iv2(v2), iv3(v3);
    if (iv1 == 0 && iv2 == 0 && iv3 == 0) {
        iv1 = 1;
    }
    n = normalize(vec3(iv1, iv2, iv3));
    return true;
}

void BoundingBox::init() {
    float lx = dx.second - dx.first, ly = dy.second - dy.first, lz = dz.second - dz.first;
    float px = lx / float(size), py = ly / float(size), pz = lz / float(size);
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
                grid[i][j][k].load(leftBottom, rightTop);
            }
        }
    }
}

void BoundingBox::initTriangles(const objset<const Triangle *> &trs) {
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
    //std::cout << min_x << " " << min_y << " " << min_z << std::endl;
    //std::cout << max_x << " " << max_y << " " << max_z << std::endl;
    float eps = 0.0001f;
    vec3 veps(eps);
    box.load(vec3(min_x, min_y, min_z) - veps, vec3(max_x, max_y, max_z) + veps);
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
    //std::cout << triangles.size() << std::endl;
    for (auto it = triangles.begin(); it != triangles.end(); it++) {
        delete *it;
    }
}

bool Island::rayIntersect(const vec3 &orig, const vec3 &dir, float &dist) {
    float oldDist = dist;
    if (intersectRayPlane(orig, dir, vec3(0.0f, -4.0f, 0.0f), normal, dist)) {
        vec3 pt = orig + dir * dist;
        float xa = (pt.x - center.x) / majorAxis, zb = (pt.z - center.z) / minorAxis;

        if (pt.x > center.x - majorAxis && pt.x < center.x + majorAxis &&
            pt.z > center.z - minorAxis && pt.z < center.z + minorAxis &&
            xa * xa + zb * zb <= 1.0f) {
                return true; 
        }
    }
    dist = oldDist;
    return false;
}

void SceneObjects::spheresIntersect(const vec3 &orig, const vec3 &dir,
        vec3 &hit, vec3 &N, Material &material, float &minDist) const {
    float spheresDist = std::numeric_limits<float>::max();
    for (auto &sphere: spheres) {
        float curDist = 0.0f;
        if (sphere->rayIntersect(orig, dir, curDist) && curDist < spheresDist && curDist < minDist) {
            spheresDist = curDist;
            hit = orig + dir * curDist;
            N = normalize(hit - sphere->center);
            material = sphere->material;
        }
    }
    minDist = std::min(spheresDist, minDist);
}

void SceneObjects::trianglesIntersect(const vec3 &orig, const vec3 &dir,
        vec3 &hit, vec3 &N, Material &material, float &minDist) const {
    float trianglesDist = std::numeric_limits<float>::max();
    for (auto &triangle : triangles) {
        float curDist = 0.0f;
        if (triangle->rayIntersect(orig, dir, curDist) && curDist < trianglesDist && curDist < minDist) {
            trianglesDist = curDist;
            hit = orig + dir * curDist;
            N = triangleNormal(triangle->v0, triangle->v1, triangle->v2);
            material = triangle->material;
        }
    }
    minDist = std::min(trianglesDist, minDist);
}

void SceneObjects::modelsIntersect(const vec3 &orig, const vec3 &dir,
        vec3 &hit, vec3 &N, Material &material, float &minDist) const {
    float modelsDist = std::numeric_limits<float>::max();
    for (auto &model : models) {
        float boxDist = 0.0f;
        if (model->boxIntersect(orig, dir, boxDist)) {
            int size = model->box.size;
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size ; j++) {
                    for (int k = 0; k < size; k++) {
                        if (model->box.grid[i][j][k].rayIntersect(orig, dir, boxDist)) {
                            float curDist = 0.0f;
                            for (const Triangle *triangle : model->box.tgrid[i][j][k]) {
                            //for (const Triangle *triangle : model->triangles) {    
                                if (triangle->rayIntersect(orig, dir, curDist) && curDist < modelsDist &&
                                    curDist < minDist) {
                                    modelsDist = curDist;
                                    hit = orig + dir * curDist;
                                    N = triangleNormal(triangle->v0, triangle->v1, triangle->v2);
                                    material = triangle->material;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    minDist = std::min(modelsDist, minDist);
}

void SceneObjects::cubesIntersect(const vec3 &orig, const vec3 &dir,
        vec3 &hit, vec3 &N, Material &material, float &minDist) const {
    float cubesDist = std::numeric_limits<float>::max();
    for (auto &cube : cubes) {
        float curDist = 0.0f;
        vec3 normal(0.0f);
        if (cube->rayIntersect(orig, dir, curDist, normal) && curDist < cubesDist &&
            curDist < minDist) {
            cubesDist = curDist;
            hit = orig + dir * curDist;
            N = normal;
            material = cube->material;
        }
    }
    minDist = std::min(cubesDist, minDist);
}

void SceneObjects::islandsIntersect(const vec3 &orig, const vec3 &dir, vec3 &hit, vec3 &N,
        Material &material, float &minDist) const {
    float islandsDist = std::numeric_limits<float>::max();
    for (auto &island : islands) {
        float curDist = 0.0f;
        if (island->rayIntersect(orig, dir, curDist) && curDist < minDist) {
            islandsDist = curDist;
            hit = orig + dir * curDist;
            N = island->normal;
            material = island->material;
            material.diffuse = (int(hit.x + hit.z)) & 1 ? island->lightGreen : island->darkGreen;
            float x = hit.x, fx = roundf(hit.x);
            float z = ceilf(hit.x + hit.z) - fx - 0.5;
            float rad = 0.15;
            if ((x - fx) * (x - fx) + (hit.z - z) * (hit.z - z) < rad * rad) {
                material.diffuse = (int(hit.x + hit.z)) & 1 ? island->darkGreen : island->lightGreen;
            }
        }
    }
    minDist = std::min(islandsDist, minDist);
}