#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <cmath>
#include <omp.h>
#include <set>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "bitmap_image.hpp"
#include "includes.h"
#include "Material.h"
#include "SceneObjects.h"
#include "Scene1.h"

using namespace glm;

bool Scene1::sceneIntersect(const vec3 &orig, const vec3 &dir,
                            const SceneObjects &sceneObjects,
                            vec3 &hit, vec3 &N, Material &material) {
    float minDist = std::numeric_limits<float>::max();
    sceneObjects.spheresIntersect(orig, dir, hit, N, material, minDist);
    sceneObjects.trianglesIntersect(orig, dir, hit, N, material, minDist);
    sceneObjects.cubesIntersect(orig, dir, hit, N, material, minDist);
    sceneObjects.modelsIntersect(orig, dir, hit, N, material, minDist);
    float checkerboardDist = std::numeric_limits<float>::max();
    if (fabs(dir.y) > fprec) {
        float d = -(orig.y + 4) / dir.y;
        vec3 pt = orig + dir * d;
        if (d > 0 && fabs(pt.x) < 10 && pt.z < -10 && pt.z > -30 && d < minDist) {
            checkerboardDist = d;
            hit = pt;
            N = vec3(0.0f, 1.0f, 0.0f);
            material.diffuse = (int(0.5f * hit.x + 1000) + int(0.5 * hit.z)) & 1 ? 
                vec3(.3f, .3f, .3f) : vec3(.3f, .2f, .1f);
        }
    }
    minDist = std::min(checkerboardDist, minDist);
    return minDist < 800.0f;
}

vec3 Scene1::traceRay(const vec3 &orig, const vec3 &dir,
                      const SceneObjects &sceneObjects, size_t depth = 0) {
    vec3 point(0.0f), N(0.0f);
    Material material;
    if (depth > refractDepth || !sceneIntersect(orig, dir, sceneObjects, point, N, material)) {
        Sphere env(vec3(0.0f), 100, Material());
        float dist(0.0f);
        env.rayIntersect(orig, dir, dist);
        vec3 p = orig + dir * dist;
        unsigned a = (atan2(p.z, p.x) / (2.0f * M_PI) + 0.5f) * (float) backgrWidth;
        unsigned b = acos(p.y / 100.0f) / M_PI * (float) backgrHeight;
        unsigned place = a + b * backgrWidth;
        return background[place % background.size()];
    }
    vec3 reflectDir = normalize(reflect(dir, N));
    vec3 refractDir = normalize(refract(dir, N, material.refractive));

    vec3 positiveOrig = point + N * fprec, negativeOrig = point - N * fprec;

    vec3 reflectOrig = isNegative(dot(reflectDir, N)) ? negativeOrig : positiveOrig;
    vec3 refractOrig = isNegative(dot(refractDir, N)) ? negativeOrig : positiveOrig;

    vec3 reflectColor = traceRay(reflectOrig, reflectDir, sceneObjects, depth + 1);
    vec3 refractColor = traceRay(refractOrig, refractDir, sceneObjects, depth + 1);

    float diffuseIntensity(0.0f), specularIntensity(0.0f);
    for (auto &light: sceneObjects.lights) {
        vec3 lightDir = normalize(light->position - point);
        float lightDistance = distance(light->position, point);
        vec3 shadowOrig = isNegative(dot(lightDir, N)) ? negativeOrig : positiveOrig;
        vec3 shadowPt(0.0f), shadowN(0.0f);
        Material empty;
        if (sceneIntersect(shadowOrig, lightDir, sceneObjects, shadowPt, shadowN, empty) &&
            distance(shadowPt, shadowOrig) < lightDistance)
            continue;

        diffuseIntensity += light->intensity * std::max(0.0f, dot(lightDir, N));
        specularIntensity += powf(std::max(0.0f, dot(-reflect(-lightDir, N), dir)),
            material.specular) * light->intensity;
    }
    return material.diffuse * diffuseIntensity * material.albedo[0] +
           vec3(1.0f) * specularIntensity * material.albedo[1] +
           reflectColor * material.albedo[2] + refractColor * material.albedo[3];
}

objset<Sphere *> Scene1::loadSpheres() {
    objset<Sphere *> spheres;
    spheres.insert(new Sphere(vec3( 7, 5, -18), 4, materials["mirror"]));
    spheres.insert(new Sphere(vec3( 0, 7, -25), 4, materials["glass"]));
    spheres.insert(new Sphere(vec3( -7, 5, -20), 4, materials["red_rubber"]));
    return spheres;
}

objset<Light *> Scene1::loadLights() {
    objset<Light *> lights;
    lights.insert(new Light(vec3(-20, 20,  20), 1.5));
    lights.insert(new Light(vec3( 30, 50, -25), 2.0)); //bad
    lights.insert(new Light(vec3( 30, 20,  30), 1.7)); // ok
    return lights;
}

objset<Triangle *> Scene1::loadTriangles() {
    objset<Triangle *> triangles;
    Material glass = materials["pyramid_glass"];
   /* triangles.insert(new Triangle(vec3(8.0f, -4.f, -17.f), vec3(8, -4, -19), vec3(7, -1, -18), glass)); //right
    triangles.insert(new Triangle(vec3(6.f, -4.f, -19.f), vec3(6, -4, -17), vec3(7, -1, -18), glass)); //left
    triangles.insert(new Triangle(vec3(6, -4, -17), vec3(8, -4, -17), vec3(7, -1, -18), glass)); //front
    triangles.insert(new Triangle(vec3(8, -4, -19), vec3(6, -4, -19), vec3(7, -1, -18), glass)); //back*/
    //Quadrangle mirror(vec3(-10, -4, -15), vec3(-10, -4, -25), vec3(-10, 4, -25), vec3(-10, 4, -15), materials["mirror"]);
    //objset<Triangle> mirTr = mirror.toTriangles();
    //triangles.insert(triangles.begin(), mirTr.begin(), mirTr.end());
    return triangles;
}

objset<Model *> Scene1::loadModels() {
    objset<Model *> models;
    vec3 bias(2.0f, -7.0f + 1.65f, -23.0f); //+1.7f
    Model *bunny = new Model("../resources/bunny.obj", 40.0f, bias, materials["gold"]);
    models.insert(bunny);
    return models;
}

objset<Cube *> Scene1::loadCubes() {
    objset<Cube *> cubes;
    cubes.insert(new Cube(vec3(-6, -4, -20), vec3(-3, -1, -23), materials["silver"]));
    return cubes;
}

void Scene1::deleteSceneObjects(SceneObjects &so) {
    for (auto sphere : so.spheres)
        delete sphere;
    for (auto triangle : so.triangles)
        delete triangle;
    for (auto light : so.lights)
        delete light;
    for (auto model : so.models)
        delete model;
    for (auto cube : so.cubes)
        delete cube;
}

void Scene1::render(const Settings &settings) {
    objset<Sphere *> spheres = loadSpheres();
    objset<Triangle *> triangles = loadTriangles();
    objset<Light *> lights = loadLights();
    objset<Model *> models = loadModels();
    objset<Cube *> cubes = loadCubes();
    
    SceneObjects sceneObjects(spheres, triangles, lights, models, cubes);

    const float fov = M_PI / 3.0f;
    std::vector<vec3> imageGrid(width * height, vec3(0));

    vec3 camera(0, 5, 0);
#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < height; ++i) {
        for (size_t j = 0; j < width; ++j) {
            float x =  (j + 0.5f) -  width / 2.0f;
            float y = -(i + 0.5f) + height / 2.0f;
            float z = -height / (2.0f * tan(fov / 2.0f));
            vec3 dir = normalize(vec3(x, y, z));
            vec3 result = traceRay(camera, dir, sceneObjects);
            imageGrid[i * width + j] = result;
        }
    }

    deleteSceneObjects(sceneObjects);

    bitmap_image image(width, height);
    for (size_t y = 0; y < image.height(); ++y) {
        for (size_t x = 0; x < image.width(); ++x) {
            vec3 &c = imageGrid[y * width + x];
            float max = std::max(c[0], std::max(c[1], c[2]));
            if (max > 1) {
                c *= 1.0 / max;
            }
            char red = 255 * std::max(0.0f, std::min(1.0f, imageGrid[y * width + x][0]));
            char green = 255 * std::max(0.0f, std::min(1.0f, imageGrid[y * width + x][1]));
            char blue = 255 * std::max(0.0f, std::min(1.0f, imageGrid[y * width + x][2]));
            rgb_t color(red, green, blue);
            image.set_pixel(x, y, color);
        }
    }

    image.save_image(settings.out);
}

int Scene1::run(const Settings &settings) {
    omp_set_num_threads(settings.threads);
    int n = -1;
    unsigned char *pixmap = stbi_load("../resources/map3.jpg", &backgrWidth, &backgrHeight, &n, 0);
    if (!pixmap || n != 3) {
        std::cout << "Error: can not load the environment map" << std::endl;
        return -1;
    }
    background = std::vector<vec3>(backgrWidth * backgrHeight);
    float colorAlpha = 1 / 255.0f;
    for (int j = backgrHeight - 1; j >= 0; j--) {
        for (int i = 0; i < backgrWidth; i++) {
            unsigned pos = i + j * backgrWidth;
            vec3 colorRGB(pixmap[pos * 3], pixmap[pos * 3 + 1], pixmap[pos * 3 + 2]);
            background[pos] = colorRGB * colorAlpha;
        }
    }
    stbi_image_free(pixmap);
    std::cout << "[Current settings]" << std::endl;
    std::cout << "Output path: " << settings.out << std::endl;
    std::cout << "Scene number: " << settings.scene << std::endl;
    std::cout << "Threads: " << settings.threads << std::endl;
    std::cout << "Real threads: " << omp_get_max_threads() << std::endl;
    std::cout << "\nRendering image..." << std::endl;
    double time0 = omp_get_wtime();
    render(settings);
    std::cout << "Ok." << std::endl;
    double time1 = omp_get_wtime();
    std::cout << "Done. Elapsed time: " << time1 - time0 << std::endl;
    return 0;
}