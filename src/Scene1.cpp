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
#include "glm/gtx/intersect.hpp"
#include "includes.h"
#include "Material.h"
#include "SceneObjects.h"
#include "Scene1.h"

using namespace glm;

vec3 refract(const vec3 &I, const vec3 &N, const float etat, const float etai = 1.0f) {
    float cosi = -std::max(-1.0f, std::min(1.0f, dot(I, N)));
    if (cosi < 0)
        return refract(I, -N, etai, etat);
    float eta = etai / etat;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? vec3(1, 0, 0) : I * eta + N * (eta * cosi - sqrtf(k));
}

bool Scene1::sceneIntersect(const vec3 &orig, const vec3 &dir,
                            const SceneObjects &sceneObjects,
                            vec3 &hit, vec3 &N, Material &material) {
    float minDist = std::numeric_limits<float>::max();
    sceneObjects.spheresIntersect(orig, dir, hit, N, material, minDist);
    sceneObjects.trianglesIntersect(orig, dir, hit, N, material, minDist);
    sceneObjects.cubesIntersect(orig, dir, hit, N, material, minDist);
    sceneObjects.modelsIntersect(orig, dir, hit, N, material, minDist);

    vec3 waterOrig(0.0f, -4.5f, 0.0f);
    vec3 waterNorm(0.0f, 1.0f, 0.0f);
    float waterDist = std::numeric_limits<float>::max();
    if (intersectRayPlane(orig, dir, waterOrig, waterNorm, waterDist) && waterDist < minDist) {
        minDist = waterDist;
        hit = orig + dir * waterDist;
        N = waterNorm;
        material = materials["water_bottom"];
    }
    minDist = std::min(waterDist, minDist);

    sceneObjects.islandsIntersect(orig, dir, hit, N, material, minDist);
    return minDist < 700.0f;
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
        unsigned a = (atan2(p.z, p.x) / (2.0f * M_PI) + 0.5f) * (float) backgroundWidth;
        unsigned b = acos(p.y / 100.0f) / M_PI * (float) backgroundHeight;
        unsigned place = a + b * backgroundWidth + 200;
        return background[place % background.size()];
    }
    
    vec3 positiveOrig = point + N * fprec, negativeOrig = point - N * fprec;

    vec3 reflectDir = normalize(reflect(dir, N));
    vec3 reflectOrig = isNegative(dot(reflectDir, N)) ? negativeOrig : positiveOrig;
    vec3 reflectColor = traceRay(reflectOrig, reflectDir, sceneObjects, depth + 1);

    vec3 refractDir = normalize(refract(dir, N, material.refractive));
    vec3 refractOrig = isNegative(dot(refractDir, N)) ? negativeOrig : positiveOrig;
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
    spheres.insert(new Sphere(vec3( -6, -2.8, -22), 1, materials["silver"]));
    spheres.insert(new Sphere(vec3( -4, -3.1, -21), 0.7, materials["glass"]));
    spheres.insert(new Sphere(vec3(-15.0f, -4.0f + 5, -50.0f), 5, materials["mirror"]));
    spheres.insert(new Sphere(vec3(28.0f, -4.0f + 5, -75.0f), 5, materials["glass"]));
    return spheres;
}

objset<Cube *> Scene1::loadCubes() {
    objset<Cube *> cubes; 
    Material wood = materials["wood"];
    /*  */
    cubes.insert(new Cube(vec3(-6, -4, -28), vec3(-5, 3, -29), wood)); //center
    cubes.insert(new Cube(vec3(-3, -4, -28), vec3(-2, 3, -29), wood)); //left
    cubes.insert(new Cube(vec3(0, -4, -28), vec3(1, 3, -29), wood)); //right
    cubes.insert(new Cube(vec3(3, -4, -28), vec3(4, 3, -29), wood)); 
    cubes.insert(new Cube(vec3(3, -4, -28), vec3(4, 3, -29), wood)); 
    cubes.insert(new Cube(vec3(-7, 1.5, -28), vec3(5, 2.5, -29), wood));
    cubes.insert(new Cube(vec3(-7, -1.5, -28), vec3(5, -0.5, -29), wood));
    /* */ 
    cubes.insert(new Cube(vec3(-8, -4, -20), vec3(-2, -3.75, -23), materials["wood"]));
    return cubes;
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
    Material box = materials["wood"];
    triangles.insert(new Triangle(vec3(-8, -3.75, -20), vec3(-8, -3.75, -23), vec3(-8, -1, -23), box));
    triangles.insert(new Triangle(vec3(-8, -3.75, -23), vec3(-2, -3.75, -23), vec3(-8, -1, -23), box));
    return triangles;
}

objset<Model *> Scene1::loadModels() {
    objset<Model *> models;
    vec3 bias(4.0f, -7.0f + 1.4f, -23.0f); //+1.7f
    Model *bunny = new Model("../resources/bunny.obj", 45.0f, bias, materials["gold"]);
    models.insert(bunny);
    return models;
}

objset<Island *> Scene1::loadIslands() {
    objset<Island *> islands;
    islands.insert(new Island(vec3(0.0f, -4.0f, -25.0f), 12, 8, materials["island"]));
    islands.insert(new Island(vec3(-15.0f, -4.0f, -50.0f), 12, 15, materials["island"]));
    islands.insert(new Island(vec3(25.0f, -4.0f, -80.0f), 15, 15, materials["island"]));
    islands.insert(new Island(vec3(-45.0f, -4.0f, -100.0f), 12, 10, materials["island"]));
    return islands;
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
    for (auto island: so.islands) {
        delete island;
    }
}

void Scene1::render(const Settings &settings) {
    objset<Sphere *> spheres = loadSpheres();
    objset<Triangle *> triangles = loadTriangles();
    objset<Light *> lights = loadLights();
    objset<Model *> models = loadModels();
    objset<Cube *> cubes = loadCubes();
    objset<Island *> islands = loadIslands();
    
    SceneObjects sceneObjects(spheres, triangles, lights, models, cubes, islands);

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

    auto getRGBColor = [](vec3 point){
        vec3 colorRGB;
        colorRGB.x = 255 * std::max(0.0f, std::min(1.0f, point.x)); //red
        colorRGB.y = 255 * std::max(0.0f, std::min(1.0f, point.y)); //green
        colorRGB.z = 255 * std::max(0.0f, std::min(1.0f, point.z)); //blue
        return colorRGB;
    };

    bitmap_image image(width, height);
    for (size_t y = 0; y < image.height(); ++y) {
        for (size_t x = 0; x < image.width(); ++x) {
            vec3 &c = imageGrid[y * width + x];
            float max = std::max(c[0], std::max(c[1], c[2]));
            if (max > 1) {
                c *= 1.0 / max;
            }
            vec3 colorRGB = getRGBColor(imageGrid[y * width + x]);
            rgb_t color(colorRGB.x, colorRGB.y, colorRGB.z);
            image.set_pixel(x, y, color);
        }
    }

    image.save_image(settings.out);
}

int Scene1::run(const Settings &settings) {
    omp_set_num_threads(settings.threads);
    int n = -1;
    unsigned char *pixmap = stbi_load("../resources/map2.jpg", &backgroundWidth, &backgroundHeight, &n, 0);
    if (!pixmap || n != 3) {
        std::cout << "Error: can not load the environment map" << std::endl;
        return -1;
    }
    background = std::vector<vec3>(backgroundWidth * backgroundHeight);
    float colorAlpha = 1 / 255.0f;
    for (int j = backgroundHeight - 1; j >= 0; j--) {
        for (int i = 0; i < backgroundWidth; i++) {
            unsigned pos = i + j * backgroundWidth;
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
    double time1 = omp_get_wtime();
    std::cout << "Done. Elapsed time: " << time1 - time0 << std::endl;
    return 0;
}