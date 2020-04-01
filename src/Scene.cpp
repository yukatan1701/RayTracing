#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <cmath>
#include <omp.h>
#include <set>

#include "bitmap_image.hpp"
#include "glm/gtx/intersect.hpp"
#include "includes.h"
#include "Material.h"
#include "SceneObjects.h"
#include "Scene.h"
#include "stb_image.h"

using namespace glm;

vec3 refract(const vec3 &I, const vec3 &N, const float etat, const float etai = 1.0f) {
    vec3 nI = normalize(I);
    float cosi = -std::max(-1.0f, std::min(1.0f, dot(I, N)));
    if (cosi < 0)
        return refract(I, -N, etai, etat);
    float eta = etai / etat;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? vec3(1, 0, 0) : I * eta + N * (eta * cosi - sqrtf(k));
}

bool Scene::sceneIntersect(const vec3 &orig, const vec3 &dir,
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
        float scale = 10.0f; // textures per 1 x 1
        float x = hit.x / scale, z = hit.z / scale;
        float xx = x - floor(x), zz = z - floor(z);
       // float x = (hit.x - floor(hit.x)), z = (hit.z - floor(hit.z));
        unsigned xpos = xx * waterWidth, zpos = zz * waterHeight;
        material.diffuse = (water[zpos * waterWidth + xpos] * 0.3f + material.diffuse) / 2.0f;
        //printf("(%.3f, %.3f, %.3f)", material.diffuse.x, material.diffuse.y, material.diffuse.z);
    }
    minDist = std::min(waterDist, minDist);

    sceneObjects.islandsIntersect(orig, dir, hit, N, material, minDist);
    return minDist < 700.0f;
}

Lightning Scene::calcLightning(const vec3 &orig, const vec3 &dir, const vec3 &point, const vec3 &N,
    const Material &mt, const SceneObjects &sceneObjects) {
    Lightning lightning;
    vec3 positiveOrig = point + N * fprec, negativeOrig = point - N * fprec;
    for (auto &light: sceneObjects.lights) {
        vec3 lightDir = normalize(light->position - point);
        vec3 shadowOrig = isNegative(dot(lightDir, N)) ? negativeOrig : positiveOrig;
        vec3 shadowPt(0.0f), shadowN(0.0f);
        Material empty;
        if (sceneIntersect(shadowOrig, lightDir, sceneObjects, shadowPt, shadowN, empty) &&
            distance(shadowPt, shadowOrig) < distance(light->position, point))
            continue;

        auto difIntDot = dot(lightDir, N);
        if (difIntDot >= 0)
            lightning.diffuse += light->intensity * difIntDot;
        auto specIntDot = dot(-reflect(-lightDir, N), dir);
        if (specIntDot >= 0)
            lightning.specular += powf(specIntDot, mt.specular) * light->intensity;
    }
    return lightning;
}

vec3 Scene::traceRay(const vec3 &orig, const vec3 &dir,
                     const SceneObjects &sceneObjects, size_t depth) {
    vec3 point(0.0f), N(0.0f);
    Material material;
    if (depth <= 0 || !sceneIntersect(orig, dir, sceneObjects, point, N, material)) {
        return sceneObjects.sky->getColorAt(dir);
    }
    
    Lightning light = calcLightning(orig, dir, point, N, material, sceneObjects);

    vec3 positiveOrig = point + N * fprec, negativeOrig = point - N * fprec;

    vec3 reflectDir = normalize(reflect(dir, N));
    vec3 reflectOrig = isNegative(dot(reflectDir, N)) ? negativeOrig : positiveOrig;
    vec3 reflectColor = traceRay(reflectOrig, reflectDir, sceneObjects, depth - 1);

    vec3 refractDir = normalize(refract(dir, N, material.refractive));
    vec3 refractOrig = isNegative(dot(refractDir, N)) ? negativeOrig : positiveOrig;
    vec3 refractColor = traceRay(refractOrig, refractDir, sceneObjects, depth - 1);
    
    return material.diffuse * light.diffuse * material.albedo[0] +
           vec3(1.0f) * light.specular * material.albedo[1] +
           reflectColor * material.albedo[2] + refractColor * material.albedo[3];
}

objset<Sphere *> Scene::loadSpheres() {
    objset<Sphere *> spheres;
    spheres.insert(new Sphere(vec3( -6, -2.8, -22), 1, materials["silver"]));
    spheres.insert(new Sphere(vec3( -4, -3.1, -21), 0.7, materials["glass"]));
    spheres.insert(new Sphere(vec3(-15.0f, -4.0f + 5, -50.0f), 5, materials["mirror"]));
    spheres.insert(new Sphere(vec3(25.0f, -4.0f + 4, -50.0f), 4, materials["glass"]));
    return spheres;
}

objset<Cube *> Scene::loadCubes() {
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

objset<Light *> Scene::loadLights() {
    objset<Light *> lights;
    lights.insert(new Light(vec3(-20, 20,  20), 1.5));
    lights.insert(new Light(vec3( 35, 50, -25), 2.0));
    lights.insert(new Light(vec3( 30, 20,  30), 1.7));
    return lights;
}

objset<Triangle *> Scene::loadTriangles() {
    objset<Triangle *> triangles;
    Material box = materials["wood"];
    triangles.insert(new Triangle(vec3(-8, -3.75, -20), vec3(-8, -3.75, -23), vec3(-8, -1, -23), box));
    triangles.insert(new Triangle(vec3(-8, -3.75, -20), vec3(-8, -1, -23), vec3(-8, -3.75, -23), box)); // back
    triangles.insert(new Triangle(vec3(-8, -3.75, -23), vec3(-2, -3.75, -23), vec3(-8, -1, -23), box));
    triangles.insert(new Triangle(vec3(-8, -3.75, -23), vec3(-8, -1, -23), vec3(-2, -3.75, -23), box)); //back
    return triangles;
}

objset<Model *> Scene::loadModels() {
    objset<Model *> models;
    vec3 bias(4.0f, -7.0f + 1.4f, -23.0f); //+1.4f
    Model *bunny = new Model("../resources/bunny.obj", 45.0f, bias, materials["gold"]);
    models.insert(bunny);
    return models;
}

objset<Island *> Scene::loadIslands() {
    objset<Island *> islands;
    islands.insert(new Island(vec3(0.0f, -4.0f, -25.0f), 12, 8, materials["island"]));
    islands.insert(new Island(vec3(-15.0f, -4.0f, -50.0f), 12, 15, materials["island"]));
    islands.insert(new Island(vec3(25.0f, -4.0f, -50.0f), 10, 12, materials["island"]));
    islands.insert(new Island(vec3(15.0f, -4.0f, -100.0f), 15, 12, materials["island"]));
    return islands;
}

void Scene::deleteSceneObjects(SceneObjects &so) {
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
    delete so.sky;
}

void Scene::render(const Settings &settings) {
    objset<Sphere *> spheres = loadSpheres();
    objset<Triangle *> triangles = loadTriangles();
    objset<Light *> lights = loadLights();
    objset<Model *> models = loadModels();
    objset<Cube *> cubes = loadCubes();
    objset<Island *> islands = loadIslands();
    
    SceneObjects sceneObjects(spheres, triangles, lights, models, cubes, islands);

    const float fov = M_PI / 3.0f;

	float angleX = 15 * M_PI / float(180);
    float angleY = 15 * M_PI / float(180);
	
    mat3 rotateX = rotate(mat4(1), -angleX, vec3(1, 0, 0));
    mat3 rotateY = rotate(mat4(1), -angleY, vec3(0, 1, 0));
    vec3 camera(-7.5, 8, -5); //(-7.5, 8, - ...)

    auto getRGBColor = [](vec3 point){
        rgb_t colorRGB;
        colorRGB.red = 255 * std::max(0.0f, std::min(1.0f, point.x));
        colorRGB.green = 255 * std::max(0.0f, std::min(1.0f, point.y));
        colorRGB.blue = 255 * std::max(0.0f, std::min(1.0f, point.z));
        return colorRGB;
    };

    bitmap_image image(width, height);
    double time0 = omp_get_wtime();
#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < height; ++i) {
        for (size_t j = 0; j < width; ++j) {
            float x =  (j + 0.5f) -  width / 2.0f;
            float y = -(i + 0.5f) + height / 2.0f;
            float z = -height / (2.0f * tan(fov / 2.0f));
            vec3 dir = normalize(rotateY * rotateX * vec3(x, y, z) + camera);
            vec3 result = traceRay(camera, dir, sceneObjects, 4);
            image.set_pixel(j, i, getRGBColor(result));
        }
    }
    double time1 = omp_get_wtime();
    std::cout << "Done. Render time: " << time1 - time0 << std::endl;
    image.save_image(settings.out);
    deleteSceneObjects(sceneObjects);
}

int Scene::run(const Settings &settings) {
    omp_set_num_threads(settings.threads);
    double time0 = omp_get_wtime();
    
    int n = -1;
    unsigned char *water_image = stbi_load("../resources/water.jpg", &waterWidth, &waterHeight, &n, 0);
    if (!water_image || n != 3) {
        std::cout << "Error: can not load water texture: " << stbi_failure_reason() <<std::endl;
        return -1;
    }
    water = std::vector<vec3>(waterWidth * waterHeight);
    for (int j = waterHeight - 1; j >= 0; j--) {
        for (int i = 0; i < waterWidth; i++) {
            unsigned pos = i + j * waterWidth;
            vec3 colorRGB(water_image[pos * 3], water_image[pos * 3 + 1], water_image[pos * 3 + 2]);
            water[pos] = colorRGB * 1.0f / 255.0f;
        }
    }
    stbi_image_free(water_image);

    std::cout << "[Current settings]" << std::endl;
    std::cout << "Output path: " << settings.out << std::endl;
    std::cout << "Scene number: " << settings.scene << std::endl;
    std::cout << "Threads: " << settings.threads << std::endl;
    std::cout << "Real threads: " << omp_get_max_threads() << std::endl;
    std::cout << "\nRendering image..." << std::endl;

    render(settings);
    double time1 = omp_get_wtime();
    std::cout << "Total time (including loading textures and models): " << time1 - time0 << std::endl;
    return 0;
}
