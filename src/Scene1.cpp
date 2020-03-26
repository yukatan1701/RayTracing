#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <cmath>
//#include <omp.h>
#include <set>
#ifndef STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#endif
#include "stb_image.h"
#include "bitmap_image.hpp"
#include "glm/vec2.hpp"
#include "glm/vec3.hpp"
#include "glm/vec4.hpp"
#include "glm/gtx/norm.hpp"
#include "glm/gtx/normal.hpp"
#include "Material.h"
#include "SceneObjects.h"
#include "Scene1.h"

const int width = 4;
const int height = 4;

double bench_t_start, bench_t_end;

int envmap_width, envmap_height;
std::vector<vec3> envmap;

using namespace glm;

bool Scene1::sceneIntersect(const vec3 &orig, const vec3 &dir, const SceneObjects &sceneObjects,
                    vec3 &hit, vec3 &N, Material &material) {
    std::cout << "Begin scene" << std::endl;
    float min_dist = std::numeric_limits<float>::max();
    float spheresDist = std::numeric_limits<float>::max();
    for (auto &sphere: sceneObjects.spheres) {
        float curDist;
        if (sphere.rayIntersect(orig, dir, curDist) && curDist < spheresDist) {
            spheresDist = curDist;
            hit = orig + dir * curDist;
            N = normalize(hit - sphere.center);
            material = sphere.material;
        }
    }
    min_dist = std::min(spheresDist, min_dist);
    float trianglesDist = std::numeric_limits<float>::max();
    for (auto &triangle : sceneObjects.triangles) {
        float curDist;
        if (triangle.rayIntersect(orig, dir, curDist) && curDist < trianglesDist && curDist < min_dist) {
            trianglesDist = curDist;
            hit = orig + dir * curDist;
            N = triangleNormal(triangle.v0, triangle.v1, triangle.v2);
            material = triangle.material;
        }
    }
    min_dist = std::min(trianglesDist, min_dist);
    float cubesDist = std::numeric_limits<float>::max();
    for (auto &cube : sceneObjects.cubes) {
        float curDist;
        vec3 normal = vec3(0.0f);
        if (cube.rayIntersect(orig, dir, curDist, normal) && curDist < cubesDist &&
            curDist < min_dist) {
            cubesDist = curDist;
            hit = orig + dir * curDist;
            N = normal;
            material = cube.material;
        }
    }
    min_dist = std::min(cubesDist, min_dist);
    std::cout << "Begin models" << std::endl;
    float modelsDist = 10000.0f;
    for (auto &model : sceneObjects.models) {
        float boxDist;
        if (model.boxIntersect(orig, dir, boxDist)) {
            int size = model.box.size;
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size ; j++) {
                    for (int k = 0; k < size; k++) {
                        if (model.box.grid[i][j][k].rayIntersect(orig, dir, boxDist)) {
                            float curDist = 0.0f;
                            for (const Triangle *triangle : model.box.tgrid[i][j][k]) {
                                if (triangle != nullptr && triangle->rayIntersect(orig, dir, curDist) && curDist < modelsDist &&
                                    curDist < min_dist) {
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
    std::cout << "End models" << std::endl;
    min_dist = std::min(modelsDist, min_dist);
    float checkerboardDist = std::numeric_limits<float>::max();
    if (fabs(dir.y) > 0.001) {
        float d = -(orig.y + 4) / dir.y;
        vec3 pt = orig + dir * d;
        if (d > 0 && fabs(pt.x) < 10 && pt.z < -10 && pt.z > -30 && d < min_dist) {
            checkerboardDist = d;
            hit = pt;
            N = vec3(0.0f, 1.0f, 0.0f);
            material.diffuse = (int(0.5f * hit.x + 1000) + int(0.5 * hit.z)) & 1 ? 
                vec3(.3f, .3f, .3f) : vec3(.3f, .2f, .1f);
        }
    }
    min_dist = std::min(checkerboardDist, min_dist);
    std::cout << "End scene" << std::endl;
    return min_dist < 1000;
}

vec3 Scene1::traceRay(const vec3 &orig, const vec3 &dir, const SceneObjects &sceneObjects,
              size_t depth = 0) {
    vec3 point, N;
    Material material;
    std::cout << "Begin trace" << std::endl;
    if (depth > 4 || !sceneIntersect(orig, dir, sceneObjects, point, N, material)) {
        std::cout << "Begin inside" << std::endl;
        Sphere env(vec3(0.0f), 100, Material());
        float dist = 0;
        env.rayIntersect(orig, dir, dist);
        vec3 p = orig + dir*dist;
        long unsigned a = (atan2(p.z, p.x) / (2 * M_PI) + 0.5f) * envmap_width;
        long unsigned b = acos(p.y / 100.0f) / M_PI * envmap_height;
        std::cout << "End inside" << a+b*envmap_width << std::endl;
        return envmap[a+b*envmap_width];
    }
    std::cout << "Begin outside" << std::endl;
    std::cout << "Begin norma" << std::endl;
    vec3 reflectDir = normalize(reflect(dir, N));
    vec3 refractDir = normalize(refract(dir, N, material.refractive));
    vec3 reflectOrig = dot(reflectDir, N) < 0.0f ? point - N * 0.001f : point + N * 0.001f;
    vec3 refractOrig = dot(refractDir, N) < 0.0f ? point - N * 0.001f : point + N * 0.001f;
    std::cout << "End norma" << std::endl;
    std::cout << "Begin reflect" << std::endl;
    vec3 reflectColor = traceRay(reflectOrig, reflectDir, sceneObjects, depth + 1);
    vec3 refractColor = traceRay(refractOrig, refractDir, sceneObjects, depth + 1);
    std::cout << "End norma" << std::endl;
    float diffuseLightIntensity = 0;
    float specularLightIntensity = 0;
    for (auto &light: sceneObjects.lights) {
        vec3 lightDir = normalize(light.position - point);
        float lightDistance = distance(light.position, point);
        vec3 shadowOrig = dot(lightDir, N) < 0 ? point - N * 0.001f : point + N * 0.001f;
        vec3 shadowPt, shadowN;
        Material tmpmaterial;
        if (sceneIntersect(shadowOrig, lightDir, sceneObjects, shadowPt, shadowN, tmpmaterial) &&
            distance(shadowPt, shadowOrig) < lightDistance)
            continue;

        diffuseLightIntensity += light.intensity * std::max(0.0f, dot(lightDir, N));
        specularLightIntensity += powf(std::max(0.0f, dot(-reflect(-lightDir, N), dir)),
            material.specular) * light.intensity;
    }
    std::cout << "End trace" << std::endl;
    return material.diffuse * diffuseLightIntensity * material.albedo[0] +
           vec3(1.0f, 1.0f, 1.0f) * specularLightIntensity * material.albedo[1] +
           reflectColor * material.albedo[2] + refractColor * material.albedo[3];
}

std::deque<Sphere> Scene1::loadSpheres(const Materials &materials) {
    std::deque<Sphere> spheres;
    spheres.push_front(Sphere(vec3( 7, 5, -18), 4, materials["mirror"]));
    spheres.push_front(Sphere(vec3( 14, 5, -25), 4, materials["glass"]));
    spheres.push_front(Sphere(vec3( -7, 5, -20), 4, materials["red_rubber"]));
    return spheres;
}

std::deque<Light> Scene1::loadLights() {
    std::deque<Light> lights;
    lights.push_front(Light(vec3(-20, 20,  20), 1.5));
    lights.push_front(Light(vec3( 30, 50, -25), 1.8));
    lights.push_front(Light(vec3( 30, 20,  30), 1.7));
    return lights;
}

std::deque<Triangle> Scene1::loadTriangles(const Materials &materials) {
    std::deque<Triangle> triangles;
    Material glass = materials["pyramid_glass"];
    triangles.push_front(Triangle(vec3(8.0f, -4.f, -17.f), vec3(8, -4, -19), vec3(7, -1, -18), glass)); //right
    triangles.push_front(Triangle(vec3(6.f, -4.f, -19.f), vec3(6, -4, -17), vec3(7, -1, -18), glass)); //left
    triangles.push_front(Triangle(vec3(6, -4, -17), vec3(8, -4, -17), vec3(7, -1, -18), glass)); //front
    triangles.push_front(Triangle(vec3(8, -4, -19), vec3(6, -4, -19), vec3(7, -1, -18), glass)); //back
    Quadrangle mirror(vec3(-10, -4, -15), vec3(-10, -4, -25), vec3(-10, 4, -25), vec3(-10, 4, -15), materials["mirror"]);
    std::deque<Triangle> mirTr = mirror.toTriangles();
    triangles.insert(triangles.begin(), mirTr.begin(), mirTr.end());
    return triangles;
}

std::deque<Model> Scene1::loadModels(const Materials &materials) {
    std::deque<Model> models;
    vec3 bias(4.0f, -7.0f + 1.7f, -25.0f);
    Model bunny("../resources/bunny.obj", 40.0f, bias, materials["gold"]);
    models.push_front(bunny);
    return models;
}

void Scene1::render(const Settings &settings) {
    const Materials materials;
    std::deque<Sphere> spheres = loadSpheres(materials);
    std::deque<Triangle> triangles = loadTriangles(materials);
    std::deque<Light> lights = loadLights();
    std::deque<Model> models = loadModels(materials);
    Cube cube(vec3(-6, -4, -20), vec3(-3, -1, -23), materials["silver"]);
    std::deque<Cube> cubes { cube };
    //cube.printFaces();
    //std::deque<Cube> cubes;
    SceneObjects sceneObjects(spheres, triangles, lights, models, cubes);

    const float fov = M_PI / 3.0f;
    std::vector<vec3> frameBuffer(width * height);

    vec3 camera(0, 5, 0);
//#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < height; ++i) {
        for (size_t j = 0; j < width; ++j) {
            float x =  (j + 0.5f) -  width / 2.0f;
            float y = -(i + 0.5f) + height / 2.0f;
            float z = -height / (2.0f * tan(fov / 2.0f));
            vec3 dir = normalize(vec3(x, y, z));
            vec3 result = traceRay(camera, dir, sceneObjects);
            frameBuffer[i * width + j] = result;
        }
    }

    bitmap_image image(width, height);
    for (size_t y = 0; y < image.height(); ++y) {
        for (size_t x = 0; x < image.width(); ++x) {
            vec3 &c = frameBuffer[y * width + x];
            float max = std::max(c[0], std::max(c[1], c[2]));
            if (max > 1) {
                c *= 1.0 / max;
            }
            char red = 255 * std::max(0.0f, std::min(1.0f, frameBuffer[y * width + x][0]));
            char green = 255 * std::max(0.0f, std::min(1.0f, frameBuffer[y * width + x][1]));
            char blue = 255 * std::max(0.0f, std::min(1.0f, frameBuffer[y * width + x][2]));
            rgb_t color(red, green, blue);
            image.set_pixel(x, y, color);
        }
    }
    image.save_image(settings.out);
}

int Scene1::run(const Settings &settings) {
    //omp_set_num_threads(settings.threads);
    int n = -1;
    unsigned char *pixmap = stbi_load("../resources/map3.jpg", &envmap_width, &envmap_height, &n, 0);
    if (!pixmap || 3!=n) {
        std::cout << "Error: can not load the environment map" << std::endl;
        return -1;
    }
    envmap = std::vector<vec3>(envmap_width * envmap_height);
    for (int j = envmap_height-1; j>=0 ; j--) {
        for (int i = 0; i<envmap_width; i++) {
            envmap[i+j*envmap_width] = vec3(pixmap[(i+j*envmap_width)*3+0],
                                            pixmap[(i+j*envmap_width)*3+1],
                                            pixmap[(i+j*envmap_width)*3+2]) * (1/255.0f);
        }
    }
    stbi_image_free(pixmap);
    std::cout << "[Current settings]" << std::endl;
    std::cout << "Output path: " << settings.out << std::endl;
    std::cout << "Scene number: " << settings.scene << std::endl;
    std::cout << "Threads: " << settings.threads << std::endl;
    //std::cout << "Real threads: " << omp_get_max_threads() << std::endl;
    std::cout << "\nRendering image..." << std::endl;
    //double time0 = omp_get_wtime();
    render(settings);
    std::cout << "Ok." << std::endl;
    //double time1 = omp_get_wtime();
    //std::cout << "Done. Elapsed time: " << time1 - time0 << std::endl;
    return 0;
}