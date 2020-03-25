#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <cmath>
#include <ctime>
#include <sys/time.h>
//#define OPENMP
#include <omp.h>
#include <set>
#include "bitmap_image.hpp"
#include "glm/vec2.hpp"
#include "glm/vec3.hpp"
#include "glm/vec4.hpp"
#include "glm/gtx/norm.hpp"
#include "glm/gtx/normal.hpp"
#include "Material.h"
#include "SceneObjects.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "Scene1.h"

const int width = 512;
const int height = 512;

double bench_t_start, bench_t_end;

int envmap_width, envmap_height;
std::vector<vec3> envmap;

using namespace glm;

static double rtclock() {
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, NULL);
    if (stat != 0)
      printf ("Error return from gettimeofday: %d", stat);
    return (Tp.tv_sec + Tp.tv_usec * 1.0e-6);
}

void bench_timer_start() {
  bench_t_start = rtclock ();
}

void bench_timer_stop() {
  bench_t_end = rtclock ();
}

vec3 reflect(const vec3 &i, const vec3 &n) {
    return i - n * 2.0f * dot(i, n);
}

vec3 refract(const vec3 &I, const vec3 &N, const float etat, const float etai = 1.0f) {
    float cosi = -std::max(-1.0f, std::min(1.0f, dot(I, N)));
    if (cosi < 0)
        return refract(I, -N, etai, etat);
    float eta = etai / etat;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? vec3(1, 0, 0) : I * eta + N * (eta * cosi - sqrtf(k));
}

bool sceneIntersect(const vec3 &orig, const vec3 &dir, const SceneObjects &sceneObjects,
                    vec3 &hit, vec3 &N, Material &material) {
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
        vec3 normal = vec3(0, 0, 1);
        if (cube.rayIntersect(orig, dir, curDist, normal) && curDist < cubesDist &&
            curDist < min_dist) {
            cubesDist = curDist;
            hit = orig + dir * curDist;
            N = normal;
            material = cube.material;
        }
    }
    min_dist = std::min(cubesDist, min_dist);
    /*
    float modelsDist = std::numeric_limits<float>::max();
    for (auto &model : sceneObjects.models) {
        float boxDist, curDist;
        vec3 n;
        if (model.boxIntersect(orig, dir, boxDist, n)) {
            vec3 hitPoint = orig + dir * boxDist;
            int size = model.box.size;
            size = 2;
            int begin = 1;
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < 1; j++) {
                    for (int k = 0; k < 1; k++) {
                        if (model.box.grid[i][j][k].rayIntersect(orig, dir, boxDist, n)) {
                            modelsDist = boxDist;
                            hit = orig + dir * boxDist;
                            N = n;
                            //material = model.box.grid[i][j][k].material;
                            /*for (const Triangle *triangle : model.box.tgrid[i][j][k]) {  
                                if (triangle->rayIntersect(orig, dir, curDist) && curDist < modelsDist &&
                                    curDist < *distances.begin()) {
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
            /*for (auto &triangle : model.triangles) {
                float curDist;
                if (triangle.rayIntersect(orig, dir, curDist) && curDist < modelsDist &&
                    curDist < *distances.begin()) {
                    modelsDist = curDist;
                    hit = orig + dir * curDist;
                    N = triangleNormal(triangle.v0, triangle.v1, triangle.v2);
                    material = triangle.material;
                }
            }
        }
    }
    distances.insert(modelsDist);*/
    float checkerboardDist = std::numeric_limits<float>::max();
    if (fabs(dir.y) > 0.001) {
        float d = -(orig.y + 4) / dir.y;
        vec3 pt = orig + dir * d;
        if (d > 0 && fabs(pt.x) < 10 && pt.z < -10 && pt.z > -30 && d < min_dist) {
            checkerboardDist = d;
            hit = pt;
            N = vec3(0, 1, 0);
            material.diffuse = (int(0.5f * hit.x + 1000) + int(0.5 * hit.z)) & 1 ? 
                vec3(.3, .3, .3) : vec3(.3, .2, .1);
        }
    }
    min_dist = std::min(checkerboardDist, min_dist);
    return min_dist < 1000;
}

vec3 traceRay(const vec3 &orig, const vec3 &dir, const SceneObjects &sceneObjects,
              size_t depth = 0) {
    vec3 point, N;
    Material material;

    if (depth > 4 || !sceneIntersect(orig, dir, sceneObjects, point, N, material)) {
        Sphere env(vec3(0,0,0), 100, Material());
        float dist = 0;
        env.rayIntersect(orig, dir, dist);
        vec3 p = orig + dir*dist;
        int a = (atan2(p.z, p.x) / (2 * M_PI) + 0.5f) * envmap_width;
        int b = acos(p.y / 100.0f) / M_PI * envmap_height;
        return envmap[a+b*envmap_width];
    }
    vec3 reflectDir = normalize(reflect(dir, N));
    vec3 refractDir = normalize(refract(dir, N, material.refractive));
    vec3 reflectOrig = dot(reflectDir, N) < 0.0f ? point - N * 0.001f : point + N * 0.001f;
    vec3 refractOrig = dot(refractDir, N) < 0.0f ? point - N * 0.001f : point + N * 0.001f;
    vec3 reflectColor = traceRay(reflectOrig, reflectDir, sceneObjects, depth + 1);
    vec3 refractColor = traceRay(refractOrig, refractDir, sceneObjects, depth + 1); 
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
    return material.diffuse * diffuseLightIntensity * material.albedo[0] +
           vec3(1.0f, 1.0f, 1.0f) * specularLightIntensity * material.albedo[1] +
           reflectColor * material.albedo[2] + refractColor * material.albedo[3];
}

std::deque<Sphere> loadSpheres(const Materials &materials) {
    std::deque<Sphere> spheres;
    spheres.push_front(Sphere(vec3( 7, 5, -18), 4, materials["mirror"]));
    spheres.push_front(Sphere(vec3( 14, 5, -25), 4, materials["glass"]));
    spheres.push_front(Sphere(vec3( -7, 5, -20), 4, materials["red_rubber"]));
    return spheres;
}

std::deque<Light> loadLights() {
    std::deque<Light> lights;
    lights.push_front(Light(vec3(-20, 20,  20), 1.5));
    lights.push_front(Light(vec3( 30, 50, -25), 1.8));
    lights.push_front(Light(vec3( 30, 20,  30), 1.7));
    return lights;
}

std::deque<Triangle> loadTriangles(const Materials &materials) {
    std::deque<Triangle> triangles;
    Material glass = materials["pyramid_glass"];
    triangles.push_front(Triangle(vec3(8, -4, -17), vec3(8, -4, -19), vec3(7, -1, -18), glass)); //right
    triangles.push_front(Triangle(vec3(6, -4, -19), vec3(6, -4, -17), vec3(7, -1, -18), glass)); //left
    triangles.push_front(Triangle(vec3(6, -4, -17), vec3(8, -4, -17), vec3(7, -1, -18), glass)); //front
    triangles.push_front(Triangle(vec3(8, -4, -19), vec3(6, -4, -19), vec3(7, -1, -18), glass)); //back
    Quadrangle mirror(vec3(-10, -4, -15), vec3(-10, -4, -25), vec3(-10, 4, -25),
                      vec3(-10, 4, -15), materials["mirror"]);
    std::deque<Triangle> mirTr = mirror.toTriangles();
    triangles.insert(triangles.begin(), mirTr.begin(), mirTr.end());
    return triangles;
}

std::deque<Model> loadModels(const Materials &materials) {
    std::deque<Model> models;
    /*vec3 bias(4.0f, -7.0f + 1.7f, -25.0f);
    Model bunny("../resources/bunny.obj", 40.0f, bias, materials["gold"]);
    models.push_front(bunny);*/
    return models;
}

void render(const Settings &settings) {
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
#ifdef OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (size_t i = 0; i < height; ++i) {
        for (size_t j = 0; j < width; ++j) {
            float x =  (j + 0.5f) -  width / 2.0f;
            float y = -(i + 0.5f) + height / 2.0f;
            float z = -height / (2.0f * tan(fov / 2.0f));
            vec3 dir = normalize(vec3(x, y, z));
            frameBuffer[i * width + j] = traceRay(camera, dir, sceneObjects);
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

int Scene1::run(const Settings &settings) const {
    #ifdef OPENMP
    omp_set_num_threads(settings.threads);
    #endif

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
    #ifdef OPENMP
    double time0 = omp_get_wtime();
    #endif
    #ifndef OPENMP
    bench_timer_start();
    #endif
    render(settings);
    #ifdef OPENMP
    double time1 = omp_get_wtime();
    std::cout << "Done. Elapsed time: " << time1 - time0 << std::endl;
    #endif
    #ifndef OPENMP
    bench_timer_stop();
    printf ("Time in seconds = %0.6lf\n", bench_t_end - bench_t_start);
    #endif
    return 0;
}