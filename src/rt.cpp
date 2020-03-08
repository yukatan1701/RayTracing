#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <cmath>
#include <omp.h>
#include "bitmap_image.hpp"
#include "glm/vec2.hpp"
#include "glm/vec3.hpp"
#include "glm/vec4.hpp"
#include "glm/gtx/intersect.hpp"
#include "glm/gtx/norm.hpp"
#include "glm/gtx/normal.hpp"
#include "Material.h"
#include "SceneObjects.h"
#include "ParseException.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

const int width = 1024;
const int height = 768;

int procsNum = 1;

int envmap_width, envmap_height;
std::vector<vec3> envmap;

using namespace glm;

struct Settings {
    std::string out;
    unsigned int scene;
    unsigned int threads;
    Settings() : out("rt.bmp"), scene(1), threads(1) {}
    Settings(std::string out, unsigned int scene, unsigned int threads) :
        out(out), scene(scene), threads(threads) {}
};

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

bool sceneIntersect(const vec3 &orig, const vec3 &dir, const std::vector<Sphere> &spheres,
                    const std::deque<Triangle> &triangles, const std::vector<Model> &models,
                    vec3 &hit, vec3 &N, Material &material) {
    float spheresDist = std::numeric_limits<float>::max();
    for (auto &sphere: spheres) {
        float curDist;
        if (sphere.rayIntersect(orig, dir, curDist) && curDist < spheresDist) {
            spheresDist = curDist;
            hit = orig + dir * curDist;
            N = normalize(hit - sphere.center);
            material = sphere.material;
        }
    }
    float trianglesDist = std::numeric_limits<float>::max();
    for (auto &triangle : triangles) {
        float curDist;
        if (triangle.rayIntersect(orig, dir, curDist) && curDist < trianglesDist && curDist < spheresDist) {
            trianglesDist = curDist;
            hit = orig + dir * curDist;
            N = triangleNormal(triangle.v0, triangle.v1, triangle.v2);
            material = triangle.material;
        }
    }
    float modelsDist = std::numeric_limits<float>::max();
    int n = 1;
    for (auto &model : models) {
        float curDist;
        if (model.cubeIntersect(orig, dir, curDist) && curDist < modelsDist &&
                    curDist < trianglesDist && curDist < spheresDist) {
            if (n == 1) {
                //std::cout << "Success!" << std::endl;
            }
            modelsDist = curDist;
            hit = orig + dir * curDist;
            N = vec3(0, -1, 0);
            material = Material();
            /*
            for (auto &triangle : model.triangles) {
                float curDist;
                if (triangle.rayIntersect(orig, dir, curDist) && curDist < modelsDist &&
                    curDist < trianglesDist && curDist < spheresDist) {
                    modelsDist = curDist;
                    hit = orig + dir * curDist;
                    N = triangleNormal(triangle.v0, triangle.v1, triangle.v2);
                    material = triangle.material;
                }
            }*/
        }
    }
    float checkerboardDist = std::numeric_limits<float>::max();
    if (fabs(dir.y) > 0.001) {
        float d = -(orig.y + 4) / dir.y;
        vec3 pt = orig + dir * d;
        if (d > 0 && fabs(pt.x) < 10 && pt.z < -10 && pt.z > -30 && d < modelsDist &&
            d < spheresDist && d < trianglesDist) {
            checkerboardDist = d;
            hit = pt;
            N = vec3(0, 1, 0);
            material.diffuse = (int(0.5f * hit.x + 1000) + int(0.5 * hit.z)) & 1 ? 
                vec3(.3, .3, .3) : vec3(.3, .2, .1);
        }
    }
    return std::min(modelsDist, std::min(trianglesDist, std::min(spheresDist, checkerboardDist))) < 1000;
}

vec3 traceRay(const vec3 &orig, const vec3 &dir, const std::vector<Sphere> &spheres,
              const std::deque<Triangle> &triangles, const std::vector<Model> &models,
              const std::vector<Light> &lights, size_t depth = 0) {
    vec3 point, N;
    Material material;

    if (depth > 4 || !sceneIntersect(orig, dir, spheres, triangles, models, point, N, material)) {
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
    vec3 reflectColor = traceRay(reflectOrig, reflectDir, spheres, triangles, models, lights, depth + 1);
    vec3 refractColor = traceRay(refractOrig, refractDir, spheres, triangles, models, lights, depth + 1); 
    float diffuseLightIntensity = 0;
    float specularLightIntensity = 0;
    for (auto &light: lights) {
        vec3 lightDir = normalize(light.position - point);
        float lightDistance = distance(light.position, point);
        vec3 shadowOrig = dot(lightDir, N) < 0 ? point - N * 0.001f : point + N * 0.001f;
        vec3 shadowPt, shadowN;
        Material tmpmaterial;
        if (sceneIntersect(shadowOrig, lightDir, spheres, triangles, models, shadowPt, shadowN, tmpmaterial) &&
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

std::vector<Sphere> loadSpheres(const Materials &materials) {
    std::vector<Sphere> spheres;
    /*spheres.push_back(Sphere(vec3(-3, 0, -16), 2, materials["ivory"]));
    spheres.push_back(Sphere(vec3(-1.0, -1.5, -12), 2, materials["glass"]));
    spheres.push_back(Sphere(vec3( 1.5, -0.5, -18), 3, materials["gold"]));/**/
    spheres.push_back(Sphere(vec3( 7, 5, -18), 4, materials["mirror"]));
    //spheres.push_back(Sphere(vec3(-5, 5, -20), 4, materials["silver"]));
    return spheres;
}

std::vector<Light> loadLights() {
    std::vector<Light> lights;
    lights.push_back(Light(vec3(-20, 20,  20), 1.5));
    lights.push_back(Light(vec3( 30, 50, -25), 1.8));
    lights.push_back(Light(vec3( 30, 20,  30), 1.7));
    return lights;
}

std::deque<Triangle> loadTriangles(const Materials &materials) {
    std::deque<Triangle> triangles;
    Material glass = materials["pyramid_glass"];
    /*triangles.push_front(Triangle(vec3(8, -4, -17), vec3(8, -4, -19), vec3(7, -1, -18), glass)); //right
    triangles.push_front(Triangle(vec3(6, -4, -19), vec3(6, -4, -17), vec3(7, -1, -18), glass)); //left
    triangles.push_front(Triangle(vec3(6, -4, -17), vec3(8, -4, -17), vec3(7, -1, -18), glass)); //front
    triangles.push_front(Triangle(vec3(8, -4, -19), vec3(6, -4, -19), vec3(7, -1, -18), glass)); //back
    Quadrangle mirror(vec3(-10, -4, -15), vec3(-10, -4, -25), vec3(-10, 4, -25),
                      vec3(-10, 4, -15), materials["mirror"]);
    std::deque<Triangle> mirTr = mirror.toTriangles();
    triangles.insert(triangles.begin(), mirTr.begin(), mirTr.end());*/
    return triangles;
}

void render(Settings &settings) {
    const Materials materials;
    std::vector<Sphere> spheres = loadSpheres(materials);
    std::deque<Triangle> triangles = loadTriangles(materials);
    std::vector<Light> lights = loadLights();
    Model bunny("../resources/bunny.obj", 100.0f, materials["ivory"]);
    std::vector<Model> models;
    models.push_back(bunny);
    //triangles.insert(triangles.begin(), bunny.triangles.begin(), bunny.triangles.end());

    const float fov = M_PI / 3.0f;
    const int samples = 8;
    std::vector<vec3> frameBuffer(width * height);

    vec3 camera(0, 5, 0);
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < height; ++i) {
        for (size_t j = 0; j < width; ++j) {
            //frameBuffer[i * width + j] = vec3(0, 0, 0);
            //for (int s = 0; s < samples; ++s) {
                float x =  (j + 0.5f) -  width / 2.0f;
                float y = -(i + 0.5f) + height / 2.0f;
                float z = -height / (2.0f * tan(fov / 2.0f));
                vec3 dir = normalize(vec3(x, y, z));
                frameBuffer[i * width + j] = traceRay(camera, dir, spheres, triangles, models, lights);
            //}
            //frameBuffer[i * width + j] /= (float) samples;
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

Settings parseArgs(int argc, char **argv) {
    int paramCount = 4;
    if (argc == 1)
        return Settings();
    if (argc < 2 * paramCount - 1) {
        throw ParseException(ParseException::NOT_ENOUGH);
    }
    Settings settings("", 0, 0);
    for (int i = 1; i < 2 * paramCount - 1; i += 2) {
        std::string arg = argv[i];
        if (arg == "-out") {
            if (settings.out.length() > 0) {
                throw ParseException(ParseException::SAME_ARGS);
            }
            settings.out = argv[i + 1];
        } else if (arg == "-scene") {
            if (settings.scene != 0) {
                throw ParseException(ParseException::SAME_ARGS);
            }
            settings.scene = atoi(argv[i + 1]);
        } else if (arg == "-threads") {
            if (settings.threads != 0) {
                throw ParseException(ParseException::SAME_ARGS);
            }
            settings.threads = atoi(argv[i + 1]);
        }
    }
    return settings;
}

int main(int argc, char **argv) {
    Settings settings;
    try {
        settings = parseArgs(argc, argv);
    } catch (ParseException &e) {
        e.printMessage();
        return -1;
    }
    procsNum = omp_get_num_procs();
    omp_set_num_threads(procsNum);
    Cube cube(vec3(0, 0, 0), vec3(1, 1, 1), Material());
    cube.printFaces();
    int n = -1;
    unsigned char *pixmap = stbi_load("../resources/map3.jpg", &envmap_width, &envmap_height, &n, 0);
    if (!pixmap || 3!=n) {
        std::cout << "Error: can not load the environment map" << std::endl;
        return -1;
    }
    envmap = std::vector<vec3>(envmap_width*envmap_height);
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
    std::cout << "\nRendering image..." << std::endl;
    double time0 = omp_get_wtime();
    render(settings);
    double time1 = omp_get_wtime();
    std::cout << "Done. Elapsed time: " << time1 - time0 << std::endl;
    return 0;
}