#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <cmath>
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

const int width = 1024;
const int height = 768;

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
                    const std::vector<Triangle> &triangles,
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
    float checkerboardDist = std::numeric_limits<float>::max();
    if (fabs(dir.y) > 0.001) {
        float d = -(orig.y + 4) / dir.y;
        vec3 pt = orig + dir * d;
        if (d > 0 && fabs(pt.x) < 10 && pt.z < -10 && pt.z > -30 && d < spheresDist && d < trianglesDist) {
            checkerboardDist = d;
            hit = pt;
            N = vec3(0, 1, 0);
            material.diffuse = (int(0.5f * hit.x + 1000) + int(0.5 * hit.z)) & 1 ? 
                vec3(.3, .3, .3) : vec3(.3, .2, .1);
        }
    }
    return std::min(trianglesDist, std::min(spheresDist, checkerboardDist)) < 1000;
}

vec3 traceRay(const vec3 &orig, const vec3 &dir, const std::vector<Sphere> &spheres,
              const std::vector<Triangle> &triangles,
              const std::vector<Light> &lights, size_t depth = 0) {
    vec3 point, N;
    Material material;

    if (depth > 4 || !sceneIntersect(orig, dir, spheres, triangles, point, N, material)) {
        return vec3(0.2, 0.7, 0.8);
    }
    vec3 reflectDir = normalize(reflect(dir, N));
    vec3 refractDir = normalize(refract(dir, N, material.refractive));
    vec3 reflectOrig = dot(reflectDir, N) < 0.0f ? point - N * 0.001f : point + N * 0.001f;
    vec3 refractOrig = dot(refractDir, N) < 0.0f ? point - N * 0.001f : point + N * 0.001f;
    vec3 reflectColor = traceRay(reflectOrig, reflectDir, spheres, triangles, lights, depth + 1);
    vec3 refractColor = traceRay(refractOrig, refractDir, spheres, triangles, lights, depth + 1); 
    float diffuseLightIntensity = 0;
    float specularLightIntensity = 0;
    for (auto &light: lights) {
        vec3 lightDir = normalize(light.position - point);
        float lightDistance = distance(light.position, point);
        vec3 shadowOrig = dot(lightDir, N) < 0 ? point - N * 0.001f : point + N * 0.001f;
        vec3 shadowPt, shadowN;
        Material tmpmaterial;
        if (sceneIntersect(shadowOrig, lightDir, spheres, triangles, shadowPt, shadowN, tmpmaterial) &&
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
    spheres.push_back(Sphere(vec3(-3, 0, -16), 2, materials["ivory"]));
    spheres.push_back(Sphere(vec3(-1.0, -1.5, -12), 2, materials["glass"]));
    spheres.push_back(Sphere(vec3( 1.5, -0.5, -18), 3, materials["gold"]));
    spheres.push_back(Sphere(vec3( 7, 5, -18), 4, materials["mirror"]));
    spheres.push_back(Sphere(vec3(-5, 5, -20), 4, materials["silver"]));
    return spheres;
}

std::vector<Light> loadLights() {
    std::vector<Light> lights;
    lights.push_back(Light(vec3(-20, 20,  20), 1.5));
    lights.push_back(Light(vec3( 30, 50, -25), 1.8));
    lights.push_back(Light(vec3( 30, 20,  30), 1.7));
    return lights;
}

std::vector<Triangle> loadTriangles(const Materials &materials) {
    std::vector<Triangle> triangles;
    Material glass = materials["pyramid_glass"];
    triangles.push_back(Triangle(vec3(8, -4, -17), vec3(8, -4, -19), vec3(7, -1, -18), glass)); //right
    triangles.push_back(Triangle(vec3(6, -4, -19), vec3(6, -4, -17), vec3(7, -1, -18), glass)); //left
    triangles.push_back(Triangle(vec3(6, -4, -17), vec3(8, -4, -17), vec3(7, -1, -18), glass)); //front
    triangles.push_back(Triangle(vec3(8, -4, -19), vec3(6, -4, -19), vec3(7, -1, -18), glass)); //back
    Quadrangle mirror(vec3(-10, -4, -15), vec3(-10, -4, -25), vec3(-10, 4, -25),
                      vec3(-10, 4, -15), materials["mirror"]);
    std::vector<Triangle> mirTr = mirror.toTriangles();
    triangles.insert(triangles.end(), mirTr.begin(), mirTr.end());
    return triangles;
}

void render(Settings &settings) {
    const Materials materials;
    std::vector<Sphere> spheres = loadSpheres(materials);
    std::vector<Triangle> triangles = loadTriangles(materials);
    std::vector<Light> lights = loadLights();

    const float fov = M_PI / 3.0f;
    std::vector<vec3> frameBuffer(width * height);
//#pragma omp parallel for
    for (size_t i = 0; i < height; ++i) {
        for (size_t j = 0; j < width; ++j) {
            float x =  (j + 0.5f) -  width / 2.0f;
            float y = -(i + 0.5f) + height / 2.0f;
            float z = -height / (2.0f * tan(fov / 2.0f));
            vec3 dir = normalize(vec3(x, y, z));
            frameBuffer[i * width + j] = traceRay(vec3(0, 0, 0), dir, spheres, triangles, lights);
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
    std::cout << "[Current settings]" << std::endl;
    std::cout << "Output path: " << settings.out << std::endl;
    std::cout << "Scene number: " << settings.scene << std::endl;
    std::cout << "Threads: " << settings.threads << std::endl;
    std::cout << "\nRendering image..." << std::endl;
    render(settings);
    std::cout << "Done." << std::endl;
    return 0;
}