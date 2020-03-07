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

vec3 refract(const vec3 &I, const vec3 &N, const float refractiveIndex) {
    float cosi = -std::max(-1.0f, std::min(1.0f, dot(I, N)));
    float etai = 1.0f, etat = refractiveIndex;
    vec3 n = N;
    if (cosi < 0) {
        cosi = -cosi;
        std::swap(etai, etat);
        n = -N;
    }
    float eta = etai / etat;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? vec3(0, 0, 0) : I * eta + n * (eta * cosi - sqrtf(k));
}

bool sceneIntersect(const vec3 &orig, const vec3 &dir, const std::vector<Sphere> &spheres,
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
    float checkerboardDist = std::numeric_limits<float>::max();
    if (fabs(dir.y) > 1e-3) {
        float d = -(orig.y + 4) / dir.y;
        vec3 pt = orig + dir * d;
        if (d > 0 && fabs(pt.x) < 10 && pt.z < -10 && pt.z > -30 && d < spheresDist) {
            checkerboardDist = d;
            hit = pt;
            N = vec3(0, 1, 0);
            material.diffuseColor = (int(0.5f * hit.x + 1000) + int(0.5 * hit.z)) & 1 ? 
                vec3(1, 1, 1) : vec3(1, 0.7f, 0.3f);
            material.diffuseColor *= 0.3f;
        }
    }
    return std::min(spheresDist, checkerboardDist) < 1000;
}

vec3 traceRay(const vec3 &orig, const vec3 &dir, const std::vector<Sphere> &spheres,
              const std::vector<Light> &lights, size_t depth = 0) {
    vec3 point, N;
    Material material;

    if (depth > 4 || !sceneIntersect(orig, dir, spheres, point, N, material)) {
        return vec3(0.2, 0.7, 0.8);
    }
    vec3 reflectDir = normalize(reflect(dir, N));
    vec3 refractDir = normalize(refract(dir, N, material.refractiveIndex));
    vec3 reflectOrig = dot(reflectDir, N) < 0 ? point - N * (float) 1e-3 : point + N * (float) 1e-3;
    vec3 refractOrig = dot(refractDir, N) < 0 ? point - N * (float) 1e-3 : point + N * (float) 1e-3;
    vec3 reflectColor = traceRay(reflectOrig, reflectDir, spheres, lights, depth + 1);
    vec3 refractColor = traceRay(refractOrig, refractDir, spheres, lights, depth + 1); 
    float diffuseLightIntensity = 0;
    float specularLightIntensity = 0;
    for (auto &light: lights) {
        vec3 lightDir = normalize(light.position - point);
        float lightDistance = distance2(light.position, point);
        vec3 shadowOrig = dot(lightDir, N) < 0 ? point - N * (float) 1e-3 : point + N * (float) 1e-3;
        vec3 shadowPt, shadowN;
        Material tmpmaterial;
        if (sceneIntersect(shadowOrig, lightDir, spheres, shadowPt, shadowN, tmpmaterial) &&
            distance2(shadowPt, shadowOrig) < lightDistance)
            continue;

        diffuseLightIntensity += light.intensity * std::max(0.0f, dot(lightDir, N));
        specularLightIntensity += powf(std::max(0.0f, dot(reflect(lightDir, N),dir)),
            material.specularExponent) * light.intensity;
    }
    return material.diffuseColor * diffuseLightIntensity * material.albedo[0] +
           vec3(1.0f, 1.0f, 1.0f) * specularLightIntensity * material.albedo[1] +
           reflectColor * material.albedo[2] + refractColor * material.albedo[3];
}

void render(Settings &settings) {
    Materials materials;
    std::vector<Sphere> spheres;
    spheres.push_back(Sphere(vec3(-3, 0, -16), 4, materials.get("ivory")));
    spheres.push_back(Sphere(vec3(-1.0, -1.5, -12), 4, materials.get("glass")));
    spheres.push_back(Sphere(vec3( 1.5, -0.5, -18), 8, materials.get("red_rubber")));
    spheres.push_back(Sphere(vec3( 7, 5, -18), 10, materials.get("mirror")));

    std::vector<Light> lights;
    lights.push_back(Light(vec3(-20, 20,  20), 1.5));
    lights.push_back(Light(vec3( 30, 50, -25), 1.8));
    lights.push_back(Light(vec3( 30, 20,  30), 1.7));

    const int fov = M_PI / 3.0f;
    std::vector<vec3> frameBuffer(width * height);
//#pragma omp parallel for
    for (size_t i = 0; i < height; ++i) {
        for (size_t j = 0; j < width; ++j) {
            float x =  (j + 0.5f) -  width / 2.0f;
            float y = -(i + 0.5f) + height / 2.0f;
            float z = -height / (2.0f * tan(fov / 2.0f));
            vec3 dir = normalize(vec3(x, y, z));
            frameBuffer[i * width + j] = traceRay(vec3(0.0f, 0.0f, 0.0f), dir, spheres, lights);
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