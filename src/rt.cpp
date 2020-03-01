#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <cmath>
#include "bitmap_image.hpp"
#include "glm/vec2.hpp"
#include "glm/vec3.hpp"
#include "glm/gtx/intersect.hpp"
#include "glm/gtx/norm.hpp"
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

struct Material {
    vec2 albedo;
    vec3 diffuseColor;
    float specularExponent;
    Material() : albedo(1, 0), diffuseColor(), specularExponent() {}
    Material(const vec2 &a, const vec3 &color, const float &spec) :
        albedo(a), diffuseColor(color), specularExponent(spec) {}
};

struct Sphere {
    vec3 center;
    float radius;
    Material material;

    Sphere(const vec3 &center, const float &radius): center(center), radius(radius) {}
    Sphere(const vec3 &c, const float &r, const Material &m): center(c), radius(r), material(m) {}

    bool rayIntersect(const vec3 &orig, const vec3 &dir, float &dist) const {
        return intersectRaySphere(orig, dir, center, radius, dist);
    }
};

struct Light {
    vec3 position;
    float intensity;
    Light(const vec3 &pos, const float &i) : position(pos), intensity(i) {}
};

vec3 reflect(const vec3 &i, const vec3 &n) {
    return i - n * 2.0f * dot(i, n);
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
    return spheresDist < 1000;
}

vec3 castRay(const vec3 &orig, const vec3 &dir, const std::vector<Sphere> &spheres,
             const std::vector<Light> &lights) {
    vec3 point, N;
    Material material;

    if (!sceneIntersect(orig, dir, spheres, point, N, material)) {
        return vec3(0.2, 0.7, 0.8);
    }
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
           vec3(1.0f, 1.0f, 1.0f) * specularLightIntensity * material.albedo[1];
}

void render(Settings &settings) {
    Material ivory(vec2(0.6, 0.3), vec3(0.4, 0.4, 0.3), 50.0);
    Material red_rubber(vec2(0.9, 0.1), vec3(0.3, 0.1, 0.1), 10.0);
    std::vector<Sphere> spheres;
    spheres.push_back(Sphere(vec3(-3,    0,   -16), 4,      ivory));
    spheres.push_back(Sphere(vec3(-1.0, -1.5, -12), 4, red_rubber));
    spheres.push_back(Sphere(vec3( 1.5, -0.5, -18), 8, red_rubber));
    spheres.push_back(Sphere(vec3( 7,    5,   -18), 10,      ivory));

    std::vector<Light> lights;
    lights.push_back(Light(vec3(-20, 20,  20), 1.5));
    lights.push_back(Light(vec3( 30, 50, -25), 1.8));
    lights.push_back(Light(vec3( 30, 20,  30), 1.7));

    const int fov = M_PI / 2.0f;
    std::vector<vec3> frameBuffer(width * height);
//#pragma omp parallel for
    for (size_t i = 0; i < height; ++i) {
        for (size_t j = 0; j < width; ++j) {
            float y = -(2.0f * (i + 0.5f) / (float) width  - 1) * tan(fov / 2.0f) * width / (float) height;
            float x = (2.0f * (j + 0.5f) / (float) height - 1) * tan(fov / 2.0f);
            vec3 dir = normalize(vec3(x, y, -1.0f));
            frameBuffer[i * width + j] = castRay(vec3(0.0f, 0.0f, 0.0f), dir, spheres, lights);
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