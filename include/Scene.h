#ifndef __SCENE__
#define __SCENE__

#include <string>
#include <iostream>
#include "Material.h"
#include "SceneObjects.h"

struct Settings {
    std::string out;
    unsigned int scene;
    unsigned int threads;
    Settings() : out("rt.bmp"), scene(1), threads(8) {}
    Settings(std::string out, unsigned int scene, unsigned int threads) :
        out(out), scene(scene), threads(threads) {}
};

class Scene {
protected:
    std::string description;
    vec3 reflect(const vec3 &i, const vec3 &n);
    vec3 refract(const vec3 &I, const vec3 &N, const float etat, const float etai = 1.0f);
public:
    Scene(const std::string &description) : description(description) {}
    std::string getDescription() const { return description; }
    virtual int run(const Settings &s) = 0;
};

#endif