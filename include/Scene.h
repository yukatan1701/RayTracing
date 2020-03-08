#ifndef __SCENE__
#define __SCENE__

#include <string>
#include <iostream>
#include "Scene1.h"
#include "Scene2.h"
#include "Scene3.h"

struct Settings {
    std::string out;
    unsigned int scene;
    unsigned int threads;
    Settings() : out("rt.bmp"), scene(1), threads(1) {}
    Settings(std::string out, unsigned int scene, unsigned int threads) :
        out(out), scene(scene), threads(threads) {}
};

class Scene {
    std::string description;
public:
    Scene(const std::string &description) : description(description) {}
    std::string getDescription() const { return description; }
    virtual void run(const Settings &s) const = 0;
};

#endif