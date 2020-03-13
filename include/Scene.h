#ifndef __SCENE__
#define __SCENE__

#include <string>
#include <iostream>

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
public:
    Scene(const std::string &description) : description(description) {}
    std::string getDescription() const { return description; }
    virtual int run(const Settings &s) const = 0;
};

#endif