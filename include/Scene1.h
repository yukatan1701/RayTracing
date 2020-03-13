#ifndef __SCENE1__
#define __SCENE1__

#include "Scene.h"

class Scene1 : public Scene {
public:
    Scene1(const std::string &description = "description") : Scene(description) {}
    int run(const Settings &s) const;
};

#endif