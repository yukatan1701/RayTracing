#ifndef __SCENE2__
#define __SCENE2__

#include "Scene.h"

class Scene2 : public Scene {
public:
    Scene2(const std::string &description = "description") : Scene(description) {}
    int run(const Settings &s);
};

#endif