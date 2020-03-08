#ifndef __SCENE2__
#define __SCENE2__

#include "Scene.h"

class Scene2 : public Scene {
public:
    Scene2(const std::string &description) : Scene(description) {}
    void run(const Settings &s) const;
};

#endif