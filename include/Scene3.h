#ifndef __SCENE3__
#define __SCENE3__

#include "Scene.h"

class Scene3 : public Scene {
public:
    Scene3(const std::string &description) : Scene(description) {}
    void run(const Settings &s) const;
};

#endif