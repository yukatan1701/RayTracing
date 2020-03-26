#ifndef __SCENE1__
#define __SCENE1__

#include "Scene.h"

class Scene1 : public Scene {
private:
    bool sceneIntersect(const vec3 &orig, const vec3 &dir,
                    const SceneObjects &sceneObjects,
                    vec3 &hit, vec3 &N, Material &material);
    vec3 traceRay(const vec3 &orig, const vec3 &dir, const SceneObjects &sceneObjects,
              size_t depth);
    std::deque<Sphere> loadSpheres(const Materials &materials);
    std::deque<Light> loadLights();
    std::deque<Triangle> loadTriangles(const Materials &materials);
    std::deque<Model> loadModels(const Materials &materials);
    void render(const Settings &settings);
public:
    Scene1(const std::string &description = "description") : Scene(description) {}
    int run(const Settings &s);
};

#endif