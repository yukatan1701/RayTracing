#ifndef __SCENE1__
#define __SCENE1__

#include "Scene.h"
#include "Material.h"
#include <cmath>

class Scene1 : public Scene {
private:
    const Materials materials;
    const int width = 800;
    const int height = 600;
    const int refractDepth = 4;
    const float fprec = 0.001f;

    int backgroundWidth, backgroundHeight;
    std::vector<vec3> background;
    bool sceneIntersect(const vec3 &orig, const vec3 &dir,
                    const SceneObjects &sceneObjects,
                    vec3 &hit, vec3 &N, Material &material);
    vec3 traceRay(const vec3 &orig, const vec3 &dir, const SceneObjects &sceneObjects,
              size_t depth);
    objset<Sphere *> loadSpheres();
    objset<Light *> loadLights();
    objset<Triangle *> loadTriangles();
    objset<Model *> loadModels();
    objset<Cube *> loadCubes();
    void deleteSceneObjects(SceneObjects &so);
    inline bool isNegative(const float &value) { return std::signbit(value); }
    void render(const Settings &settings);
public:
    Scene1(const std::string &description = "description") : Scene(description) {}
    int run(const Settings &s);
};

#endif