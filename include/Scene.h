#ifndef __SCENE__
#define __SCENE__
#include "includes.h"
#include "Material.h"
#include "SceneObjects.h"
#include <cmath>

struct Settings {
    std::string out;
    unsigned int scene;
    unsigned int threads;
    Settings() : out("rt.bmp"), scene(1), threads(8) {}
    Settings(std::string out, unsigned int scene, unsigned int threads) :
        out(out), scene(scene), threads(threads) {}
};

class Scene {
private:
    const Materials materials;
    const int width = 800;
    const int height = 600;
    const int refractDepth = 4;
    const float fprec = 0.001f;

    int waterWidth, waterHeight;
    std::vector<vec3> background, water;
    Lightning calcLightning(const vec3 &orig, const vec3 &dir, const vec3 &point, const vec3 &N,
        const Material &mt, const SceneObjects &sceneObjects);
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
    objset<Island *> loadIslands();
    void deleteSceneObjects(SceneObjects &so);
    inline bool isNegative(const float &value) const { return std::signbit(value); }
    void render(const Settings &settings);
public:
    Scene() {}
    int run(const Settings &s);
};

#endif