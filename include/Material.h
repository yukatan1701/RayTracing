#ifndef __MATERIAL__
#define __MATERIAL__

#include "includes.h"
#include <vector>
#include <string>

struct Material {
    glm::vec4 albedo; // brightness, smoothness, reflectivity, refractivity
    glm::vec3 diffuse;
    float specular;
    float refractive;
    std::string name;
    Material() : refractive(1), albedo(1, 0, 0, 0), diffuse(), specular(), name("none") {}
    Material(const float &r, const glm::vec4 &a, const glm::vec3 &color, const float &spec,
             const std::string &name) : refractive(r), albedo(a), diffuse(color),
             specular(spec), name(name) {}
};

class Materials {
private:
    std::vector<Material> materials;
public:
    Materials();
    Material get(const std::string &name) const;
    Material operator [](const std::string &name) const { return get(name); }
};

#endif