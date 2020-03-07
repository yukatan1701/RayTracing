#ifndef __MATERIAL__
#define __MATERIAL__

#include "glm/vec3.hpp"
#include "glm/vec4.hpp"
#include <vector>
#include <string>

struct Material {
    glm::vec4 albedo;
    glm::vec3 diffuseColor;
    float refractiveIndex;
    float specularExponent;
    std::string name;
    Material() : refractiveIndex(1), albedo(1, 0, 0, 0), diffuseColor(), specularExponent(), name("none") {}
    Material(const float &r, const glm::vec4 &a, const glm::vec3 &color, const float &spec,
             const std::string &name) : refractiveIndex(r), albedo(a), diffuseColor(color),
             specularExponent(spec), name(name) {}
};

class Materials {
private:
    std::vector<Material> materials;
public:
    Materials();
    const Material& get(const std::string &name);
};

#endif