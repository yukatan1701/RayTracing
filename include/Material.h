#ifndef __MATERIAL__
#define __MATERIAL__

#include "glm/vec3.hpp"
#include "glm/vec4.hpp"
#include <vector>
#include <string>

struct Material {
    float refractive;
    glm::vec4 albedo; //( яркость, матовость, отражающая способность, прозрачность)
    glm::vec3 diffuse;
    float specular;
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
};

#endif