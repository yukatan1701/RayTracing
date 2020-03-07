#include "Material.h"
#include <iostream>

using namespace glm;
/*
Materials::Materials() {
    materials.push_back(Material());
    materials.push_back(Material(1.0, vec4(0.6,  0.3, 0.1, 0.0), vec3(0.4, 0.4, 0.3),   50., "ivory"));
    materials.push_back(Material(1.5, vec4(0.0,  0.5, 0.1, 0.8), vec3(0.6, 0.7, 0.8),  125., "glass"));
    materials.push_back(Material(1.0, vec4(0.9,  0.1, 0.0, 0.0), vec3(0.3, 0.1, 0.1),   10., "red_rubber"));
    materials.push_back(Material(1.0, vec4(0.0, 10.0, 0.8, 0.0), vec3(1.0, 1.0, 1.0), 1425., "mirror"));
}

const Material& Materials::get(const std::string &name) {
    for (auto &item : materials) {
        if (item.name == name)
            return item;
    }
    std::cerr << "Failed to find material of type '" << name <<
                    "'. Load standard material." << std::endl; 
    return materials[0];
}*/