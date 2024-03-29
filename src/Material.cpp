
#include "Material.h"
#include <iostream>

using namespace glm;

Materials::Materials() {
    materials.push_back(Material());
    materials.push_back(Material(1.5, vec4(0.0,  0.5, 0.1, 0.8), vec3(0.6, 0.7, 0.8),  125., "glass"));
    materials.push_back(Material(1.0, vec4(0.0, 10.0, 0.8, 0.0), vec3(1.0, 1.0, 1.0), 1425., "mirror"));
    materials.push_back(Material(1.0, vec4(0.5, 3.0, 0.2, 0.0), vec3(120/255.0f, 85/255.0f, 15/255.0f), 125, "gold"));
    materials.push_back(Material(1.0, vec4(0.5, 3.0, 0.2, 0.0), vec3(0.3, 0.3, 0.3), 125, "silver"));
    materials.push_back(Material(1.0, vec4(0.9,  0.1, 0.2, 0.0), vec3(0, 60, 179) * 0.3f / 255.0f, 10., "water_bottom"));
    materials.push_back(Material(1.0, vec4(0.9,  0.1, 0.2, 0.5), vec3(0, 43, 128) * 0.3f / 255.0f, 125., "water_top"));
    materials.push_back(Material(1.0, vec4(0.9,  0.1, 0.0, 0.0), vec3(136, 204, 0) * 0.3f / 255.0f, 10., "island"));
    materials.push_back(Material(1.0, vec4(0.9,  0.1, 0.0, 0.0), vec3(102, 51, 0) * 0.5f / 255.0f, 10., "wood"));
}

Material Materials::get(const std::string &name) const {
    for (auto &item : materials) {
        if (item.name == name)
            return item;
    }
    std::cerr << "Failed to find material of type '" << name <<
                    "'. Load standard material." << std::endl; 
    return Material();
}