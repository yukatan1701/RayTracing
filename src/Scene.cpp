#include "Scene.h"

vec3 Scene::reflect(const vec3 &i, const vec3 &n) {
    return i - n * 2.0f * dot(i, n);
}

vec3 Scene::refract(const vec3 &I, const vec3 &N, const float etat, const float etai) {
    float cosi = -std::max(-1.0f, std::min(1.0f, dot(I, N)));
    if (cosi < 0)
        return refract(I, -N, etai, etat);
    float eta = etai / etat;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? vec3(1.0f) : I * eta + N * (eta * cosi - sqrtf(k));
}