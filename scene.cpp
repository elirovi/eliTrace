#include "scene.h"
#include "scene_types.h"
#include <string.h>
#include <algorithm>

Object *initSphere(point3 center, float radius, Material mat) {
    Object *ret;
    ret = (Object *)malloc(sizeof(Object));
    ret->geom.type = SPHERE;
    ret->geom.sphere.center = center;
    ret->geom.sphere.radius = radius;
    memcpy(&(ret->mat), &mat, sizeof(Material));
    return ret;
}

Object *initPlane(vec3 normal, float d, Material mat) {
    Object *ret;
    ret = (Object *)malloc(sizeof(Object));
    ret->geom.type = PLANE;
    ret->geom.plane.normal = normalize(normal);
    ret->geom.plane.dist = d;
    memcpy(&(ret->mat), &mat, sizeof(Material));
    return ret;
}

Object* initCylinder(vec3 direction, point3 center, float radius, float length, Material mat){
    Object *ret;
    ret=(Object *)malloc(sizeof(Object));
    ret->geom.type=CYLINDER;
    ret->geom.cylinder.direction= normalize(direction);
    ret->geom.cylinder.center= center;
    ret->geom.cylinder.radius= radius;
    ret->geom.cylinder.length= length;
    memcpy(&(ret->mat), &mat, sizeof(Material));
    return ret;
}

Object* initCone(point3 top, point3 base, float alpha, Material mat){
    Object *ret;
    ret=(Object*)malloc(sizeof(Object));
    ret->geom.type=CONE;
    ret->geom.cone.top=top;
    ret->geom.cone.base=base;
    ret->geom.cone.alpha=alpha;
    memcpy(&(ret->mat),&mat,sizeof(Material));
    return ret;
}

void freeObject(Object *obj) {
    free(obj);
}

Light *initLight(point3 position, color3 color) {
    Light *light = (Light*)malloc(sizeof(Light));
    light->position = position;
    light->color = color;
    return light;
}

void freeLight(Light *light) {
    free(light);
}

Scene * initScene() {
    return new Scene;
}

void freeScene(Scene *scene) {
    std::for_each(scene->objects.begin(), scene->objects.end(), freeObject);
    std::for_each(scene->lights.begin(), scene->lights.end(), freeLight);
    delete scene;
}

void setCamera(Scene *scene, point3 position, point3 at, vec3 up, float fov, float aspect) {
    scene->cam.fov = fov;
    scene->cam.aspect = aspect;
    scene->cam.position = position;
    scene->cam.zdir = normalize(at-position);
    scene->cam.xdir = normalize(cross(up, scene->cam.zdir));
    scene->cam.ydir = normalize(cross(scene->cam.zdir, scene->cam.xdir));
    scene->cam.center = 1.f / tanf ((scene->cam.fov * M_PI / 180.f) * 0.5f) * scene->cam.zdir;
}

void addObject(Scene *scene, Object *obj) {
    scene->objects.push_back(obj);
}

void addLight(Scene *scene, Light *light) {
    scene->lights.push_back(light);
}

void setSkyColor(Scene *scene, color3 c) {
    scene->skyColor = c;
}
