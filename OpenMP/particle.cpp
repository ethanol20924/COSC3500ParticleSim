#include "particle.hpp"

#include <iostream>

Particle::Particle() {
    this->x = 0.0f;
    this->y = 0.0f;
    this->dx = 0.0f;
    this->dy = 0.0f;
    this->radius = 0.0f;
    this->mass = 0.0f;
}

Particle::Particle(float x, float y, float radius, float mass) {
    this->x = x;
    this->y = y;
    this->dx = 0.0f;
    this->dy = 0.0f;
    this->radius = radius;
    this->mass = mass;
}

Particle::Particle(float x, float y, float dx, float dy, float radius, float mass) {
    this->x = x;
    this->y = y;
    this->dx = dx;
    this->dy = dy;
    this->radius = radius;
    this->mass = mass;
}

void Particle::update(float timeStep) {
    x = x + dx * timeStep;
    y = y + dy * timeStep;
}