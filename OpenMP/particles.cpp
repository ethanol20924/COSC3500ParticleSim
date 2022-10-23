#include "particles.hpp"

#include <iostream>

Particles::Particles() {
    this->timeStep = 0.0f;
}

/**
 * @brief Construct a new Particles:: Particles object. Will generate Particle::Particle objects at random locations without overlapping, with a random initial velocity.
 */
Particles::Particles(uint number, float width, float height, float particleSize, float particleMass, float maxSpeed, float timeStep) {
    this->timeStep = timeStep;
    this->width = width;
    this->height = height;
    
    for (uint i = 0; i < number; i++) {
        bool intersecting = true;
        while (intersecting) {
            float x = particleSize + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(width - 2 * particleSize)));
            float y = particleSize + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(height - 2 * particleSize)));
            float dx = -maxSpeed + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2 * maxSpeed)));
            float dy = -maxSpeed + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2 * maxSpeed)));
            Particle *newParticle = new Particle(x, y, dx, dy, particleSize, particleMass);

            if (particles.empty()) {
                intersecting = false;
                particles.push_back(*newParticle);
            } else {
                bool isIntersecting = false;
                vector<Particle>::iterator it;
                for (it = particles.begin(); it != particles.end(); it++) {
                    if (checkCollision(newParticle, &(*it))) {
                        isIntersecting = true;
                    }
                }
                if (!isIntersecting) {
                    intersecting = false;
                    particles.push_back(*newParticle);
                } else {
                    delete newParticle;
                }
            }
        }
    }
}

/**
 * @brief Construct a new Particles:: Particles object
 * 
 * @param config Config structure
 */
Particles::Particles(SimConfig_t *config) {
    this->timeStep = config->timeStep;
    this->width = config->simWidth;
    this->height = config->simHeight;

    uint numParticles = config->numParticles;
    float particleSize = config->particleSize;
    float maxSpeed = config->maxSpeed;

    for (uint i = 0; i < numParticles; i++) {
        bool intersecting = true;
        while (intersecting) {
            float x = particleSize + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(width - 2 * particleSize)));
            float y = particleSize + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(height - 2 * particleSize)));
            float dx = -maxSpeed + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2 * maxSpeed)));
            float dy = -maxSpeed + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2 * maxSpeed)));
            Particle *newParticle = new Particle(x, y, dx, dy, particleSize, particleMass);

            if (particles.empty()) {
                intersecting = false;
                particles.push_back(*newParticle);
            } else {
                bool isIntersecting = false;
                vector<Particle>::iterator it;
                for (it = particles.begin(); it != particles.end(); it++) {
                    if (checkCollision(newParticle, &(*it))) {
                        isIntersecting = true;
                    }
                }
                if (!isIntersecting) {
                    intersecting = false;
                    particles.push_back(*newParticle);
                } else {
                    delete newParticle;
                }
            }
        }
    }
}

bool Particles::checkCollision(Particle *p1, Particle *p2) {
    // AABB collision check
    if (p1->get_x() + p1->get_radius() + p2->get_radius() > p2->get_x()
    && p1->get_x() < p2->get_x() + p1->get_radius() + p2->get_radius()
    && p1->get_y() + p1->get_radius() + p2->get_radius() > p2->get_y()
    && p1->get_y() < p2->get_y() + p1->get_radius() + p2->get_radius()) {
        // Actual collision check
        float distance = sqrtf(
            ((p1->get_x() - p2->get_x()) * (p1->get_x() - p2->get_x()))
            + ((p1->get_y() - p2->get_y()) * (p1->get_y() - p2->get_y()))
        );
        if (distance < p1->get_radius() + p2->get_radius()) {
            return true;
        }
    }
    return false;
}