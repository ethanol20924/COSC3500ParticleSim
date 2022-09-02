#include "particles.hpp"

#include <iostream>

Particles::Particles() {
    this->timeStep = 0.0f;
}

Particles::Particles(uint number, float width, float height, float particleSize, float particleMass, float maxSpeed, float timeStep) {
    this->timeStep = timeStep;
    
    for (uint i = 0; i < number; i++) {
        bool intersecting = true;
        while (intersecting) {
            float x = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/width));
            float y = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/height));
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

void Particles::updateCollisions() {
    // https://gamedevelopment.tutsplus.com/tutorials/when-worlds-collide-simulating-circle-circle-collisions--gamedev-769

    collisions.clear();  // Clear the array, don't want to deallocate memory for the particles

    // Iterate over everything lol
    vector<Particle>::iterator it1;
    vector<Particle>::iterator it2;
    for (it1 = particles.begin(); it1 != particles.end(); it1++) {
        for (it2 = particles.begin(); it2 != particles.end(); it2++) {
            // Check that particles are different
            if (it1 != it2) {
                if (checkCollision(&(*it1), &(*it2))) {
                    it1->set_dx(it1->get_dx() * (it1->get_mass() - it2->get_mass()) + (2 * it2->get_mass() * it2->get_dx()) / (it1->get_mass() + it2->get_mass()));
                    it1->set_dy(it1->get_dy() * (it1->get_mass() - it2->get_mass()) + (2 * it2->get_mass() * it2->get_dy()) / (it1->get_mass() + it2->get_mass()));
                }
            }
        }
    }
}

void Particles::updateMovements() {
    // https://gamedevelopment.tutsplus.com/tutorials/when-worlds-collide-simulating-circle-circle-collisions--gamedev-769

    vector<Particle>::iterator it1;
    for (it1 = particles.begin(); it1 != particles.end(); it1++) {
        it1->update(timeStep);
    }
}

void Particles::updateTime() {
    currentTime += timeStep;
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

CollideEvent::CollideEvent(Particle *affecting, Particle *affected) {
    affectingParticle = affecting;
    affectedParticle = affected;
}