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
                    // cout << "Before: " << it1->get_dx() << " | " << it2->get_dx() << "\n";

                    float theta1 = atan2f(it1->get_dy(), it1->get_dx());
                    float theta2 = atan2f(it2->get_dy(), it2->get_dx());
                    float phi = atan2f(it1->get_dy() - it2->get_dy(), it1->get_dx() - it2->get_dx());
                    float speed1 = it1->get_speed();
                    float speed2 = it2->get_speed();

                    float newXVel1 = ((speed1 * cosf(theta1 - phi) * (it1->get_mass() - it2->get_mass()) + 2 * it2->get_mass() * speed2 * cosf(theta2 - phi)) / (it1->get_mass() + it2->get_mass())) * cosf(phi)
                                    + speed1 * sinf(theta1 - phi) * cosf(phi + M_PI / 2);
                    float newYVel1 = ((speed1 * cosf(theta1 - phi) * (it1->get_mass() - it2->get_mass()) + 2 * it2->get_mass() * speed2 * cosf(theta2 - phi)) / (it1->get_mass() + it2->get_mass())) * sinf(phi)
                                    + speed1 * sinf(theta1 - phi) * sinf(phi + M_PI / 2);
                    float newXVel2 = ((speed2 * cosf(theta2 - phi) * (it2->get_mass() - it1->get_mass()) + 2 * it1->get_mass() * speed1 * cosf(theta1 - phi)) / (it1->get_mass() + it2->get_mass())) * cosf(phi)
                                    + speed2 * sinf(theta2 - phi) * cosf(phi + M_PI / 2);
                    float newYVel2 = ((speed2 * cosf(theta2 - phi) * (it2->get_mass() - it1->get_mass()) + 2 * it1->get_mass() * speed1 * cosf(theta1 - phi)) / (it1->get_mass() + it2->get_mass())) * sinf(phi)
                                    + speed2 * sinf(theta2 - phi) * sinf(phi + M_PI / 2);

                    it1->set_dx(newXVel1);
                    it1->set_dy(newYVel1);
                    it2->set_dx(newXVel2);
                    it2->set_dy(newYVel2);

                    it1->update(timeStep);
                    it2->update(timeStep);

                    // cout << "After: " << it1->get_dx() << " | " << it2->get_dx() << "\n";
                }
            }
        }
    }
}

void Particles::updateMovements() {
    // https://gamedevelopment.tutsplus.com/tutorials/when-worlds-collide-simulating-circle-circle-collisions--gamedev-769

    vector<Particle>::iterator p_it;
    for (p_it = particles.begin(); p_it != particles.end(); p_it++) {
        p_it->update(timeStep);
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