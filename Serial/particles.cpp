#include "particles.hpp"

#include <iostream>

Particles::Particles() {
    this->timeStep = 0.0f;
}

/**
 * @brief Construct a new Particles:: Particles object. Will generate Particle:: Particle objects at random locations without overlapping, with a random initial velocity.
 * 
 * @param number Number of particles to generate
 * @param width Width of the simulation space
 * @param height Height of the simulation space
 * @param particleSize Size of each particle
 * @param particleMass Mass of each particle
 * @param maxSpeed The maximum initial speed
 * @param timeStep The time step of the simulation
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
 * @brief Updates the movements of particles that are colliding.
 * 
 */
void Particles::updateCollisions() {
    vector<Particle>::iterator it1;
    vector<Particle>::iterator it2;
    for (it1 = particles.begin(); it1 != particles.end(); it1++) {
        for (it2 = particles.begin(); it2 != particles.end(); it2++) {
            // Check that particles are different
            if (it1 != it2) {
                if (checkCollision(&(*it1), &(*it2))) {
                    // https://gamedevelopment.tutsplus.com/tutorials/when-worlds-collide-simulating-circle-circle-collisions--gamedev-769
                    float newXVel1 = (it1->get_dx() * (it1->get_mass() - it2->get_mass()) + (2 * it2->get_mass() * it2->get_dx())) / (it1->get_mass() + it2->get_mass());
                    float newXVel2 = (it2->get_dx() * (it2->get_mass() - it1->get_mass()) + (2 * it1->get_mass() * it1->get_dx())) / (it1->get_mass() + it2->get_mass());
                    float newYVel1 = (it1->get_dy() * (it1->get_mass() - it2->get_mass()) + (2 * it2->get_mass() * it2->get_dy())) / (it1->get_mass() + it2->get_mass());
                    float newYVel2 = (it2->get_dy() * (it2->get_mass() - it1->get_mass()) + (2 * it1->get_mass() * it1->get_dy())) / (it1->get_mass() + it2->get_mass());

                    it1->set_dx(newXVel1);
                    it2->set_dx(newXVel2);
                    it1->set_dy(newYVel1);
                    it2->set_dy(newYVel2);

                    // Update here to ensure particles don't double update
                    it1->update(timeStep);
                    it2->update(timeStep);
                }
            }
        }
    }
}

/**
 * @brief Updates the positions of particles. Also checks if particles should bounce off a wall.
 * 
 */
void Particles::updateMovements() {
    vector<Particle>::iterator it;
    for (it = particles.begin(); it != particles.end(); it++) {
        float x = it->get_x();
        float y = it->get_y();
        float radius = it->get_radius();

        if (x < radius || x > width - radius) {
            it->set_dx(-1 * it->get_dx());
        }
        if (y < radius || y > width - radius) {
            it->set_dy(-1 * it->get_dy());
        }

        it->update(timeStep);
    }
}

/**
 * @brief Increments the simulation by updating the current time.
 * 
 */
void Particles::updateTime() {
    currentTime += timeStep;
}

/**
 * @brief Checks for collisions between two particles
 * 
 * @param p1 Pointer to first particle
 * @param p2 Pointer to second particle
 * @return true 
 * @return false 
 */
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