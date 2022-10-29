#ifndef PARTICLE_H
#define PARTICLE_H

#include <cmath>

/**
 * @brief Particle class to represent a single particle.
 * 
 */
class Particle {
    private:
        float x;  // m
        float y;  // m
        float dx;  // m/s
        float dy;  // m/s
        float radius;  // m
        float mass;  // kg

    public:
        Particle();
        Particle(float x, float y, float radius, float mass);  // Construct with initial position
        Particle(float x, float y, float dx, float dy, float radius, float mass);  // Construct with initial position and velocity

        int particleNum;
        bool hasCollided = false;

        float get_x() {return x;}
        float get_y() {return y;}
        float get_dx() {return dx;}
        float get_dy() {return dy;}
        float get_mass() { return mass;}
        float get_radius() { return radius;}

        void set_dx(float dx) { this->dx = dx; }
        void set_dy(float dy) { this->dy = dy; }

        void update(float timeStep);
        
        template<class Archive>
        void serialize(Archive &archive) { archive(x, y, radius); }
};

#endif