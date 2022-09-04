#ifndef PARTICLE_H
#define PARTICLE_H

#include <cmath>
#define _USE_MATH_DEFINES

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

        float get_x() {return x;}
        float get_y() {return y;}
        float get_dx() {return dx;}
        float get_dy() {return dy;}
        float get_mass() { return mass;}
        float get_radius() { return radius;}
        float get_speed() { return sqrtf(powf(dx, 2.0) + powf(dy, 2.0)); }

        void set_dx(float dx) { this->dx = dx; }
        void set_dy(float dy) { this->dy = dy; }

        void update(float timeStep);
        
        template<class Archive>
        void serialize(Archive &archive) { archive(x, y, radius); }
};

#endif