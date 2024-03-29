#ifndef PARTICLES_H
#define PARTICLES_H

#include "particle.hpp"
#include "libs/cereal/types/vector.hpp"

#include <vector>
#include <random>
#include <fstream>
using namespace std;

/**
 * @brief Particles class to represent a collection of particles. Handles the calcuations and updating of each particle.
 * 
 */
class Particles {
    private:
        vector<Particle> particles;

        float timeStep;
        float currentTime = 0;

        float width;
        float height;

    public:
        Particles();
        Particles(uint number, float width, float height, float particleSize, float particleMass, float maxSpeed, float timeStep);
        
        void updateCollisions();
        void updateMovements();
        void updateTime();

        bool checkCollision(Particle *p1, Particle *p2);

        vector<Particle> getParticles() { return particles; }
        
        template <class Archive>
        void serialize(Archive &archive) {
            archive(CEREAL_NVP(currentTime), CEREAL_NVP(particles));
        }
};

#endif // PARTICLES_H