#ifndef PARTICLES_H
#define PARTICLES_H

#include "particle.hpp"
#include "libs/cereal/types/vector.hpp"

#include <vector>
#include <random>
#include <fstream>
using namespace std;

#define MAX_PARTICLES_PER_CELL 4
#define NUM_THREADS 1
#define DATA_PER_PARTICLE 4

#define DEBUG_GRID false
#define DEBUG_COLLISION false
#define DEBUG_MOVEMENT false

/**
 * @brief Config structure for simulation initialisation
 * 
 */
struct SimConfig_t {
    uint numParticles;
    float simWidth;
    float simHeight;
    float particleSize;
    float particleMass;
    float minSpeed;
    float maxSpeed;
    float timeStep;
};

/**
 * @brief Particles class to represent a collection of particles. Handles the calcuations and updating of each particle.
 * 
 */
class Particles {
    private:
        uint *blankCounters;  // Host side blank counter
        uint *gridCounters;  // array storing the number of particles in a cell

        uint *blankCells;  // Host side cells
        uint *gridCells;  // array storing the indexes of the particles in each cell, max 4 particles per cell

        vector<Particle> particles;
        float *temp_particles;  // same as below
        float *d_particles;  // format: x, y, dx, dy

        float *newVels;  // array for holding new velocities
        float *blankVels;  // blanks

        bool *inGridFlags;  // used for allocating particle to grid
        bool *blankFlags;

        uint numParticles;
        float particleMass;
        float particleSize;

        float timeStep;
        float currentTime = 0;

        uint numRows;
        uint numCols;

        float width;
        float height;

        float gridSize;

    public:
        Particles();
        Particles(SimConfig_t *config);

        ~Particles();

        void updateTime();
        void updateGrid();
        void updateCollisions();
        void updateMovements();

        static inline bool checkAABBCircle(float p1x, float p1y, float p1r, float p2x, float p2y, float p2r);
        static inline bool checkAABBRect(float p1x, float p1y, float p1w, float p1h, float p2x, float p2y, float p2w, float p2h);
        bool checkCollision(Particle *p1, Particle *p2);

        template <class Archive>
        void serialize(Archive &archive) {
            archive(CEREAL_NVP(currentTime), CEREAL_NVP(particles));
        }
};

#endif // PARTICLES_H