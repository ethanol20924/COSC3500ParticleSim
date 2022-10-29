#ifndef PARTICLES_H
#define PARTICLES_H

#include "particle.hpp"
#include "libs/cereal/types/vector.hpp"

#include <vector>
#include <random>
#include <fstream>
using namespace std;

/**
 * @brief Config structure for simulation initialisation
 * 
 */
struct SimConfig_t {
    uint numParticles;
    float simWidth;
    float simHeight;
    uint gridWidth;  // not actually width, it's a count
    uint gridHeight;
    float particleSize;
    float particleMass;
    float maxSpeed;
    float timeStep;
};

struct Cell {
    uint cellNum;
    vector<uint> elements;  // Vector of particles overlapping the cell. Contains the index of the particle in the 'particles' vector
};

/**
 * @brief Structure representing each grid row
 * 
 */
struct GridRow {
    uint rowNum;
    vector<Cell> cells;  // Stores all the cells in said row
};

/**
 * @brief Particles class to represent a collection of particles. Handles the calcuations and updating of each particle.
 * 
 */
class Particles {
    private:
        vector<Particle> particles;  // Vector of particles. We will use this to prevent doubled up collisions
        vector<GridRow> rows;  // Vector of rows
        vector<GridRow> buffer;  // Dirty buffer to allow for parallelism

        uint numParticles;

        float timeStep;
        float currentTime = 0;

        float width;
        float height;

        uint numRows;
        uint numCols;

        float gridWidth;
        float gridHeight;

    public:
        Particles();
        Particles(SimConfig_t *config);

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