#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;

#include "particles.hpp"
#include "libs/cereal/archives/json.hpp"

// SIM SETUP PARAMS
#define NUM_PARTICLES 100
#define BOX_WIDTH 1  // m
#define BOX_HEIGHT 1  // m
#define TIMESTEP 0.01  // s

#define SIM_TIME 5  // s

#define PARTICLE_SIZE 0.01  // m
#define PARTICLE_MASS 0.1  // kg
#define MAX_PARTICLE_SPEED 0.25  // m/s

#define OUTPUT_FILENAME "../out/test1.json"

int main() {
    remove(OUTPUT_FILENAME);

    Particles *particles = new Particles(NUM_PARTICLES, BOX_WIDTH, BOX_HEIGHT, PARTICLE_SIZE, PARTICLE_MASS, MAX_PARTICLE_SPEED, TIMESTEP);

    {    
        ofstream file;
        file.open(OUTPUT_FILENAME, fstream::out | fstream::app | ios::binary);
        cereal::JSONOutputArchive oarchive(file);

        for (float time = 0.0f; time < SIM_TIME; time += TIMESTEP) {
            particles->updateCollisions();
            particles->updateMovements();
            particles->updateTime();
            oarchive(*particles);
        }

        file.close();
    }
    
    return 0;
}