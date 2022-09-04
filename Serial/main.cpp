#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>
using namespace std;

#include "particles.hpp"
#include "libs/cereal/archives/json.hpp"

// SIM SETUP PARAMS
#define NUM_PARTICLES 150
#define BOX_WIDTH 1  // m
#define BOX_HEIGHT 1  // m
#define TIMESTEP 0.01  // s

#define SIM_TIME 5  // s

#define PARTICLE_SIZE 0.01  // m
#define PARTICLE_MASS 0.1  // kg
#define MAX_PARTICLE_SPEED 0.25  // m/s

#define OUTPUT_FILENAME "../out/test2.json"

int main() {
    remove(OUTPUT_FILENAME);

    auto progStart = chrono::steady_clock::now();
    auto cumulativeSimTime = chrono::duration<double, milli>(0);
    uint frames = 0;

    Particles *particles = new Particles(NUM_PARTICLES, BOX_WIDTH, BOX_HEIGHT, PARTICLE_SIZE, PARTICLE_MASS, MAX_PARTICLE_SPEED, TIMESTEP);

    {    
        ofstream file;
        file.open(OUTPUT_FILENAME, fstream::out | fstream::app | ios::binary);
        cereal::JSONOutputArchive oarchive(file);

        for (float time = 0.0f; time < SIM_TIME; time += TIMESTEP) {
            auto frameStart = chrono::steady_clock::now();

            particles->updateCollisions();
            particles->updateMovements();
            particles->updateTime();

            auto frameEnd = chrono::steady_clock::now();
            auto frameTime = frameEnd - frameStart;
            cumulativeSimTime += chrono::duration<double, milli>(frameTime);
            frames++;

            oarchive(*particles);
        }

        file.close();
    }

    auto progEnd = chrono::steady_clock::now();
    auto totalTime = progEnd - progStart;

    cout << "===== Performance Profiling =====" << endl;
    cout << "----- Simulation Settings -----" << endl;
    cout << "Number of particles: " << NUM_PARTICLES << endl;
    cout << "Width: " << BOX_WIDTH << "m" << endl;
    cout << "Height: " << BOX_HEIGHT << "m" << endl;
    cout << "Timestep: " << TIMESTEP << "s" << endl;
    cout << "Simulation time: " << SIM_TIME << "s" << endl;
    cout << "Total frames: " << frames << endl;
    cout << "----- Timings -----" << endl;
    cout << "Total execution time: " << chrono::duration<double, milli>(totalTime).count() << "ms" << endl;
    cout << "Simulation computation time: " << cumulativeSimTime.count() << "ms" << endl;
    cout << "Average frame time: " << cumulativeSimTime.count() / frames << "ms" << endl;
    
    return 0;
}