#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>
using namespace std;

#include "particles.hpp"
#include "libs/cereal/archives/json.hpp"

// INITIAL SIM SETUP PARAMS
#define NUM_PARTICLES 100
#define BOX_WIDTH 2  // m
#define BOX_HEIGHT 2  // m
#define NUM_ROWS 100  // m
#define TIMESTEP 0.01  // s

#define SIM_TIME 5  // s

#define PARTICLE_SIZE 0.001  // m
#define PARTICLE_MASS 0.1  // kg
#define MAX_PARTICLE_SPEED 0.25  // m/s

#define OUTPUT_ENABLED true  // Enable or disable JSON serialisation
#define OUTPUT_FILENAME "../out/test.json"

#define COUT_TO_FILE false  // Note: we can specify outputs with SLURM script
#define PROFILE_FILENAME "../out/profile.txt"

// REPEATED SIM SETUPS
#define MULTIPLE_SIMS false  // Enable batch running multiple simulations one after another
#define NUM_NEW_PARTICLES 500  // Number of particles to add after each sim
#define MAX_PARTICLES 1000  // Maximum number of particles allowed
#define REPEAT_SIM 5  // Number of time to run the same sim for averaging

int main() {
    #if OUTPUT_ENABLED
    remove(OUTPUT_FILENAME);
    #endif

    #if COUT_TO_FILE
    remove(PROFILE_FILENAME);
    #endif

    uint max_particles = MULTIPLE_SIMS ? MAX_PARTICLES : NUM_PARTICLES;

    #if MULTIPLE_SIMS
    cout << "<<<<< BEGIN BATCH SIMULATION >>>>>" << endl;
    cout << "Initial particle count: " << NUM_PARTICLES << endl;
    cout << "Maximum particle count: " << MAX_PARTICLES << endl;
    cout << "Samples per particle count: " << REPEAT_SIM << endl;
	cout << endl;
    #endif

    for (uint num_particles = NUM_PARTICLES; num_particles <= max_particles; num_particles += NUM_NEW_PARTICLES) {

		auto averageTotalTime = chrono::duration<double, milli>(0);
		auto averageSimTime = chrono::duration<double, milli>(0);

		for (uint i = 0; i < REPEAT_SIM; i++) {

			auto progStart = chrono::steady_clock::now();
			auto cumulativeSimTime = chrono::duration<double, milli>(0);
			uint frames = 0;

			auto simConfig = new SimConfig_t;
			simConfig->numParticles = num_particles;
			simConfig->simWidth = BOX_WIDTH;
			simConfig->simHeight = BOX_HEIGHT;
			simConfig->particleSize = PARTICLE_SIZE;
			simConfig->particleMass = PARTICLE_MASS;
			simConfig->gridWidth = NUM_ROWS;
			simConfig->gridHeight = NUM_ROWS;
			simConfig->maxSpeed = MAX_PARTICLE_SPEED;
			simConfig->timeStep = TIMESTEP;

			Particles *particles = new Particles(simConfig);

			{
			#if OUTPUT_ENABLED
			ofstream file;
			file.open(OUTPUT_FILENAME, fstream::out | fstream::app | ios::binary);
			cereal::JSONOutputArchive oarchive(file);
			#endif

			for (float time = 0.0f; time < SIM_TIME; time += TIMESTEP) {
				auto frameStart = chrono::steady_clock::now();

				particles->updateGrid();
				particles->updateCollisions();
				particles->updateMovements();
				particles->updateTime();

				// Timing is done here as we don't really want to time the serialisation part
				auto frameEnd = chrono::steady_clock::now();
				auto frameTime = frameEnd - frameStart;
				cumulativeSimTime += chrono::duration<double, milli>(frameTime);
				frames++;

				#if OUTPUT_ENABLED
				oarchive(*particles);
				#endif
			}
			#if OUTPUT_ENABLED
			file.close();
			#endif
			}

			auto progEnd = chrono::steady_clock::now();
			auto totalTime = progEnd - progStart;

			averageTotalTime += chrono::duration<double, milli>(totalTime);
			averageSimTime += cumulativeSimTime;

		}

		averageTotalTime /= REPEAT_SIM;
		averageSimTime /= REPEAT_SIM;

        #if COUT_TO_FILE
        freopen(PROFILE_FILENAME,"a",stdout);
        #endif

        cout << "===== Performance Profiling =====" << endl;
        cout << "----- Simulation Settings -----" << endl;
        cout << "Number of particles: " << num_particles << endl;
        cout << "Width: " << BOX_WIDTH << "m" << endl;
        cout << "Height: " << BOX_HEIGHT << "m" << endl;
        cout << "Grid Width: " << static_cast<float>(BOX_WIDTH) / static_cast<float>(NUM_ROWS) << "m" << endl;
        cout << "Grid Height: " << static_cast<float>(BOX_HEIGHT) / static_cast<float>(NUM_ROWS) << "m" << endl;
        cout << "Timestep: " << TIMESTEP << "s" << endl;
        cout << "Simulation time: " << SIM_TIME << "s" << endl;
        cout << "Frames per run: " << (SIM_TIME / TIMESTEP) << endl;
        cout << "----- Timings -----" << endl;
        cout << "Average execution time: " << averageTotalTime.count() << "ms" << endl;
        cout << "Average simulation computation time: " << averageSimTime.count() << "ms" << endl;
        cout << "Average frame time: " << averageSimTime.count() / (SIM_TIME / TIMESTEP) << "ms" << endl;
        cout << "===== End Profiling =====" << endl;
        cout << endl;

    }

    #if MULTIPLE_SIMS
    cout << "<<<<< END BATCH SIMULATION >>>>>" << endl;
    #endif
    
    return 0;
}
