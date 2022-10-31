#include "particles.hpp"

#include <iostream>
#include <algorithm>
#include <immintrin.h>
#include <cuda.h>

Particles::Particles() {
    this->timeStep = 0.0f;
}

/**
 * @brief Construct a new Particles:: Particles object
 * 
 * @param config Config structure
 */
Particles::Particles(SimConfig_t *config) {
    this->timeStep = config->timeStep;
    this->width = config->simWidth;
    this->height = config->simHeight;

    this->numParticles = config->numParticles;
    this->particleSize = config->particleSize;
    this->particleMass = config->particleMass;
    float maxSpeed = config->maxSpeed;

    temp_particles = new float[this->numParticles * 4]();

    for (uint i = 0; i < this->numParticles; i++) {
        bool intersecting = true;
        while (intersecting) {
            float x = this->particleSize + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(width - 2 * this->particleSize)));
            float y = this->particleSize + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(height - 2 * this->particleSize)));
            float dx = -maxSpeed + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2 * maxSpeed)));
            float dy = -maxSpeed + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2 * maxSpeed)));
            Particle *newParticle = new Particle(x, y, dx, dy, this->particleSize, this->particleMass);
            newParticle->particleNum = i;

            if (particles.empty()) {
                intersecting = false;
                particles.push_back(*newParticle);

                temp_particles[i * 4] = newParticle->get_x();
                temp_particles[i * 4 + 1] = newParticle->get_y();
                temp_particles[i * 4 + 2] = newParticle->get_dx();
                temp_particles[i * 4 + 3] = newParticle->get_dy();
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

                    temp_particles[i * 4] = newParticle->get_x();
                    temp_particles[i * 4 + 1] = newParticle->get_y();
                    temp_particles[i * 4 + 2] = newParticle->get_dx();
                    temp_particles[i * 4 + 3] = newParticle->get_dy();
                } else {
                    delete newParticle;
                }
            }
        }
    }

    // Allocate device particle storage
    cudaMalloc((void **)&d_particles, sizeof(float) * this->numParticles * 4); 
    cudaMemcpy(d_particles, temp_particles, sizeof(float) * this->numParticles * 4, cudaMemcpyHostToDevice);

    this->gridSize = config->particleSize * 2;

    this->numRows = this->height / this->gridSize;
    this->numCols = this->width / this->gridSize;

    // std::cout << this->numRows << std::endl;

    // Allocate grid
    cudaMalloc((void **)&gridCounters, sizeof(uint) * this->numRows * this->numCols);
    cudaMalloc((void **)&gridCells, sizeof(uint) * this->numRows * this->numCols * MAX_PARTICLES_PER_CELL);

    // Allocate blank grid host side
    blankCounters = new uint[this->numRows * this->numCols]();
    blankCells = new uint[this->numRows * this->numCols * MAX_PARTICLES_PER_CELL]();

    // Allocate velocity buffers
    cudaMalloc((void **)&newVels, sizeof(float) * this->numParticles * 2);
    blankVels = new float[this->numParticles * 2]();
}

/**
 * @brief Increments the simulation by updating the current time.
 * 
 */
void Particles::updateTime() {
    this->currentTime += this->timeStep;
}

__global__ void d_updateGrid(float *particles, uint *counters, uint *cells, float size, uint cols) {
    uint particle = blockIdx.x * blockDim.x + threadIdx.x;

    // Figure out which grid we're in
    uint gridX = floor(particles[particle * 4] / size);
    uint gridY = floor(particles[particle * 4 + 1] / size);

    // Add to counter and cell
    atomicInc(&counters[gridY * cols + gridX]);
    for (uint i = 0; i < 4; i++) {
        __threadfence();
        if (cells[gridY * cols * MAX_PARTICLES_PER_CELL + gridX * MAX_PARTICLES_PER_CELL + i] == 0) {
            cells[gridY * cols * MAX_PARTICLES_PER_CELL + gridX * MAX_PARTICLES_PER_CELL + i] = particle + 1;
            break;
        }
    }
}

void Particles::updateGrid() {
    cudaMemcpy(gridCounters, blankCounters, sizeof(uint) * this->numRows * this->numCols, cudaMemcpyHostToDevice);
    cudaMemcpy(gridCells, blankCells, sizeof(uint) * this->numRows * this->numCols * MAX_PARTICLES_PER_CELL, cudaMemcpyHostToDevice);

    d_updateGrid<<<(this->numParticles + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>>(d_particles, gridCounters, gridCells, this->particleSize, this->numCols);
}

__global__ void d_updateCollisions(float *particles, uint *counters, uint *cells, float *vels, float size, float mass) {
    uint particle = blockIdx.x * blockDim.x + threadIdx.x;

    // Check current cell
    float x1 = particles[particle * 4];
    float y1 = particles[particle * 4 + 1];
    uint gridX = floor(x1 / size);
    uint gridY = floor(y1 / size);

    // Variables to determine which quadrant of the cell the particle is in
    int addX = particles[particle * 4] >= size * (gridX + 0.5) ? 1 : -1;
    int addY = particles[particle * 4 + 1] >= size * (gridY + 0.5) ? 1 : -1;

    // Variables for the other particle
    uint otherParticle;
    float x2, y2;

    bool collided = false;

    __threadfence();
    if (counters[gridY * cols + gridX] > 1 && !collided) {
        for (uint i = 0; i < 4; i++) {
            otherParticle = cells[gridY * cols * MAX_PARTICLES_PER_CELL + gridX * MAX_PARTICLES_PER_CELL + i] - 1
            if (otherParticle == 0) {
                break;
            } else {
                x2 = particles[otherParticle * 4];
                y2 = particles[otherParticle * 4 + 1];
                collided = d_checkCollision(x1, y1, x2, y2, size / 2);
            }

            if (collided) {
                break;
            }
        }
    } 
    
    if (counters[gridY * cols + gridX + addX] > 1 && !collided) {
        for (uint i = 0; i < 4; i++) {
            otherParticle = cells[gridY * cols * MAX_PARTICLES_PER_CELL + (gridX + addX) * MAX_PARTICLES_PER_CELL + i] - 1
            if (otherParticle == 0) {
                break;
            } else {
                x2 = particles[otherParticle * 4];
                y2 = particles[otherParticle * 4 + 1];
                collided = d_checkCollision(x1, y1, x2, y2, size / 2);
            }

            if (collided) {
                break;
            }
        }
    }
    
    if (counters[(gridY + addY) * cols + gridX] > 1 && !collided) {
        for (uint i = 0; i < 4; i++) {
            otherParticle = cells[(gridY + addY) * cols * MAX_PARTICLES_PER_CELL + gridX * MAX_PARTICLES_PER_CELL + i] - 1
            if (otherParticle == 0) {
                break;
            } else {
                x2 = particles[otherParticle * 4];
                y2 = particles[otherParticle * 4 + 1];
                collided = d_checkCollision(x1, y1, x2, y2, size / 2);
            }

            if (collided) {
                break;
            }
        }
    }
    
    if (counters[(gridY + addY) * cols + gridX + addX] > 1 && !collided) {
        for (uint i = 0; i < 4; i++) {
            otherParticle = cells[(gridY + addY) * cols * MAX_PARTICLES_PER_CELL + (gridX + addX) * MAX_PARTICLES_PER_CELL + i] - 1
            if (otherParticle == 0) {
                break;
            } else {
                x2 = particles[otherParticle * 4];
                y2 = particles[otherParticle * 4 + 1];
                collided = d_checkCollision(x1, y1, x2, y2, size / 2);
            }

            if (collided) {
                break;
            }
        }
    }

    // Resolve collision
    if (collided) {
        float dx2 = particles[otherParticle * 4 + 2];
        float dy2 = particles[otherParticle * 4 + 3];

        vels[particle * 2] = (2 * mass * dx2) / (2 * mass);
        vels[particle * 2 + 1] = (2 * mass * dy2) / (2 * mass);
    }
}

__device__ bool d_checkCollision(float x1, float y1, float x2, float y2, float r) {
    return x1 + 2 * r > x2
        && x1 < x2 + 2 * r
        && y1 + 2 * r > y2
        && y1 < y2 + 2 * r;
}

void Particles::updateCollisions() {
    cudaMemcpy(newVels, blankVels, sizeof(float) * this->numParticles * 2, cudaMemcpyHostToDevice);

    d_updateCollisions<<<(this->numParticles + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>>(d_particles, gridCounters, gridCells, newVels, this->particleSize, this->particleMass);
}

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
        it->hasCollided = false;
    }
}

inline bool Particles::checkAABBCircle(float p1x, float p1y, float p1r, float p2x, float p2y, float p2r) {
    return p1x + p1r + p2r > p2x
        && p1x < p2x + p1r + p2r
        && p1y + p1r + p2r > p2y
        && p1y < p2y + p1r + p2r;
}

inline bool Particles::checkAABBRect(float p1x, float p1y, float p1w, float p1h, float p2x, float p2y, float p2w, float p2h) {
    return p1x + p1w + p2w > p2x
        && p1x < p2x + p1w + p2w
        && p1y + p1h + p2h > p2y
        && p1y < p2y + p1h + p2h;
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
    if (checkAABBCircle(p1->get_x(), p1->get_y(), p1->get_radius(), p2->get_x(), p2->get_y(), p2->get_radius())) {
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