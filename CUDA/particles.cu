#include "particles.hpp"

#include <iostream>
#include <algorithm>
#include <immintrin.h>
#include <cuda.h>

// ================================================================== //

// https://leimao.github.io/blog/Proper-CUDA-Error-Checking/
#define CHECK_CUDA_ERROR(val) check((val), #val, __FILE__, __LINE__)
template <typename T>
void check(T err, const char* const func, const char* const file,
           const int line)
{
    if (err != cudaSuccess)
    {
        std::cerr << "CUDA Runtime Error at: " << file << ":" << line
                  << std::endl;
        std::cerr << cudaGetErrorString(err) << " " << func << std::endl;
        // We don't exit when we encounter CUDA errors in this example.
        // std::exit(EXIT_FAILURE);
    }
}

#define CHECK_LAST_CUDA_ERROR() checkLast(__FILE__, __LINE__)
void checkLast(const char* const file, const int line)
{
    cudaError_t err{cudaGetLastError()};
    if (err != cudaSuccess)
    {
        std::cerr << "CUDA Runtime Error at: " << file << ":" << line
                  << std::endl;
        std::cerr << cudaGetErrorString(err) << std::endl;
        // We don't exit when we encounter CUDA errors in this example.
        // std::exit(EXIT_FAILURE);
    }
}

// ================================================================== //

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

    temp_particles = new float[this->numParticles * DATA_PER_PARTICLE]();

    for (uint i = 0; i < this->numParticles; i++) {
        bool intersecting = true;
        while (intersecting) {
            float x = this->particleSize + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(this->width - 2 * this->particleSize)));
            float y = this->particleSize + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(this->height - 2 * this->particleSize)));
            float dx = -maxSpeed + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2 * maxSpeed)));
            float dy = -maxSpeed + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2 * maxSpeed)));
            Particle *newParticle = new Particle(x, y, dx, dy, this->particleSize, this->particleMass);
            newParticle->particleNum = i;

            if (particles.empty()) {
                intersecting = false;
                particles.push_back(*newParticle);

                temp_particles[i * DATA_PER_PARTICLE] = newParticle->get_x();
                temp_particles[i * DATA_PER_PARTICLE + 1] = newParticle->get_y();
                temp_particles[i * DATA_PER_PARTICLE + 2] = newParticle->get_dx();
                temp_particles[i * DATA_PER_PARTICLE + 3] = newParticle->get_dy();
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

                    temp_particles[i * DATA_PER_PARTICLE] = newParticle->get_x();
                    temp_particles[i * DATA_PER_PARTICLE + 1] = newParticle->get_y();
                    temp_particles[i * DATA_PER_PARTICLE + 2] = newParticle->get_dx();
                    temp_particles[i * DATA_PER_PARTICLE + 3] = newParticle->get_dy();
                } else {
                    delete newParticle;
                }
            }
        }
    }

    // Allocate device particle storage
    CHECK_CUDA_ERROR(cudaMalloc((void **)&d_particles, sizeof(float) * this->numParticles * DATA_PER_PARTICLE));
    CHECK_CUDA_ERROR(cudaMemcpy(d_particles, temp_particles, sizeof(float) * this->numParticles * DATA_PER_PARTICLE, cudaMemcpyHostToDevice));

    this->gridSize = config->particleSize * 2;

    this->numRows = this->height / this->gridSize;
    this->numCols = this->width / this->gridSize;

    // std::cout << this->numRows << std::endl;

    // Allocate grid
    CHECK_CUDA_ERROR(cudaMalloc((void **)&gridCounters, sizeof(uint) * this->numRows * this->numCols));
    CHECK_CUDA_ERROR(cudaMalloc((void **)&gridCells, sizeof(uint) * this->numRows * this->numCols * MAX_PARTICLES_PER_CELL));

    // Allocate blank grid host side
    blankCounters = new uint[this->numRows * this->numCols]();
    blankCells = new uint[this->numRows * this->numCols * MAX_PARTICLES_PER_CELL]();

    // Allocate velocity buffers
    CHECK_CUDA_ERROR(cudaMalloc((void **)&newVels, sizeof(float) * this->numParticles * 2));
    blankVels = new float[this->numParticles * 2]();
}

Particles::~Particles() {
    CHECK_CUDA_ERROR(cudaFree(gridCounters));
    CHECK_CUDA_ERROR(cudaFree(gridCells));
    CHECK_CUDA_ERROR(cudaFree(newVels));

    delete [] blankCounters;
    delete [] blankCells;
    delete [] blankVels;
}

/**
 * @brief Increments the simulation by updating the current time.
 * 
 */
void Particles::updateTime() {
    this->currentTime += this->timeStep;
}

__global__ void d_updateGrid(float *particles, uint *counters, uint *cells, float size, uint cols, uint numParticles) {
    uint particle = blockIdx.x * blockDim.x + threadIdx.x;
    // printf("%d\n", blockDim.x);

    if (particle < numParticles) {
        // printf("s: %d\n", particle);
        #if DEBUG_GRID
        printf("Start update grid for particle %d\n", particle);
        printf("Particle %d; X: %f; Y: %f; dX: %f; dY: %f\n", particle, particles[particle * DATA_PER_PARTICLE], particles[particle * DATA_PER_PARTICLE + 1], particles[particle * 4 + 2], particles[particle * 4 + 3]);
        #endif

        // Figure out which grid we're in
        uint gridX = floor(particles[particle * DATA_PER_PARTICLE] / size);
        uint gridY = floor(particles[particle * DATA_PER_PARTICLE + 1] / size);

        #if DEBUG_GRID
        printf("Particle %d at; X: %d; Y: %d\n", particle, gridX, gridY);
        #endif

        // Add to counter and cell
        int idx = atomicAdd(&counters[gridY * cols + gridX], 1);
        uint *cell = &cells[gridY * cols * MAX_PARTICLES_PER_CELL + gridX * MAX_PARTICLES_PER_CELL + idx];
        *cell = particle + 1;

        #if DEBUG_GRID
        if (particle == 0) {
            for (uint x = 0; x < cols; x++) {
                for (uint y = 0; y < cols; y++) {
                    printf("%d particle(s) at %d, %d: ", counters[y * cols + x], x, y);
                    for (uint z = 0; z < MAX_PARTICLES_PER_CELL; z++) {
                        int e = cells[y * cols * MAX_PARTICLES_PER_CELL + x * MAX_PARTICLES_PER_CELL + z] - 1;
                        if (e == -1) {
                            break;
                        } else {
                            printf("%d; ", e);
                        }
                    }
                    printf("\n");
                }
            }
        }

        printf("End update grid for particle %d\n", particle);
        #endif
        // printf("e: %d\n", particle);
    }
}

void Particles::updateGrid() {
    CHECK_CUDA_ERROR(cudaMemcpy(gridCounters, blankCounters, sizeof(uint) * this->numRows * this->numCols, cudaMemcpyHostToDevice));
    CHECK_CUDA_ERROR(cudaMemcpy(gridCells, blankCells, sizeof(uint) * this->numRows * this->numCols * MAX_PARTICLES_PER_CELL, cudaMemcpyHostToDevice));

    d_updateGrid<<<(this->numParticles + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>>(d_particles, gridCounters, gridCells, this->gridSize, this->numCols, this->numParticles);
    CHECK_LAST_CUDA_ERROR();
}

__device__ bool d_checkCollision(float x1, float y1, float x2, float y2, float r) {
    if (x1 + 2 * r > x2 && x1 < x2 + 2 * r && y1 + 2 * r > y2 && y1 < y2 + 2 * r) {
        float distance = sqrtf((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
        return distance < 2 * r;
    }
    return false;
}

__global__ void d_updateCollisions(float *particles, uint *counters, uint *cells, float *vels, float size, float mass, uint cols, uint numParticles) {
    int particle = blockIdx.x * blockDim.x + threadIdx.x;

    if (particle < numParticles) {
        #if DEBUG_COLLISION
        printf("Start update collision for particle %d\n", particle);
        #endif
        // Check current cell
        float x1 = particles[particle * DATA_PER_PARTICLE];
        float y1 = particles[particle * DATA_PER_PARTICLE + 1];
        uint gridX = floor(x1 / size);
        uint gridY = floor(y1 / size);

        // Variables for the other particle
        int otherParticle;
        float x2, y2;

        bool collided = false;

        if (!collided) {
            #if DEBUG_COLLISION
            printf("Particle %d checking cell %d, %d\n", particle, gridX, gridY);
            #endif
            for (uint i = 0; i < DATA_PER_PARTICLE; i++) {
                otherParticle = cells[gridY * cols * MAX_PARTICLES_PER_CELL + gridX * MAX_PARTICLES_PER_CELL + i] - 1;
                if (otherParticle == -1) {
                    break;
                } else if (particle == otherParticle) {
                    continue;
                } else {
                    x2 = particles[otherParticle * DATA_PER_PARTICLE];
                    y2 = particles[otherParticle * DATA_PER_PARTICLE + 1];
                    collided = d_checkCollision(x1, y1, x2, y2, size / 2);
                }
                #if DEBUG_COLLISION
                printf("Particle %d checking %d\n", particle, otherParticle);
                #endif
                if (collided) {
                    break;
                }
            }
        } 
        
        if (!collided && gridX + 1 < cols) {
            #if DEBUG_COLLISION
            printf("Particle %d checking cell %d, %d\n", particle, gridX + 1, gridY);
            #endif
            for (uint i = 0; i < 4; i++) {
                otherParticle = cells[gridY * cols * MAX_PARTICLES_PER_CELL + (gridX + 1) * MAX_PARTICLES_PER_CELL + i] - 1;
                if (otherParticle == -1) {
                    break;
                } else if (particle == otherParticle) {
                    continue;
                } else {
                    x2 = particles[otherParticle * DATA_PER_PARTICLE];
                    y2 = particles[otherParticle * DATA_PER_PARTICLE + 1];
                    collided = d_checkCollision(x1, y1, x2, y2, size / 2);
                }
                #if DEBUG_COLLISION
                printf("Particle %d checking %d\n", particle, otherParticle);
                #endif
                if (collided) {
                    break;
                }
            }
        }
        
        if (!collided && gridX - 1 < cols) {
            #if DEBUG_COLLISION
            printf("Particle %d checking cell %d, %d\n", particle, gridX - 1, gridY);
            #endif
            for (uint i = 0; i < 4; i++) {
                otherParticle = cells[gridY * cols * MAX_PARTICLES_PER_CELL + (gridX - 1) * MAX_PARTICLES_PER_CELL + i] - 1;
                if (otherParticle == -1) {
                    break;
                } else if (particle == otherParticle) {
                    continue;
                } else {
                    x2 = particles[otherParticle * DATA_PER_PARTICLE];
                    y2 = particles[otherParticle * DATA_PER_PARTICLE + 1];
                    collided = d_checkCollision(x1, y1, x2, y2, size / 2);
                }
                #if DEBUG_COLLISION
                printf("Particle %d checking %d\n", particle, otherParticle);
                #endif
                if (collided) {
                    break;
                }
            }
        }

        if (!collided && gridY + 1 < cols) {
            #if DEBUG_COLLISION
            printf("Particle %d checking cell %d, %d\n", particle, gridX, gridY + 1);
            #endif
            for (uint i = 0; i < 4; i++) {
                otherParticle = cells[(gridY + 1) * cols * MAX_PARTICLES_PER_CELL + gridX * MAX_PARTICLES_PER_CELL + i] - 1;
                if (otherParticle == -1) {
                    break;
                } else if (particle == otherParticle) {
                    continue;
                } else {
                    x2 = particles[otherParticle * DATA_PER_PARTICLE];
                    y2 = particles[otherParticle * DATA_PER_PARTICLE + 1];
                    collided = d_checkCollision(x1, y1, x2, y2, size / 2);
                }
                #if DEBUG_COLLISION
                printf("Particle %d checking %d\n", particle, otherParticle);
                #endif
                if (collided) {
                    break;
                }
            }
        }

        if (!collided && gridY - 1 < cols) {
            #if DEBUG_COLLISION
            printf("Particle %d checking cell %d, %d\n", particle, gridX, gridY - 1);
            #endif
            for (uint i = 0; i < 4; i++) {
                otherParticle = cells[(gridY - 1) * cols * MAX_PARTICLES_PER_CELL + gridX * MAX_PARTICLES_PER_CELL + i] - 1;
                if (otherParticle == -1) {
                    break;
                } else if (particle == otherParticle) {
                    continue;
                } else {
                    x2 = particles[otherParticle * DATA_PER_PARTICLE];
                    y2 = particles[otherParticle * DATA_PER_PARTICLE + 1];
                    collided = d_checkCollision(x1, y1, x2, y2, size / 2);
                }
                #if DEBUG_COLLISION
                printf("Particle %d checking %d\n", particle, otherParticle);
                #endif
                if (collided) {
                    break;
                }
            }
        }
        
        if (!collided && gridX + 1 < cols && gridY + 1 < cols) {
            #if DEBUG_COLLISION
            printf("Particle %d checking cell %d, %d\n", particle, gridX + 1, gridY + 1);
            #endif
            for (uint i = 0; i < 4; i++) {
                otherParticle = cells[(gridY + 1) * cols * MAX_PARTICLES_PER_CELL + (gridX + 1) * MAX_PARTICLES_PER_CELL + i] - 1;
                if (otherParticle == -1) {
                    break;
                } else if (particle == otherParticle) {
                    continue;
                } else {
                    x2 = particles[otherParticle * DATA_PER_PARTICLE];
                    y2 = particles[otherParticle * DATA_PER_PARTICLE + 1];
                    collided = d_checkCollision(x1, y1, x2, y2, size / 2);
                }
                #if DEBUG_COLLISION
                printf("Particle %d checking %d\n", particle, otherParticle);
                #endif
                if (collided) {
                    break;
                }
            }
        }

        if (!collided && gridX - 1 < cols && gridY + 1 < cols) {
            #if DEBUG_COLLISION
            printf("Particle %d checking cell %d, %d\n", particle, gridX - 1, gridY + 1);
            #endif
            for (uint i = 0; i < 4; i++) {
                otherParticle = cells[(gridY + 1) * cols * MAX_PARTICLES_PER_CELL + (gridX - 1) * MAX_PARTICLES_PER_CELL + i] - 1;
                if (otherParticle == -1) {
                    break;
                } else if (particle == otherParticle) {
                    continue;
                } else {
                    x2 = particles[otherParticle * DATA_PER_PARTICLE];
                    y2 = particles[otherParticle * DATA_PER_PARTICLE + 1];
                    collided = d_checkCollision(x1, y1, x2, y2, size / 2);
                }
                #if DEBUG_COLLISION
                printf("Particle %d checking %d\n", particle, otherParticle);
                #endif
                if (collided) {
                    break;
                }
            }
        }

        if (!collided && gridX + 1 < cols && gridY - 1 < cols) {
            #if DEBUG_COLLISION
            printf("Particle %d checking cell %d, %d\n", particle, gridX + 1, gridY - 1);
            #endif
            for (uint i = 0; i < 4; i++) {
                otherParticle = cells[(gridY - 1) * cols * MAX_PARTICLES_PER_CELL + (gridX + 1) * MAX_PARTICLES_PER_CELL + i] - 1;
                if (otherParticle == -1) {
                    break;
                } else if (particle == otherParticle) {
                    continue;
                } else {
                    x2 = particles[otherParticle * DATA_PER_PARTICLE];
                    y2 = particles[otherParticle * DATA_PER_PARTICLE + 1];
                    collided = d_checkCollision(x1, y1, x2, y2, size / 2);
                }
                #if DEBUG_COLLISION
                printf("Particle %d checking %d\n", particle, otherParticle);
                #endif
                if (collided) {
                    break;
                }
            }
        }

        if (!collided && gridX - 1 < cols && gridY - 1 < cols) {
            #if DEBUG_COLLISION
            printf("Particle %d checking cell %d, %d\n", particle, gridX - 1, gridY - 1);
            #endif
            for (uint i = 0; i < 4; i++) {
                otherParticle = cells[(gridY - 1) * cols * MAX_PARTICLES_PER_CELL + (gridX - 1) * MAX_PARTICLES_PER_CELL + i] - 1;
                if (otherParticle == -1) {
                    break;
                } else if (particle == otherParticle) {
                    continue;
                } else {
                    x2 = particles[otherParticle * DATA_PER_PARTICLE];
                    y2 = particles[otherParticle * DATA_PER_PARTICLE + 1];
                    collided = d_checkCollision(x1, y1, x2, y2, size / 2);
                }
                #if DEBUG_COLLISION
                printf("Particle %d checking %d\n", particle, otherParticle);
                #endif
                if (collided) {
                    break;
                }
            }
        }

        // Resolve collision
        if (collided) {
            #if DEBUG_COLLISION
            printf("Particle %d collided with %d\n", particle, otherParticle);
            #endif

            float dx2 = particles[otherParticle * DATA_PER_PARTICLE + 2];
            float dy2 = particles[otherParticle * DATA_PER_PARTICLE + 3];

            vels[particle * 2] = (2 * mass * dx2) / (2 * mass);
            vels[particle * 2 + 1] = (2 * mass * dy2) / (2 * mass);
        } else {
            vels[particle * 2] = particles[particle * DATA_PER_PARTICLE + 2];
            vels[particle * 2 + 1] = particles[particle * DATA_PER_PARTICLE + 3];
        }

        #if DEBUG_COLLISION
        printf("Particle %d set to; dX: %f; dY: %f\n", particle, vels[particle * 2], vels[particle * 2 + 1]);
        printf("End update collision for particle %d\n", particle);
        #endif
    }
}

void Particles::updateCollisions() {
    CHECK_CUDA_ERROR(cudaMemcpy(newVels, blankVels, sizeof(float) * this->numParticles * 2, cudaMemcpyHostToDevice));
    
    d_updateCollisions<<<(this->numParticles + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>>(d_particles, gridCounters, gridCells, newVels, this->gridSize, this->particleMass, this->numCols, this->numParticles);
    CHECK_LAST_CUDA_ERROR();
}

__global__ void d_updateMovements(float *particles, float *vels, float timeStep, float width, float radius, uint numParticles) {
    uint particle = blockIdx.x * blockDim.x + threadIdx.x;

    if (particle < numParticles) {
        #if DEBUG_MOVEMENT
        printf("Start update movements for particle %d\n", particle);
        #endif

        float dx = vels[particle * 2];
        float dy = vels[particle * 2 + 1];

        #if DEBUG_MOVEMENT
        printf("Particle %d initially at; dX: %f; dY: %f\n", particle, dx, dy);
        #endif

        float x = particles[particle * DATA_PER_PARTICLE];
        float y = particles[particle * DATA_PER_PARTICLE + 1];

        if (x < radius || x > width - radius) {
            dx *= -1;
            #if DEBUG_MOVEMENT
            printf("Particle %d bounce on X at %f\n", particle, x);
            #endif
        }
        if (y < radius || y > width - radius) {
            dy *= -1;
            #if DEBUG_MOVEMENT
            printf("Particle %d bounce on Y at %f\n", particle, y);
            #endif
        }

        float newX = particles[particle * DATA_PER_PARTICLE] + dx * timeStep;
        float newY = particles[particle * DATA_PER_PARTICLE + 1] + dy * timeStep;

        #if DEBUG_MOVEMENT
        printf("Particle %d moved from %f, %f to %f, %f\n", particle, particles[particle * DATA_PER_PARTICLE], particles[particle * DATA_PER_PARTICLE + 1], newX, newY);
        #endif

        particles[particle * DATA_PER_PARTICLE] = newX;
        particles[particle * DATA_PER_PARTICLE + 1] = newY;
        particles[particle * DATA_PER_PARTICLE + 2] = dx;
        particles[particle * DATA_PER_PARTICLE + 3] = dy;

        #if DEBUG_MOVEMENT
        printf("End update movements for particle %d\n", particle);
        #endif
    }
}

void Particles::updateMovements() {
    d_updateMovements<<<(this->numParticles + NUM_THREADS - 1) / NUM_THREADS, NUM_THREADS>>>(d_particles, newVels, this->timeStep, this->width, this->particleSize, this->numParticles);
    CHECK_LAST_CUDA_ERROR();

    CHECK_CUDA_ERROR(cudaMemcpy(temp_particles, d_particles, sizeof(float) * this->numParticles * DATA_PER_PARTICLE, cudaMemcpyDeviceToHost));

    // Can we omp this lmao
    for (uint i = 0; i < this->numParticles; i++) {
        Particle *temp = &(particles.at(i));

        float x = temp_particles[i * DATA_PER_PARTICLE];
        float y = temp_particles[i * DATA_PER_PARTICLE + 1];

        temp->set_x(x);
        temp->set_y(y);

        // temp->update(this->timeStep);
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