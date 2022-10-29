#include "particles.hpp"

#include <iostream>

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
    float particleSize = config->particleSize;
    float particleMass = config->particleMass;
    float maxSpeed = config->maxSpeed;

    for (uint i = 0; i < this->numParticles; i++) {
        bool intersecting = true;
        while (intersecting) {
            float x = particleSize + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(width - 2 * particleSize)));
            float y = particleSize + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(height - 2 * particleSize)));
            float dx = -maxSpeed + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2 * maxSpeed)));
            float dy = -maxSpeed + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2 * maxSpeed)));
            Particle *newParticle = new Particle(x, y, dx, dy, particleSize, particleMass);
            newParticle->particleNum = i;

            if (particles.empty()) {
                intersecting = false;
                particles.push_back(*newParticle);
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
                } else {
                    delete newParticle;
                }
            }
        }
    }

    this->numRows = config->gridHeight;
    this->numCols = config->gridWidth;

    this->gridHeight = height / numRows;
    this->gridWidth = width / numCols;

    // std::cout << this->numRows << std::endl;

    // Generate Grid
    for (uint i = 0; i < this->numRows; i++) {
        GridRow newRow;
        newRow.rowNum = i;
        for (uint j = 0; j < this->numCols; j++) {
            Cell newCell;
            newCell.cellNum = j;
            newRow.cells.push_back(newCell);
        }
        this->rows.push_back(newRow);
    }

    // Generate Buffer
    for (uint i = 0; i < this->numRows; i++) {
        GridRow newRow;
        newRow.rowNum = i;
        for (uint j = 0; j < this->numCols; j++) {
            Cell newCell;
            newCell.cellNum = j;
            newRow.cells.push_back(newCell);
        }
        this->buffer.push_back(newRow);
    }
}

/**
 * @brief Increments the simulation by updating the current time.
 * 
 */
void Particles::updateTime() {
    this->currentTime += this->timeStep;
}

void Particles::updateGrid() {
    // TODO: write changes into a "dirty" object then copy over

    // Empty grid and add particles
    uint i, j, k;
    #pragma omp parallel private(i, j) num_threads(64)
    {
        #pragma omp for
        for (i = 0; i < numRows; i++) {
            auto row = rows.at(i);
            auto bufRow = &(buffer.at(i));
            float gridY = gridHeight * (i + 0.5);
            for (j = 0; j < numCols; j++) {
                auto cell = row.cells.at(j);
                auto bufCell = &(bufRow->cells.at(j));
                float gridX = gridWidth * (cell.cellNum + 0.5);
                bufCell->elements.clear();
                for (k = 0; k < numParticles; k++) {
                    auto p = particles.at(k);
                    // Check if AABB intersects with grid
                    if (checkAABBRect(p.get_x(), p.get_y(), p.get_radius(), p.get_radius(), gridX, gridY, 0.5 * gridWidth, 0.5 * gridHeight)) {
                        bufCell->elements.push_back(p.particleNum);
                    }
                }
            }
        }
    }

    for (i = 0; i < numRows; i++) {
        auto row = &(rows.at(i));
        auto bufRow = buffer.at(i);
        for (j = 0; j < numCols; j++) {
            auto cell = &(row->cells.at(j));
            auto bufCell = bufRow.cells.at(j);
            cell->elements = bufCell.elements;
        }
    }
}

void Particles::updateCollisions() {
    uint i, j;
    #pragma omp parallel for private(i, j) num_threads(64)
    for (i = 0; i < numRows; i++) {
        auto row = rows.at(i);
        for (j = 0; j < numRows; j++) {
            auto cell = row.cells.at(j);
            if (cell.elements.size() > 1) {
                for (int i1 : cell.elements) {
                    auto p1 = &(particles.at(i1));
                    for (int i2 : cell.elements) {
                        auto p2 = &(particles.at(i2));
                        if (i1 != i2 && !p1->hasCollided && !p2->hasCollided) {
                            if (checkCollision(p1, p2)) {
                                // https://gamedevelopment.tutsplus.com/tutorials/when-worlds-collide-simulating-circle-circle-collisions--gamedev-769
                                float newXVel1 = (p1->get_dx() * (p1->get_mass() - p2->get_mass()) + (2 * p2->get_mass() * p2->get_dx())) / (p1->get_mass() + p2->get_mass());
                                float newXVel2 = (p2->get_dx() * (p2->get_mass() - p1->get_mass()) + (2 * p1->get_mass() * p1->get_dx())) / (p1->get_mass() + p2->get_mass());
                                float newYVel1 = (p1->get_dy() * (p1->get_mass() - p2->get_mass()) + (2 * p2->get_mass() * p2->get_dy())) / (p1->get_mass() + p2->get_mass());
                                float newYVel2 = (p2->get_dy() * (p2->get_mass() - p1->get_mass()) + (2 * p1->get_mass() * p1->get_dy())) / (p1->get_mass() + p2->get_mass());

                                p1->set_dx(newXVel1);
                                p2->set_dx(newXVel2);
                                p1->set_dy(newYVel1);
                                p2->set_dy(newYVel2);

                                p1->hasCollided = true;
                                p2->hasCollided = true;
                            }
                        }
                    }
                }
            }
        }
    }
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