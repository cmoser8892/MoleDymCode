//
// Created by cm on 06.05.21.
//
#include <iostream>
#include "../Headerfiles/verlet.h"

/**
 * @fn void verletStep1(Positions_t &positions, Velocities_t &velocities, const Forces_t &forces ,double timestep)
 * @brief calculates the first step with a constant mass of 1
 * @param positions
 * @param velocities
 * @param forces
 * @param timestep
 */
void verletStep1(Positions_t &positions, Velocities_t &velocities, const Forces_t &forces ,double timestep)
{
    double mass = 1.0;
    //velocity
    velocities += 0.5 * forces * timestep/mass;
    //postions
    positions += velocities * timestep;
}

void verletStep1Atoms(Atoms &atoms,double timestep) {
    //velocity
    for(int i = 0; i < atoms.nb_atoms(); i++) {
        atoms.velocities.col(i) += 0.5 * atoms.forces.col(i) * timestep/atoms.mass(i);
    }
    //position
    atoms.positions += atoms.velocities*timestep;
}

/**
 * @fn void verletStep2(Velocities_t &velocities, Forces_t &forces, double timestep)
 * @brief calculates the second step with a constant mass of 1
 * @param velocities
 * @param forces
 * @param timestep
 */
void verletStep2(Velocities_t &velocities, Forces_t &forces, double timestep)
{
    double mass = 1.0;
    velocities += 0.5 * forces * timestep/mass;
}

void verletStep2Atoms(Atoms &atoms,double timestep) {
    for(int i = 0; i < atoms.nb_atoms(); i++) {
        atoms.velocities.col(i) += 0.5 *atoms.forces.col(i) * timestep/atoms.mass(i);
    }
}

/**
 * @fn void verletIntegratorConstantForce(Positions_t &positions, Velocities_t &velocities, Forces_t &forces, double timestep, unsigned int nbSteps)
 * @brief Propagates positions and velocities forward without the influence of any force
 * @param positions
 * @param velocities
 * @param forces
 * @param timestep
 * @param nbSteps
 */
void verletIntegratorConstantForce(Positions_t &positions, Velocities_t &velocities, Forces_t &forces
        , double timestep, unsigned int nbSteps)
{
    double mass = 1;
    // For loop
    for(int i = 0; i < nbSteps; ++i) {
        //std::cout << "Step: " << i << std::endl;
        verletStep1(positions,velocities,forces,timestep);
        verletStep2(velocities,forces,timestep);
    }
}

/**
 * @fn
 * @brief
 * @param atoms
 * @param timestep
 * @param nbSteps
 */
void verletIntegratorConstantForceAtoms(Atoms &atoms, double timestep, unsigned int nbSteps) {
    // For loop
    for(int i = 0; i < nbSteps; ++i) {
        //std::cout << "Step: " << i << std::endl;
        verletStep1Atoms(atoms,timestep);
        verletStep2Atoms(atoms,timestep);
    }
}