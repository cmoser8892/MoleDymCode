//
// Created by cm on 06.05.21.
//
#include <iostream>
#include "../Headerfiles/verlet.h"

/**
 * @fn void verletStep1(Positions_t &positions, Velocities_t &velocities, const Forces_t &forces ,double timestep)
 * @brief calculates the first step
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
/**
 * @fn void verletStep2(Velocities_t &velocities, Forces_t &forces, double timestep)
 * @brief
 * @param velocities
 * @param forces
 * @param timestep
 */
void verletStep2(Velocities_t &velocities, Forces_t &forces, double timestep)
{
    double mass = 1.0;
    velocities += 0.5 * forces * timestep/mass;
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