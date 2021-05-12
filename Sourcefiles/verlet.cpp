//
// Created by cm on 06.05.21.
//
#include <iostream>
#include "../Headerfiles/verlet.h"

void verletStep1(Positions_t &positions, Velocities_t &velocities, const Forces_t &forces ,double timestep)
{
    double mass = 1.0;
    //velocity
    velocities += 0.5 * forces * timestep/mass;
    //postions
    positions += velocities * timestep;
}

void verletStep2(Velocities_t &velocities, Forces_t &forces, double timestep)
{
    double mass = 1.0;
    velocities += 0.5 * forces * timestep/mass;
}

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