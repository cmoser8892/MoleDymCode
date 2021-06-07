//
// Created by cm on 06.05.21.
//

#ifndef __VERLET_H_
#define __VERLET_H_

#include "../Headerfiles/types.h"
#include "../Headerfiles/atoms.h"

void verletStep1(Positions_t &positions, Velocities_t &velocities, const Forces_t &forces ,double timestep);
void verletStep1Atoms(Atoms &atoms,double timestep);

void verletStep2(Velocities_t &velocities, Forces_t &forces, double timestep);
void verletStep2Atoms(Atoms &atoms,double timestep);

//
void verletIntegratorConstantForce(Positions_t &positions, Velocities_t &velocities, Forces_t &forces, double timestep, unsigned int nbSteps);
void verletIntegratorConstantForceAtoms(Atoms &atoms, double timestep, unsigned int nbSteps);
#endif //__VERLET_H_
