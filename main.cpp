#include <iostream>

#include "Headerfiles/verlet.h"
#include "Headerfiles/types.h"
#include "Headerfiles//lenardJonesDirectionSummation.h"

int main() {
    std::cout << "Molecular Dynamics Project" << std::endl;
    //initalization
    int nbAtoms = 2;
    Positions_t  positions(3,nbAtoms);
    Velocities_t  velocities(3, nbAtoms);
    Forces_t forces(3,nbAtoms);
    //
    positions.setZero();
    velocities.setZero();
    forces.setZero();
    forces.row(0) = 1.0;
    double timestep = 1.0;
    //main loop
    verletIntegratorConstantForce(positions,velocities,forces,timestep,10);
    //std::cout << positions << std::endl;
    //TEST add for compilation
    // go to the minimum to check the forces
    positions.col(0) = 0;
    positions(0,1) = 3;
    Atoms a(positions);
    lendardJonesDirectSummation(a);
    return 0;
}
