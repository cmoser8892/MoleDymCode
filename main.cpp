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
    /** Milestone 4 */
    //TEST add for compilation
    // go to the minimum to check the forces they should be near 0 otherwise i have a mistake
    double testSigma = 1;
    Atoms atoms(nbAtoms);
    //use the real minimum as it is calculated
    atoms.positions(0,1) =  pow(2.0,(1.0/6.0)) * testSigma;
    lendardJonesDirectSummation(atoms, 1,testSigma);
    std::cout << "R ~ 1.12sigma, Forces should be near 0" << std::endl;
    std::cout << atoms.forces << std::endl;
    //test if th forces filp the sign
    std::cout << "Filptest: Depending on sigma the forces should once be postive and once negative" << std::endl;
    atoms.positions(0,1) = 0.5 *testSigma;
    lendardJonesDirectSummation(atoms, 1,testSigma);
    std::cout << "Repulsion:" << std::endl;
    std::cout << atoms.forces << std::endl;
    atoms.positions(0,1) = 2 *testSigma;
    lendardJonesDirectSummation(atoms, 1,testSigma);
    std::cout << "Attraction:" << std::endl;
    std::cout << atoms.forces << std::endl;
    return 0;
}
