#include <iostream>

#include "Headerfiles/verlet.h"
#include "Headerfiles/types.h"
#include "Headerfiles//lenardJonesDirectionSummation.h"
#include "Headerfiles/xyz.h"
#include<fstream>

using namespace std;

int main() {
    /** Vars */
    double epsilon = 1;
    double sigma = 1;
    double mass = 1;
    double timeStep = 0.001 * sqrt((mass*sigma*sigma)/epsilon);
    double totalTime = 0.001  * sqrt((mass*sigma*sigma)/epsilon); //TODO change back
    double safeDumpTime = 1  * sqrt((mass*sigma*sigma)/epsilon);
    double currentTime = 0;
    double energy = 0;
    /** Init */
    std::cout << "Molecular Dynamics Project" << std::endl;
    auto [names, positions, velocities]{read_xyz_with_velocities("../AJupyter/lj54.xyz")};
    Atoms atoms(names,positions,velocities);

    /** Initial State */
    energy = lendardJonesDirectSummation(atoms,epsilon,sigma);

    int i = 0;
    /** Loop */
    do {
        i++;
        std::cout << "Step:" << i << std::endl;
        //verlet1
        verletStep1(atoms.positions,atoms.velocities,atoms.forces,timeStep);
        //forces
        energy = lendardJonesDirectSummation(atoms,epsilon,sigma);
        energy += calculateKineticEnergy(atoms); //TODO: implementation
        //TODO: dumps
        //verlet2
        verletStep2(atoms.velocities,atoms.forces,timeStep);
        //update current time
        currentTime += timeStep;
        //safe data dumps
        std::cout << "Writing Dump at:" << currentTime << std::endl;
        write_xyz("testTrajectory.xyz",atoms);
    }while(currentTime < totalTime);

    return 0;
}
