#include <iostream>

#include "Headerfiles/verlet.h"
#include "Headerfiles/types.h"
#include "Headerfiles//lenardJonesDirectionSummation.h"
#include "Headerfiles/xyz.h"
#include "Headerfiles/helperfunctions.h"

using namespace std;

int main() {
    /** Vars */
    double epsilon = 1;
    double sigma = 1;
    double mass = 1;
    /** Times */
    double timeStep = 0.001 * sqrt((mass*sigma*sigma)/epsilon);
    double totalTime = 3  * sqrt((mass*sigma*sigma)/epsilon); //TODO change back
    double safeDumpTime = 1  * sqrt((mass*sigma*sigma)/epsilon);
    int safeAtStep = safeDumpTime/timeStep;
    double currentTime = 0;
    /** global */
    double energy = 0;
    /** Init */
    std::cout << "Molecular Dynamics Project" << std::endl;
    auto [names, positions, velocities]{read_xyz_with_velocities("../AJupyter/lj54.xyz")};
    Atoms atoms(names,positions,velocities);

    /** Initial State */
    energy = lendardJonesDirectSummation(atoms,epsilon,sigma);
    dumpData(atoms,"/home/cm/CLionProjects/MoleDymCode/cmake-build-debug/TrajectoryDumps","Trajectory",
             (unsigned int) (totalTime/timeStep),0);
    int i = 0;
    /** Loop */
    while (currentTime < totalTime) {
        i++;
        //verlet1
        verletStep1(atoms.positions, atoms.velocities, atoms.forces, timeStep);
        //forces
        energy = lendardJonesDirectSummation(atoms, epsilon, sigma);
        energy += calculateKineticEnergy(atoms); //TODO: implementation
        //TODO: dump energy somewhere
        //verlet2
        verletStep2(atoms.velocities, atoms.forces, timeStep);
        //update current time
        currentTime += timeStep;
        //safe data dumps
        if((i%safeAtStep) == 0)
        {
            std::cout << "Writing Dump at:" << currentTime << std::endl;
            dumpData(atoms,"/home/cm/CLionProjects/MoleDymCode/cmake-build-debug/TrajectoryDumps","Trajectory",
                     (unsigned int) (totalTime/timeStep),i/safeAtStep);
        }
        //std::cout << "Step:" << i << std::endl;
        //std::cout << currentTime << std::endl;
    }

    return 0;
}
