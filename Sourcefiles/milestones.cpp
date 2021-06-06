//
// Created by cm on 06.06.21.
//
/** Include Everything */
#include "../Headerfiles/verlet.h"
#include "../Headerfiles/types.h"
#include "../Headerfiles//lenardJonesDirectionSummation.h"
#include "../Headerfiles/xyz.h"
#include "../Headerfiles/helperfunctions.h"
#include "../Headerfiles/milestones.h"

int milestone4Code() {
    /** Vars */
    double epsilon = 1;
    double sigma = 1;
    double mass = 1;
    /** Times */
    double timeStep = 0.01 * sqrt((mass*sigma*sigma)/epsilon);
    double totalTime = 100  * sqrt((mass*sigma*sigma)/epsilon);
    double safeDumpTime = 1  * sqrt((mass*sigma*sigma)/epsilon);
    int safeAtStep = safeDumpTime/timeStep;
    double currentTime = 0;
    /** global */
    double energy = 0;
    std::vector<double> energyStorage(totalTime/timeStep);
    /** Init */
    auto [names, positions, velocities]{read_xyz_with_velocities("../AJupyter/lj54.xyz")};
    Atoms atoms(names,positions,velocities);
    /** Initial State */
    if(0) { //this means the first image is not the initial state but only kinda
        energy = lendardJonesDirectSummation(atoms,epsilon,sigma);
        std::cout << "Writing Dump at:" << currentTime << std::endl;
        dumpData(atoms,"/home/cm/CLionProjects/MoleDymCode/cmake-build-debug/TrajectoryDumps","Trajectory",
                 (unsigned int) (totalTime/timeStep),0); }
    int i = 0;
    /** Loop */
    while (currentTime <= totalTime) {
        //verlet1
        verletStep1(atoms.positions, atoms.velocities, atoms.forces, timeStep);
        //forces and energy
        energy = lendardJonesDirectSummation(atoms, epsilon, sigma);
        //verlet2
        verletStep2(atoms.velocities, atoms.forces, timeStep);
        //safe data dumps
        energy += calculateKineticEnergy(atoms);
        energyStorage[i] = energy;
        if((i%safeAtStep) == 0)
        {
            std::cout << "Writing Dump at:" << currentTime << std::endl;
            std::cout << energyStorage[i] << std::endl;
            dumpData(atoms,"/home/cm/CLionProjects/MoleDymCode/cmake-build-debug/TrajectoryDumps","Trajectory",1000,(unsigned int) i/safeAtStep);
        }
        //update time and counter
        currentTime += timeStep;
        i++;
        //std::cout << "Step:" << i << std::endl;
        //std::cout << currentTime << std::endl;
    }
    //energy dump for ploting
    dumpEnergy(energyStorage,"/home/cm/CLionProjects/MoleDymCode/AJupyter","energy");
    return 0;
}