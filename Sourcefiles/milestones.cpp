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
#include "../Headerfiles/berendsenThermostat.h"

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

int milestone5Code(int argc, char *argv[]) {
    int returnValue = 0;
    /** Vars */
    double epsilon = 1;
    double sigma = 1;
    double mass = 12*atomicUnit; // 12C6
    unsigned int nbAtoms = 4;
    double targetTemperatur = 275; //about roomtemp
    /** Times */
    double timeStep = 0.01 * sqrt((mass * sigma * sigma) / epsilon); //around 10e-15
    double totalTime = 10000 *timeStep;
    double safeDumpTime = 100 * timeStep;
    double relaxationTimeFactor = 80.0;
    double relaxationTime = relaxationTimeFactor*timeStep;
    /** SafeLocations */
    std::string trajectorySafeLocation = "/home/cm/CLionProjects/MoleDymCode/cmake-build-debug/TrajectoryDumps";
    std::string energyDataSafeLocation = "/home/cm/CLionProjects/MoleDymCode/AJupyter";
    std::string trajectoryBaseName = "Trajectory";
    std::string energyName = "energy";
    /** Argument stuff */
    if(argc == 1) {
        /** Do nothing */
        std::cout << "No arguments given" << std::endl;
    }
    else if(argc >= 2){
        /** rewrite varibles */
        //atoms in the simulation
        nbAtoms = std::atoi(argv[1]); //arg2
        std::cout << nbAtoms << std::endl;
        //new safe locations given
        if(argc >= 3) {
            //relaxation time rescaling
            relaxationTimeFactor = std::atof(argv[2]); //arg3
            std::cout << relaxationTimeFactor << std::endl;
            relaxationTime = relaxationTimeFactor*timeStep;
            if(argc >= 4) {
                trajectorySafeLocation = argv[3]; //arg4
                energyDataSafeLocation = argv[4]; //arg5
            }
        }
    }
    /** Simulation related Variables */
    int safeAtStep = safeDumpTime/timeStep;
    double currentTime = 0;
    /** global */
    double energy = 0;
    double kineticEnergy = 0;
    std::vector<double> energyStorage(totalTime/timeStep);
    /** Init */
    Positions_t  p = createLatticeCube(nbAtoms,sigma);
    Atoms atoms(p,mass);
    setANameInAtoms(atoms, 'X');
    /** Initial State */
    int i = 0;
    /** Loop */
    while (currentTime <= totalTime) {
        //verlet1
        verletStep1Atoms(atoms, timeStep);
        //forces and energy
        energy = lendardJonesDirectSummation(atoms, epsilon, sigma);
        //verlet2
        verletStep2Atoms(atoms, timeStep);
        //velocity rescaling
        berendsenThermostat(atoms,targetTemperatur,timeStep,relaxationTime);
        //energy
        kineticEnergy = calculateKineticEnergy(atoms);
        energy += kineticEnergy;
        energyStorage[i] = energy;
        //check for the temperatur to increase the relaxation and keep the sim happy
        //btw this increased the relaxationtime slowly to infinity
        if(abs(calculateCurrentTemperatur(atoms)-targetTemperatur) < 10.) {
            //std::cout << "Increase the relaxation Time " << relaxationTime << std::endl;
            relaxationTime *= 100000;
        }
        //safe
        if ((i % safeAtStep) == 0) {
            //std::cout << "Writing Dump at:" << currentTime << " with " << i/safeAtStep << std::endl;
            //std::cout << energyStorage[i] << std::endl;
            //std::cout << kineticEnergy << " " << energyStorage[i]-kineticEnergy << " " << calculateCurrentTemperatur(atoms) << std::endl;
            dumpData(atoms, trajectorySafeLocation, trajectoryBaseName,
                     1000, (unsigned int) i / safeAtStep);
            if(checkMoleculeTrajectories(atoms,2* pow(2.0, 1.0/6.0)) == false) {
                std::cerr << "Cube Exploded at: " << i << std::endl;
                returnValue = -1;
                break;
            }
        }
        //update time and counter
        currentTime += timeStep;
        i++;
    }
    //energy dump for ploting
    //std::cout << "Dumping the energy" << std::endl;
    dumpEnergy(energyStorage, energyDataSafeLocation, energyName);
    return returnValue;
}