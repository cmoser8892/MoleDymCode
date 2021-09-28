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
#include "../Headerfiles/neighbors.h"
#include "../Headerfiles/embeddedAtomPotential.h"

/**
 * @fn int milestone4Code()
 * @brief function that contains all of the Milestone4 code
 * @param argc
 * @param argv
 * @return not relevant
 */
int milestone4Code(int argc, char *argv[]) {
    std::string filename = "energy";
    std::string location = "/home/cm/CLionProjects/MDCode/AData";
    /** Vars */
    double epsilon = 1;
    double sigma = 1;
    double mass = 1;
    //prefactor for the timestep with a bad name
    double norre = 0.1; //variablename by tschuli (dont ask ur friends for names lol)
    /** Value Read */
    if(argc == 1) {
        /** Do nothing */
        std::cout << "No arguments given" << std::endl;
    }
    else if(argc == 2){
        std::cout << "Additional Arguments" << std::endl;
        /** rewrite varibles */
        norre = atof(argv[1]);
        std::string number(argv[1]);
        std::cout << number << std::endl;
        filename += number;
    }
    /** Times */
    double preFactor = sqrt((mass*sigma*sigma)/epsilon);
    double timeStep = norre * preFactor;
    double totalTime = 10000  * timeStep;
    int safeAtStep = 100;
    double currentTime = 0;
    /** global */
    double energy = 0;
    std::vector<double> energyStorage(totalTime/timeStep);
    /** Init */
    auto [names, positions, velocities]{read_xyz_with_velocities("../AData/lj54.xyz")};
    Atoms atoms(names,positions,velocities);
    /** Initial State */
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
            //std::cout << "Writing Dump at:" << currentTime << std::endl;
            //std::cout << energyStorage[i] << std::endl;
            //dumpData(atoms,"/home/cm/CLionProjects/MoleDymCode/cmake-build-debug/TrajectoryDumps","Trajectory",1000,(unsigned int) i/safeAtStep);
        }
        //update time and counter
        currentTime += timeStep;
        i++;
        //std::cout << "Step:" << i << std::endl;
        //std::cout << currentTime << std::endl;
    }
    //energy dump for ploting
    dumpVectorData(energyStorage, location, filename);
    return 0;
}

/**
 * @fn int milestone5Code(int argc, char *argv[])
 * @brief function that contains all of the Milestone5 code
 * @param argc
 * @param argv
 * @return not relvant
 */
int milestone5Code(int argc, char *argv[]) {
    int returnValue = 0;
    /** Vars */
    double epsilon = 1;
    double sigma = 1;
    double mass = 12*atomicUnit; // 12C6 arbitrary kinda more stable for timerecoirding
    unsigned int nbAtoms = 160;
    double targetTemperatur = 275;
    /** Times */
    double timeStep = 0.01 * sqrt((mass * sigma * sigma) / epsilon); //around 10e-15
    double totalTime = 10000 * timeStep;
    double safeDumpTime = 100 * timeStep;
    double relaxationTimeFactor = 50.0;
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
    /** global */
    int safeAtStep = safeDumpTime/timeStep;
    double currentTime = 0;
    double energy = 0;
    double kineticEnergy = 0;
    std::vector<double> energyStorage(totalTime/timeStep);
    /** Init */
    Positions_t  p = createLatticeCube(nbAtoms,sigma);
    //Positions_t  p = createLatticesLongRod(nbAtoms,3,sigma+0.000001); //experimented with different shapes
    Atoms atoms(p,mass);
    setANameInAtoms(atoms, "X");
    Positions_t controlCube = generateCapsel(atoms, 100);
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
            if(checkMoleculeTrajectories(atoms,controlCube) == false) {
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
    dumpVectorData(energyStorage, energyDataSafeLocation, energyName);
    return returnValue;
}

/**
 * @fn int milestone6Code(int argc, char *argv[])
 * @brief function that contains all of the Milestone6 code
 * @param argc
 * @param argv
 * @return not relvant
 */
int milestone6Code(int argc, char *argv[]) {
    int returnValue = 0;
    bool once = false;
    /** Variables for the simulation */
    /** Atoms variables */
    double epsilon = 1;
    double sigma = 1; //*pow(2.0, 1.0/6.0); //stable in lattice (more or less)
    double mass = 12*atomicUnit; // 12C6 arbitrary kinda more stable for timerecoirding
    unsigned int nbAtoms = 60;
    bool thermostatUsed = true;
    double targetTemperatur = 275; //about 0
    double cutoffRange = 2.5 * sigma;

    /** Times */
    double timeStep = 0.01 * sqrt((mass * sigma * sigma) / epsilon); //around 10e-15
    double totalTime = 10000 *timeStep;
    double safeDumpTime = 100 * timeStep;
    double relaxationTimeFactor = 50.0;
    double relaxationTime = relaxationTimeFactor*timeStep;
    int safeAtStep = safeDumpTime/timeStep; //bad casting lol

    /** SafeLocations */
    std::string trajectorySafeLocation = "/home/cm/CLionProjects/MoleDymCode/cmake-build-debug/TrajectoryDumps";
    std::string energyDataSafeLocation = "/home/cm/CLionProjects/MoleDymCode/AJupyter";
    std::string trajectoryBaseName = "Trajectory";
    std::string energyName = "energy";
    std::vector<double> energyStorage(totalTime/timeStep);
    double energy = 0;
    double kineticEnergy = 0;


    /** Argument Processing */
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
    /** set up atoms */
    //init of the atoms
    Positions_t  p = createLatticeCube(nbAtoms,sigma+0.000001);
    //Positions_t  p = createLatticesLongRod(nbAtoms,3,sigma+0.000001); //changed for test was not used in data aquisation
    Atoms atoms(p,mass);
    setANameInAtoms(atoms, "X");
    // neighbor list
    NeighborList neighborList(cutoffRange);
    neighborList.update(atoms);
    Positions_t controlCube = generateCapsel(atoms, 100);
    /** Loop */
    int i = 0;
    double currentTime = 0;
    //
    while (currentTime <= totalTime) {
        //verlet1
        verletStep1Atoms(atoms, timeStep);
        //forces and energy
        energy = lenardJonesDirectSummationWithCutoff(atoms,cutoffRange,epsilon,sigma);
        //verlet2
        verletStep2Atoms(atoms, timeStep);
        if(thermostatUsed == true) {
            //velocity rescaling
            berendsenThermostat(atoms, targetTemperatur, timeStep, relaxationTime);
        }
        //energy
        kineticEnergy = calculateKineticEnergy(atoms);
        energy += kineticEnergy;
        energyStorage[i] = energy;
        //relaxation time
        if(thermostatUsed == true) {
            if (abs(calculateCurrentTemperatur(atoms) - targetTemperatur) < 100.) {
                if (once == false) {
                    //std::cout << "Increase the relaxation Time " << relaxationTime << std::endl;
                    relaxationTime *= 1000000;
                    once = true;
                }
            }
        }

        //Dumping the data and checking for an explosion
        if ((i % safeAtStep) == 0) {
            //std::cout << "Writing Dump at:" << currentTime << " with " << i/safeAtStep << std::endl;
            //std::cout << energyStorage[i] << std::endl;
            //std::cout << kineticEnergy << " " << energyStorage[i]-kineticEnergy << " " << calculateCurrentTemperatur(atoms) << std::endl;
            dumpData(atoms, trajectorySafeLocation, trajectoryBaseName,
                     1000, (unsigned int) i / safeAtStep);
            if(thermostatUsed == true) {
                if (checkMoleculeTrajectories(atoms, controlCube) == false) {
                    std::cerr << "Cube Exploded at: " << i << std::endl;
                    returnValue = -1;
                    break;
                }
            }
        }
        //update time and counter
        currentTime += timeStep;
        i++;
    }
    //
    dumpVectorData(energyStorage, energyDataSafeLocation, energyName);
    return returnValue;
}

///global storage :-)
std::vector<double> kineticEnergyStorage;
std::vector<double> potentialEnergyStorage;
std::vector<double> setEnergyStorage;
std::vector<double> setTemperaturStorage;

/**
 * @fn int milestone7Code(int argc, char *argv[])
 * @brief function that contains all of the Milestone7 code
 * @param argc
 * @param argv
 * @return not relvant
 */
int milestone7Code(int argc, char *argv[]) {
    int returnValue = 0;
    bool once = false;
    /** Gold Params
        double cutoff = 2.0; //A
        double A = 0.2061; //eV
        double xi = 1.790; //eV
        double p = 10.229; //
        double q = 4.036; //
        double re = 4.079 / sqrt(2)); //distance
    */
    /** Cluster generation */
    Names_t names;
    Positions_t positions;
    /** File basenames*/
    std::string dataLocation = "/home/cm/CLionProjects/MoleDymCode/AData";
    std::string kineticEnergyFile =     "kineticEnergy";
    std::string potentialEnergyFile =   "potentialEnergy";
    std::string temperaturFile =        "temperatur";
    std::string energyFile =            "energy";
    /** Argument Processing */
    if(argc == 1) {
        /** Do nothing */
        std::cout << "No arguments given" << std::endl;
        std::string location =  "/home/cm/CLionProjects/MoleDymCode/AData/Clusters";
        std::string filename = location + "/cluster" + "3" +".xyz";
        //auto [tupleNames, tuplePositions]{read_xyz("../AData/cluster_923.xyz")};
        auto [tupleNames, tuplePositions]{read_xyz(filename)};
        names = tupleNames;
        positions = tuplePositions;
    }
    else {
        /** Create cluster based on layers given */
        int layers = atoi(argv[1]);
        std::string number(argv[1]);
        //generate a cluster based on the layers
        std::string location =  "/home/cm/CLionProjects/MoleDymCode/AData/Clusters";
        generateClusterHull(layers, location); //stupid shit does not work
        //get from the preprocessed files
        std::string filename = location + "/cluster" + number +".xyz";
        auto [tupleNames, tuplePositions]{read_xyz(filename)};
        names = tupleNames;
        positions = tuplePositions;
        //make the data identifyable
        dataLocation = location; //change Data location
        kineticEnergyFile += number;
        potentialEnergyFile += number;
        temperaturFile += number;
        energyFile += number;
        std::cout << "LayerID" << number << std::endl;
    }
    /** Set up atoms */
    double atomicMassAu = 196.97; // 197Au79
    double mass = atomicMassAu/massCorrectionFactor; //mass is in u convert it to a correct mass for gupta
    Atoms atoms(names,positions,mass);
    /** Data */
    std::vector<double> meanEnergyStorage;
    std::vector<double> meanTemperaturStorage;
    /** set up data */
    //write all the data in the data stucture
    SimulationData_t data;
    ///
    data.simulationID = 0;
    data.maxTrajectoryNumber = 100000;  //used for numberconvention in the trajectory name
    data.trajectorySafeLocation = "/home/cm/CLionProjects/MoleDymCode/cmake-build-debug/TrajectoryDumps";
    data.trajectoryBaseName = "Trajectory";
    data.doDumping = false;
    data.totalEnergyRecording = true;
    ///
    data.controlCube = generateCapsel(atoms,400); //always has to be generated otherwise crash
    data.timeStep = 1; //in fs
    data.simulationTime = 10 * data.timeStep;
    data.relaxationTime = 1000 * data.timeStep;
    data.cutoffDistance = 10.0;
    data.targetTemperatur = 300;
    /** Main simulation */
    std::cout << "Starting Simulation" << std::endl;
    //preheating
    int runs = 50; //50
    /// simulation
    ///increase till room temp
    for (int i = 0; i < runs; ++i) {
        for( int j = 0; j < 100; ++j ) {
            auto[energy, temperatur]{simulationBuildStone(data, atoms)};
            ///
            printData(i,energy ,temperatur);
            ///
        }
    }
    ///relax a bit so the temp is in all cases stable
    data.relaxationTime *= 1e3333; // basically thermostat has now no effect: to infinity and beyond
    runs = 100; //10ps
    for ( int i = 0; i < runs; ++i) {
        auto[energy, temperatur]{simulationBuildStone(data, atoms)};
        ///
        printData(i,energy,temperatur);
        ///
    }
    ///data get for plot
    runs = 120; //60
    //10 - 100 ps
    for(int i = 0; i < runs; ++i) {
        int subruns = 600;
        ///
        //std::vector<double> setEnergyStorage;
        //std::vector<double> setTemperaturStorage;
        ///clear global vars could fit this into a class to make it pretty lol
        setEnergyStorage.clear();
        setTemperaturStorage.clear();
        ///
        depositRescaledHeat(5e-3*atoms.nb_atoms(),atoms);
        for(int j = 0; j < subruns/2; ++j) {
            auto[totalEnergy, temperatur]{simulationBuildStone(data, atoms)};
            ///
            printData(i,totalEnergy,temperatur);
            ///
        }
        ///clear global
        setEnergyStorage.clear();
        setTemperaturStorage.clear();
        for(int j = 0; j < (subruns/2); ++j) {
            auto[totalEnergy, temperatur]{simulationBuildStone(data, atoms)};
            ///
            printData(i,totalEnergy,temperatur);
            ///
        }
        //calculateMean
        meanEnergyStorage.push_back(averageVector(setEnergyStorage));
        meanTemperaturStorage.push_back(averageVector(setTemperaturStorage));
    }
    //safe the information
    dumpVectorData(kineticEnergyStorage,dataLocation,kineticEnergyFile);
    dumpVectorData(potentialEnergyStorage,dataLocation,potentialEnergyFile);
    dumpVectorData(meanEnergyStorage,dataLocation,energyFile);
    dumpVectorData(meanTemperaturStorage,dataLocation,temperaturFile);
    ////
    return returnValue;
}

/**
 * @fn std::tuple<double energy, double temperature> simulationBuildStone(SimulationData_t data,Atoms &atoms) )
 * @brief runs one simulation(-part) and saves the trajectory (allways with
 * @param data
 * @param atoms
 * @return the total energy (first) and then the temperatur at that step
 */
std::tuple<double, double> simulationBuildStone(SimulationData_t &data, Atoms &atoms) {
    double returnEnergy = 0;
    double kineticEnergy = 0;
    double potentialEnergy = 0;
    double returnTemperatur = calculateCurrentTemperaturEV(atoms);
    // simulation initialization
    int i = 0;
    double currentTime = 0;
    NeighborList list(data.cutoffDistance);
    //run the simulation for some time
    while (currentTime < data.simulationTime) {
        /// computation
        verletStep1Atoms(atoms,data.timeStep);
        //update list before
        list.update(atoms);
        potentialEnergy = gupta(atoms,list);
        returnEnergy = potentialEnergy;
        //last step
        verletStep2Atoms(atoms,data.timeStep);
        //thermostat
        berendsenThermostatEV(atoms, data.targetTemperatur, data.timeStep, data.relaxationTime);
        //add the kinetic energy to the potential
        kineticEnergy = calculateKineticEnergy(atoms);
        returnEnergy += kineticEnergy;
        //calculate the temperature for each step (debugging)
        returnTemperatur = calculateCurrentTemperaturEV(atoms);
        setEnergyStorage.push_back(returnEnergy);
        setTemperaturStorage.push_back(returnTemperatur);
        ////
        if(data.totalEnergyRecording == true) {
            kineticEnergyStorage.push_back(kineticEnergy);
            potentialEnergyStorage.push_back(potentialEnergy);
        }
        ///basic loop stuff last
        currentTime += data.timeStep;
        i++;
    }
    //save the trajectory
    if(data.doDumping == true) {
        dumpData(atoms, data.trajectorySafeLocation, data.trajectoryBaseName, data.maxTrajectoryNumber, data.simulationID);
    }
    //check if valid
    if (checkMoleculeTrajectories(atoms,  data.controlCube) == false) {
        std::cerr << "Cube Exploded at: " << data.simulationID << std::endl;
        //exit the program
        exit(EXIT_FAILURE);
    }
    //return stuff
    return {returnEnergy, returnTemperatur};
}

/**
 * @fn int writeClusters()
 * @brief little mockup to create the files for the Clusters, couse it refuses to work in the normal code
 * @return nothing
 */
int writeClusters() {
    int layers = 12;
    //generate a cluster based on the layers
    std::string location =  "/home/cm/CLionProjects/MoleDymCode/AData/Clusters";
    generateClusterHull(layers, location);
    return 0;
}