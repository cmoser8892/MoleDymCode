//
// Created by cm on 31.05.21.
//

#include <iomanip>
#include <iostream>
#include <fstream>

#include "../Headerfiles/helperfunctions.h"
#include "../Headerfiles/xyz.h"
/**
 * @fn void dumpData(Atoms &atoms,std::string location, std::string namingScheme,unsigned int expectedNumberOfDumps ,unsigned int number)
 * @brief Creates a file in a series of files in .xyz format to be used for the simulation
 * @param atoms
 * @param location
 * @param namingScheme
 * @param expectedNumberOfDumps
 * @param number
 */
void dumpData(Atoms &atoms,std::string location, std::string namingScheme,unsigned int expectedNumberOfDumps ,unsigned int number) {
    //create whole location string
    unsigned int condensedNumber = expectedNumberOfDumps * 10 + number;
    std::ostringstream convert;
    convert << condensedNumber;
    std::string condensedNum = convert.str().erase(0,1);
    std::string total = location +"/" + namingScheme +condensedNum+".xyz";
    //
    write_xyz(total,atoms);
}

/**
 * @fn void dumpEnergy( std::vector<double> data,std::string location, std::string name)
 * @brief dumps the date from the vector data to a file
 * @param data
 * @param location
 * @param name
 */
void dumpEnergy( std::vector<double> data,
                 std::string location, std::string name) {
    std::string total = location + "/" + name +".txt";
    std::ofstream file(total);
    for(int i = 0; i < data.size();++i) {
     file << data[i] << std::endl;
    }
    file.close();
}

/**
 * @fn void setANameInAtoms(Atoms &atoms, char name)
 * @brief sets the name of an atom to X so data can be dumped
 * @param atoms
 * @param name
 */
void setANameInAtoms(Atoms &atoms, std::string name) {
    for (int i = 0; i < atoms.nb_atoms(); ++i) {
        atoms.names[i] = name;
    }
}

/**
 * @fn Positions_t createLatticeCube(unsigned int numberOfAtoms, double latticeConstant)
 * @brief creates positions in an arranged cube kinda inefficient
 * @param numberOfAtoms
 * @param latticeConstant
 * @return Positions of the atoms
 */
Positions_t createLatticeCube(unsigned int numberOfAtoms, double latticeConstant) {
    /**
     * 1. First determin the Side lenght of the Cube
     * 2. Fill it up
     * */
    Positions_t returnValue(3, numberOfAtoms);
    int cubeSideLenght = ceil(pow(numberOfAtoms,1./3.));
    //
    int xCounter = 0;
    int yCounter = 0;
    int zCounter = 0;
    int atomNumberCounter = 0;
    //
    while(zCounter < cubeSideLenght) {
        while(yCounter < cubeSideLenght) {
            while (xCounter < cubeSideLenght) {
                //cancel
                if(atomNumberCounter == numberOfAtoms) {
                    break;
                }
                // Set value all the other stuff just fancy
                Vector_t set(xCounter*latticeConstant,yCounter*latticeConstant,zCounter*latticeConstant);
                returnValue.col(atomNumberCounter) = set;
                atomNumberCounter++;
                //
                xCounter++;
            }
            xCounter = 0;
            yCounter++;
            if(atomNumberCounter == numberOfAtoms) {
                break;
            }
        }
        yCounter = 0;
        zCounter++;
        if(atomNumberCounter == numberOfAtoms) {
            break;
        }
    }
    return returnValue;
}

/**
 * @fn Positions_t createLatticesLongRod(unsigned int numberOfAtoms, unsigned int baseSideLength, double latticeConstant)
 * @brief this creates Positons for a rod with a baseSideLength till the atoms are used up; modified generateLatticeCube Code
 * @param numberOfAtoms
 * @param baseSideLength
 * @param latticeConstant
 * @return Positions of the atoms
 */
Positions_t createLatticesLongRod(unsigned int numberOfAtoms, unsigned int baseSideLength, double latticeConstant) {
    Positions_t returnValue(3, numberOfAtoms);
    //
    int xCounter = 0;
    int yCounter = 0;
    int zCounter = 0;
    int atomNumberCounter = 0;
    // modified cube code !!
    while(1) {
        while(yCounter < baseSideLength) {
            while (xCounter < baseSideLength) {
                //cancel only breakout here
                if(atomNumberCounter == numberOfAtoms) {
                    break;
                }
                // Set value all the other stuff just fancy
                Vector_t set(xCounter*latticeConstant,yCounter*latticeConstant,zCounter*latticeConstant);
                returnValue.col(atomNumberCounter) = set;
                atomNumberCounter++;
                //
                xCounter++;
            }
            xCounter = 0;
            yCounter++;
            if(atomNumberCounter == numberOfAtoms) {
                break;
            }
        }
        yCounter = 0;
        zCounter++;
        if(atomNumberCounter == numberOfAtoms) {
            break;
        }
    }
    //
    return returnValue;
}

/**
 * @fn double calculateKineticEnergy(Atoms &atoms)
 * @brief calculates teh Energy of the Atoms; only really works in kartesian coordinates mass just set to 1
 * @param atoms
 * @return the kinetic Energy in the System
 */
double calculateKineticEnergy(Atoms &atoms) {
    double energyReturn = 0;
    for ( int i =0; i < atoms.nb_atoms(); ++i) {
        double thisEnergy = 0.5 * atoms.mass(i) * (pow(atoms.velocities(0,i),2)
                                       +pow(atoms.velocities(1,i),2)
                                       +pow(atoms.velocities(2,i),2));
        energyReturn += thisEnergy;
    }
    return energyReturn;
}

/**
 * @fn double calculateCurrentTemperatur(Atoms &atoms)
 * @brief calculates the current Temperatur of the simulation using boltzmann
 * @param atoms
 * @return the temperatur of the system
 */
double calculateCurrentTemperatur(Atoms &atoms) {
    double totalKineticEnergy = calculateKineticEnergy(atoms);
    //totalKineticEnergy = calculateEnergyWithQuadradicMeanVelocity(atoms);
    double temperatur = 2./3. * (totalKineticEnergy/boltzmannConstant);
    return temperatur;
}

/**
 * @fn double temperaturDampening(double initalTemperatur, double targetTemperatur, double relaxationTime, double timestep)
 * @brief exponential dampening of the temperatur as a function
 * @param initalTemperatur
 * @param targetTemperatur
 * @param relaxationTime
 * @param timestep
 * @return temperatur that the system should have
 */
double temperaturDampening(double initalTemperatur, double targetTemperatur, double relaxationTime, double timestep) {
    double currentTemperatur = targetTemperatur + (initalTemperatur - targetTemperatur)* exp(-timestep/relaxationTime);
    return currentTemperatur;
}

/**
 * @fn
 * @brief
 * @param atoms
 * @return
 */
double calculateEnergyWithQuadradicMeanVelocity(Atoms &atoms) {
    /** same as other method btw */
    Vector_t totalVelocity(0,0,0);
    double totalMass = 0;
    double kineticEnergy = 0;
    for(int i = 0; i <atoms.nb_atoms(); ++i) {
        totalVelocity += pow(atoms.velocities.col(i),2);
        totalMass += atoms.mass(i);
    }
    totalVelocity = totalVelocity/atoms.nb_atoms();
    kineticEnergy = 0.5 * totalMass *(totalVelocity(0) +totalVelocity(1) +totalVelocity(2));
    return kineticEnergy;
}

/**
 * @fn bool checkMoleculeTrajectories(Atoms &atoms, double scaling)
 * @brief creates a capsul defined by two points min and max and applies a sclaing to the capsul
 * @param atoms
 * @param cubeFactor
 * @return true or false if ALL the atoms are inside the capsule
 */
bool checkMoleculeTrajectories(Atoms &atoms, double scaling) {
    /** Basic idea, check weather or not the atoms fling themself outside of another cube thats a bit bigger than the original one */
    Positions_t controlPostitions = generateCapsel(atoms, 2);
    bool returnValue = true;

    for(int i = 0; i < atoms.nb_atoms();++i) {
        /** check for every position if it still is inside of the cube
         * as it is alinged just need to check for xmin and xmax
         * */
        //minum Position
        if(compareVectorsBigSmall(controlPostitions.col(0),atoms.positions.col(i)) != true) {
            returnValue = false;
            break;
        }
        //maximum Position
        if(compareVectorsBigSmall(atoms.positions.col(i),controlPostitions.col(1)) != true) {
            returnValue = false;
            break;
        }

    }
    return returnValue;
}

/**
 * @fn Positions_t generateCapsel(Atoms &atoms, double cubeFactor)
 * @brief creates two points of capsel around the structure so that the atoms can be checked to not escape
 * @param atoms
 * @param cubeFactor
 * @return
 */
Positions_t generateCapsel(Atoms &atoms, double scaling) {
    Positions_t returnValue(3,2);
    returnValue.setZero();
    Vector_t min{atoms.positions.row(0).minCoeff(),atoms.positions.row(1).minCoeff(),atoms.positions.row(2).minCoeff()};
    Vector_t max{atoms.positions.row(0).maxCoeff(),atoms.positions.row(1).maxCoeff(),atoms.positions.row(2).maxCoeff()};
    Vector_t middlePoint = 0.5*(max -min);
    //scaling is done from the middlepoint
    Vector_t middleToMin = middlePoint + scaling * (min- middlePoint);
    Vector_t middleToMax = middlePoint + scaling * (max- middlePoint) ;
    //write back
    returnValue.col(0) = middleToMin;
    returnValue.col(1) = middleToMax;
    return returnValue;
}

/**
 * @fn bool compareVectorsBigSmall(Vector_t v1, Vector_t v2)
 * @brief checks weather v1 is bigger than v2 in all directions
 * @param v1
 * @param v2
 * @return true or false
 */
bool compareVectorsBigSmall(Vector_t v1, Vector_t v2) {
    bool returnValue = true;
    if(v1(0) < v2(0)) {
        returnValue = false;
    }
    if(v1(1) < v2(1)) {
        returnValue = false;
    }
    if(v1(2) < v2(2)) {
        returnValue = false;
    }
    return returnValue;
}

/**
 * @fn double calculateDistanceBetweenVectors(Vector_t distanceVector)
 * @brief calculates the length of the Vector
 * @param distanceVector
 * @return double the length of the Vektor
 */
double calculateDistanceBetweenVectors(Vector_t distanceVector) {
    //pythagoras
    double dist = distanceVector(0)* distanceVector(0) +
                  distanceVector(1)* distanceVector(1) +
                  distanceVector(2)* distanceVector(2);
    return sqrt(dist);
}