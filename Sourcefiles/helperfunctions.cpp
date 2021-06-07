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
void setANameInAtoms(Atoms &atoms, char name) {
    for (int i = 0; i < atoms.nb_atoms(); ++i) {
        atoms.names[i] = name;
    }
}

/**
 * @fn Positions_t createLatticeCube(unsigned int numberOfAtoms, double latticeConstant)
 * @brief creates positions in an arranged cube kinda inefficient
 * @param numberOfAtoms
 * @param latticeConstant
 * @return
 */
Positions_t createLatticeCube(unsigned int numberOfAtoms, double latticeConstant) {
    /**
     * 1. First determin the Side lenght of the Cube
     * 2. Fill it up
     * */
    Positions_t returnValue(3, numberOfAtoms);
    int cubeSideLenght = pow(numberOfAtoms,1./3.) +1;
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
        }
        yCounter = 0;
        zCounter++;
    }
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


bool checkMoleculeTrajectories(Atoms &atoms) {

}