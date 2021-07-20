//
// Created by cm on 31.05.21.
//

#ifndef MYPROJECT_HELPERFUNCTIONS_H
#define MYPROJECT_HELPERFUNCTIONS_H

#include "../Headerfiles/atoms.h"

/** Data dump*/
void dumpData( Atoms &atoms, std::string location, std::string namingScheme, unsigned int expectedNumberOfDumps ,unsigned int number);
void dumpVectorData(std::vector<double> data, std::string location, std::string name);
void setANameInAtoms(Atoms &atoms, std::string name);

/**  calculation Functions */
Positions_t createLatticeCube(unsigned int numberOfAtoms, double latticeConstant = 1.0);
Positions_t createLatticesLongRod(unsigned int numberOfAtoms, unsigned int baseSideLength = 3, double latticeConstant = 1.0);
double calculateKineticEnergy(Atoms &atoms);
double calculateCurrentTemperatur(Atoms &atoms);
double calculateCurrentTemperaturEV(Atoms &atoms);

/** Other helpers*/
double temperaturDampening(double initalTemperatur, double targetTemperatur, double relaxationTime, double timestep);
double calculateEnergyWithQuadradicMeanVelocity(Atoms &atoms);
bool checkMoleculeTrajectories(Atoms &atoms, Positions_t controlPositions);
Positions_t generateCapsel(Atoms &atoms, double scaling);
bool compareVectorsBigSmall(Vector_t v1, Vector_t v2);
double calculateDistanceBetweenVectors(Vector_t distanceVector);
void depositHeat(double heat, Atoms &atoms);
void depositRescaledHeat(double heat, Atoms &atoms);
void printAtomsVelocitiesAndPositions(Atoms &atoms);
double averageVector(std::vector<double> values);
void printData(int step, double energy, double temperatur);
void generateClusterHull(unsigned int layers, std::string location);

#endif //MYPROJECT_HELPERFUNCTIONS_H
