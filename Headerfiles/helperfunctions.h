//
// Created by cm on 31.05.21.
//

#ifndef MYPROJECT_HELPERFUNCTIONS_H
#define MYPROJECT_HELPERFUNCTIONS_H

#include "../Headerfiles/atoms.h"

/** Data dump*/
void dumpData( Atoms &atoms, std::string location, std::string namingScheme, unsigned int expectedNumberOfDumps ,unsigned int number);
void dumpEnergy( std::vector<double> data, std::string location, std::string name);
void setANameInAtoms(Atoms &atoms, std::string name);

/**  calculation Functions */
Positions_t createLatticeCube(unsigned int numberOfAtoms, double latticeConstant = 1.0);
Positions_t generateLatticesLongRod(unsigned int numberOfAtoms, unsigned int baseSideLength = 3, double latticeConstant = 1.0);
double calculateKineticEnergy(Atoms &atoms);
double calculateCurrentTemperatur(Atoms &atoms);

/** Other helpers*/
double temperaturDampening(double initalTemperatur, double targetTemperatur, double relaxationTime, double timestep);
double calculateEnergyWithQuadradicMeanVelocity(Atoms &atoms);
bool checkMoleculeTrajectories(Atoms &atoms, double cubeFactor);
Positions_t generateCube(Atoms &atoms, double cubeFactor);
Positions_t generateCapsle(Atoms &atoms, double cubeFactor);

#endif //MYPROJECT_HELPERFUNCTIONS_H
