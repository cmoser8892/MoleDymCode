//
// Created by cm on 31.05.21.
//

#ifndef MYPROJECT_HELPERFUNCTIONS_H
#define MYPROJECT_HELPERFUNCTIONS_H

#include "../Headerfiles/atoms.h"

void dumpData( Atoms &atoms,
        std::string location, std::string namingScheme, unsigned int expectedNumberOfDumps ,unsigned int number);

void dumpEnergy( std::vector<double> data,
                 std::string location, std::string name);

Positions_t createLatticeCube(unsigned int numberOfAtoms, double latticeConstant);
#endif //MYPROJECT_HELPERFUNCTIONS_H
