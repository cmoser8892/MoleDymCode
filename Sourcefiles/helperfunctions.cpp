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
 * @fn
 * @brief
 * @param numberOfAtoms
 * @param latticeConstant
 * @return
 */
Positions_t createLatticeCube(unsigned int numberOfAtoms, double latticeConstant) {
    /**
     * 1. First determin the Side lenght of the Cube
     * 2. Fill it up
     * */
}