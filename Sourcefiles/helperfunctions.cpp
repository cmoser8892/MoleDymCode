//
// Created by cm on 31.05.21.
//

#include <iomanip>
#include <iostream>
#include <fstream>

#include "../Headerfiles/helperfunctions.h"
#include "../Headerfiles/xyz.h"

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

void dumpEnergy( std::vector<double> data,
                 std::string location, std::string name) {
    std::string total = location + "/" + name +".txt";
    std::ofstream file(total);
    for(int i = 0; i < data.size();++i) {
     file << data[i] << std::endl;
    }
    file.close();
}