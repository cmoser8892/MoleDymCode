#include <iostream>


#include "Headerfiles/milestones.h"

using namespace std;

/**
 * Notes:
 * TEST FIRST for the thermostat
 *
 */
int main() {
    std::cout << "Molecular Dynamics Project" << std::endl;
    /** Variables */
    unsigned int nbAtoms = 10;
    double sigma = 1.0;
    double epsilon = 1.0;
    double latticeConstant = 1* sigma;
    /** Calls */
    Positions_t p = createLatticeCube(nbAtoms, latticeConstant);

    std::cout << p << std::endl;


    return 0;
    //return milestone4Code();
}
