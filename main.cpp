#include <iostream>

#include "Headerfiles/verlet.h"
#include "Headerfiles/types.h"

int main() {
    std::cout << "Molecular Dynamics Project" << std::endl;
    //
    double xNow = 0;
    double yNow = 0;
    double zNow = 0;
    double vxNow = 0;
    double vyNow = 0;
    double vzNow = 0;
    //
    double forceX = 1.0;
    double forceY = 0.0;
    double forceZ = 0.0;
    double timestep = 1.0;
    //
    Positions_t x;
    int nb_atoms = 10;
    Positions_t potions(3,nb_atoms);
    potions(2,1) = 1.0;
    auto pos2{potions.col(1)};
    Positions_t distance_vector{potions.col(0) - potions.col(1)};
    /** main loop */
    verletIntegratorConstantForce(xNow,yNow,zNow,
                                  vxNow,vyNow,vzNow,
                                  10,forceX,forceY,forceZ,timestep);
    return 0;
}
