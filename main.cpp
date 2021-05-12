#include <iostream>

#include "Headerfiles/verlet.h"
#include "Eigen/Dense"
using Positions_t = Eigen::Array3Xd;


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
    /** main loop */
    verletIntegratorConstantForce(xNow,yNow,zNow,
                                  vxNow,vyNow,vzNow,
                                  10,forceX,forceY,forceZ,timestep);
    return 0;
}
