//
// Created by cm on 12.05.21.
//

#ifndef MYPROJECT_TYPES_H
#define MYPROJECT_TYPES_H

/** make some type changes to the Eigennames so that they are more readable */
#include <iostream>
#include <string>
#include <vector>

#include "Eigen/Core"
/** Redefinitions*/
using Positions_t = Eigen::Array3Xd;
using Velocities_t = Eigen::Array3Xd;
using Forces_t = Eigen::Array3Xd;

using Vector_t = Eigen::Array3d;

using Names_t = std::vector<std::string>;

/** Constants */
static double boltzmannConstant = 1.380649e-23;


#endif //MYPROJECT_TYPES_H
