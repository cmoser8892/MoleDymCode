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

/** Build stuff */
//#define NOPRINTING

/** Redefinitions*/
using Positions_t = Eigen::Array3Xd;
using Velocities_t = Eigen::Array3Xd;
using Forces_t = Eigen::Array3Xd;
using Mass_t = Eigen::ArrayXd;

using Vector_t = Eigen::Array3d;

using Names_t = std::vector<std::string>;

using Energies_t = Eigen::ArrayXd;

/** Constants */
static const double boltzmannConstant = 1.380649e-23; ///J/K
static const double atomicUnit = 1.660539e-27; ///kg
static const double electronVolt = 1.602176634e-19; ///J
static const double massCorrectionFactor = 1.6e-29/atomicUnit; ///kg given if the timestep is 1fs

static const double boltzmannElectronVolt = boltzmannConstant/electronVolt; ///eV/K

#endif //MYPROJECT_TYPES_H
