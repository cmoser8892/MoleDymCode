//
// Created by cm on 02.06.21.
//

#include "../Headerfiles/berendsenThermostat.h"
#include "../Headerfiles/helperfunctions.h"
#include <iostream>


/**
 * @fn void berendsenThermostat(Atoms &atoms, double targetTemperature, double timestep, double relaxationTime)
 * @brief performs a velocity rescaling using the brendsen Equations
 * @param atoms
 * @param targetTemperature
 * @param timestep
 * @param relaxationTime
 */
void berendsenThermostat(Atoms &atoms, double targetTemperature, double timestep, double relaxationTime) {
    //initialization
    double currentTemperatur = calculateCurrentTemperatur(atoms);
    unsigned int numberOfAtoms = atoms.nb_atoms();
    //calculate the rescaling factor
    double rescalingFactor = sqrt(1 + (targetTemperature/currentTemperatur -1) * (timestep/relaxationTime));

    //rescale the velocities
    atoms.velocities = rescalingFactor * atoms.velocities;
}
