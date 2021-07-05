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
    //calculate the rescaling factor
    double rescalingFactor = sqrt(1 + (targetTemperature/currentTemperatur -1) * (timestep/relaxationTime));
    //rescale the velocities
    atoms.velocities = rescalingFactor * atoms.velocities;
}

/**
 * @fn void berendsenThermostatEV(Atoms &atoms, double targetTemperature, double timestep, double relaxationTime)
 * @brief performs a velocity rescaling using the brendsen Equations with another variable
 * @param atoms
 * @param targetTemperature
 * @param timestep
 * @param relaxationTime
 */
void berendsenThermostatEV(Atoms &atoms, double targetTemperature, double timestep, double relaxationTime) {
    //calculate the rescaling factor
    double rescalingFactor = sqrt(1 + (targetTemperature/calculateCurrentTemperaturEV(atoms) -1) * (timestep/relaxationTime));
    //rescale the velocities
    atoms.velocities = rescalingFactor * atoms.velocities;
}