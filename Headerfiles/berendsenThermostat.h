//
// Created by cm on 02.06.21.
//

#ifndef MYPROJECT_BERENDSENTHERMOSTAT_H
#define MYPROJECT_BERENDSENTHERMOSTAT_H

#include "../Headerfiles/atoms.h"

void berendsenThermostat(Atoms &atoms, double targetTemperature, double timestep, double relaxationTime);
void berendsenThermostatEV(Atoms &atoms, double targetTemperature, double timestep, double relaxationTime);
#endif //MYPROJECT_BERENDSENTHERMOSTAT_H
