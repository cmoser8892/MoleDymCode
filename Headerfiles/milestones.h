//
// Created by cm on 06.06.21.
//

#ifndef MYPROJECT_MILESTONES_H
#define MYPROJECT_MILESTONES_H

/** Include Everything */
#include "../Headerfiles/verlet.h"
#include "../Headerfiles/types.h"
#include "../Headerfiles//lenardJonesDirectionSummation.h"
#include "../Headerfiles/xyz.h"
#include "../Headerfiles/helperfunctions.h"
#include "../Headerfiles/milestones.h"

int milestone4Code();
int milestone5Code(int argc = 0, char *argv[] = NULL);
int milestone6Code(int argc = 0, char *argv[] = NULL);
int milestone7Code(int argc = 0, char *argv[] = NULL);

/** Big Helper */
typedef struct SimulationData {
    unsigned int simulationID;

    ///basic simulation stuff
    Atoms &atoms;
    double timeStep;
    double simulationTime;
    double relaxationTime;
    double cutoffDistance;
    double targetTemperatur;

    ///data related stuff
    unsigned int maxTrajectoryNumber;
    std::string trajectorySafeLocation;
    std::string trajectoryBaseName;
}SimulationData_t;

std::tuple<double energy, double temperature> simulationBuildStone(SimulationData_t data);

#endif //MYPROJECT_MILESTONES_H
