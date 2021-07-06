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
//all in a struct so it is still readable and you dont confuse the stuff as easily
typedef struct SimulationData {
    unsigned int simulationID = 0;

    ///basic simulation stuff
    Positions_t controlCube;
    double timeStep = 1e-15;
    double simulationTime = 0;
    double relaxationTime = 0;
    double cutoffDistance = 0;
    double targetTemperatur = 0;

    ///data related stuff
    unsigned int maxTrajectoryNumber = 0;
    std::string trajectorySafeLocation = "";
    std::string trajectoryBaseName = "";
}SimulationData_t;

std::tuple<double, double > simulationBuildStone(SimulationData_t &data, Atoms &atoms);

#endif //MYPROJECT_MILESTONES_H
