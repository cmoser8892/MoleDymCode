//
// Created by cm on 07.06.21.
//

#include <gtest/gtest.h>
#include "../Headerfiles/berendsenThermostat.h"
#include "../Headerfiles/helperfunctions.h"

/** Test */
TEST(BerendsenExponentialDamp, ThermostatTest) {
    /** Init */
    unsigned int nbAtoms = 50;
    double distance = 1;
    double timestep = 1;
    double relaxationTime = 1000000;
    double simulationTime = 1000;
    double targetTemperatur = 0;
    /**  set  */
    Positions_t p = createLatticeCube(nbAtoms,distance);
    Atoms atoms(p,atomicUnit);
    atoms.velocities.setOnes();
    double initialTemperatur = calculateCurrentTemperatur(atoms);
    /** run a small simulation */
    double currentTime = 0;
    while(currentTime < simulationTime) {
        double testTemperatur = temperaturDampening(initialTemperatur,targetTemperatur,relaxationTime,timestep);
        berendsenThermostat(atoms,targetTemperatur,timestep,relaxationTime);
        double currentTemperatur = calculateCurrentTemperatur(atoms);
        EXPECT_NEAR(testTemperatur,currentTemperatur,1e-5);
        currentTime += timestep;
    }
}

TEST(BerendsenExponentialDamp, ThermostatTest1Atom) {
    /** Init */
    unsigned int nbAtoms = 1;
    double distance = 1;
    double timestep = 1;
    double relaxationTime = 100000;
    double simulationTime = 1000;
    double targetTemperatur = 0;
    /**  set  */
    Positions_t p = createLatticeCube(nbAtoms,distance);
    Atoms atoms(p,atomicUnit);
    atoms.velocities.setOnes();
    double initialTemperatur = calculateCurrentTemperatur(atoms);
    /** run a small simulation */
    double currentTime = 0;
    while(currentTime < simulationTime) {
        double testTemperatur = temperaturDampening(initialTemperatur,targetTemperatur,relaxationTime,timestep);
        berendsenThermostat(atoms,targetTemperatur,timestep,relaxationTime);
        double currentTemperatur = calculateCurrentTemperatur(atoms);
        EXPECT_NEAR(testTemperatur,currentTemperatur,1e-5);
        currentTime += timestep;
    }
}

TEST(BerendsenExponentialDamp, NilTest) {
    // a test to see that the temperatur is not corrupting data with /0
    unsigned int nbAtoms = 8;
    Positions_t p = createLatticeCube(nbAtoms);
    Atoms atoms(p,5);
    double timeStep = 1e-15;
    double relaxationTime = 100* timeStep;
    double targetTemperatur = 275;
    berendsenThermostat(atoms,targetTemperatur,timeStep,relaxationTime);
    for(int i = 0; i < nbAtoms; ++i) {
        for( int j = 0; j < 3; ++j) {
            EXPECT_EQ(atoms.velocities(j,i),0);
        }
    }
}

