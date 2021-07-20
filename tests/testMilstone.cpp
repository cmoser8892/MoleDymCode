//
// Created by cm on 06.07.21.
//

#include <gtest/gtest.h>
#include <filesystem>

#include "../Headerfiles/milestones.h"

TEST(TestMilestone, TestSimulationBlock) {
    //test for the dumping first
    //then just check weather the initial temperatur is smaller than later
    double atomicMassAu = 196.97; // 197Au79
    double mass = atomicMassAu / massCorrectionFactor;
    Positions_t positions = createLatticeCube(8);
    Atoms atoms(positions,mass);
    SimulationData_t data;
    ///
    data.simulationID = 1;
    data.maxTrajectoryNumber = 100;
    data.trajectorySafeLocation = "/home/cm/CLionProjects/MoleDymCode/cmake-build-debug/tests/TestDumps";
    data.trajectoryBaseName = "Trajectory";
    ///
    data.controlCube = generateCapsel(atoms,10); //allways has to be generated otherwise crash
    data.timeStep = 1;
    data.simulationTime = 10 * data.timeStep;
    data.relaxationTime = data.simulationTime;
    data.cutoffDistance = 3.0;
    data.targetTemperatur = 296;
    //check if the file is constructed correctly
    //first delete stuff
    std::filesystem::remove_all("TestDumps");
    std::filesystem::create_directory("TestDumps");
    {
        data.simulationID = 1;
        data.maxTrajectoryNumber = 100;
        std::filesystem::path path{"TestDumps/Trajectory001.xyz"};
        auto[energy,temperatur]{simulationBuildStone(data,atoms)};
        ASSERT_FALSE(std::filesystem::exists(path));
    }
    {
        data.simulationID = 55;
        data.maxTrajectoryNumber = 100;
        data.doDumping = true;
        std::filesystem::path path{"TestDumps/Trajectory055.xyz"};
        auto[energy,temperatur]{simulationBuildStone(data,atoms)};
        EXPECT_TRUE(std::filesystem::exists(path));
    }
    {
        data.simulationTime = 2*data.timeStep;
        auto[energy,temperatur]{simulationBuildStone(data,atoms)};
    }
}
