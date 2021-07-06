//
// Created by cm on 06.07.21.
//

#include <gtest/gtest.h>
#include "../Headerfiles/milestones.h"
#include "../Headerfiles/atoms.h"

TEST(TestMilestone, TestSimulationBlock) {
    //test for the dumping first
    //then just check weather the initial temperatur is smaller than later
    //maybe do some death tests with bad variables ?! dont know yet
    Positions_t positions = createLatticeCube(8);
    Atoms atoms(positions,10);
    SimulationData_t data;
    data.controlCube = generateCapsel(atoms,2); //allways has to be generated otherwise crash
    //
    auto[energy,temperatur]{simulationBuildStone(data,atoms)};

    ASSERT_TRUE(false) << "Implement me!!";
}
