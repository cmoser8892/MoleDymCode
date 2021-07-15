//
// Created by cm on 24.06.21.
//

#include <gtest/gtest.h>
//
#include "../Headerfiles/atoms.h"
#include "../Headerfiles/helperfunctions.h"
#include "../Headerfiles/xyz.h"

TEST(TestHelperfunctions, cubicEncapsulation) {
    Positions_t p = createLatticeCube(8,10);
    Atoms atoms(p);
    Vector_t min{atoms.positions.row(0).minCoeff(),atoms.positions.row(1).minCoeff(),atoms.positions.row(2).minCoeff()};
    Vector_t max{atoms.positions.row(0).maxCoeff(),atoms.positions.row(1).maxCoeff(),atoms.positions.row(2).maxCoeff()};
    Positions_t controlPostitions = generateCapsel(atoms, 2);

    //check weather all atoms are inside
    for(int i = 0; i < atoms.nb_atoms(); ++i) {
        //should be easier too just compare the values of the Vectors
        //the positions should be smaller than the max
        EXPECT_EQ(compareVectorsBigSmall(controlPostitions.col(1),atoms.positions.col(i)),true);
        //the positions should be smaller than the min
        EXPECT_EQ(compareVectorsBigSmall(atoms.positions.col(i),controlPostitions.col(0)),true);
    }
}

TEST(TestHelperfunctions, encapsulatingFunction) {
    Positions_t p = createLatticesLongRod(8,2,10);
    Atoms atoms(p);
    Vector_t min{atoms.positions.row(0).minCoeff(),atoms.positions.row(1).minCoeff(),atoms.positions.row(2).minCoeff()};
    Vector_t max{atoms.positions.row(0).maxCoeff(),atoms.positions.row(1).maxCoeff(),atoms.positions.row(2).maxCoeff()};
    Positions_t controlPostitions = generateCapsel(atoms, 2);

    //check weather all atoms are inside
    for(int i = 0; i < atoms.nb_atoms(); ++i) {
        //should be easier too just compare the values of the Vectors
        //the positions should be smaller than the max
        EXPECT_EQ(compareVectorsBigSmall(controlPostitions.col(1),atoms.positions.col(i)),true);
        //the positions should be smaller than the min
        EXPECT_EQ(compareVectorsBigSmall(atoms.positions.col(i),controlPostitions.col(0)),true);
    }
}

TEST(TestHelperfunctions, encapsulatingFunctionWithHexagon) {
    double mass = 196.97*atomicUnit; // 197Au79
    auto [names, positions]{read_xyz("../../AData/cluster_923.xyz")};
    Atoms atoms(names,positions,mass);
    Vector_t min{atoms.positions.row(0).minCoeff(),atoms.positions.row(1).minCoeff(),atoms.positions.row(2).minCoeff()};
    Vector_t max{atoms.positions.row(0).maxCoeff(),atoms.positions.row(1).maxCoeff(),atoms.positions.row(2).maxCoeff()};
    Positions_t controlPostitions = generateCapsel(atoms, 1.1);

    //check weather all atoms are inside
    for(int i = 0; i < atoms.nb_atoms(); ++i) {
        //should be easier too just compare the values of the Vectors
        //the positions should be smaller than the max
        EXPECT_EQ(compareVectorsBigSmall(controlPostitions.col(1),atoms.positions.col(i)),true);
        //the positions should be smaller than the min
        EXPECT_EQ(compareVectorsBigSmall(atoms.positions.col(i),controlPostitions.col(0)),true);
    }
}


TEST(TestHelperfunctions, vectorCompare) {
    {
        //standard cases
        Vector_t v1(0, 0, 0);
        Vector_t v2(1, 1, 1);
        ASSERT_EQ(compareVectorsBigSmall(v2, v1), true);
        ASSERT_EQ(compareVectorsBigSmall(v1, v2), false);
    }
    {
        //minus
        Vector_t v1(0, 0, 0);
        Vector_t v2(-1, -1, -1);
        ASSERT_EQ(compareVectorsBigSmall(v1, v2), true);
        ASSERT_EQ(compareVectorsBigSmall(v2, v1), false);
    }
    {
        //minus bild test
        Vector_t v1(0, 0, 0);
        Vector_t v2(-1, 1, -1);
        ASSERT_EQ(compareVectorsBigSmall(v1, v2), false);
        ASSERT_EQ(compareVectorsBigSmall(v2, v1), false);
    }
    {
        //edge test
        Vector_t v1(1, 1, 1);
        Vector_t v2(1, 1, 1);
        ASSERT_EQ(compareVectorsBigSmall(v1, v2), true);
        ASSERT_EQ(compareVectorsBigSmall(v2, v1), true);
    }
    {
        //just x bigger
        Vector_t v1(1, 1, 1);
        Vector_t v2(2, 1, 1);
        ASSERT_EQ(compareVectorsBigSmall(v1, v2), false);
        ASSERT_EQ(compareVectorsBigSmall(v2, v1), true);
    }
    {
        //just y bigger
        Vector_t v1(1, 1, 1);
        Vector_t v2(1, 2, 1);
        ASSERT_EQ(compareVectorsBigSmall(v1, v2), false);
        ASSERT_EQ(compareVectorsBigSmall(v2, v1), true);
    }
    {
        //just y bigger
        Vector_t v1(1, 1, 1);
        Vector_t v2(1, 1, 2);
        ASSERT_EQ(compareVectorsBigSmall(v1, v2), false);
        ASSERT_EQ(compareVectorsBigSmall(v2, v1), true);
    }
}

TEST(TestHelperfunctions, testcheckTrajectories) {
    //create postions
    Atoms atoms(createLatticeCube(8));
    //create the capsul
    Positions_t controlCube = generateCapsel(atoms, 1.5);
    //check
    EXPECT_EQ(checkMoleculeTrajectories(atoms,controlCube),true);
    //move on atom
    atoms.positions.col(0) = 10;
    EXPECT_EQ(checkMoleculeTrajectories(atoms,controlCube),false);
    //move on atom
    atoms.positions.col(0) = -10;
    EXPECT_EQ(checkMoleculeTrajectories(atoms,controlCube),false);
    //move back
    atoms.positions.col(0) = 0;
    EXPECT_EQ(checkMoleculeTrajectories(atoms,controlCube),true);
    //create new control
    controlCube = generateCapsel(atoms, 1.0);
    EXPECT_EQ(checkMoleculeTrajectories(atoms,controlCube),true);
    //make too small
    controlCube = generateCapsel(atoms, 0.9);
    EXPECT_EQ(checkMoleculeTrajectories(atoms,controlCube),false);
    //make way to big
    controlCube = generateCapsel(atoms, 100);
    //std::cout << controlCube << std::endl;
    atoms.positions.col(0) = 50;
    EXPECT_EQ(checkMoleculeTrajectories(atoms,controlCube), true);
}

TEST(TestHelperfunctions, heatdepostitTest) {
    double energy = 10;
    Positions_t lattice = createLatticeCube(8);
    //constant mass
    Atoms atoms(lattice,1);
    depositHeat(energy,atoms);
    EXPECT_EQ(calculateCurrentTemperaturEV(atoms),2./3. * (energy/boltzmannElectronVolt));
}

TEST(TestHelperfunctions, basicHeatEnergyTest) {
    double energy = 10;
    Positions_t lattice = createLatticeCube(8);
    //constant mass
    Atoms atoms(lattice,1);
    atoms.velocities.setOnes();
    double originalValue = calculateKineticEnergy(atoms);
    depositRescaledHeat(energy,atoms);
    double value = calculateKineticEnergy(atoms);
    EXPECT_NEAR(originalValue+energy,value,1e-5);
    //second Part energy is zero
    atoms.velocities.setZero();
    originalValue = calculateKineticEnergy(atoms);
    depositRescaledHeat(energy,atoms);
    value = calculateKineticEnergy(atoms);
    EXPECT_NEAR(originalValue+energy,value,1e-5);
}

TEST(TestHelperfunctions, heatEnergyTest) {
    double energy = 1;
    Positions_t lattice = createLatticeCube(8);
    //constant mass
    Atoms atoms(lattice,1);
    atoms.velocities.setRandom();
    double originalValue = calculateKineticEnergy(atoms);
    depositRescaledHeat(energy,atoms);
    double value = calculateKineticEnergy(atoms);
    EXPECT_NEAR(originalValue+energy,value,1e-5);
}