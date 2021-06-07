//
// Created by cm on 06.05.21.
//
#include <gtest/gtest.h>
#include "Headerfiles/verlet.h"

/**
 * Info:
 * Use EXPECT_NEAR or ASSERT_NEAR bc double!!
 */

TEST(VerletTest, IntegratorCheckConstantForceXYZDirectionSingleAtom)
{
    //single Atom test
    int nbAtoms = 1;
    Positions_t  positions(3,nbAtoms);
    Velocities_t  velocities(3, nbAtoms);
    Forces_t forces(3,nbAtoms);
    //
    positions.setZero();
    velocities.setZero();
    forces.setZero();
    forces.row(0) = 1.0;
    forces.row(0) = 1.0;
    double timestep = 1.0;
    //run fkt
    verletIntegratorConstantForce(positions,velocities,forces, timestep, 10);
    //check
    ASSERT_NEAR(positions(0),50,1e-6);
    ASSERT_NEAR(velocities(0),10,1e-6);
}

TEST(VerletTest, IntegratorCheckConstantForceXYZDirectionMultipleAtom)
{
    //single Atom test
    int nbAtoms = 10;
    Positions_t  positions(3,nbAtoms);
    Velocities_t  velocities(3, nbAtoms);
    Forces_t forces(3,nbAtoms);
    //
    positions.setZero();
    velocities.setZero();
    forces.setZero();
    forces.row(0) = 1.0;
    forces.row(1) = 1.0;
    forces.row(2) = 1.0;
    double timestep = 1.0;
    //run fkt
    verletIntegratorConstantForce(positions,velocities,forces, timestep, 10);
    //check
    for(int i = 0; i < nbAtoms; ++i)
    {
        //check positions
        ASSERT_NEAR(positions(0,i),50,1e-6);
        ASSERT_NEAR(positions(1,i),50,1e-6);
        ASSERT_NEAR(positions(2,i),50,1e-6);
        //check velocities
        ASSERT_NEAR(velocities(0,i),10,1e-6);
        ASSERT_NEAR(velocities(0,i),10,1e-6);
        ASSERT_NEAR(velocities(0,i),10,1e-6);
    }
}

TEST(VerletTest, IntegratorCheckConstantForceXYZDirectionMultipleAtomNewMethod)
{
    //single Atom test
    int nbAtoms = 10;
    Atoms atoms(10);
    atoms.forces.row(0) = 1.0;
    atoms.forces.row(1) = 1.0;
    atoms.forces.row(2) = 1.0;
    double timestep = 1.0;
    //run fkt
    verletIntegratorConstantForceAtoms(atoms, timestep, 10);
    //check
    for(int i = 0; i < nbAtoms; ++i)
    {
        //check positions
        ASSERT_NEAR(atoms.positions(0,i),50,1e-6);
        ASSERT_NEAR(atoms.positions(1,i),50,1e-6);
        ASSERT_NEAR(atoms.positions(2,i),50,1e-6);
        //check velocities
        ASSERT_NEAR(atoms.velocities(0,i),10,1e-6);
        ASSERT_NEAR(atoms.velocities(0,i),10,1e-6);
        ASSERT_NEAR(atoms.velocities(0,i),10,1e-6);
    }
}
