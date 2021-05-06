//
// Created by cm on 06.05.21.
//
#include <gtest/gtest.h>
#include "verlet.h"

/**
 * Info:
 * Use EXPECT_NEAR or ASSERT_NEAR bc double!!
 */

TEST(VerletTest, IntegratorCheckConstantForceXYZDirection)
{
    //
    double xNow = 0;
    double yNow = 0;
    double zNow = 0;
    double vxNow = 0;
    double vyNow = 0;
    double vzNow = 0;
    //
    double forceX = 1.0;
    double forceY = 1.0;
    double forceZ = 1.0;
    double timestep = 1.0;
    unsigned int numberOfSteps = 10;
    //
    verletIntegratorConstantForce(xNow,yNow,zNow,
                                  vxNow,vyNow,vzNow,
                                  numberOfSteps,forceX,forceY,forceZ,timestep);
    /** Check position */
    ASSERT_NEAR(xNow,50,1e-6);
    ASSERT_NEAR(yNow,50,1e-6);
    ASSERT_NEAR(zNow,50,1e-6);
    /** Check velocities */
    ASSERT_NEAR(vxNow,10,1e-6);
    ASSERT_NEAR(vyNow,10,1e-6);
    ASSERT_NEAR(vzNow,10,1e-6);
}

TEST(VerletTest, IntegratorCheckConstantForceXYZDirectionInitialVariables)
{
    //
    double xNow = 10;
    double yNow = 10;
    double zNow = 10;
    double vxNow = 10;
    double vyNow = 10;
    double vzNow = 10;
    //
    double forceX = 1.0;
    double forceY = 1.0;
    double forceZ = 1.0;
    double timestep = 1.0;
    unsigned int numberOfSteps = 10;
    //
    verletIntegratorConstantForce(xNow,yNow,zNow,
                                  vxNow,vyNow,vzNow,
                                  numberOfSteps,forceX,forceY,forceZ,timestep);
    /** Check position */
    ASSERT_NEAR(xNow,160,1e-6);
    ASSERT_NEAR(yNow,160,1e-6);
    ASSERT_NEAR(zNow,160,1e-6);
    /** Check velocities */
    ASSERT_NEAR(vxNow,20,1e-6);
    ASSERT_NEAR(vyNow,20,1e-6);
    ASSERT_NEAR(vzNow,20,1e-6);
}

TEST(VerletTest, IntegratorCheckConstantForceXYZDirectionNegativeForce)
{
    //
    double xNow = 0;
    double yNow = 0;
    double zNow = 0;
    double vxNow = 0;
    double vyNow = 0;
    double vzNow = 0;
    //
    double forceX = -1.0;
    double forceY = -1.0;
    double forceZ = -1.0;
    double timestep = 1.0;
    unsigned int numberOfSteps = 10;
    //
    verletIntegratorConstantForce(xNow,yNow,zNow,
                                  vxNow,vyNow,vzNow,
                                  numberOfSteps,forceX,forceY,forceZ,timestep);
    /** Check position */
    ASSERT_NEAR(xNow,-50,1e-6);
    ASSERT_NEAR(yNow,-50,1e-6);
    ASSERT_NEAR(zNow,-50,1e-6);
    /** Check velocities */
    ASSERT_NEAR(vxNow,-10,1e-6);
    ASSERT_NEAR(vyNow,-10,1e-6);
    ASSERT_NEAR(vzNow,-10,1e-6);
}

/**
TEST(Test Suit, Test name)
{
    ASSERT_TRUE(false) << "Implement me!!";
}
*/