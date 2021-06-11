/*
* Copyright 2021 Lars Pastewka
*
* ### MIT license
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in
* all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*/

#include <gtest/gtest.h>

#include "../Headerfiles/lenardJonesDirectionSummation.h"
#include "../Headerfiles/helperfunctions.h"

TEST(LJDirectSummationTest, Forces) {
    constexpr int nb_atoms = 10;
    constexpr double epsilon = 0.7;  // choose different to 1 to pick up missing factors
    constexpr double sigma = 0.3;
    constexpr double delta = 0.00001;  // difference used for numerical (finite difference) computation of forces

    Atoms atoms(nb_atoms);
    atoms.positions.setRandom();  // random numbers between -1 and 1

    // compute and store energy of the indisturbed configuration
    double e0{lendardJonesDirectSummation(atoms, epsilon, sigma)};
    Forces_t forces0{atoms.forces};

    // loop over all atoms and compute forces from a finite differences approximation
    Forces_t dummy_forces(3, nb_atoms);  // we don't actually need these
    for (int i{0}; i < nb_atoms; ++i) {
        // loop over all Cartesian directions
        for (int j{0}; j < 3; ++j) {
            // move atom to the right
            atoms.positions(j, i) += delta;
            double eplus{lendardJonesDirectSummation(atoms, epsilon, sigma)};
            // move atom to the left
            atoms.positions(j, i) -= 2 * delta;
            double eminus{lendardJonesDirectSummation(atoms, epsilon, sigma)};
            // move atom back to original position
            atoms.positions(j, i) += delta;

            // finite-differences forces
            double fd_force{-(eplus - eminus) / (2 * delta)};

            // check whether finite-difference and analytic forces agree
            if (abs(forces0(j, i)) > 1e-10) {
                EXPECT_NEAR(abs(fd_force - forces0(j, i)) / forces0(j, i), 0, 1e-5);

            } else {
                EXPECT_NEAR(fd_force, forces0(j, i), 1e-10);
            }
        }
    }
}

//own Tests
TEST(LJDirectSummationTest,ForceMinimumTwoAtoms) {
    //bad implementation i know, but too lazy to write functions
    /** at the point sigma = 2**1/6 the force needs to be near zero; This test performs the test in all the cardinal directions first
     * and then with a random weight */
    unsigned int nbAtoms = 2;
    Atoms atoms(nbAtoms);
    double testSigma = 1;
    double testEpsilon = 0.3;
    double minimumDistance = pow(2.0, 1.0/6.0) * testSigma;

    /**  three cardinal directions */
    /** x-Direction */
    //calc
    atoms.positions.setZero();
    atoms.forces.setZero();
    atoms.positions(0,1) = minimumDistance;
    lendardJonesDirectSummation(atoms,testEpsilon,testSigma);
    //check
    EXPECT_NEAR(atoms.forces(0,0),0,1e-5);
    EXPECT_NE(atoms.forces(0,0),0);
    EXPECT_EQ(atoms.forces(1,0),0);
    EXPECT_EQ(atoms.forces(2,0),0);

    /** y-Direction */
    //calc
    atoms.positions.setZero();
    atoms.forces.setZero();
    atoms.positions(1,1) = minimumDistance;
    lendardJonesDirectSummation(atoms,testEpsilon,testSigma);
    //check
    EXPECT_NEAR(atoms.forces(1,0),0,1e-5);
    EXPECT_NE(atoms.forces(1,0),0);
    EXPECT_EQ(atoms.forces(0,0),0);
    EXPECT_EQ(atoms.forces(2,0),0);

    /** z-Direction */
    //calc
    atoms.positions.setZero();
    atoms.forces.setZero();
    atoms.positions(2,1) = minimumDistance;
    lendardJonesDirectSummation(atoms,testEpsilon,testSigma);
    //check
    EXPECT_NEAR(atoms.forces(2,0),0,1e-5);
    EXPECT_NE(atoms.forces(2,0),0);
    EXPECT_EQ(atoms.forces(0,0),0);
    EXPECT_EQ(atoms.forces(1,0),0);

    /** General Direction */
    /** direction in (1,1,1)*/
    {
        //calc
        atoms.positions.setZero();
        atoms.forces.setZero();
        Vector_t direction(1, 1, 1);
        atoms.positions.col(1) = direction * minimumDistance / sqrt(3);
        lendardJonesDirectSummation(atoms, testEpsilon, testSigma);
        //check
        EXPECT_NEAR(atoms.forces(0, 0), 0, 1e-5);
        EXPECT_NE(atoms.forces(0, 0), 0);
        EXPECT_NEAR(atoms.forces(1, 0), 0, 1e-5);
        EXPECT_NE(atoms.forces(1, 0), 0);
        EXPECT_NEAR(atoms.forces(2, 0), 0, 1e-5);
        EXPECT_NE(atoms.forces(2, 0), 0);
    }

    /** direction (1,1,0) */
    {
        //calc
        atoms.positions.setZero();
        atoms.forces.setZero();
        Vector_t direction(1, 1, 0);
        atoms.positions.col(1) = direction * minimumDistance / sqrt(2);
        lendardJonesDirectSummation(atoms, testEpsilon, testSigma);
        //check
        EXPECT_NEAR(atoms.forces(0, 0), 0, 1e-5);
        EXPECT_NE(atoms.forces(0, 0), 0);
        EXPECT_NEAR(atoms.forces(1, 0), 0, 1e-5);
        EXPECT_NE(atoms.forces(1, 0), 0);
        EXPECT_EQ(atoms.forces(2, 0), 0);
    }

    /** direction (-1,-1,-1) */
    {
        //calc
        atoms.positions.setZero();
        atoms.forces.setZero();
        Vector_t direction(-1, -1, -1);
        atoms.positions.col(1) = direction * minimumDistance / sqrt(3);
        lendardJonesDirectSummation(atoms, testEpsilon, testSigma);
        //check
        EXPECT_NEAR(atoms.forces(0, 0), 0, 1e-5);
        EXPECT_NE(atoms.forces(0, 0), 0);
        EXPECT_NEAR(atoms.forces(1, 0), 0, 1e-5);
        EXPECT_NE(atoms.forces(1, 0), 0);
        EXPECT_NEAR(atoms.forces(2, 0), 0, 1e-5);
        EXPECT_NE(atoms.forces(2, 0), 0);
    }

    /** direction (-1,-1,0) */
    {
        //calc
        atoms.positions.setZero();
        atoms.forces.setZero();
        Vector_t direction(1, 1, 0);
        atoms.positions.col(1) = direction * minimumDistance / sqrt(2);
        lendardJonesDirectSummation(atoms, testEpsilon, testSigma);
        //check
        EXPECT_NEAR(atoms.forces(0, 0), 0, 1e-5);
        EXPECT_NE(atoms.forces(0, 0), 0);
        EXPECT_NEAR(atoms.forces(1, 0), 0, 1e-5);
        EXPECT_NE(atoms.forces(1, 0), 0);
        EXPECT_EQ(atoms.forces(2, 0), 0);
    }
}

TEST(LJDirectSummationTest,ComparistionBetweenOldAndNew) {
    // should not place too high of a value in this test the one below is the important one
    int nbAtoms = 4;
    double sigma = 1;
    Positions_t  p = createLatticeCube(nbAtoms);
    Atoms atoms(p);
    EXPECT_NEAR(lendardJonesDirectSummation(atoms), lenardJonesDirectSummationWithCutoff(atoms,8*sigma),1e-3);
}

TEST(LJDirectSummationTest, ModifiedForcesTest) {
    constexpr int nb_atoms = 10;
    constexpr double epsilon = 0.7;  // choose different to 1 to pick up missing factors
    constexpr double sigma = 0.3;
    constexpr double delta = 0.00001;  // difference used for numerical (finite difference) computation of forces
    constexpr double cutoff = 2.5*sigma;

    Atoms atoms(nb_atoms);
    atoms.positions.setRandom();  // random numbers between -1 and 1

    // compute and store energy of the indisturbed configuration
    double e0{lenardJonesDirectSummationWithCutoff(atoms,cutoff, epsilon, sigma)};
    Forces_t forces0{atoms.forces};

    // loop over all atoms and compute forces from a finite differences approximation
    Forces_t dummy_forces(3, nb_atoms);  // we don't actually need these
    for (int i{0}; i < nb_atoms; ++i) {
        // loop over all Cartesian directions
        for (int j{0}; j < 3; ++j) {
            // move atom to the right
            atoms.positions(j, i) += delta;
            double eplus{lenardJonesDirectSummationWithCutoff(atoms,cutoff, epsilon, sigma)};
            // move atom to the left
            atoms.positions(j, i) -= 2 * delta;
            double eminus{lenardJonesDirectSummationWithCutoff(atoms,cutoff, epsilon, sigma)};
            // move atom back to original position
            atoms.positions(j, i) += delta;

            // finite-differences forces
            double fd_force{-(eplus - eminus) / (2 * delta)};

            // check whether finite-difference and analytic forces agree
            if (abs(forces0(j, i)) > 1e-10) {
                EXPECT_NEAR(abs(fd_force - forces0(j, i)) / forces0(j, i), 0, 1e-5);

            } else {
                EXPECT_NEAR(fd_force, forces0(j, i), 1e-10);
            }
        }
    }
}