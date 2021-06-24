//
// Created by cm on 24.06.21.
//

#include <gtest/gtest.h>
//
#include "../Headerfiles/atoms.h"
#include "../Headerfiles/helperfunctions.h"

TEST(TestHelperfunctions, cubicEncapsulation)
{
    Positions_t p = createLatticeCube(8,10);
    Atoms atoms(p);
    Vector_t min{atoms.positions.row(0).minCoeff(),atoms.positions.row(1).minCoeff(),atoms.positions.row(2).minCoeff()};
    Vector_t max{atoms.positions.row(0).maxCoeff(),atoms.positions.row(1).maxCoeff(),atoms.positions.row(2).maxCoeff()};
    Positions_t controlPostitions = generateCapsel(atoms, 2);

    //check weather all atoms are inside
    for(int i = 0; i < atoms.nb_atoms(); ++i) {
        //should be easier too just compare the values of the Vectors
        //the positions should be smaller than the min
        EXPECT_EQ(compareVectorsBigSmall(atoms.positions.col(i),controlPostitions.col(0)),true);
        //the positions should be smaller than the max
        EXPECT_EQ(compareVectorsBigSmall(atoms.positions.col(i),controlPostitions.col(0)),true);
    }
}

TEST(TestHelperfunctions, encapsulatingFunction)
{
    //ASSERT_TRUE(false) << "Implement me!!";
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