//
// Created by pastewka on 21.04.21.
//

#include <gtest/gtest.h>

// Demonstrate some basic assertions.
TEST(HelloTest, BasicAssertions) {
  // Expect two strings not to be equal.
  EXPECT_STRNE("hello", "world");
  //ASSERT_TRUE(false);
  // Expect equality.
  EXPECT_EQ(7 * 6, 42);
}

/**
TEST(Test Suit, Test name)
{
    ASSERT_TRUE(false) << "Implement me!!";
}
*/