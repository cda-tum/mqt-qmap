#include "na/zonar/DLEBMF.hpp"

#include <gtest/gtest.h>

TEST(DLEBMF, Get) {
  const na::DLEBMF matrix({{false, false}, {false, true}});
  EXPECT_FALSE(matrix.get(0, 0));
  EXPECT_FALSE(matrix.get(0, 1));
  EXPECT_FALSE(matrix.get(1, 0));
  EXPECT_TRUE(matrix.get(1, 1));
}

TEST(DLEBMF, String) {
  const na::DLEBMF matrix({{true, false, true, false}, {false, false, false, false}, {false, false, true, false}});
  EXPECT_EQ(matrix.toString(), "1 0 1 0\n0 0 0 0\n0 0 1 0");
}