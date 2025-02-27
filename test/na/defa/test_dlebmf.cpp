#include "na/defa/DLEBMF.hpp"

#include <gtest/gtest.h>

TEST(DLEBMF, GetSparse) {
  const auto matrix = na::DLEBMF::fromSparseMatrix(2, 2, {{}, {1}});
  EXPECT_FALSE(matrix.get(0, 0));
  EXPECT_FALSE(matrix.get(0, 1));
  EXPECT_FALSE(matrix.get(1, 0));
  EXPECT_TRUE(matrix.get(1, 1));
}

TEST(DLEBMF, GetDense) {
  const auto matrix = na::DLEBMF::fromDenseMatrix({{false, false}, {false, true}});
  EXPECT_FALSE(matrix.get(0, 0));
  EXPECT_FALSE(matrix.get(0, 1));
  EXPECT_FALSE(matrix.get(1, 0));
  EXPECT_TRUE(matrix.get(1, 1));
}

TEST(DLEBMF, StringSparse) {
  const auto matrix = na::DLEBMF::fromSparseMatrix(3, 4, {{0, 2}, {}, {2}});
  EXPECT_EQ(matrix.toString(), "1 0 1 0\n0 0 0 0\n0 0 1 0");
}

TEST(DLEBMF, StringDense) {
  const auto matrix = na::DLEBMF::fromDenseMatrix({{true, false, true, false}, {false, false, false, false}, {false, false, true, false}});
  EXPECT_EQ(matrix.toString(), "1 0 1 0\n0 0 0 0\n0 0 1 0");
}