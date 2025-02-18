#include "na/zonar/DLEBMF.hpp"

#include <algorithm>
#include <cstddef>
#include <initializer_list>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace na {
DLEBMF::DLEBMF(const size_t rows, const size_t cols) : rows(rows), cols(cols) {
  createEmptyColumns();
}

auto DLEBMF::fromDenseMatrix(
    const std::initializer_list<std::initializer_list<bool>>& entries)
    -> DLEBMF {
  if (checkRectangularShape(entries)) {
    throw std::invalid_argument("All columns must have the same length.");
  }
  DLEBMF instance(entries.size(), entries.begin() != entries.end()
                                      ? entries.begin()->size()
                                      : 0);

  size_t r = 0;
  for (const auto& row : entries) {
    size_t c = 0;
    Column* currentCol = instance.matrix.get();
    Cell* lastInRow = nullptr;
    for (const auto val : row) {
      if (val) {
        auto currentCell = std::make_unique<Cell>();
        currentCell->row = r;
        currentCell->col = c;
        currentCell->left = lastInRow;
        currentCell->up = currentCol->bottom;
        if (lastInRow != nullptr) {
          lastInRow->right = currentCell.get();
        }
        lastInRow = currentCell.get();
        currentCol->size++;
        if (currentCol->down == nullptr) {
          currentCol->down = std::move(currentCell);
          currentCol->bottom = currentCol->down.get();
        } else {
          currentCol->bottom->down = std::move(currentCell);
          currentCol->bottom = currentCol->bottom->down.get();
        }
      }
      currentCol = currentCol->right.get();
      ++c;
    }
    ++r;
  }
  return instance;
}
auto DLEBMF::fromSparseMatrix(
    const std::size_t rows, const std::size_t cols,
    const std::initializer_list<std::initializer_list<std::size_t>>& entries)
    -> DLEBMF {
  if (entries.size() != rows) {
    throw std::invalid_argument(
        "Number of rows does not match the number of rows in the entries.");
  }
  // Get the maximum column index
  const auto maxColId = std::accumulate(
      entries.begin(), entries.end(), 0,
      [](const std::size_t acc, const std::initializer_list<std::size_t>& col) {
        return col.size() == 0
                   ? acc
                   : std::max(acc, *std::max_element(col.begin(), col.end()));
      });
  if (maxColId >= cols) {
    throw std::invalid_argument("The maximum column index in the entries "
                                "exceeds the number of columns.");
  }
  // Check for duplicate indices in the same row
  if (!checkUniqueIndices(entries)) {
    throw std::invalid_argument(
        "Duplicate indices in the same row are not allowed.");
  }

  DLEBMF instance(rows, cols);

  size_t r = 0;
  for (const auto& row : entries) {
    Column* currentCol = instance.matrix.get();
    Cell* lastInRow = nullptr;
    std::vector colIdxs(row.begin(), row.end());
    std::sort(colIdxs.begin(), colIdxs.end());
    for (const auto c : colIdxs) {
      while (currentCol->col < c) {
        currentCol = currentCol->right.get();
      }
      auto currentCell = std::make_unique<Cell>();
      currentCell->row = r;
      currentCell->col = c;
      currentCell->left = lastInRow;
      currentCell->up = currentCol->bottom;
      if (lastInRow != nullptr) {
        lastInRow->right = currentCell.get();
      }
      lastInRow = currentCell.get();
      currentCol->size++;
      if (currentCol->down == nullptr) {
        currentCol->down = std::move(currentCell);
        currentCol->bottom = currentCol->down.get();
      } else {
        currentCol->bottom->down = std::move(currentCell);
        currentCol->bottom = currentCol->bottom->down.get();
      }
    }
    ++r;
  }
  return instance;
}

auto DLEBMF::createEmptyColumns() -> void {
  if (cols == 0) {
    matrix = nullptr;
  } else {
    matrix = std::make_unique<Column>();
    matrix->col = 0;
    Column* lastCol = matrix.get();
    for (size_t i = 1; i < cols; ++i) {
      auto currentCol = std::make_unique<Column>();
      currentCol->col = i;
      currentCol->left = lastCol;
      lastCol->right = std::move(currentCol);
      lastCol = lastCol->right.get();
    }
  }
}

auto DLEBMF::checkRectangularShape(
    const std::initializer_list<std::initializer_list<bool>>& entries) -> bool {
  if (entries.size() == 0) {
    return true;
  }
  const auto* firstCol = entries.begin();
  const size_t nCols = firstCol->size();
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  return std::any_of(firstCol + 1, entries.end(),
                     [nCols](const auto& col) { return col.size() != nCols; });
}

auto DLEBMF::checkUniqueIndices(
    const std::initializer_list<std::initializer_list<std::size_t>>& entries)
    -> bool {
  for (const auto& row : entries) {
    std::vector colIdxs(row.begin(), row.end());
    std::sort(colIdxs.begin(), colIdxs.end());
    if (std::adjacent_find(colIdxs.begin(), colIdxs.end()) != colIdxs.end()) {
      return false;
    }
  }
  return true;
}

auto DLEBMF::get(const size_t row, const size_t col) const -> bool {
  if (row >= rows || col >= cols) {
    throw std::out_of_range("Row or column index out of range.");
  }
  const auto* currentCol = matrix.get();
  for (size_t i = 0; i < col; ++i) {
    currentCol = currentCol->right.get();
  }
  if (currentCol->size == 0) {
    return false;
  }
  const auto* currentCell = currentCol->down.get();
  while (currentCell->row < row && currentCell->down != nullptr) {
    currentCell = currentCell->down.get();
  }
  return currentCell->row == row;
}

auto DLEBMF::factorize() -> const std::vector<Factor>& {
  // 1. remove all empty columns
  // 2. remove all duplicate columns
  // 3. select a column; choose the one with the fewest 1s
  // 4. select a subset of rows with 1s in the selected column; choose all
  // subsets including the first row with a one and all subsequent rows with a 1
  // starting with the full subset; backtrack
  // 5. remove the 1s in the selected rows in all columns where there is a 1
  // in all of the selected rows
  // 6. backtrack and repeat from step 1 until there are no more columns
  // 7. restore all duplicate columns
  // 8. restore all empty columns
}

auto DLEBMF::Factor::toString() const -> std::string {
  std::stringstream ss;
  ss << "rows: [";
  for (size_t i = 0; i < rows.size(); ++i) {
    ss << rows[i];
    if (i < rows.size() - 1) { // not for the last row
      ss << ", ";              // separate indices by a comma
    }
  }
  ss << "]\ncols: [";
  for (size_t i = 0; i < cols.size(); ++i) {
    ss << cols[i];
    if (i < cols.size() - 1) { // not for the last column
      ss << ", ";              // separate indices by a comma
    }
  }
  ss << "]";
  return ss.str();
}

auto DLEBMF::toString() const -> std::string {
  std::stringstream ss;
  for (size_t r = 0; r < rows; ++r) {
    for (size_t c = 0; c < cols; ++c) {
      ss << (get(r, c) ? "1" : "0");
      if (c < cols - 1) { // not in the last column
        ss << ' ';
      }
    }
    if (r < rows - 1) { // not in the last row
      ss << '\n';
    }
  }
  return ss.str();
}

} // namespace na
