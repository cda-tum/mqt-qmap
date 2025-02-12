#pragma once
#include <cstddef>
#include <initializer_list>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

namespace na {

/**
 * @brief A class to solve the exact binary matrix factorization problem with
 * dancing links.
 * @details This class is both, the representation of the boolean-valued matrix
 * and allows to compute the exact binary matrix factorization using the
 * dancing links technique.
 *
 * @par
 * The representation of the matrix only stores the @c true entries. The entries
 * are stored as a doubly linked list, i.e., each entray has a pointer to every
 * surrounding @c true entry. Additionally, there is one control row on top of
 * the matrix that holds the first cell in every column. Additionally, the
 * control head of each column stores the number of @c true entries in this
 * column.
 */
class DLEBMF {
private:
  size_t rows = 0; ///< Number of rows in the matrix.
  size_t cols = 0; ///< Number of columns in the matrix.

  /**
   * @brief A struct to represent a @c true entry in the matrix.
   * @details The cells are created as unique pointers. The owner of the
   * respective cell is always the cell above or the control head of the column
   * if the cell is the topmost in this column.
   */
  struct Cell {
    size_t row = 0; ///< Row index of the cell in the original matrix.
    size_t col = 0; ///< Column index of the cell in the original matrix.
    /// Pointer to the next cell in the row to the right.
    Cell* right = nullptr;
    /**
     * @brief Pointer to the next cell in the column below.
     * @details This is the owner of the cell below if it exists, otherwise the
     * pointer is @c nullptr.
     */
    std::unique_ptr<Cell> down = nullptr;
    /// Pointer to the previous cell in the row to the left.
    Cell* left = nullptr;
    Cell* up = nullptr; ///< Pointer to the previous cell in the column above.
  };

  /**
   * @brief A struct to represent a column in the matrix.
   * @details This is an additional node on top of each column and does not
   * represent any entry. Additionally, the size, i.e., the number of @c true
   * entries in that column is stored here.
   */
  struct Column {
    size_t col = 0;  ///< Column index of the column in the original matrix.
    size_t size = 0; ///< Number of @c true entries in this column.
    std::unique_ptr<Column> right = nullptr; ///< Pointer to the next column.
    /// Pointer to the first, topmost cell in the column.
    std::unique_ptr<Cell> down = nullptr;
    /// Pointer to the previous column to the left.
    Column* left = nullptr;
    /// Pointer to the last, bottommost cell in the column.
    Cell* bottom = nullptr;
  };

  /// Pointer to the first column in the matrix that holds the rest of the
  /// matrix.
  std::unique_ptr<Column> matrix;

  /**
   * @brief Create an empty matrix with empty columns.
   * @details This method creates the control heads for @ref cols many columns.
   * All columns will be empty, i.e., the resulting matrix does not contain
   * any @c true entries.
   */
  auto createEmptyColumns() -> void;
  /**
   * @brief Check if the given list of lists is representing a matrix, i.e., has
   * a rectangular shape.
   * @param entries is the list of lists with boolean values to check.
   * @return @c true if the list of lists is rectangular, @c false otherwise.
   */
  static auto checkRectangularShape(
      const std::initializer_list<std::initializer_list<bool>>& entries)
      -> bool;

public:
  /**
   * @brief Create an empty matrix with no rows and columns.
   */
  DLEBMF() = default;
  /**
   * @brief Create a new matrix with the given number of rows and columns. All
   * entries in this matrix are initialized with @c false.
   */
  DLEBMF(size_t rows, size_t cols);
  /**
   * @brief Create a new matrix with the given entries.
   * @param entries is a list of lists with boolean values. The list of lists
   * must represent a matrix, i.e., all inner lists must be of the same size.
   * @throws std::invalid_argument if the given list of lists is not
   * rectangular.
   * @note The arguments of this constructor can easily be written as a
   * brace-init list. For example, the matrix
   * \f[
   * \begin{pmatrix}
   * 1 & 0 & 1 & 0 \\
   * 0 & 0 & 0 & 0 \\
   * 0 & 0 & 1 & 0
   * \end{pmatrix}
   * \f]
   * can be created with
   * @code
   * na::DLEBMF matrix({{true, false, true, false}, {false, false, false,
   * false}, {false, false, true, false}});
   * @endcode
   */
  DLEBMF(const std::initializer_list<std::initializer_list<bool>>& entries);

  /**
   * @brief Get the value of the entry at the given row and column, this can
   * either be @c true or @c false.
   * @param row is the row index of the entry.
   * @param col is the column index of the entry.
   * @return @c true if the entry is @c true, @c false otherwise.
   */
  [[nodiscard]] auto get(size_t row, size_t col) const -> bool;

  /**
   * @brief A class to represent a factor of the exact binary matrix
   * factorization.
   * @details A factor is a submatrix of the original matrix that has a
   * rectangular shape. This means that if the entries (i, j) and (k, l) are @c
   * true in the factor, then all entries (i, l), (k, j), (i, j), and (k, l) are
   * @c true in the factor as well. This kind of submatrix can be efficiently
   * stored as a list of row and column indices of the @c true entries.
   */
  class Factor {
  private:
    std::vector<size_t> rows; ///< Row indices of the @c true entries.
    std::vector<size_t> cols; ///< Column indices of the @c true entries.
    Factor() = default;

  public:
    /**
     * @brief Get the row indices of the @c true entries in this factor.
     * @return A vector with the row indices of the @c true entries in this
     * factor.
     */
    [[nodiscard]] auto getRows() const -> const std::vector<size_t>& {
      return rows;
    }
    /**
     * @brief Get the column indices of the @c true entries in this factor.
     * @return A vector with the column indices of the @c true entries in this
     * factor.
     */
    [[nodiscard]] auto getCols() const -> const std::vector<size_t>& {
      return cols;
    }
    /**
     * @brief Returns a string representation of the factor.
     * @details The string contains the vectors of row and column indices
     * in a human-readable format.
     * @return A string representation of the factor.
     */
    [[nodiscard]] auto toString() const -> std::string;
    /**
     * @brief Output operator to print the factor to an output stream.
     * @param os is the output stream to write to.
     * @param factor is the factor to write to the output stream.
     * @return The output stream after writing the factor to it.
     */
    friend auto operator<<(std::ostream& os, const Factor& factor)
        -> std::ostream& {
      os << factor.toString();
      return os;
    }
  };

  [[nodiscard]] auto factorize() -> const std::vector<Factor>&;

  /**
   * @brief Returns a string representation of the matrix.
   * @details The matrix @code na::DLEBMF({{true, false, true, false}, {false,
   * false, false, false}, {false, false, true, false}}) @endcode will be
   * represented as: @verbatim 1 0 1 0
   * 0 0 0 0
   * 0 0 1 0@endverbatim
   * @return A string representation of the matrix where @c true entries are
   * represented by @c 1 and @c false entries are represented by @c 0. Rows are
   * separated by newline characters and columns are separated by spaces.
   */
  [[nodiscard]] auto toString() const -> std::string;
  /**
   * @brief Output operator to print the matrix to an output stream.
   * @param os is the output stream to write to.
   * @param obj is the matrix to write to the output stream.
   * @return The output stream after writing the matrix to it.
   */
  friend auto operator<<(std::ostream& os, const DLEBMF& obj) -> std::ostream& {
    os << obj.toString();
    return os;
  }
};

} // namespace na