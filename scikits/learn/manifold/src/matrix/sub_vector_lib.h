/**
 * \file sub_vector_lib.h
 * This files contains the description of in-matrix vectors, column, line, diagonal or subdiagonal
 */

#ifndef IN_VECTOR_LIB
#define IN_VECTOR_LIB

#include <matrix/matrix_lib.h>
#include <matrix/sub_matrix_lib.h>

namespace Matrix
{
  /// The sub diagonal vector class
  template<class Matrix_Type>
  struct SubDiagonalVector
  {
    /// Typedef for <<
    typedef typename std::ostream ostream;

    typedef Matrix_Type MatrixType;
    typedef typename MatrixType::Data_Type DataType;
    /** @name Constructors and destructors (subvector)*/
    /// @{
    /**
     * Constructor from a matrix
     * @param mat is the matrix whose diagonal components will be accessible with this vector
     */
    SubDiagonalVector(MatrixType& mat)
      :matrix(mat)
    {
    }
    ///@}

    /** @name Other operations (subvector)*/
    /// @{
    /**
     * Returns the length of a vector
     * @return the length of the vector
     */
    unsigned long length() const
    {
      return std::min(matrix.height(), matrix.width());
    }

    unsigned long width() const
    {
      return 1U;
    }

    unsigned long height() const
    {
      return length();
    }

    /**
     * Returns the size of a vector
     * @return the size of the vector
     */
    unsigned long size() const
    {
      return std::min(matrix.height(), matrix.width());
    }

    /**
     * Subscript operator
     * @param i is the line selected
     * @return a reference to an element in the diagonal
     */
    DataType& operator()(unsigned long i)
    {
      assert(i >= 0);
      assert(i < matrix.height());
      assert(i < matrix.width());
      return matrix(i, i);
    }

    /**
     * Subscript operator
     * @param i is the line selected
     * @return a reference to an element in the diagonal
     */
    DataType operator()(unsigned long i) const
    {
      assert(i >= 0);
      assert(i < matrix.height());
      assert(i < matrix.width());
      return matrix(i, i);
    }

    /**
     * Subscript operator
     * @param i is the line selected
     * @return a reference to an element in the diagonal
     */
    DataType operator()(unsigned long i, unsigned long j) const
    {
      assert(i >= 0);
      assert(i < matrix.height());
      assert(i < matrix.width());
      assert(j == 0);
      return matrix(i, i);
    }
    ///@}

    /**
     * The matrix inside the vector class
     */
    MatrixType& matrix;
  };

  /// The sub diagonal vector class
  template<class Matrix_Type>
  struct SubDownDiagonalVector
  {
    typedef Matrix_Type MatrixType;
    typedef typename MatrixType::Data_Type DataType;
    /** @name Constructors and destructors (subvector)*/
    /// @{
    /**
     * Constructor from a matrix
     * @param mat is the matrix whose sub diagonal components will be accessible with this vector
     */
    SubDownDiagonalVector(MatrixType& mat)
      :matrix(mat)
    {
    }
    ///@}

    /** @name Other operations (subvector)*/
    /// @{
    /**
     * Returns the length of a vector
     * @return the length of the vector
     */
    unsigned long length() const
    {
      return std::min(matrix.height() - 1, matrix.width());
    }

    /**
     * Returns the size of a vector
     * @return the size of the vector
     */
    unsigned long size() const
    {
      return std::min(matrix.height() - 1, matrix.width());
    }

    /**
     * Subscript operator
     * @param i is the line selected
     * @return a reference to an element in the subdiagonal
     */
    DataType& operator()(unsigned long i)
    {
      assert(i >= 0);
      assert(i < matrix.height() - 1);
      assert(i < matrix.width());
      return matrix(i + 1, i);
    }

    /**
     * Subscript operator
     * @param i is the line selected
     * @return a reference to an element in the subdiagonal
     */
    DataType operator()(unsigned long i) const
    {
      assert(i >= 0);
      assert(i < matrix.height() - 1);
      assert(i < matrix.width());
      return matrix(i + 1, i);
    }
    ///@}

    /**
     * The matrix inside the vector class
     */
    MatrixType& matrix;
  };

  /// The upper diagonal vector class
  template<class Matrix_Type>
  struct SubUpDiagonalVector
  {
    typedef Matrix_Type MatrixType;
    typedef typename MatrixType::Data_Type DataType;
    /** @name Constructors and destructors (subvector)*/
    /// @{
    /**
     * Constructor from a matrix
     * @param mat is the matrix whose sub diagonal components will be accessible with this vector
     */
    SubUpDiagonalVector(MatrixType& mat)
      :matrix(mat)
    {
    }
    ///@}

    /** @name Other operations (subvector)*/
    /// @{
    /**
     * Returns the length of a vector
     * @return the length of the vector
     */
    unsigned long length() const
    {
      return std::min(matrix.height(), matrix.width() - 1);
    }

    /**
     * Returns the size of a vector
     * @return the size of the vector
     */
    unsigned long size() const
    {
      return std::min(matrix.height(), matrix.width() - 1);
    }

    /**
     * Subscript operator
     * @param i is the line selected
     * @return a reference to an element in the subdiagonal
     */
    DataType& operator()(unsigned long i)
    {
      assert(i >= 0);
      assert(i < matrix.height());
      assert(i < matrix.width() - 1);
      return matrix(i, i + 1);
    }

    /**
     * Subscript operator
     * @param i is the line selected
     * @return a reference to an element in the subdiagonal
     */
    DataType operator()(unsigned long i) const
    {
      assert(i >= 0);
      assert(i < matrix.height());
      assert(i < matrix.width() - 1);
      return matrix(i, i + 1);
    }
    ///@}

    /**
     * The matrix inside the vector class
     */
    MatrixType& matrix;
  };
}

#endif

