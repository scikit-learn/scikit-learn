/**
 * \file diagonal_matrix_lib.h
 * Defines the class of diagonal matrixes, based on vectors
 */

#ifndef DIAGONALMATRIXLIB
#define DIAGONALMATRIXLIB

#include <matrix/diagonal_matrix_iterator.h>

namespace Matrix
{
  /// A diagonal matrix is dependent on a Vector type that will indicate its size and inner data
  template<class Vector>
  struct DiagonalMatrix
  {
    /// Static dimensions
    enum{staticHeight = Vector::staticHeight * Vector::staticWidth, staticWidth = Vector::staticHeight * Vector::staticWidth} ;

    /// The type of data inside the matrix
    typedef typename Vector::Data_Type Data_Type;
    /// The type of self
    typedef DiagonalMatrix Self;
    /// Typedef of operation result
    typedef Matrix<Data_Type, staticHeight, staticWidth> Result;
    /// Typedef of comparison result
    typedef Matrix<bool, staticHeight, staticWidth> ComparisonResult;
    /// The type of the transpose
    typedef Self Transpose;
    /// Typedef of iterator to an element
    typedef DiagonalIterator<Self> iterator;
    /// Typedef of const iterator to an element
    typedef DiagonalIterator<const Self> const_iterator;

    /// Typedef for <<
    typedef typename std::ostream ostream;
    /// Typedef for >>
    typedef typename std::istream istream;

    /// Typedef for data traits
    typedef DataTypeTraits<Data_Type> DataTraits;

    /**
     * Constructor
     * @param vector is the diagonal vector that is used
     */
    DiagonalMatrix(const Vector& vector)
      :vector(vector)
    {
    }

    /** @name Other operations*/
    /// @{
    /**
     * Calculates the height of the matrix
     * @return the height of the matrix
     */
    unsigned long height() const
    {
      return vector.height() * vector.width();
    }

    /**
     * Calculates the width of the matrix
     * @return the width of the matrix
     */
    unsigned long width() const
    {
      return vector.height() * vector.width();
    }

    /**
     * Copy operator
     * @param mat is the submatrix that will be copied
     * @return the copied submatrix
     */
    template<class OtherMatrixType>
    Self& operator=(const OtherMatrixType mat)
    {
      assert(height() == mat.height());
      assert(width() == mat.width());
      for (unsigned long j = 0; j < width(); ++j)
        (*this)(j, j) = mat(j, j);
      return *this;
    }

    /**
     * Assignement operator
     * @param rhs is the value that will be copied
     * @return the modified submatrix
     */
    Self& operator=(const Data_Type rhs)
    {
      for (unsigned long j = 0; j < width(); ++j)
        (*this)(j, j) = rhs;

      return (*this);
    }

    /**
     * Subscript operator
     * @param i is the line selected
     * @param j is the column selected
     * @return the element in the matrix at the emplacement [i,j]
     * \pre i==j
     */
    Data_Type& operator()(unsigned long i, unsigned long j)
    {
      assert(i < height());
      assert(j < width());
      assert(i==j);
      return vector(i);
    }

    /**
     * Const subscript operator
     * @param i is the line selected
     * @param j is the column selected
     * @return the element in the matrix at the emplacement [i,j]
     */
    const Data_Type operator()(unsigned long i, unsigned long j) const
    {
      assert(i < height());
      assert(j < width());
      return (i==j ? vector(i) : DataTraits::zero(vector(0)));
    }

    /**
     * Return an iterator to the beginning
     * @return an iterator to the beginning
     */
    iterator begin()
    {
      return DiagonalIterator<Self>(this, 0);
    }

    /**
     * Return an iterator to the beginning
     * @return an iterator to the beginning
     */
    const_iterator begin() const
    {
      return DiagonalIterator<const Self>(this, 0);
    }

    /**
     * Return an iterator to the end
     * @return an iterator to the end
     */
    const_iterator end() const
    {
      return DiagonalIterator<const Self>(this, height() * width());
    }

    /**
     * Return an iterator to the end
     * @return an iterator to the end
     */
    iterator end()
    {
      return DiagonalIterator<Self>(this, height() * width());
    }
    ///@}

  private:
    /// The real matrix
    Vector vector;
  };
}
#endif
