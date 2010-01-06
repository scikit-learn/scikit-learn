/**
 * \file matrix_lib.h
 * The main header file which should be included
 * \todo - check the speed of the multiplication implementation -> row vs columns !
 * \todo - pseudoinverse
 * \todo - other multiplication algorithm
 */

#ifndef MATRIX_LIB
#define MATRIX_LIB

#include <iosfwd>
#include <cassert>
#include <string>
#include <algorithm>
#include <functional>
#include <numeric>

#include <boost/mpl/not.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>

#include <matrix/matrix_traits.h>
#include <matrix/static_container_matrix.h>
#include <matrix/dynamic_container_matrix.h>
#include <matrix/matrix_functions.h>

namespace Matrix
{
  /// The general matrix class
  template<class DataType, unsigned long Height, unsigned long Width>
  struct Matrix : public MatrixContainer<DataType, Height, Width>
  {
  public:
    /// Only to have a link on the type used in the matrix
    typedef DataType Data_Type;
    /// Typedef of parent
    typedef MatrixContainer<DataType, Height, Width> Parent;
    /// Typedef of self
    typedef Matrix Self;
    /// Typedef of operation result
    typedef Matrix<DataType, Height, Width> Result;
    /// Typedef of comparison result
    typedef Matrix<bool, Height, Width> ComparisonResult;
    /// Typedef of transposition
    typedef Matrix<DataType, Width, Height> Transpose;
    /// Typedef of pointer to an element
    typedef typename Parent::pointer pointer;
    /// Typedef of iterator to an element
    typedef typename Parent::iterator iterator;
    /// Typedef of const iterator to an element
    typedef typename Parent::const_iterator const_iterator;
    /// Typedef of column selection
    typedef typename ColumnVector<DataType, Height>::Type Column;
    /// Typedef of line selection
    typedef typename LineVector<DataType, Width>::Type Line;
    /// Typedef for <<
    typedef typename std::ostream ostream;
    /// Typedef for >>
    typedef typename std::istream istream;

    /// Typedef for data traits
    typedef DataTypeTraits<DataType> DataTraits;
    /// Typedef for own traits
    typedef DataTypeTraits<Self> OwnTraits;

    /// Static dimensions
    enum{staticHeight = Height, staticWidth = Width} ;

    using Parent::mMatrix;
    using Parent::height;
    using Parent::width;
    using Parent::begin;
    using Parent::end;
    using Parent::operator();

    /** @name Constructors and destructors */
    /// @{
    /**
     * The constructor of the matrix
     * @param height is the number of lines
     * @param width is the number of columns
     * @param init_val is the initial value that has to be put in the matrix
     */
    Matrix(unsigned long height, unsigned long width, DataType init_val = DataTraits::zero(0))
      :Parent(height, width, init_val)
    {
    }

    /**
     * The constructor of the matrix with functor
     * @param height is the number of lines
     * @param width is the number of columns
     * @param functor is the functor class used to fill the matrix
     */
    template<class Functor>
    Matrix(Functor& functor, unsigned long height, unsigned long width)
      :Parent(functor, height, width)
    {
    }

    /**
     * Copy constructor
     * @param other is the copied matrix
     */
    template<class MatrixType>
    explicit Matrix(const MatrixType& other)
      :Parent(other, true)
    {
      BOOST_STATIC_ASSERT(Height == 0U || MatrixType::staticHeight == 0 || Height == MatrixType::staticHeight);
      BOOST_STATIC_ASSERT(Width == 0U || MatrixType::staticWidth == 0 || Width == MatrixType::staticWidth);
    }

    /**
     * Constructor from a table of elements
     * @param height is the height of the new matrix
     * @param width is the width of the new matrix
     * @param elements is a pointer to the table of DataType elements
     */
    Matrix(unsigned long height, unsigned long width, const DataType* elements)
      :Parent(height, width, elements)
    {
    }

    /**
     * Bogus constructor
     */
    explicit Matrix(void)
      :Parent()
    {
    }
    ///@}

    /** @name Other operations */
    /// @{
    /**
     * Gets a column in the matrix
     * @param j is the column selected
     * @return the column in the matrix at the emplacement [:,j]
     */
    const Column operator[](unsigned long j) const
    {
      assert(j < width());
      Column temp_vector(height(), 1U, DataTraits::zero((*this)(0, 0)));

      for(unsigned long i = 0; i < height(); ++i)
        temp_vector(i, 0) = (*this)(i, j);

      return temp_vector;
    }

    /**
     * Transforms a matrix to a bool, useful for implict conversion
     */
    operator const bool() const
    {
      BOOST_STATIC_ASSERT((boost::is_same<bool, DataType>::value));
      bool result = true;
      for(unsigned long j = 0; j < width(); ++j)
      {
        for(unsigned long i = 0; i < height(); ++i)
        {
          result &= static_cast<bool>((*this)(i, j));
        }
      }
      return result;
    }

    /**
     * Assignement operator : assigning a matrix
     * @param rhs is the copied matrix
     * @return the new matrix
     */
    template<class Assignable>
    Self& operator=(const Assignable rhs)
    {
      Parent::operator=(rhs);
      return *this;
    }

    /**
     * Subscript operator
     * @param i is the case selected
     * @return the element in the matrix at the emplacement [i % height() , i / height()]
     */
    DataType& operator()(unsigned long i)
    {
      return Parent::mMatrix[i];
    }

    /**
     * Subscript operator
     * @param i is the case selected
     * @return the element in the matrix at the emplacement [i % height() , i / height()]
     */
    const DataType operator()(unsigned long i) const
    {
      return Parent::mMatrix[i];
    }
    ///@}

    /** @name Shift */
    /// @{
    /**
     * <<= operator : left shift by a constant
     * @param rhs is the shift asked
     * @return the modified matrix
     */
    Self& operator<<=(const DataType rhs)
    {
      for(unsigned long j = 0; j < width(); ++j)
      {
        for(unsigned long i = 0; i < height(); ++i)
        {
          (*this)(i,j) <<= rhs;
        }
      }
      return (*this);
    }

    /**
     * <<= operator : left shift by a matrix
     * @param rhs is the matrix shift asked
     * @return the modified matrix
     */
    template<class OtherMatrixType>
    typename boost::enable_if<boost::mpl::not_<boost::is_arithmetic<OtherMatrixType> >, Self&>::type operator<<=(const OtherMatrixType rhs)
    {
      BOOST_STATIC_ASSERT(Height == 0U || OtherMatrixType::staticHeight == 0 || Height == OtherMatrixType::staticHeight);
      BOOST_STATIC_ASSERT(Width == 0U || OtherMatrixType::staticWidth == 0 || Width == OtherMatrixType::staticWidth);
      assert(height() == rhs.height());
      assert(width() == rhs.width());
      for(unsigned long j = 0; j < width(); ++j)
      {
        for(unsigned long i = 0; i < height(); ++i)
        {
          (*this)(i,j) <<= rhs(i, j);
        }
      }
      return (*this);
    }

    /**
     * >>= operator
     * @param rhs is the shift asked
     * @return the modified matrix
     */
    Self& operator>>=(const DataType rhs)
    {
      for(unsigned long j = 0; j < width(); ++j)
      {
        for(unsigned long i = 0; i < height(); ++i)
        {
          (*this)(i,j) >>= rhs;
        }
      }
      return (*this);
    }

    /**
     * >>= operator : right shift by a matrix
     * @param rhs is the matrix shift asked
     * @return the modified matrix
     */
    template<class OtherMatrixType>
    typename boost::enable_if<boost::mpl::not_<boost::is_arithmetic<OtherMatrixType> >, Self&>::type operator>>=(const OtherMatrixType rhs)
    {
      BOOST_STATIC_ASSERT(Height == 0U || OtherMatrixType::staticHeight == 0 || Height == OtherMatrixType::staticHeight);
      BOOST_STATIC_ASSERT(Width == 0U || OtherMatrixType::staticWidth == 0 || Width == OtherMatrixType::staticWidth);
      assert(height() == rhs.height());
      assert(width() == rhs.width());
      for(unsigned long j = 0; j < width(); ++j)
      {
        for(unsigned long i = 0; i < height(); ++i)
        {
          (*this)(i,j) >>= rhs(i, j);
        }
      }
      return (*this);
    }
    ///@}

    /** @name Usual operations */
    /// @{
    /**
     * Tests if the matrix is square or not
     * @return true if the matrix is square
     */
    bool isSquare()
    {
      return (height() == width());
    }

    /**
     * Finds the absolute max in the matrix and returns it
     * @param low is the upper coordinate
     * @param left is the left coordinate
     * @param high is the lower coordinate
     * @param right is the right coordinate
     * @param absisse is the returned absisse of the max
     * @param ordinate is the returned ordinate of the max
     * @return the max in the submatrix
     */
    DataType max(unsigned long low = 0, unsigned long left = 0, unsigned long high = std::numeric_limits<unsigned long>::max(), unsigned long right = std::numeric_limits<unsigned long>::max(), unsigned long* absisse = NULL, unsigned long* ordinate = NULL)
    {
      if(high == std::numeric_limits<unsigned long>::max())
        high = height();
      if(right == std::numeric_limits<unsigned long>::max())
        right = width();
      unsigned long absi = left;
      unsigned long ord = low;
      DataType tmp_max = (*this)(ord, absi);
      for(unsigned long j = left; j < right; ++j)
      {
        for(unsigned long i = low; i < high; ++i)
        {
          if(tmp_max < DataTraits::absolute((*this)(i, j)))
          {
            tmp_max = DataTraits::absolute((*this)(i, j));
            absi = j;
            ord = i;
          }
        }
      }
      if(absisse != NULL)
      {
       *absisse = absi;
      }
      if(ordinate != NULL)
      {
       *ordinate = ord;
      }
      return (*this)(ord, absi);
    }

    /**
     * Calculates the kronecker product of a matrix
     * also called Direct Product, or Tensor Product
     * @return the tensor product of the matrix
     */
    template<unsigned long otherHeight, unsigned long otherWidth>
    const Matrix<DataType, Height * otherHeight, Width * otherWidth> kronecker(const Matrix<DataType, otherHeight, otherWidth> rhs)
    {
      Matrix<DataType, Height * otherHeight, Width * otherWidth> res(width() * rhs.width(), height() * rhs.height());

      unsigned long r_width = rhs.width();
      unsigned long r_height = rhs.height();
      SubMatrix<Matrix<DataType, Height * otherHeight, Width * otherWidth> > sub(res, 0, 0, r_width, r_height);
      for(unsigned long j = 0; j < width(); ++j)
      {
        sub.setLineOffset(j * r_width);
        for(unsigned long i = 0; i < height(); ++i)
        {
          sub.setColumnOffset(i * r_height);
          sub = (*this)(i,j) * rhs;
        }
      }
      return res;
    }

    /**
     * Aliases for kronecker product
     */

    template<unsigned long otherHeight, unsigned long otherWidth>
    const Matrix<DataType, Height * otherHeight, Width * otherWidth> kroneckerProduct(const Matrix<DataType, otherHeight, otherWidth> rhs){ return kronecker(rhs); }
    template<unsigned long otherHeight, unsigned long otherWidth>
    const Matrix<DataType, Height * otherHeight, Width * otherWidth> tensorProduct(const Matrix<DataType, otherHeight, otherWidth> rhs){ return kronecker(rhs); }
    template<unsigned long otherHeight, unsigned long otherWidth>
    const Matrix<DataType, Height * otherHeight, Width * otherWidth> directProduct(const Matrix<DataType, otherHeight, otherWidth> rhs){ return kronecker(rhs); }

    ///@}

    /** @name Factory */
    /// @{
    /**
     * Creates the left matrix needed for the inversion of a square matrix
     * @param matrix is the matrix from which the left matrix is created
     * @param height is the height of the matrix being inverted
     * @param width is the width of the matrix being inverted
     * @param abs_pivot is the absisse of the pivot used in the matrix
     * @param ord_pivot is the ordinate of the pivot used in the matrix
     * @param inverse is the inverse of the pivot
     */
    static Self leftPivotForSquareMatrix(const Self& matrix, const unsigned long height, const unsigned long width, const unsigned long abs_pivot, const unsigned long ord_pivot, DataType inverse)
    {
      assert(abs_pivot < matrix.width());
      assert(ord_pivot < matrix.height());

      unsigned long m_width = matrix.width();
      unsigned long m_height = matrix.height();

      Self pivot(m_height - 1, m_width, DataTraits::zero(matrix(0, 0)));

      for(unsigned long i = 0, j = 0; i < m_height - 1; ++i, ++j)
      {
        // if we are at the same line as the ordonate
        if(i == ord_pivot)
          ++j;
        for(unsigned long k = 0, l = 0; k < m_width; ++k, ++l)
        {
          //std::cout << i << j << k << l << "\n";
          if(k != ord_pivot)
          {
            if(l != i)
              pivot(i, k) = DataTraits::zero(matrix(0, 0));
            else
              pivot(i, k) = DataTraits::one(matrix(0, 0));
          }
          else
          {
            pivot(i, k) = - matrix(j, abs_pivot) * inverse;
            --l;
          }
        }
      }
      return pivot;
    }

    /**
     * Creates the right matrix needed for the inversion of a square matrix
     * @param matrix is the matrix from which the right matrix is created
     * @param height is the height of the matrix being inverted
     * @param width is the width of the matrix being inverted
     * @param abs_pivot is the absisse of the pivot used in the matrix
     * @param ord_pivot is the ordinate of the pivot used in the matrix
     * @param inverse is the inverse of the pivot
     */
    static Self rightPivotForSquareMatrix(const Self& matrix, const unsigned long height, const unsigned long width, const unsigned long abs_pivot, const unsigned long ord_pivot, DataType inverse)
    {
      assert(abs_pivot < matrix.width());
      assert(ord_pivot < matrix.height());

      unsigned long m_width = matrix.width();
      unsigned long m_height = matrix.height();

      Self pivot(m_height, m_width - 1, DataTraits::zero(matrix(0, 0)));

      for(unsigned long k = 0, l = 0; k < m_width - 1; ++k, ++l)
      {
        // if we are at the same column as the absisse
        if(k == abs_pivot)
          ++l;
        for(unsigned long i = 0, j = 0; i < m_height; ++i, ++j)
        {
          //std::cout << i << j << k << l << "\n";
          if(i != abs_pivot)
          {
            if(k != j)
              pivot(i, k) = DataTraits::zero(matrix(0, 0));
            else
              pivot(i, k) = DataTraits::one(matrix(0, 0));
          }
          else
          {
            pivot(i, k) = (- inverse) * matrix(ord_pivot, l);
            --j;
          }
        }
      }
      return pivot;
    }
 
	/**
     * Creates the identity matrix
     * @param size is the size of the identity
     * @param exemplar is a sample of the new matrix, useful when dealing with matrix of matrixes
     */
    static Self Identity(unsigned long size, DataType exemplar = DataTraits::zero(0))
    {
      Self m(size, size, DataTraits::zero(exemplar));
      for (unsigned long i = 0; i < size; i++)
        m(i,i) = DataTraits::one(exemplar);
      return m;
    }
    ///@}
  };

  typedef Matrix<int> Matrixi;
  typedef Matrix<short> Matrixs;
  typedef Matrix<long> Matrixl;
  typedef Matrix<float> Matrixf;
  typedef Matrix<double> Matrixd;
}
#endif
