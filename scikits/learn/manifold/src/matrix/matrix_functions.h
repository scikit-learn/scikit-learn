/**
 * \file matrix_functions.h
 * Contains functions that can be used on every
 */

#ifndef MATRIX_FUNCTIONS
#define MATRIX_FUNCTIONS

#include <cassert>
#include <ios>
#include <iosfwd>
#include <functional>

#include <boost/mpl/if.hpp>
#include <boost/mpl/equal_to.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/not_equal_to.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/has_xxx.hpp>

#include <boost/utility/enable_if.hpp>

#include <matrix/matrix_traits.h>

namespace Matrix
{
  template<class Vector>
  struct DiagonalMatrix;

  /** @name Metaprogramming helpers */
  /// @{
  BOOST_MPL_HAS_XXX_TRAIT_DEF(const_iterator)
  BOOST_MPL_HAS_XXX_TRAIT_DEF(Type)

  /// Struct used to choose the size
  template<bool cond, unsigned long a, unsigned long b>
  struct choiceBetweenSizes
  {
    enum {
      Size = (cond) ? a:b
    };
  };

  /// Struct to maximize static returns
  template<class MatrixClass1, class MatrixClass2>
  struct returnTypeStaticSize
  {
    typedef boost::mpl::int_<choiceBetweenSizes<MatrixClass1::staticHeight == 0, MatrixClass2::staticHeight, MatrixClass1::staticHeight>::Size > Height;
    typedef boost::mpl::int_<choiceBetweenSizes<MatrixClass1::staticWidth == 0, MatrixClass2::staticWidth, MatrixClass1::staticWidth>::Size > Width;

    typedef Matrix<typename MatrixClass1::Data_Type, Height::value, Width::value> Type;
  };

  /// Struct to maximize static returns from bool functions
  template<class MatrixClass1, class MatrixClass2>
  struct returnBoolTypeStaticSize
  {
    typedef boost::mpl::int_<choiceBetweenSizes<MatrixClass1::staticHeight == 0, MatrixClass2::staticHeight, MatrixClass1::staticHeight>::Size > Height;
    typedef boost::mpl::int_<choiceBetweenSizes<MatrixClass1::staticWidth == 0, MatrixClass2::staticWidth, MatrixClass1::staticWidth>::Size > Width;

    typedef Matrix<bool, Height::value, Width::value> Type;
  };

  /// Struct for multiplication returns
  template<class MatrixClass1, class MatrixClass2>
  struct returnTypeMultiplicationSize
  {
    typedef boost::mpl::int_<MatrixClass1::staticHeight> Height;
    typedef boost::mpl::int_<MatrixClass2::staticWidth> Width;

    typedef Matrix<typename MatrixClass1::Data_Type, Height::value, Width::value> Type;
  };

  /// Struct to authorize self operations
  template<class MatrixClass1, class MatrixClass2>
  struct selfReturn
  {
    typedef boost::mpl::int_<MatrixClass1::staticHeight> Height1;
    typedef boost::mpl::int_<MatrixClass2::staticHeight> Height2;
    typedef boost::mpl::int_<MatrixClass1::staticWidth> Width1;
    typedef boost::mpl::int_<MatrixClass2::staticWidth> Width2;

    typedef boost::mpl::and_<
          boost::mpl::or_<
            boost::mpl::or_<
              boost::mpl::equal_to<Height1, boost::mpl::int_<0> >,
              boost::mpl::equal_to<Height2, boost::mpl::int_<0> >
              >,
            boost::mpl::equal_to<Height1, Height2>
            >,
          boost::mpl::or_<
            boost::mpl::or_<
              boost::mpl::equal_to<Width1, boost::mpl::int_<0> >,
              boost::mpl::equal_to<Width2, boost::mpl::int_<0> >
              >,
            boost::mpl::equal_to<Width1, Width2>
            >
          > Result;
  };
  ///@}
#ifdef USE_ITERATOR_FUNCTIONS
}
# include <matrix/iterator_based_functions.h>
namespace Matrix
{
#endif

#ifndef USE_ITERATOR_FUNCTIONS
#ifndef USE_EXPRESSION_TEMPLATES
  /** @name Addition */
  /// @{
  /**
   * + operator : adding something
   * @param lhs is the left member
   * @param rhs is the right member
   * @return the new matrix
   */
  template<class MatrixType1, class MatrixType2>
  const typename returnTypeStaticSize<typename MatrixType1::Result, typename MatrixType2::Result>::Type operator+(const MatrixType1& lhs, const MatrixType2& rhs)
  {
    typename returnTypeStaticSize<typename MatrixType1::Result, typename MatrixType2::Result>::Type tempMat(lhs);
    tempMat += rhs;
    return (tempMat);
  }

  /**
   * + operator : adding a constant
   * @param lhs is the left member
   * @param rhs is the right member
   * @return the new matrix
   */
  template<class MatrixType>
  const typename MatrixType::Result operator+(const typename MatrixType::Data_Type lhs, const MatrixType& rhs)
  {
    typename MatrixType::Result tempMat(rhs);
    tempMat += lhs;
    return (tempMat);
  }

  /**
   * + operator : adding a constant
   * @param lhs is the left member
   * @param rhs is the right member
   * @return the new matrix
   */
  template<class MatrixType>
  const typename MatrixType::Result operator+(const MatrixType& lhs, const typename MatrixType::Data_Type rhs)
  {
    return (rhs + lhs);
  }

  /**
   * += operator : adding a constant
   * @param lhs is the matrix
   * @param rhs is the value to add
   * @return the modified matrix
   */
  template<class MatrixType>
  MatrixType& operator+=(MatrixType& lhs, const typename MatrixType::Data_Type rhs)
  {
    //std::for_each(begin(), end(), boost::bind(DataTraits::addAssign, _1, rhs));
    for(unsigned long j = 0; j < lhs.width(); ++j)
    {
      for(unsigned long i = 0; i < lhs.height(); ++i)
      {
        lhs(i, j) += rhs;
      }
    }
    return lhs;
  }

  /**
   * += operator : adding another matrix
   * @param lhs is the matrix
   * @param rhs is the matrix to add
   * @return the modified matrix
   */
  template<class MatrixType1, class MatrixType2>
  typename boost::enable_if<typename selfReturn<typename MatrixType1::Self, typename MatrixType2::Self>::Result, MatrixType1&>::type operator+=( MatrixType1& lhs, const MatrixType2& rhs)
  {
    assert(lhs.height() == rhs.height());
    assert(lhs.width() == rhs.width());

    for(unsigned long j = 0; j < lhs.width(); ++j)
    {
      for(unsigned long i = 0; i < lhs.height(); ++i)
      {
        lhs(i,j) += rhs(i,j);
      }
    }
    return lhs;
  }
  ///@}

  /** @name Unary operations */
  /// @{
  /**
   * - operator : opposes matrix
   * @return the opposed matrix
   */
  template<class MatrixType>
  const typename MatrixType::Result operator-(const MatrixType& matrix)
  {
    typename MatrixType::Result temp_mat(matrix.height(), matrix.width(), matrix(0, 0));
    std::transform(matrix.begin(), matrix.end(), temp_mat.begin(), std::negate<typename MatrixType::Data_Type>());
    return (temp_mat);
  }
  ///@}

  /** @name Substraction */
  /// @{
  /**
   * - operator : substracting something
   * @param lhs is the left member
   * @param rhs is the right member
   * @return the new matrix
   */
  template<class MatrixType1, class MatrixType2>
  const typename returnTypeStaticSize<typename MatrixType1::Result, typename MatrixType2::Result>::Type operator-(const MatrixType1& lhs, const MatrixType2& rhs)
  {
    typename returnTypeStaticSize<typename MatrixType1::Result, typename MatrixType2::Result>::Type tempMat(lhs);
    tempMat -= rhs;
    return (tempMat);
  }

  /**
   * - operator : substracting to a constant
   * @param lhs is the left member
   * @param rhs is the right member
   * @return the new matrix
   */
  template<class MatrixType>
  const typename MatrixType::Result operator-(const typename MatrixType::Data_Type lhs, const MatrixType& rhs)
  {
    typename MatrixType::Result tempMat(rhs);
    tempMat -= lhs;
    return -(tempMat);
  }

  /**
   * - operator : substracting a constant
   * @param lhs is the left member
   * @param rhs is the right member
   * @return the new matrix
   */
  template<class MatrixType>
  const typename MatrixType::Result operator-(const MatrixType& lhs, const typename MatrixType::Data_Type rhs)
  {
    typename MatrixType::Result tempMat(lhs);
    tempMat -= rhs;
    return (tempMat);
  }

  /**
   * -= operator : substracting a constant
   * @param lhs is the matrix
   * @param rhs is the value to add
   * @return the modified matrix
   */
  template<class MatrixType>
  MatrixType& operator-=(MatrixType& lhs, const typename MatrixType::Data_Type rhs)
  {
    //std::for_each(begin(), end(), boost::bind(DataTraits::addAssign, _1, rhs));
    for(unsigned long j = 0; j < lhs.width(); ++j)
    {
      for(unsigned long i = 0; i < lhs.height(); ++i)
      {
        lhs(i, j) -= rhs;
      }
    }
    return lhs;
  }

  /**
   * -= operator : substracting another matrix
   * @param lhs is the matrix
   * @param rhs is the matrix to add
   * @return the modified matrix
   */
  template<class MatrixType1, class MatrixType2>
  typename boost::enable_if<typename selfReturn<typename MatrixType1::Self, typename MatrixType2::Self>::Result, MatrixType1&>::type operator-=(MatrixType1& lhs, const MatrixType2& rhs)
  {
    assert(lhs.height() == rhs.height());
    assert(lhs.width() == rhs.width());

    for(unsigned long j = 0; j < lhs.width(); ++j)
    {
      for(unsigned long i = 0; i < lhs.height(); ++i)
      {
        lhs(i,j) -= rhs(i,j);
      }
    }
    return lhs;
  }
  ///@}

  /** @name Multiplication */
  /// @{
  /**
   * * operator : multiply by something
   * @param lhs is the left member
   * @param rhs is the right member
   * @return the new matrix
   */
  template<class MatrixType>
  const typename MatrixType::Result operator*(const MatrixType& lhs, const typename MatrixType::Data_Type rhs)
  {
    typename MatrixType::Result tempMat(lhs);
    tempMat *= rhs;
    return (tempMat);
  }

  /**
   * *= operator : multiply by something
   * @param lhs is the matrix
   * @param rhs is the value to multiply
   * @return the modified matrix
   */
  template<class MatrixType>
  MatrixType& operator*=(MatrixType& lhs, const typename MatrixType::Data_Type rhs)
  {
    for(unsigned long j = 0; j < lhs.width(); ++j)
    {
      for(unsigned long i = 0; i < lhs.height(); ++i)
      {
        lhs(i, j) *= rhs;
      }
    }
    return lhs;
  }

  /**
   * * operator : multiply the matrix by another one
   * @param lhs is the left member
   * @param rhs is the right member
   * @return the modified matrix
   */
  template<class MatrixType1, class MatrixType2>
  const typename returnTypeMultiplicationSize<typename MatrixType1::Result, typename MatrixType2::Result>::Type operator*(const MatrixType1& lhs, const MatrixType2& rhs)
  {
    assert(lhs.width() == rhs.height());

    typedef typename returnTypeMultiplicationSize<typename MatrixType1::Result, typename MatrixType2::Result>::Type returnType;
    returnType result(lhs.height(), rhs.width(), returnType::DataTraits::zero(lhs(0, 0) * rhs(0, 0)));

    for(unsigned long j = 0; j < rhs.width(); ++j)
    {
      for(unsigned long k = 0; k < lhs.width(); ++k)
      {
        for(unsigned long i = 0; i < lhs.height(); ++i)
        {
          result(i, j) += lhs(i,k) * rhs(k,j);
        }
      }
    }
    return result;
  }

  /**
   * * operator : useful for non-commutative DataType
   * @param lhs is a DataType data
   * @param rhs is the matrix being left-multiplied
   * @return lhs * rhs
   */
  template<class MatrixType>
  const typename MatrixType::Result operator*(const typename MatrixType::Data_Type lhs, const MatrixType& rhs)
  {
    typename MatrixType::Result tempMat(rhs);
    for(unsigned long j = 0; j < tempMat.width(); ++j)
    {
      for(unsigned long i = 0; i < tempMat.height(); ++i)
      {
        tempMat(i, j) = lhs * rhs(i, j);
      }
    }
    return (tempMat);
  }

  /**
   * * operator : multiply a diagonal matrix by a matrix
   * @param lhs is the left member
   * @param rhs is the right member
   * @return the modified matrix
   */
  template<class MatrixType1, class MatrixType2>
  const typename returnTypeMultiplicationSize<typename DiagonalMatrix<MatrixType1>::Result, typename MatrixType2::Result>::Type operator*(const DiagonalMatrix<MatrixType1>& lhs, const MatrixType2& rhs)
  {
    assert(lhs.width() == rhs.height());

    typedef typename returnTypeMultiplicationSize<typename DiagonalMatrix<MatrixType1>::Result, typename MatrixType2::Result>::Type returnType;
    returnType result(lhs.height(), rhs.width());

    for(unsigned long j = 0; j < rhs.width(); ++j)
    {
      for(unsigned long i = 0; i < lhs.height(); ++i)
      {
        result(i, j) = lhs(i,i) * rhs(i,j);
      }
    }
    return result;
  }

  /**
   * * operator : multiply the matrix by a diagonal matrix
   * @param lhs is the left member
   * @param rhs is the right member
   * @return the modified matrix
   */
#ifdef __GNUC__
  template<class MatrixType1, class MatrixType2>
  const typename returnTypeMultiplicationSize<typename MatrixType1::Result, typename DiagonalMatrix<MatrixType2>::Result>::Type operator*(const MatrixType1& lhs, const DiagonalMatrix<MatrixType2>& rhs)
  {
    assert(lhs.width() == rhs.height());

    typedef typename returnTypeMultiplicationSize<typename MatrixType1::Result, typename DiagonalMatrix<MatrixType2>::Result>::Type returnType;
    returnType result(lhs.height(), rhs.width());

    for(unsigned long j = 0; j < rhs.width(); ++j)
    {
      for(unsigned long i = 0; i < lhs.height(); ++i)
      {
        result(i, j) = lhs(i,j) * rhs(j,j);
      }
    }
    return result;
  }
#endif

  /**
   * * operator : multiply a diagonal matrix by a diagonal matrix
   * @param lhs is the left member
   * @param rhs is the right member
   * @return the modified matrix
   */
  template<class MatrixType1, class MatrixType2>
  const DiagonalMatrix<typename returnTypeStaticSize<typename MatrixType1::Result, typename MatrixType2::Result>::Type> operator*(const DiagonalMatrix<MatrixType1>& lhs, const DiagonalMatrix<MatrixType2>& rhs)
  {
    assert(lhs.height() == rhs.height());

    typedef DiagonalMatrix<typename returnTypeStaticSize<typename MatrixType1::Result, typename MatrixType2::Result>::Type> returnType;
    typename returnTypeStaticSize<typename MatrixType1::Result, typename MatrixType2::Result>::Type vector(lhs.height());

    for(unsigned long j = 0; j < rhs.width(); ++j)
    {
      vector(j) = lhs(j,j) * rhs(j,j);
    }
    return returnType(vector);
  }
  ///@}
#endif
#endif

  /**
   * * operator : multiply the matrix by another one
   * @param lhs is the left member
   * @param rhs is the right member
   * @return the modified matrix
   */
  template<class MatrixType1, class MatrixType2>
  typename boost::enable_if<has_Type<returnTypeMultiplicationSize<typename MatrixType1::Result, typename MatrixType2::Result> >, MatrixType1&>::type operator*=(MatrixType1& lhs, const MatrixType2& rhs)
  {
    lhs = lhs * rhs;
    return lhs;
  }

  /** @name Shift */
  /// @{
  /**
   * << operator : left shift by something
   * @param lhs is the left member
   * @param rhs is the right member
   * @return the new matrix
   */
  template<class MatrixType, class Something>
  const typename MatrixType::Result operator<<(const MatrixType& lhs, const Something rhs)
  {
    typename MatrixType::Result tempMat(lhs);
    tempMat <<= rhs;
    return (tempMat);
  }

  /**
   * >> operator : right shift by something
   * @param lhs is the left member
   * @param rhs is the right member
   * @return the new matrix
   */
  template<class MatrixType, class Something>
  const typename MatrixType::Result operator>>(const MatrixType& lhs, const Something rhs)
  {
    typename MatrixType::Result tempMat(lhs);
    tempMat >>= rhs;
    return (tempMat);
  }
  ///@}

  /** @name Usual operations */
  /// @{

  /**
   * Transposes the matrix
   * @return the transposed matrix
   */
  template<class MatrixType>
  const typename MatrixType::Transpose transpose(const MatrixType& mat)
  {
    typename MatrixType::Transpose temp_mat(mat.width(), mat.height(), MatrixType::DataTraits::conjuge(mat(0, 0)));
    for(unsigned long j = 0; j < mat.height(); ++j)
    {
      for(unsigned long i = 0; i < mat.width(); ++i)
      {
        temp_mat(i, j) = MatrixType::DataTraits::conjuge(mat(j,i));
      }
    }
    return (temp_mat);
  }

  /**
   * Calculates the inner product of a matrix
   * @return the inner product of the matrix
   */
  template<class MatrixType>
  const Matrix<typename MatrixType::Data_Type, MatrixType::staticWidth, MatrixType::staticWidth> inner(const MatrixType& mat)
  {
    return transpose(mat) * mat;
  }

  /**
   * Calculates the outer product of a matrix
   * @return the outer product of the matrix
   */
  template<class MatrixType>
  const Matrix<typename MatrixType::Data_Type, MatrixType::staticHeight, MatrixType::staticHeight> outer(const MatrixType& mat)
  {
    return mat * transpose(mat);
  }

  /**
   * Calculates the sum of the elements of a matrix
   * @return the sum of all elements in the matrix
   */
  template<class MatrixType>
  typename MatrixType::Data_Type sum(const MatrixType& me)
  {
    typename MatrixType::Data_Type sumTemp(MatrixType::DataTraits::zero(me(0, 0)));
    for(unsigned long j = 0; j < me.width(); ++j)
    {
      for(unsigned long i = 0; i < me.height(); ++i)
      {
        sumTemp += me(i,j);
      }
    }
    return sumTemp;
  }

  /**
   * Calculates the sum of the columns of a matrix
   * @return a line vector containing the sum of each column of the matrix
   */
  template<class MatrixType>
  const typename MatrixType::Line columnSum(const MatrixType& me)
  {
    typename MatrixType::Line sumTemp(1U, me.width(), MatrixType::DataTraits::zero(me(0, 0)));
    for(unsigned long j = 0; j < me.width(); ++j)
    {
      for(unsigned long i = 0; i < me.height(); ++i)
      {
        sumTemp(j) += me(i,j);
      }
    }
    return sumTemp;
  }

  /**
   * Calculates the sum of the lines of the a matrix
   * @return a column vector containing the sum of each line of the matrix
   */
  template<class MatrixType>
  const typename MatrixType::Column lineSum(const MatrixType& me)
  {
    typename MatrixType::Column sumTemp(me.height(), 1U, MatrixType::DataTraits::zero(me(0, 0)));
    for(unsigned long j = 0; j < me.width(); ++j)
    {
      for(unsigned long i = 0; i < me.height(); ++i)
      {
        sumTemp(i) += me(i,j);
      }
    }
    return sumTemp;
  }

  /**
   * Calculates the euclidian norm of the matrix
   * @return square the norm2 of the matrix
   */
  template<class MatrixType>
  typename MatrixType::Data_Type norm2(const MatrixType& me)
  {
    typename MatrixType::Data_Type norm_temp = MatrixType::DataTraits::zero(me(0, 0));
    for(unsigned long j = 0; j < me.width(); ++j)
    {
      for(unsigned long i = 0; i < me.height(); ++i)
      {
        norm_temp += MatrixType::DataTraits::conjuge(me(i,j)) * me(i,j);
      }
    }
    return norm_temp;
  }

  /**
   * Calculates the euclidian norm of the columns of matrix
   * @return a line vector containing the square of the norm of each column of the matrix
   */
  template<class MatrixType>
  const typename MatrixType::Line columnNorm2(const MatrixType& me)
  {
    typename MatrixType::Line norm_temp(1U, me.width(), MatrixType::DataTraits::zero(me(0, 0)));
    for(unsigned long j = 0; j < me.width(); ++j)
    {
      for(unsigned long i = 0; i < me.height(); ++i)
      {
        norm_temp(j) += MatrixType::DataTraits::conjuge(me(i,j)) * me(i,j);
      }
    }
    return norm_temp;
  }

  /**
   * Calculates the euclidian norm of the lines of the matrix
   * @return a column vector containing the square of the norm of each line of the matrix
   */
  template<class MatrixType>
  const typename MatrixType::Column lineNorm2(const MatrixType& me)
  {
    typename MatrixType::Column norm_temp(me.height(), 1U, MatrixType::DataTraits::zero(me(0, 0)));
    for(unsigned long j = 0; j < me.width(); ++j)
    {
      for(unsigned long i = 0; i < me.height(); ++i)
      {
        norm_temp(i) += MatrixType::DataTraits::conjuge(me(i,j)) * me(i,j);
      }
    }
    return norm_temp;
  }

  /**
   * Calculates the absolute norm of the matrix
   * @return the norm1 of the matrix
   */
  template<class MatrixType>
  typename MatrixType::Data_Type norm1(const MatrixType& me)
  {
    typename MatrixType::Data_Type norm_temp = MatrixType::DataTraits::zero(me(0, 0));
    for(unsigned long j = 0; j < me.width(); ++j)
    {
      for(unsigned long i = 0; i < me.height(); ++i)
      {
        norm_temp += MatrixType::DataTraits::absolute(me(i,j));
      }
    }
    return norm_temp;
  }

  /**
   * Calculates the absolute norm of the matrix
   * @return a line vector containing the absolute norm of each column of the matrix
   */
  template<class MatrixType>
  const typename MatrixType::Line columnNorm1(const MatrixType& me)
  {
    typename MatrixType::Line norm_temp(1U, me.width(), MatrixType::DataTraits::zero(me(0, 0)));
    for(unsigned long j = 0; j < me.width(); ++j)
    {
      for(unsigned long i = 0; i < me.height(); ++i)
      {
        norm_temp(j) += MatrixType::DataTraits::absolute(me(i,j));
      }
    }
    return norm_temp;
  }

  /**
   * Calculates the absolute norm of the matrix
   * @return a column vector containing the absolute norm of each line of the matrix
   */
  template<class MatrixType>
  const typename MatrixType::Column lineNorm1(const MatrixType& me)
  {
    typename MatrixType::Column norm_temp(me.height(), 1U, MatrixType::DataTraits::zero(me(0, 0)));
    for(unsigned long j = 0; j < me.width(); ++j)
    {
      for(unsigned long i = 0; i < me.height(); ++i)
      {
        norm_temp(i) += MatrixType::DataTraits::absolute(me(i,j));
      }
    }
    return norm_temp;
  }

  /**
   * Calculates the trace of a matrix
   * @return the trace of the matrix
   */
  template<class MatrixType>
  typename MatrixType::Data_Type trace(const MatrixType& matrix)
  {
    BOOST_STATIC_ASSERT(MatrixType::staticHeight == MatrixType::staticWidth);
    typename MatrixType::Data_Type count = MatrixType::DataTraits::zero(matrix(0, 0));
    for(unsigned long i = 0; i < matrix.height(); ++i)
    {
      count += matrix(i,i);
    }
    return count;
  }

  /**
   * Equality test
   * @param lhs is the left member
   * @param rhs is the right member
   * @return true if the matrix are equal
   */
  template<class MatrixType1, class MatrixType2>
  const typename returnBoolTypeStaticSize<typename MatrixType1::Result, typename MatrixType2::Result>::Type operator==(const MatrixType1& lhs, const MatrixType2& rhs)
  {
    typename returnBoolTypeStaticSize<typename MatrixType1::Result, typename MatrixType2::Result>::Type result(lhs.height(), lhs.width(), false);

    unsigned long width = lhs.width();
    unsigned long height = lhs.height();
    assert(lhs.height() == rhs.height());
    assert(lhs.width() == rhs.width());
    for(unsigned long j = 0; j < width; ++j)
    {
      for(unsigned long i = 0; i < height; ++i)
      {
        result(i,j) = (DataTypeTraits<typename MatrixType1::Data_Type>::absolute(lhs(i,j) - rhs(i,j)) <= DataTypeTraits<typename MatrixType1::Data_Type>::epsilon(lhs(i,j) - rhs(i,j)));
      }
    }
    return result;
  }

  /**
   * Equality test
   * @param lhs is the left member
   * @param rhs is the right member
   * @return true if the matrix are equal
   */
  template<class MatrixType>
  const typename MatrixType::ComparisonResult operator==(const MatrixType& lhs, const typename MatrixType::Data_Type rhs)
  {
    typename MatrixType::ComparisonResult result(lhs.height(), lhs.width(), false);

    unsigned long width = lhs.width();
    unsigned long height = lhs.height();
    for(unsigned long j = 0; j < width; ++j)
    {
      for(unsigned long i = 0; i < height; ++i)
      {
        result(i,j) = (DataTypeTraits<typename MatrixType::Data_Type>::absolute(lhs(i,j) - rhs) <= DataTypeTraits<typename MatrixType::Data_Type>::epsilon(lhs(i,j) - rhs));
      }
    }
    return result;
  }

  template<class MatrixType>
  const typename MatrixType::ComparisonResult operator==(const typename MatrixType::Data_Type lhs, const MatrixType& rhs)
  {
    return rhs == lhs;
  }

  /**
   * Different test
   * @param lhs is the left member
   * @param rhs is the right member
   * @return true if the matrix are different
   */
  template<class MatrixType1, class MatrixType2>
  const typename returnBoolTypeStaticSize<typename MatrixType1::Result, typename MatrixType2::Result>::Type operator!=(const MatrixType1& lhs, const MatrixType2& rhs)
  {
    typename returnBoolTypeStaticSize<typename MatrixType1::Result, typename MatrixType2::Result>::Type result(lhs.height(), lhs.width(), false);
    unsigned long width = lhs.width();
    unsigned long height = lhs.height();
    assert(lhs.height() == rhs.height());
    assert(lhs.width() == rhs.width());
    for(unsigned long j = 0; j < width; ++j)
    {
      for(unsigned long i = 0; i < height; ++i)
      {
        result(i,j) = (DataTypeTraits<typename MatrixType1::Data_Type>::absolute(lhs(i,j) - rhs(i,j)) >= DataTypeTraits<typename MatrixType1::Data_Type>::epsilon(lhs(i,j) - rhs(i,j)));
      }
    }
    return result;
  }

  /**
   * Different test
   * @param lhs is the left member
   * @param rhs is the right member
   * @return true if the matrix are different
   */
  template<class MatrixType>
  const typename MatrixType::ComparisonResult operator!=(const MatrixType& lhs, const typename MatrixType::Data_Type rhs)
  {
    typename MatrixType::ComparisonResult result(lhs.height(), lhs.width(), false);

    unsigned long width = lhs.width();
    unsigned long height = lhs.height();
    for(unsigned long j = 0; j < width; ++j)
    {
      for(unsigned long i = 0; i < height; ++i)
      {
        result(i,j) = (DataTypeTraits<typename MatrixType::Data_Type>::absolute(lhs(i,j) - rhs) >= DataTypeTraits<typename MatrixType::Data_Type>::epsilon(lhs(i,j) - rhs));
      }
    }
    return result;
  }

  template<class MatrixType>
  const typename MatrixType::ComparisonResult operator!=(const typename MatrixType::Data_Type lhs, const MatrixType& rhs)
  {
    return rhs != lhs;
  }

  /**
   * Inferiority test
   * @param lhs is the left member
   * @param rhs is the right member
   * @return true if lhs is inferior to rhs
   */
  template<class MatrixType1, class MatrixType2>
  const typename returnBoolTypeStaticSize<typename MatrixType1::Result, typename MatrixType2::Result>::Type operator<=(const MatrixType1& lhs, const MatrixType2& rhs)
  {
    typename returnBoolTypeStaticSize<typename MatrixType1::Result, typename MatrixType2::Result>::Type result(lhs.height(), lhs.width(), false);
    unsigned long width = lhs.width();
    unsigned long height = lhs.height();
    assert(lhs.height() == rhs.height());
    assert(lhs.width() == rhs.width());
    for(unsigned long j = 0; j < width; ++j)
    {
      for(unsigned long i = 0; i < height; ++i)
      {
        result(i,j) = (lhs(i,j) <= rhs(i,j));
      }
    }
    return result;
  }

  /**
   * Inferiority test
   * @param lhs is the left member
   * @param rhs is the right member
   * @return true if lhs is inferior to rhs
   */
  template<class MatrixType>
  const typename MatrixType::ComparisonResult operator<=(const MatrixType& lhs, const typename MatrixType::Data_Type rhs)
  {
    typename MatrixType::ComparisonResult result(lhs.height(), lhs.width(), false);

    unsigned long width = lhs.width();
    unsigned long height = lhs.height();
    for(unsigned long j = 0; j < width; ++j)
    {
      for(unsigned long i = 0; i < height; ++i)
      {
        result(i,j) = (lhs(i,j) <= rhs);
      }
    }
    return result;
  }

  template<class MatrixType>
  const typename MatrixType::ComparisonResult operator<=(const typename MatrixType::Data_Type lhs, const MatrixType& rhs)
  {
    return rhs >= lhs;
  }

  /**
   * Inferiority test - strict
   * @param lhs is the left member
   * @param rhs is the right member
   * @return true if lhs is inferior to rhs
   */
  template<class MatrixType1, class MatrixType2>
  const typename returnBoolTypeStaticSize<typename MatrixType1::Result, typename MatrixType2::Result>::Type operator<(const MatrixType1& lhs, const MatrixType2& rhs)
  {
    typename returnBoolTypeStaticSize<typename MatrixType1::Result, typename MatrixType2::Result>::Type result(lhs.height(), lhs.width(), false);
    unsigned long width = lhs.width();
    unsigned long height = lhs.height();
    assert(lhs.height() == rhs.height());
    assert(lhs.width() == rhs.width());
    for(unsigned long j = 0; j < width; ++j)
    {
      for(unsigned long i = 0; i < height; ++i)
      {
        result(i,j) = (lhs(i,j) < rhs(i,j));
      }
    }
    return result;
  }

  /**
   * Inferiority test - strict
   * @param lhs is the left member
   * @param rhs is the right member
   * @return true if lhs is inferior to rhs
   */
  template<class MatrixType>
  const typename MatrixType::ComparisonResult operator<(const MatrixType& lhs, const typename MatrixType::Data_Type rhs)
  {
    typename MatrixType::ComparisonResult result(lhs.height(), lhs.width(), false);

    unsigned long width = lhs.width();
    unsigned long height = lhs.height();
    for(unsigned long j = 0; j < width; ++j)
    {
      for(unsigned long i = 0; i < height; ++i)
      {
        result(i,j) = (lhs(i,j) < rhs);
      }
    }
    return result;
  }

  template<class MatrixType>
  const typename MatrixType::ComparisonResult operator<(const typename MatrixType::Data_Type lhs, const MatrixType& rhs)
  {
    return rhs > lhs;
  }

  /**
   * Superiority test
   * @param lhs is the left member
   * @param rhs is the right member
   * @return true if lhs is superior to rhs
   */
  template<class MatrixType1, class MatrixType2>
  const typename returnBoolTypeStaticSize<typename MatrixType1::Result, typename MatrixType2::Result>::Type operator>=(const MatrixType1& lhs, const MatrixType2& rhs)
  {
    typename returnBoolTypeStaticSize<typename MatrixType1::Result, typename MatrixType2::Result>::Type result(lhs.height(), lhs.width(), false);
    unsigned long width = lhs.width();
    unsigned long height = lhs.height();
    assert(lhs.height() == rhs.height());
    assert(lhs.width() == rhs.width());
    for(unsigned long j = 0; j < width; ++j)
    {
      for(unsigned long i = 0; i < height; ++i)
      {
        result(i,j) = (lhs(i,j) >= rhs(i,j));
      }
    }
    return result;
  }

  /**
   * Superiority test
   * @param lhs is the left member
   * @param rhs is the right member
   * @return true if lhs is superior to rhs
   */
  template<class MatrixType>
  const typename MatrixType::ComparisonResult operator>=(const MatrixType& lhs, const typename MatrixType::Data_Type rhs)
  {
    typename MatrixType::ComparisonResult result(lhs.height(), lhs.width(), false);

    unsigned long width = lhs.width();
    unsigned long height = lhs.height();
    for(unsigned long j = 0; j < width; ++j)
    {
      for(unsigned long i = 0; i < height; ++i)
      {
        result(i,j) = (lhs(i,j) >= rhs);
      }
    }
    return result;
  }

  template<class MatrixType>
  const typename MatrixType::ComparisonResult operator>=(const typename MatrixType::Data_Type lhs, const MatrixType& rhs)
  {
    return rhs <= lhs;
  }

  /**
   * Superiority test - strict
   * @param lhs is the left member
   * @param rhs is the right member
   * @return true if lhs is superior to rhs
   */
  template<class MatrixType1, class MatrixType2>
  const typename returnBoolTypeStaticSize<typename MatrixType1::Result, typename MatrixType2::Result>::Type operator>(const MatrixType1& lhs, const MatrixType2& rhs)
  {
    typename returnBoolTypeStaticSize<typename MatrixType1::Result, typename MatrixType2::Result>::Type result(lhs.height(), lhs.width(), false);
    unsigned long width = lhs.width();
    unsigned long height = lhs.height();
    assert(lhs.height() == rhs.height());
    assert(lhs.width() == rhs.width());
    for(unsigned long j = 0; j < width; ++j)
    {
      for(unsigned long i = 0; i < height; ++i)
      {
        result(i,j) = (lhs(i,j) > rhs(i,j));
      }
    }
    return result;
  }

  /**
   * Superiority test - strict
   * @param lhs is the left member
   * @param rhs is the right member
   * @return true if lhs is superior to rhs
   */
  template<class MatrixType>
  const typename MatrixType::ComparisonResult operator>(const MatrixType& lhs, const typename MatrixType::Data_Type rhs)
  {
    typename MatrixType::ComparisonResult result(lhs.height(), lhs.width(), false);

    unsigned long width = lhs.width();
    unsigned long height = lhs.height();
    for(unsigned long j = 0; j < width; ++j)
    {
      for(unsigned long i = 0; i < height; ++i)
      {
        result(i,j) = (lhs(i,j) > rhs);
      }
    }
    return result;
  }

  template<class MatrixType>
  const typename MatrixType::ComparisonResult operator>(const typename MatrixType::Data_Type lhs, const MatrixType& rhs)
  {
    return rhs < lhs;
  }
  ///@}

  /**
   * This << operator dumps the matrix to an ostream
   * @param stream is the stream on which the matrix will be dumped
   * @param matrix is the matrix being dumped
   * @return the stream with the dumped matrix
   */
  template<class MatrixType>
  typename MatrixType::ostream& operator<<(typename MatrixType::ostream &stream, const MatrixType& matrix)
  {
    std::ios_base::fmtflags flag = stream.flags();
    stream << "Size : " << matrix.height() << " * " << matrix.width()<< "\n";
    stream.flags(std::ios_base::scientific | stream.flags());
    unsigned long width = matrix.width();
    unsigned long height = matrix.height();
    for(unsigned long i = 0; i < height; ++i)
    {
      for(unsigned long j = 0; j < width; ++j)
      {
        stream << matrix(i, j) << ",\t";
      }
      stream << "\n";
    }
    stream.flags(flag);
    return stream;
  }

  /**
   * This >> operator reads a matrix from a istream
   * @param stream is the stream from which the matrix will be read
   * @param matrix is the matrix being read
   * @return the stream without the read matrix
   */
  template<class MatrixType>
  typename MatrixType::istream& operator>>(typename MatrixType::istream &stream, MatrixType& matrix)
  {
    std::string tmp;
    double height;
    double width;
    stream >> tmp >> tmp >> height >> tmp >> width;
    if(!stream.good())
      return stream;
    matrix = MatrixType(static_cast<unsigned long>(height+ 0.5), static_cast<unsigned long>(width + 0.5));
    for(unsigned long i = 0; i < height; ++i)
    {
      for(unsigned long j = 0; j < width; ++j)
      {
        stream >> matrix(i, j) >> tmp;
      }
    }
    return stream;
  }

  /**
   * This xml operator dumps the matrix to an ostream
   * @param stream is the stream on which the matrix will be dumped
   * @param matrix is the matrix being dumped
   * @return the stream with the dumped matrix
   */
  template<class MatrixType>
  typename MatrixType::ostream& toXML(typename MatrixType::ostream &stream, const MatrixType& matrix)
  {
    std::ios_base::fmtflags flag = stream.flags();
    stream << "<matrix><size width=\"" << matrix.width() << "\" height=\"" << matrix.height()<< "\" />\n";
    stream.flags(std::ios_base::scientific | stream.flags());
    unsigned long width = matrix.width();
    unsigned long height = matrix.height();
    for(unsigned long i = 0; i < height; ++i)
    {
      stream << "  <line>\n";
      for(unsigned long j = 0; j < width; ++j)
      {
        stream << matrix(i, j) << "\t";
      }
      stream << "  </line>\n";
    }
    stream.flags(flag);
    stream << "</matrix>\n";
    return stream;
  }

  /**
   * Tests if matrix is hermitian
   * @param matrix is the matrix to test
   * @return true if matrix is hermitian
   */
  template<class MatrixType>
  bool isHermitian(const MatrixType& matrix)
  {
    bool flag = true;
    for(unsigned long j = 0; j < matrix.width(); ++j)
    {
      for(unsigned long i = j + 1; i < matrix.height(); ++i)
      {
        flag &= (matrix(i,j) == MatrixType::DataTraits::conjuge(matrix(i, j)));
      }
    }

    return flag;
  }
}

#ifdef USE_BLAS
#include "cblas.h"
#endif
#endif
