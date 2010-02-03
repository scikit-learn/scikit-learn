/**
* \file iterator_based_functions.h
* Usual functions based on iterators
*/

#ifndef ITERATORBASEDFUNTIONS
#define ITERATORBASEDFUNTIONS

#include <matrix/matrix_functions.h>

namespace Matrix
{
  /** @name Addition (iterators) */
  /// @{
  /**
   * Iteration-based addition with a constant
   * @param matrix is the add matrix
   * @param data is the constant to add
   * @return a new matrix
   */
  template<class MatrixType>
  const typename MatrixType::Result itAdd(const MatrixType& matrix, typename MatrixType::Data_Type data)
  {
    typename MatrixType::Result result(matrix.height(), matrix.width());

    typename MatrixType::const_iterator it = matrix.begin(), itend = matrix.end();
    typename MatrixType::Result::iterator itresult = result.begin();
    for(; it != itend; ++it, ++itresult)
      *itresult = *it + data;

    return result;
  }

  /**
   * Iteration-based self-addition with a constant
   * @param matrix is the add matrix
   * @param data is the constant to add
   * @return a new matrix
   */
  template<class MatrixType>
  MatrixType& itSelfAdd(MatrixType& matrix, typename MatrixType::Data_Type data)
  {
    typename MatrixType::iterator it = matrix.begin(), itend = matrix.end();
    for(; it != itend; ++it)
      *it += data;

    return matrix;
  }

  /**
   * Iteration-based addition with a constant
   * @param matrix is the add matrix
   * @param data is the constant to add
   * @return a new matrix
   */
  template<class MatrixType>
  const typename MatrixType::Result itAdd(typename MatrixType::Data_Type data, const MatrixType& matrix)
  {
    typename MatrixType::Result result(matrix.height(), matrix.width());

    typename MatrixType::const_iterator it = matrix.begin(), itend = matrix.end();
    typename MatrixType::Result::iterator itresult = result.begin();
    for(; it != itend; ++it, ++itresult)
      *itresult = data + *it;

    return result;
  }

  /**
   * Iteration-based addition
   * @param matrix1 is the first term
   * @param matrix2 is the second term
   * @return a new matrix
   */
  template<class MatrixType1, class MatrixType2>
  const typename returnTypeStaticSize<typename MatrixType1::Result, typename MatrixType2::Result>::Type itAdd(const MatrixType1& matrix1, const MatrixType2& matrix2)
  {
    assert(matrix1.height() == matrix2.height());
    assert(matrix1.width() == matrix2.width());
    typename returnTypeStaticSize<typename MatrixType1::Result, typename MatrixType2::Result>::Type result(matrix1.height(), matrix1.width());

    std::transform(matrix1.begin(), matrix1.end(), matrix2.begin(), result.begin(), std::plus<typename MatrixType1::Data_Type>());
    return result;
  }

  /**
   * Iteration-based self addition
   * @param matrix1 is the matrix
   * @param matrix2 is the matrix to add
   * @return the modified matrix
   */
  template<class MatrixType1, class MatrixType2>
  typename boost::enable_if<typename selfReturn<typename MatrixType1::Self, typename MatrixType2::Self>::Result, MatrixType1&>::type itSelfAdd(MatrixType1& matrix1, const MatrixType2& matrix2)
  {
    assert(matrix1.height() == matrix2.height());
    assert(matrix1.width() == matrix2.width());

    typename MatrixType1::iterator it = matrix1.begin(), itend = matrix1.end();
    typename MatrixType2::const_iterator it2 = matrix2.begin();
    for(; it != itend; ++it)
    {
      *it += *it2;
      ++it2;
    }

    return matrix1;
  }
  ///@}

  /** @name Substraction (iterators) */
  /// @{
  /**
   * Iteration-based substraction with a constant
   * @param matrix is the matrix to substract
   * @param data is the substracted constant
   * @return a new matrix
   */
  template<class MatrixType>
  const typename MatrixType::Result itSub(const MatrixType& matrix, typename MatrixType::Data_Type data)
  {
    typename MatrixType::Result result(matrix.height(), matrix.width());

    typename MatrixType::const_iterator it = matrix.begin(), itend = matrix.end();
    typename MatrixType::Result::iterator itresult = result.begin();
    for(; it != itend; ++it, ++itresult)
      *itresult = *it - data;

    return result;
  }

  /**
   * Iteration-based self-substraction with a constant
   * @param matrix is the substracted matrix
   * @param data is the constant to substract with
   * @return the modified matrix
   */
  template<class MatrixType>
  MatrixType& itSelfSub(MatrixType& matrix, typename MatrixType::Data_Type data)
  {
    typename MatrixType::iterator it = matrix.begin(), itend = matrix.end();
    for(; it != itend; ++it)
      *it -= data;

    return matrix;
  }

  /**
   * Iteration-based substraction from a constant
   * @param matrix is the substracted matrix
   * @param data is the constant to substract
   * @return a new matrix
   */
  template<class MatrixType>
  const typename MatrixType::Result itSub(typename MatrixType::Data_Type data, const MatrixType& matrix)
  {
    typename MatrixType::Result result(matrix.height(), matrix.width());

    typename MatrixType::const_iterator it = matrix.begin(), itend = matrix.end();
    typename MatrixType::Result::iterator itresult = result.begin();
    for(; it != itend; ++it, ++itresult)
      *itresult = data - *it;

    return result;
  }

  /**
   * Iteration-based substraction
   * @param matrix1 is the first term
   * @param matrix2 is the second term
   * @return a new matrix
   */
  template<class MatrixType1, class MatrixType2>
  const typename returnTypeStaticSize<typename MatrixType1::Result, typename MatrixType2::Result>::Type itSub(const MatrixType1& matrix1, const MatrixType2& matrix2)
  {
    assert(matrix1.height() == matrix2.height());
    assert(matrix1.width() == matrix2.width());
    typename returnTypeStaticSize<typename MatrixType1::Result, typename MatrixType2::Result>::Type result(matrix1.height(), matrix1.width());

    std::transform(matrix1.begin(), matrix1.end(), matrix2.begin(), result.begin(), std::minus<typename MatrixType1::Data_Type>());
    return result;
  }

  /**
   * Iteration-based self substraction
   * @param matrix1 is the matrix
   * @param matrix2 is the matrix to substract
   * @return the modified matrix
   */
  template<class MatrixType1, class MatrixType2>
  typename boost::enable_if<typename selfReturn<typename MatrixType1::Self, typename MatrixType2::Self>::Result, MatrixType1&>::type itSelfSub(MatrixType1& matrix1, const MatrixType2& matrix2)
  {
    assert(matrix1.height() == matrix2.height());
    assert(matrix1.width() == matrix2.width());

    typename MatrixType1::iterator it = matrix1.begin(), itend = matrix1.end();
    typename MatrixType2::const_iterator it2 = matrix2.begin();
    for(; it != itend; ++it)
    {
      *it -= *it2;
      ++it2;
    }

    return matrix1;
  }
  ///@}

  /** @name Unary operations (iterators) */
  /// @{
  /**
   * Iteration-based oppose
   * @param matrix is the matrix to opppose
   * @return a new matrix
   */
  template<class MatrixType>
  const typename MatrixType::Result itOppose(const MatrixType& matrix)
  {
    typename MatrixType::Result result(matrix.height(), matrix.width());

    typename MatrixType::const_iterator it = matrix.begin(), itend = matrix.end();
    typename MatrixType::Result::iterator itresult = result.begin();
    for(; it != itend; ++it, ++itresult)
      *itresult = - *it;

    return result;
  }
  ///@}

  /** @name Multiplication (iterators) */
  /// @{
  /**
   * Iteration-based multiplication with a constant
   * @param matrix is the matrix to multiply
   * @param data is the multiplied constant
   * @return a new matrix
   */
  template<class MatrixType>
  const typename MatrixType::Result itMult(const MatrixType& matrix, typename MatrixType::Data_Type data)
  {
    typename MatrixType::Result result(matrix.height(), matrix.width());

    typename MatrixType::const_iterator it = matrix.begin(), itend = matrix.end();
    typename MatrixType::Result::iterator itresult = result.begin();
    for(; it != itend; ++it, ++itresult)
      *itresult = *it * data;

    return result;
  }

  /**
   * Iteration-based self-multiplication with a constant
   * @param matrix is the multiplied matrix
   * @param data is the constant to multiply with
   * @return the modified matrix
   */
  template<class MatrixType>
  MatrixType& itSelfMult(MatrixType& matrix, typename MatrixType::Data_Type data)
  {
    typename MatrixType::iterator it = matrix.begin(), itend = matrix.end();
    for(; it != itend; ++it)
      *it *= data;

    return matrix;
  }

  /**
   * Iteration-based multiplication with from a constant
   * @param matrix is the multiplied matrix
   * @param data is the constant to multiply
   * @return a new matrix
   */
  template<class MatrixType>
  const typename MatrixType::Result itMult(typename MatrixType::Data_Type data, const MatrixType& matrix)
  {
    typename MatrixType::Result result(matrix.height(), matrix.width());

    typename MatrixType::const_iterator it = matrix.begin(), itend = matrix.end();
    typename MatrixType::Result::iterator itresult = result.begin();
    for(; it != itend; ++it, ++itresult)
      *itresult = data * *it;

    return result;
  }

  /**
   * Iteration-based multiplication
   * @param matrix1 is the first term
   * @param matrix2 is the second term
   * @return a new matrix
   */
  template<class MatrixType1, class MatrixType2>
  const typename returnTypeMultiplicationSize<typename MatrixType1::Result, typename MatrixType2::Result>::Type itMult(const MatrixType1& matrix1, const MatrixType2& matrix2)
  {
    assert(matrix1.width() == matrix2.height());
    typedef typename returnTypeMultiplicationSize<typename MatrixType1::Result, typename MatrixType2::Result>::Type ReturnType;
    ReturnType result(matrix1.height(), matrix2.width(), DataTypeTraits<typename ReturnType::Data_Type>::zero(matrix1(0, 0) * matrix2(0, 0)));

/*    typename MatrixType1::const_iterator it1 = matrix1.begin();
    typename MatrixType2::const_iterator it2 = matrix2.begin();
    for(typename ReturnType::iterator it = result.begin(); it != result.end(); ++it)
    {
      typename MatrixType1::const_iterator itLine = it1.lineBegin(), itEndLine = it1.lineEnd();
      typename MatrixType2::const_iterator itColumn = it2.columnBegin(), itEndColumn = it2.columnEnd();
      for(; itLine != itEndLine; ++itLine, ++itColumn)
      {
        *it += *itLine * *itColumn;
      }
      ++it1;++it2;
    }*/
    typename MatrixType1::const_iterator it1 = matrix1.begin();
    typename MatrixType2::const_iterator it2 = matrix2.begin();
    typename ReturnType::iterator itresult = result.begin();

    for(typename MatrixType2::const_iterator it = matrix2.begin(); it != matrix2.end(); ++it)
    {
      typename ReturnType::iterator itResultColumn = itresult.columnBegin();
      typename ReturnType::const_iterator itEndColumn = itresult.columnEnd();
      for(; itResultColumn != itEndColumn; ++it1, ++itResultColumn)
      {
        *itResultColumn += *it1 * *it;
      }
      if(!(it1 != matrix1.end()))
      {
        it1 = matrix1.begin();
        std::advance(itresult, result.height());
      }
    }
    return result;
  }

  /**
   * Iteration-based multiplication with a diagonal matrix
   * @param lhs is the left member
   * @param rhs is the right member
   * @return the modified matrix
   */
  template<class MatrixType1, class MatrixType2>
  const typename returnTypeMultiplicationSize<typename DiagonalMatrix<MatrixType1>::Result, typename MatrixType2::Result>::Type itMult(const DiagonalMatrix<MatrixType1>& lhs, const MatrixType2& rhs)
  {
    assert(lhs.width() == rhs.height());

    typedef typename returnTypeMultiplicationSize<typename DiagonalMatrix<MatrixType1>::Result, typename MatrixType2::Result>::Type ReturnType;
    ReturnType result(lhs.height(), rhs.width());

    typename DiagonalMatrix<MatrixType1>::const_iterator it1 = lhs.begin();
    typename MatrixType2::const_iterator it2 = rhs.begin();
    typename ReturnType::iterator itresult = result.begin();
    typename ReturnType::iterator itresultColumn = itresult.columnBegin();
    for(typename MatrixType2::const_iterator it2Column = it2.columnBegin(); it2Column != it2.columnEnd(); ++it2Column)
    {
      typename ReturnType::iterator itresultLine = itresultColumn.lineBegin();
      for(typename MatrixType2::const_iterator it2Line = it2Column.lineBegin(); it2Line != it2Column.lineEnd(); ++it2Line)
      {
        *itresultLine = *it1 * *it2Line;
        ++itresultLine;
      }
      std::advance(it1, lhs.height() + 1);
      ++itresultColumn;
    }
    return result;
  }

  /**
   * Iteration-based multiplication with a diagonal matrix
   * @param lhs is the left member
   * @param rhs is the right member
   * @return the modified matrix
   */
#ifdef __GNUC__
  template<class MatrixType1, class MatrixType2>
  const typename returnTypeMultiplicationSize<typename MatrixType1::Result, typename DiagonalMatrix<MatrixType2>::Result>::Type itMult(const MatrixType1& lhs, const DiagonalMatrix<MatrixType2>& rhs)
  {
    assert(lhs.width() == rhs.height());

    typedef typename returnTypeMultiplicationSize<typename MatrixType1::Result, typename DiagonalMatrix<MatrixType2>::Result>::Type ReturnType;
    ReturnType result(lhs.height(), rhs.width());

    typename MatrixType1::const_iterator it1 = lhs.begin();
    typename DiagonalMatrix<MatrixType2>::const_iterator it2 = rhs.begin();
    typename ReturnType::iterator itresult = result.begin();
    typename ReturnType::iterator itresultColumn = itresult.columnBegin();
    for(typename MatrixType1::const_iterator it1Line = it1.lineBegin(); it1Line != it1.lineEnd(); ++it1Line)
    {
      for(typename MatrixType1::const_iterator it1Column = it1Line.columnBegin(); it1Column != it1Line.columnEnd(); ++it1Column)
      {
        *itresult = *it1Column * *it2;
        ++itresult;
      }
      std::advance(it2, lhs.height() + 1);
    }
    return result;
  }
#endif

  /**
   * Iteration-based multiplication with a diagonal matrix
   * @param lhs is the left member
   * @param rhs is the right member
   * @return the modified matrix
   */
  template<class MatrixType1, class MatrixType2>
  const DiagonalMatrix<typename returnTypeStaticSize<typename MatrixType1::Result, typename MatrixType2::Result>::Type> itMult(const DiagonalMatrix<MatrixType1>& lhs, const DiagonalMatrix<MatrixType2>& rhs)
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

#ifdef USE_ITERATOR_FUNCTIONS
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
    return itAdd(lhs, rhs);
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
    return itAdd(lhs, rhs);
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
    return itAdd(lhs, rhs);
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
    return itSelfAdd(lhs, rhs);
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
    return itSelfAdd(lhs, rhs);
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
    return itOppose(matrix);
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
    return itSub(lhs, rhs);
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
    return itSub(lhs, rhs);
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
    return itSub(lhs, rhs);
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
    return itSelfSub(lhs, rhs);
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
    return itSelfSub(lhs, rhs);
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
    return itMult(lhs, rhs);
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
    return itSelfMult(lhs, rhs);
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
    return itMult(lhs, rhs);
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
    return itMult(lhs, rhs);
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
    return itMult(lhs, rhs);
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
    return itMult(lhs, rhs);
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
    return itMult(lhs, rhs);
  }
  ///@}
#endif
}

#endif
