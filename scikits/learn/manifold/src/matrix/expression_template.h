/**
 * \file expression_template.h
 * Defines a set of templates that will be used for standard operations
 */

#ifndef EXPRESSIONTEMPLATE
#define EXPRESSIONTEMPLATE

#include <matrix/expression_template_iterators.h>

namespace Matrix
{
  /// Basic unary template
  template<class MatrixType, unsigned char op>
  struct UnaryExpression
  {
  };

  enum
  {Trans='t', Minus='-'};

  /// The unary transpose expression
  template<class MatrixType>
  struct UnaryExpression<MatrixType, 't'>
  {
    /// Static dimensions
    enum{staticHeight = MatrixType::staticWidth, staticWidth = MatrixType::staticHeight} ;
    /// The type of data inside the matrix
    typedef typename MatrixType::Data_Type Data_Type;
    /// The type of self
    typedef UnaryExpression<MatrixType, 't'> Self;
    /// Typedef of operation result
    typedef Matrix<Data_Type, staticHeight, staticWidth> Result;
    /// Typedef of comparison result
    typedef Matrix<bool, staticHeight, staticWidth> ComparisonResult;
    /// Typedef of column selection
    typedef typename ColumnVector<Data_Type, staticHeight>::Type Column;
    /// Typedef of line selection
    typedef typename LineVector<Data_Type, staticWidth>::Type Line;
    /// Typedef for <<
    typedef typename std::ostream ostream;
    /// Typedef for >>
    typedef typename std::istream istream;

    /// Typedef for data traits
    typedef DataTypeTraits<Data_Type> DataTraits;

    /**
     * Default constructor
     * @param matrix is the matrix that will be transposed
     * \warning verify the reference thing
     */
    UnaryExpression<MatrixType, 't'>(const MatrixType& matrix)
      :matrix(matrix)
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
      return matrix.width();
    }

    /**
     * Calculates the width of the matrix
     * @return the width of the matrix
     */
    unsigned long width() const
    {
      return matrix.height();
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
      return matrix(j, i);
    }
    ///@}

  private:
    /// A const ref to the matrix that will be transposed
    const MatrixType& matrix;
  };

  /**
   * Encapsulates the creation of a transposed expression
   * @param matrix is the expressiont emplate to transform
   * @return the transposed expression
   */
  template<class MatrixType>
  const UnaryExpression<MatrixType, 't'> makeTranspose(const MatrixType& matrix)
  {
    return UnaryExpression<MatrixType, 't'>(matrix);
  }

  /// The unary oppose expression
  template<class MatrixType>
  struct UnaryExpression<MatrixType, '-'>
  {
    /// Static dimensions
    enum{staticHeight = MatrixType::staticHeight, staticWidth = MatrixType::staticWidth} ;
    /// The type of data inside the matrix
    typedef typename MatrixType::Data_Type Data_Type;
    /// The type of self
    typedef UnaryExpression<MatrixType, '-'> Self;
    /// Typedef of operation result
    typedef Matrix<Data_Type, staticHeight, staticWidth> Result;
    /// Typedef of comparison result
    typedef Matrix<bool, staticHeight, staticWidth> ComparisonResult;
    /// Typedef of column selection
    typedef typename ColumnVector<Data_Type, staticHeight>::Type Column;
    /// Typedef of line selection
    typedef typename LineVector<Data_Type, staticWidth>::Type Line;
    /// Typedef for <<
    typedef typename std::ostream ostream;
    /// Typedef for >>
    typedef typename std::istream istream;

    /// Typedef for data traits
    typedef DataTypeTraits<Data_Type> DataTraits;

    /// Const iterator
    typedef OpposeExpressionIterator<const MatrixType> const_iterator;

    /**
     * Default constructor
     * @param matrix is the matrix that will be transposed
     * \warning verify the reference thing
     */
    UnaryExpression<MatrixType, '-'>(const MatrixType& matrix)
      :matrix(matrix)
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
      return matrix.height();
    }

    /**
     * Calculates the width of the matrix
     * @return the width of the matrix
     */
    unsigned long width() const
    {
      return matrix.width();
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
      return - matrix(i, j);
    }
    ///@}

    /**
     * Iterator to the beginning of the matrix
     * @return a const iterator on the beginning of the matrix/expression
     */
    const_iterator begin() const
    {
      return const_iterator(matrix.begin());
    }

    /**
     * Iterator to the end of the matrix
     * @return a const iterator on the end of the matrix/expression
     */
    const_iterator end() const
    {
      return const_iterator(matrix.end());
    }
  private:
    /// A const ref to the matrix that will be opposed
    const MatrixType& matrix;
  };

  /**
   * Encapsulates the creation of a opposed expression
   * @param matrix is the expressiont emplate to transform
   * @return the opposed expression
   */
  template<class MatrixType>
  const UnaryExpression<MatrixType, '-'> makeOppose(const MatrixType& matrix)
  {
    return UnaryExpression<MatrixType, '-'>(matrix);
  }

  /// Basic binary template
  template<class MatrixType1, class MatrixType2, unsigned char op>
  struct BinaryExpression
  {
  };

  enum
  {Add = '+', Sub = '-', Mult = '*'};

  /// The binary addition expression
  template<class MatrixType1, class MatrixType2>
  struct BinaryExpression<MatrixType1, MatrixType2, '+'>
  {
    /// The type of data inside the matrix
    typedef typename MatrixType1::Data_Type Data_Type;
    /// The type of self
    typedef BinaryExpression<MatrixType1, MatrixType2, '+'> Self;
    /// Typedef of operation result
    typedef typename returnTypeStaticSize<MatrixType1, MatrixType2>::Type Result;
    /// Static dimensions
    enum{staticHeight = Result::staticHeight, staticWidth = Result::staticWidth} ;
    /// Typedef of comparison result
    typedef Matrix<bool, staticHeight, staticWidth> ComparisonResult;
    /// Typedef of column selection
    typedef typename ColumnVector<Data_Type, staticHeight>::Type Column;
    /// Typedef of line selection
    typedef typename LineVector<Data_Type, staticWidth>::Type Line;
    /// Typedef for <<
    typedef typename std::ostream ostream;
    /// Typedef for >>
    typedef typename std::istream istream;

    /// Typedef for data traits
    typedef DataTypeTraits<Data_Type> DataTraits;

    /// Const iterator
    typedef AdditionExpressionIterator<const MatrixType1, const MatrixType2> const_iterator;

    /**
     * Default constructor
     * @param matrix1 is the first matrix that will be added
     * @param matrix2 is the second matrix that will be added
     * \warning verify the reference thing
     */
    BinaryExpression<MatrixType1, MatrixType2, '+'>(const MatrixType1& matrix1, const MatrixType2& matrix2)
      :matrix1(matrix1), matrix2(matrix2)
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
      return matrix1.height();
    }

    /**
     * Calculates the width of the matrix
     * @return the width of the matrix
     */
    unsigned long width() const
    {
      return matrix1.width();
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
      return matrix1(i, j) + matrix2(i, j);
    }
    ///@}

    /**
     * Iterator to the beginning of the matrix
     * @return a const iterator on the beginning of the matrix/expression
     */
    const_iterator begin() const
    {
      return const_iterator(matrix1.begin(), matrix2.begin());
    }

    /**
     * Iterator to the end of the matrix
     * @return a const iterator on the end of the matrix/expression
     */
    const_iterator end() const
    {
      return const_iterator(matrix1.end(), matrix2.end());
    }
  private:
    /// A const ref to the first matrix that will be added
    const MatrixType1& matrix1;
    /// A const ref to the second matrix that will be added
    const MatrixType2& matrix2;
  };

  /// The binary addition expression
  template<class MatrixType>
  struct BinaryExpression<MatrixType, typename MatrixType::Data_Type, '+'>
  {
    /// The type of data inside the matrix
    typedef typename MatrixType::Data_Type Data_Type;
    /// The type of self
    typedef BinaryExpression<MatrixType, typename MatrixType::Data_Type, '+'> Self;
    /// Typedef of operation result
    typedef typename MatrixType::Result Result;
    /// Static dimensions
    enum{staticHeight = Result::staticHeight, staticWidth = Result::staticWidth} ;
    /// Typedef of comparison result
    typedef Matrix<bool, staticHeight, staticWidth> ComparisonResult;
    /// Typedef of column selection
    typedef typename ColumnVector<Data_Type, staticHeight>::Type Column;
    /// Typedef of line selection
    typedef typename LineVector<Data_Type, staticWidth>::Type Line;
    /// Typedef for <<
    typedef typename std::ostream ostream;
    /// Typedef for >>
    typedef typename std::istream istream;

    /// Typedef for data traits
    typedef DataTypeTraits<Data_Type> DataTraits;

    /// Const iterator
    typedef AdditionConstantExpressionIterator<const MatrixType> const_iterator;
    /**
     * Default constructor
     * @param matrix is the matrix that will be added
     * @param element is the element that will be added
     * \warning verify the reference thing
     */
    BinaryExpression<MatrixType, typename MatrixType::Data_Type, '+'>(const MatrixType& matrix, const Data_Type& element)
      :matrix(matrix), element(element)
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
      return matrix.height();
    }

    /**
     * Calculates the width of the matrix
     * @return the width of the matrix
     */
    unsigned long width() const
    {
      return matrix.width();
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
      return matrix(i, j) + element;
    }
    ///@}

    /**
     * Iterator to the beginning of the matrix
     * @return a const iterator on the beginning of the matrix/expression
     */
    const_iterator begin() const
    {
      return const_iterator(matrix.begin(), element);
    }

    /**
     * Iterator to the end of the matrix
     * @return a const iterator on the end of the matrix/expression
     */
    const_iterator end() const
    {
      return const_iterator(matrix.end(), element);
    }
  private:
    /// A const ref to the first matrix that will be added
    const MatrixType& matrix;
    /// A const ref to the second matrix that will be added
    const Data_Type& element;
  };

  /// The binary addition expression
  template<class MatrixType>
  struct BinaryExpression<typename MatrixType::Data_Type, MatrixType, '+'>
  {
    /// The type of data inside the matrix
    typedef typename MatrixType::Data_Type Data_Type;
    /// The type of self
    typedef BinaryExpression<MatrixType, typename MatrixType::Data_Type, '+'> Self;
    /// Typedef of operation result
    typedef typename MatrixType::Result Result;
    /// Static dimensions
    enum{staticHeight = Result::staticHeight, staticWidth = Result::staticWidth} ;
    /// Typedef of comparison result
    typedef Matrix<bool, staticHeight, staticWidth> ComparisonResult;
    /// Typedef of column selection
    typedef typename ColumnVector<Data_Type, staticHeight>::Type Column;
    /// Typedef of line selection
    typedef typename LineVector<Data_Type, staticWidth>::Type Line;
    /// Typedef for <<
    typedef typename std::ostream ostream;
    /// Typedef for >>
    typedef typename std::istream istream;

    /// Typedef for data traits
    typedef DataTypeTraits<Data_Type> DataTraits;

    /// Const iterator
    typedef ConstantAdditionExpressionIterator<const MatrixType> const_iterator;

    /**
     * Default constructor
     * @param matrix is the matrix that will be added
     * @param element is the element that will be added
     * \warning verify the reference thing
     */
    BinaryExpression<typename MatrixType::Data_Type, MatrixType, '+'>(const Data_Type& element, const MatrixType& matrix)
      :matrix(matrix), element(element)
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
      return matrix.height();
    }

    /**
     * Calculates the width of the matrix
     * @return the width of the matrix
     */
    unsigned long width() const
    {
      return matrix.width();
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
      return element + matrix(i, j);
    }
    ///@}

    /**
     * Iterator to the beginning of the matrix
     * @return a const iterator on the beginning of the matrix/expression
     */
    const_iterator begin() const
    {
      return const_iterator(matrix.begin(), element);
    }

    /**
     * Iterator to the end of the matrix
     * @return a const iterator on the end of the matrix/expression
     */
    const_iterator end() const
    {
      return const_iterator(matrix.end(), element);
    }
  private:
    /// A const ref to the first matrix that will be added
    const MatrixType& matrix;
    /// A const ref to the second matrix that will be added
    const Data_Type& element;
  };

  /**
   * Encapsulates the creation of a addition expression
   * @param matrix is the expression template to transform
   * @return the addition expression
   */
  template<class MatrixType1, class MatrixType2>
  const BinaryExpression<MatrixType1, MatrixType2, '+'> makeAdd(const MatrixType1& matrix1, const MatrixType2& matrix2)
  {
    return BinaryExpression<MatrixType1, MatrixType2, '+'>(matrix1, matrix2);
  }

  /// The binary addition expression
  template<class MatrixType1, class MatrixType2>
  struct BinaryExpression<MatrixType1, MatrixType2, '-'>
  {
    /// The type of data inside the matrix
    typedef typename MatrixType1::Data_Type Data_Type;
    /// The type of self
    typedef BinaryExpression<MatrixType1, MatrixType2, '-'> Self;
    /// Typedef of operation result
    typedef typename returnTypeStaticSize<MatrixType1, MatrixType2>::Type Result;
    /// Static dimensions
    enum{staticHeight = Result::staticHeight, staticWidth = Result::staticWidth} ;
    /// Typedef of comparison result
    typedef Matrix<bool, staticHeight, staticWidth> ComparisonResult;
    /// Typedef of column selection
    typedef typename ColumnVector<Data_Type, staticHeight>::Type Column;
    /// Typedef of line selection
    typedef typename LineVector<Data_Type, staticWidth>::Type Line;
    /// Typedef for <<
    typedef typename std::ostream ostream;
    /// Typedef for >>
    typedef typename std::istream istream;

    /// Typedef for data traits
    typedef DataTypeTraits<Data_Type> DataTraits;

    /// Const iterator
    typedef SubstractionExpressionIterator<const MatrixType1, const MatrixType2> const_iterator;

    /**
     * Default constructor
     * @param matrix is the matrix that will be transposed
     * \warning verify the reference thing
     */
    BinaryExpression<MatrixType1, MatrixType2, '-'>(const MatrixType1& matrix1, const MatrixType2& matrix2)
      :matrix1(matrix1), matrix2(matrix2)
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
      return matrix1.height();
    }

    /**
     * Calculates the width of the matrix
     * @return the width of the matrix
     */
    unsigned long width() const
    {
      return matrix1.width();
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
      return matrix1(i, j) - matrix2(i, j);
    }

    /**
     * Iterator to the beginning of the matrix
     * @return a const iterator on the beginning of the matrix/expression
     */
    const_iterator begin() const
    {
      return const_iterator(matrix1.begin(), matrix2.begin());
    }

    /**
     * Iterator to the end of the matrix
     * @return a const iterator on the end of the matrix/expression
     */
    const_iterator end() const
    {
      return const_iterator(matrix1.end(), matrix2.end());
    }
    ///@}

  private:
    /// A const ref to the first matrix that will be substracted
    const MatrixType1& matrix1;
    /// A const ref to the second matrix that will be substracted
    const MatrixType2& matrix2;
  };

  /// The binary substraction expression
  template<class MatrixType>
  struct BinaryExpression<MatrixType, typename MatrixType::Data_Type, '-'>
  {
    /// The type of data inside the matrix
    typedef typename MatrixType::Data_Type Data_Type;
    /// The type of self
    typedef BinaryExpression<MatrixType, typename MatrixType::Data_Type, '-'> Self;
    /// Typedef of operation result
    typedef typename MatrixType::Result Result;
    /// Static dimensions
    enum{staticHeight = Result::staticHeight, staticWidth = Result::staticWidth} ;
    /// Typedef of comparison result
    typedef Matrix<bool, staticHeight, staticWidth> ComparisonResult;
    /// Typedef of column selection
    typedef typename ColumnVector<Data_Type, staticHeight>::Type Column;
    /// Typedef of line selection
    typedef typename LineVector<Data_Type, staticWidth>::Type Line;
    /// Typedef for <<
    typedef typename std::ostream ostream;
    /// Typedef for >>
    typedef typename std::istream istream;

    /// Typedef for data traits
    typedef DataTypeTraits<Data_Type> DataTraits;

    /// Const iterator
    typedef SubstractionConstantExpressionIterator<const MatrixType> const_iterator;

    /**
     * Default constructor
     * @param matrix is the matrix that will be substracted
     * @param element is the element that will be substracted
     * \warning verify the reference thing
     */
    BinaryExpression<MatrixType, typename MatrixType::Data_Type, '-'>(const MatrixType& matrix, const Data_Type& element)
      :matrix(matrix), element(element)
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
      return matrix.height();
    }

    /**
     * Calculates the width of the matrix
     * @return the width of the matrix
     */
    unsigned long width() const
    {
      return matrix.width();
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
      return matrix(i, j) - element;
    }

    /**
     * Iterator to the beginning of the matrix
     * @return a const iterator on the beginning of the matrix/expression
     */
    const_iterator begin() const
    {
      return const_iterator(matrix.begin(), element);
    }

    /**
     * Iterator to the end of the matrix
     * @return a const iterator on the end of the matrix/expression
     */
    const_iterator end() const
    {
      return const_iterator(matrix.end(), element);
    }
    ///@}

  private:
    /// A const ref to the first matrix that will be substracted
    const MatrixType& matrix;
    /// A const ref to the second matrix that will be substracted
    const Data_Type& element;
  };

  /// The binary substraction expression
  template<class MatrixType>
  struct BinaryExpression<typename MatrixType::Data_Type, MatrixType, '-'>
  {
    /// The type of data inside the matrix
    typedef typename MatrixType::Data_Type Data_Type;
    /// The type of self
    typedef BinaryExpression<MatrixType, typename MatrixType::Data_Type, '-'> Self;
    /// Typedef of operation result
    typedef typename MatrixType::Result Result;
    /// Static dimensions
    enum{staticHeight = Result::staticHeight, staticWidth = Result::staticWidth} ;
    /// Typedef of comparison result
    typedef Matrix<bool, staticHeight, staticWidth> ComparisonResult;
    /// Typedef of column selection
    typedef typename ColumnVector<Data_Type, staticHeight>::Type Column;
    /// Typedef of line selection
    typedef typename LineVector<Data_Type, staticWidth>::Type Line;
    /// Typedef for <<
    typedef typename std::ostream ostream;
    /// Typedef for >>
    typedef typename std::istream istream;

    /// Typedef for data traits
    typedef DataTypeTraits<Data_Type> DataTraits;

    /// Const iterator
    typedef ConstantSubstractionExpressionIterator<const MatrixType> const_iterator;

    /**
     * Default constructor
     * @param matrix is the matrix that will be substracted
     * @param element is the element that will be substracted
     * \warning verify the reference thing
     */
    BinaryExpression<typename MatrixType::Data_Type, MatrixType, '-'>(const Data_Type& element, const MatrixType& matrix)
      :matrix(matrix), element(element)
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
      return matrix.height();
    }

    /**
     * Calculates the width of the matrix
     * @return the width of the matrix
     */
    unsigned long width() const
    {
      return matrix.width();
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
      return element - matrix(i, j);
    }

    /**
     * Iterator to the beginning of the matrix
     * @return a const iterator on the beginning of the matrix/expression
     */
    const_iterator begin() const
    {
      return const_iterator(matrix.begin(), element);
    }

    /**
     * Iterator to the end of the matrix
     * @return a const iterator on the end of the matrix/expression
     */
    const_iterator end() const
    {
      return const_iterator(matrix.end(), element);
    }
    ///@}

  private:
    /// A const ref to the first matrix that will be substracted
    const MatrixType& matrix;
    /// A const ref to the second matrix that will be substracted
    const Data_Type& element;
  };

  /**
   * Encapsulates the creation of a substraction expression
   * @param matrix is the expression template to transform
   * @return the substraction expression
   */
  template<class MatrixType1, class MatrixType2>
  const BinaryExpression<MatrixType1, MatrixType2, '-'> makeSub(const MatrixType1& matrix1, const MatrixType2& matrix2)
  {
    return BinaryExpression<MatrixType1, MatrixType2, '-'>(matrix1, matrix2);
  }

  /// The binary multiplication expression
  template<class MatrixType1, class MatrixType2>
  struct BinaryExpression<MatrixType1, MatrixType2, '*'>
  {
    /// The type of data inside the matrix
    typedef typename MatrixType1::Data_Type Data_Type;
    /// The type of self
    typedef BinaryExpression<MatrixType1, MatrixType2, '*'> Self;
    /// Typedef of operation result
    typedef typename returnTypeMultiplicationSize<MatrixType1, MatrixType2>::Type Result;
    /// Static dimensions
    enum{staticHeight = Result::staticHeight, staticWidth = Result::staticWidth} ;
    /// Typedef of comparison result
    typedef Matrix<bool, staticHeight, staticWidth> ComparisonResult;
    /// Typedef of column selection
    typedef typename ColumnVector<Data_Type, staticHeight>::Type Column;
    /// Typedef of line selection
    typedef typename LineVector<Data_Type, staticWidth>::Type Line;
    /// Typedef for <<
    typedef typename std::ostream ostream;
    /// Typedef for >>
    typedef typename std::istream istream;

    /// Typedef for data traits
    typedef DataTypeTraits<Data_Type> DataTraits;

    /// Const iterator
    typedef MultiplicationExpressionIterator<const MatrixType1, const MatrixType2> const_iterator;
    /**
     * Default constructor
     * @param matrix is the matrix that will be transposed
     * \warning verify the reference thing
     */
    BinaryExpression<MatrixType1, MatrixType2, '*'>(const MatrixType1& matrix1, const MatrixType2& matrix2)
      :matrix1(matrix1), matrix2(matrix2)
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
      return matrix1.height();
    }

    /**
     * Calculates the width of the matrix
     * @return the width of the matrix
     */
    unsigned long width() const
    {
      return matrix2.width();
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
      Data_Type result = DataTypeTraits<Data_Type>::zero(matrix1(i, 0) * matrix2(0, j));
      for(unsigned long k = 0; k < matrix1.width(); ++k)
      {
        result += matrix1(i,k) * matrix2(k,j);
      }
      return result;
    }

    /**
     * Iterator to the beginning of the matrix
     * @return a const iterator on the beginning of the matrix/expression
     */
    const_iterator begin() const
    {
      return const_iterator(matrix1.begin(), matrix2.begin());
    }

    /**
     * Iterator to the end of the matrix
     * @return a const iterator on the end of the matrix/expression
     */
    const_iterator end() const
    {
      return const_iterator(matrix1.end(), matrix2.end());
    }
    ///@}

  private:
    /// A const ref to the first matrix that will be multiplied
    const MatrixType1& matrix1;
    /// A const ref to the second matrix that will be multiplied
    const MatrixType2& matrix2;
  };

  /// The binary multiplication expression
  template<class MatrixType1, class MatrixType2>
  struct BinaryExpression<DiagonalMatrix<MatrixType1>, MatrixType2, '*'>
  {
    /// The type of data inside the matrix
    typedef typename MatrixType1::Data_Type Data_Type;
    /// The type of self
    typedef BinaryExpression<DiagonalMatrix<MatrixType1>, MatrixType2, '*'> Self;
    /// Typedef of operation result
    typedef typename returnTypeMultiplicationSize<DiagonalMatrix<MatrixType1>, MatrixType2>::Type Result;
    /// Static dimensions
    enum{staticHeight = Result::staticHeight, staticWidth = Result::staticWidth} ;
    /// Typedef of comparison result
    typedef Matrix<bool, staticHeight, staticWidth> ComparisonResult;
    /// Typedef of column selection
    typedef typename ColumnVector<Data_Type, staticHeight>::Type Column;
    /// Typedef of line selection
    typedef typename LineVector<Data_Type, staticWidth>::Type Line;
    /// Typedef for <<
    typedef typename std::ostream ostream;
    /// Typedef for >>
    typedef typename std::istream istream;

    /// Typedef for data traits
    typedef DataTypeTraits<Data_Type> DataTraits;

    /**
     * Default constructor
     * @param matrix is the matrix that will be transposed
     * \warning verify the reference thing
     */
    BinaryExpression<DiagonalMatrix<MatrixType1>, MatrixType2, '*'>(const DiagonalMatrix<MatrixType1>& matrix1, const MatrixType2& matrix2)
      :matrix1(matrix1), matrix2(matrix2)
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
      return matrix1.height();
    }

    /**
     * Calculates the width of the matrix
     * @return the width of the matrix
     */
    unsigned long width() const
    {
      return matrix2.width();
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
      return matrix1(i,i) * matrix2(i,j);
    }
    ///@}

  private:
    /// A const ref to the first matrix that will be multiplied
    const DiagonalMatrix<MatrixType1>& matrix1;
    /// A const ref to the second matrix that will be multiplied
    const MatrixType2& matrix2;
  };

  /// The binary multiplication expression
  template<class MatrixType1, class MatrixType2>
  struct BinaryExpression<MatrixType1, DiagonalMatrix<MatrixType2>, '*'>
  {
    /// The type of data inside the matrix
    typedef typename MatrixType1::Data_Type Data_Type;
    /// The type of self
    typedef BinaryExpression<MatrixType1, DiagonalMatrix<MatrixType2>, '*'> Self;
    /// Typedef of operation result
    typedef typename returnTypeMultiplicationSize<MatrixType1, DiagonalMatrix<MatrixType2> >::Type Result;
    /// Static dimensions
    enum{staticHeight = Result::staticHeight, staticWidth = Result::staticWidth} ;
    /// Typedef of comparison result
    typedef Matrix<bool, staticHeight, staticWidth> ComparisonResult;
    /// Typedef of column selection
    typedef typename ColumnVector<Data_Type, staticHeight>::Type Column;
    /// Typedef of line selection
    typedef typename LineVector<Data_Type, staticWidth>::Type Line;
    /// Typedef for <<
    typedef typename std::ostream ostream;
    /// Typedef for >>
    typedef typename std::istream istream;

    /// Typedef for data traits
    typedef DataTypeTraits<Data_Type> DataTraits;

    /**
     * Default constructor
     * @param matrix is the matrix that will be transposed
     * \warning verify the reference thing
     */
    BinaryExpression<MatrixType1, DiagonalMatrix<MatrixType2>, '*'>(const MatrixType1& matrix1, const DiagonalMatrix<MatrixType2>& matrix2)
      :matrix1(matrix1), matrix2(matrix2)
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
      return matrix1.height();
    }

    /**
     * Calculates the width of the matrix
     * @return the width of the matrix
     */
    unsigned long width() const
    {
      return matrix2.width();
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
      return matrix1(i,j) * matrix2(j,j);
    }
    ///@}

  private:
    /// A const ref to the first matrix that will be multiplied
    const MatrixType1& matrix1;
    /// A const ref to the second matrix that will be multiplied
    const DiagonalMatrix<MatrixType2>& matrix2;
  };

  /// The binary multiplication expression
  template<class MatrixType>
  struct BinaryExpression<MatrixType, typename MatrixType::Data_Type, '*'>
  {
    /// The type of data inside the matrix
    typedef typename MatrixType::Data_Type Data_Type;
    /// The type of self
    typedef BinaryExpression<MatrixType, typename MatrixType::Data_Type, '*'> Self;
    /// Typedef of operation result
    typedef MatrixType Result;
    /// Static dimensions
    enum{staticHeight = Result::staticHeight, staticWidth = Result::staticWidth} ;
    /// Typedef of comparison result
    typedef Matrix<bool, staticHeight, staticWidth> ComparisonResult;
    /// Typedef of column selection
    typedef typename ColumnVector<Data_Type, staticHeight>::Type Column;
    /// Typedef of line selection
    typedef typename LineVector<Data_Type, staticWidth>::Type Line;
    /// Typedef for <<
    typedef typename std::ostream ostream;
    /// Typedef for >>
    typedef typename std::istream istream;

    /// Typedef for data traits
    typedef DataTypeTraits<Data_Type> DataTraits;

    /// Const iterator
    typedef ConstantMultiplicationExpressionIterator<const MatrixType> const_iterator;

    /**
     * Default constructor
     * @param matrix is the matrix that will be substracted
     * @param element is the element that will be substracted
     * \warning verify the reference thing
     */
    BinaryExpression<MatrixType, typename MatrixType::Data_Type, '*'>(const MatrixType& matrix, const Data_Type& element)
      :matrix(matrix), element(element)
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
      return matrix.height();
    }

    /**
     * Calculates the width of the matrix
     * @return the width of the matrix
     */
    unsigned long width() const
    {
      return matrix.width();
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
      return matrix(i, j) * element;
    }

    /**
     * Iterator to the beginning of the matrix
     * @return a const iterator on the beginning of the matrix/expression
     */
    const_iterator begin() const
    {
      return const_iterator(matrix.begin(), element);
    }

    /**
     * Iterator to the end of the matrix
     * @return a const iterator on the end of the matrix/expression
     */
    const_iterator end() const
    {
      return const_iterator(matrix.end(), element);
    }
    ///@}

  private:
    /// A const ref to the first matrix that will be multiplied
    const MatrixType& matrix;
    /// A const ref to the second matrix that will be multiplied
    const Data_Type& element;
  };

  /// The binary multiplication expression
  template<class MatrixType>
  struct BinaryExpression<typename MatrixType::Data_Type, MatrixType, '*'>
  {
    /// The type of data inside the matrix
    typedef typename MatrixType::Data_Type Data_Type;
    /// The type of self
    typedef BinaryExpression<MatrixType, typename MatrixType::Data_Type, '*'> Self;
    /// Typedef of operation result
    typedef MatrixType Result;
    /// Static dimensions
    enum{staticHeight = Result::staticHeight, staticWidth = Result::staticWidth} ;
    /// Typedef of comparison result
    typedef Matrix<bool, staticHeight, staticWidth> ComparisonResult;
    /// Typedef of column selection
    typedef typename ColumnVector<Data_Type, staticHeight>::Type Column;
    /// Typedef of line selection
    typedef typename LineVector<Data_Type, staticWidth>::Type Line;
    /// Typedef for <<
    typedef typename std::ostream ostream;
    /// Typedef for >>
    typedef typename std::istream istream;

    /// Typedef for data traits
    typedef DataTypeTraits<Data_Type> DataTraits;

    /// Const iterator
    typedef MultiplicationConstantExpressionIterator<const MatrixType> const_iterator;

    /**
     * Default constructor
     * @param matrix is the matrix that will be multiplied
     * @param element is the element that will be multiplied
     * \warning verify the reference thing
     */
    BinaryExpression<typename MatrixType::Data_Type, MatrixType, '*'>(const Data_Type& element, const MatrixType& matrix)
      :matrix(matrix), element(element)
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
      return matrix.height();
    }

    /**
     * Calculates the width of the matrix
     * @return the width of the matrix
     */
    unsigned long width() const
    {
      return matrix.width();
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
      return element * matrix(i, j);
    }

    /**
     * Iterator to the beginning of the matrix
     * @return a const iterator on the beginning of the matrix/expression
     */
    const_iterator begin() const
    {
      return const_iterator(matrix.begin(), element);
    }

    /**
     * Iterator to the end of the matrix
     * @return a const iterator on the end of the matrix/expression
     */
    const_iterator end() const
    {
      return const_iterator(matrix.end(), element);
    }
    ///@}

  private:
    /// A const ref to the first matrix that will be multiplied
    const MatrixType& matrix;
    /// A const ref to the second matrix that will be multiplied
    const Data_Type& element;
  };

  /**
   * Encapsulates the creation of a multiplication expression
   * @param matrix is the expression template to transform
   * @return the multiplication expression
   */
  template<class MatrixType1, class MatrixType2>
  const BinaryExpression<MatrixType1, MatrixType2, '*'> makeMult(const MatrixType1& matrix1, const MatrixType2& matrix2)
  {
    return BinaryExpression<MatrixType1, MatrixType2, '*'>(matrix1, matrix2);
  }
}
#endif
