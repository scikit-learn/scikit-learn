/**
 * \file sub_vector_matrix_lib.h
 * Describes a view of a matrix by means of vector matrix
 */

#ifndef SUBVECTORMATRIX
#define SUBVECTORMATRIX

#include <matrix/sub_vector_matrix_iterator.h>

namespace Matrix
{
  /// A view of a matrix with some lines, but with every column
  template<class Matrix_Type, class ColumnVectorForLine>
  class SubLineMatrix
  {
  public:
    typedef Matrix_Type MatrixType;
    /// Static dimensions
    enum{staticHeight = ColumnVectorForLine::staticHeight, staticWidth = MatrixType::staticWidth} ;

    typedef ColumnVectorForLine IterateLineVector;
    /// Only to have a link on the type used in the matrix
    typedef typename MatrixType::Data_Type Data_Type;
    /// Typedef of self
    typedef SubLineMatrix Self;
    /// Typedef of operation result
    typedef Matrix<Data_Type, staticHeight, staticWidth> Result;
    /// Typedef of comparison result
    typedef Matrix<bool, staticHeight, staticWidth> ComparisonResult;
    /// Typedef of transposition
    typedef Matrix<Data_Type, staticWidth, staticHeight> Transpose;
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
    /// Typedef of iterator to an element
    typedef SubLineMatrixIterator<Self> iterator;
    /// Typedef of const iterator to an element
    typedef SubLineMatrixIterator<const Self> const_iterator;

    /**
     * Constructor of a submatrix with column view
     * @param mat is the matrix from which a submatrix is used
     * @param lines is a column matrix of the used lines
     */
    SubLineMatrix(MatrixType& mat, const ColumnVectorForLine& lines)
      :matrix(&mat), lines(lines)
    {}

    /** @name Other operations (submatrix)*/
    /// @{
    /**
     * Returns the height of the matrix
     * @return the height of the matrix
     */
    unsigned long realHeight() const
    {
      return matrix->height();
    }
    /**
     * Returns the width of the matrix
     * @return the width of the matrix
     */
    unsigned long realWidth() const
    {
      return matrix->width();
    }
    /**
     * Returns the height of the matrix
     * @return the height of the matrix
     */
    unsigned long height() const
    {
      return lines.height();
    }
    /**
     * Returns the width of the matrix
     * @return the width of the matrix
     */
    unsigned long width() const
    {
      return matrix->width();
    }
    /**
     * Returns the number of columns of the submatrix
     * @return the number of columns of the submatrix
     */
    unsigned long nbColumns() const
    {
      return width();
    }
    /**
     * Returns the number of lines of the submatrix
     * @return the number of lines of the submatrix
     */
    unsigned long nbLines() const
    {
      return height();
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
        for (unsigned long i = 0; i < height(); ++i)
          (*this)(i,j) = mat(i,j);
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
        for (unsigned long i = 0; i < height(); ++i)
          (*this)(i,j) = rhs;

      return (*this);
    }

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
     * Return an iterator to the beginning
     * @return an iterator to the beginning
     */
    iterator begin()
    {
      typename MatrixType::iterator index = matrix->begin();
      std::advance(index, lines(0));
      return iterator(this, index, lines.begin(), lines.begin(), matrix);
    }

    /**
     * Return an iterator to the beginning
     * @return an iterator to the beginning
     */
    const_iterator begin() const
    {
      typename MatrixType::const_iterator index = matrix->cbegin();
      std::advance(index, lines(0));
      return const_iterator(this, index, lines.begin(), lines.begin(), matrix);
    }

    /**
     * Return an iterator to the end
     * @return an iterator to the end
     */
    const_iterator end() const
    {
      typename MatrixType::const_iterator index = matrix->cbegin();
      std::advance(index, lines(0) + matrix->width() * matrix->height());
      return const_iterator(this, index, lines.begin(), lines.begin(), matrix);
    }

    /**
     * Return an iterator to the end
     * @return an iterator to the end
     */
    iterator end()
    {
      typename MatrixType::iterator index = matrix->begin();
      std::advance(index, lines(0) + matrix->width() * matrix->height());
      return iterator(this, index, lines.begin(), lines.begin(), matrix);
    }

    /**
     * Subscript operator
     * @param i is the line selected
     * @param j is the column selected
     * @return the element in the matrix at the emplacement [i,j]
     */
    Data_Type& operator()(unsigned long i, unsigned long j)
    {
      assert(i < height());
      assert(j < width());
      return (*matrix)(lines(i), j);
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
      return (*matrix)(lines(i), j);
    }
    ///@}

  protected:
    /// The matrix used by the submatrix
    MatrixType* matrix;
    /// The lines used by this sub matrix
    ColumnVectorForLine lines;
  };

  /// A view of a matrix with every line, but with some columns
  template<class Matrix_Type, class ColumnVectorForColumn>
  class SubColumnMatrix
  {
  public:
    typedef Matrix_Type MatrixType;
    /// Static dimensions
    enum{staticHeight = MatrixType::staticHeight, staticWidth = ColumnVectorForColumn::staticHeight} ;

    typedef ColumnVectorForColumn IterateColumnVector;
    /// Only to have a link on the type used in the matrix
    typedef typename MatrixType::Data_Type Data_Type;
    /// Typedef of self
    typedef SubColumnMatrix Self;
    /// Typedef of operation result
    typedef Matrix<Data_Type, staticHeight, staticWidth> Result;
    /// Typedef of comparison result
    typedef Matrix<bool, staticHeight, staticWidth> ComparisonResult;
    /// Typedef of transposition
    typedef Matrix<Data_Type, staticWidth, staticHeight> Transpose;
    /// Typedef of column selection
    typedef typename ColumnVector<Data_Type, staticHeight>::Type Column;
    /// Typedef of line selection
    typedef typename LineVector<Data_Type, staticWidth>::Type Line;
    /// Typedef for <<
    typedef typename std::ostream ostream;

    /// Typedef for data traits
    typedef DataTypeTraits<Data_Type> DataTraits;

    /// Typedef of iterator to an element
     typedef SubColumnMatrixIterator<Self> iterator;
    /// Typedef of const iterator to an element
     typedef SubColumnMatrixIterator<const Self> const_iterator;

    /**
     * Constructor of a submatrix with column view
     * @param mat is the matrix from which a submatrix is used
     * @param columns is a column matrix of the used columns
     */
    SubColumnMatrix(MatrixType& mat, const ColumnVectorForColumn& columns)
      :matrix(&mat), columns(columns)
    {}

    /** @name Other operations (submatrix)*/
    /// @{
    /**
     * Returns the height of the matrix
     * @return the height of the matrix
     */
    unsigned long realHeight() const
    {
      return matrix->height();
    }
    /**
     * Returns the width of the matrix
     * @return the width of the matrix
     */
    unsigned long realWidth() const
    {
      return matrix->width();
    }
    /**
     * Returns the height of the matrix
     * @return the height of the matrix
     */
    unsigned long height() const
    {
      return matrix->height();
    }
    /**
     * Returns the width of the matrix
     * @return the width of the matrix
     */
    unsigned long width() const
    {
      return columns.height();
    }
    /**
     * Returns the number of columns of the submatrix
     * @return the number of columns of the submatrix
     */
    unsigned long nbColumns() const
    {
      return columns.height();
    }
    /**
     * Returns the number of lines of the submatrix
     * @return the number of lines of the submatrix
     */
    unsigned long nbLines() const
    {
      return matrix->height();
    }

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
     * Return an iterator to the beginning
     * @return an iterator to the beginning
     */
    iterator begin()
    {
      typename MatrixType::iterator index = matrix->begin();
      std::advance(index, columns(0) * matrix->height());
      return iterator(this, index, columns.begin(), columns.begin(), columns.end());
    }

    /**
     * Return an iterator to the beginning
     * @return an iterator to the beginning
     */
    const_iterator begin() const
    {
      typename MatrixType::const_iterator index = matrix->cbegin();
      std::advance(index, columns(0) * matrix->height());
      return const_iterator(this, index, columns.begin(), columns.begin(), columns.end());
    }

    /**
     * Return an iterator to the end
     * @return an iterator to the end
     */
    const_iterator end() const
    {
      typename MatrixType::const_iterator index = matrix->cbegin();
      std::advance(index, columns(width()) * matrix->height());
      return const_iterator(this, index, columns.end(), columns.begin(), columns.end());
    }

    /**
     * Return an iterator to the end
     * @return an iterator to the end
     */
    iterator end()
    {
      typename MatrixType::iterator index = matrix->begin();
      std::advance(index, columns(width()) * matrix->height());
      return iterator(this, index, columns.end(), columns.begin(), columns.end());
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
        for (unsigned long i = 0; i < height(); ++i)
          (*this)(i,j) = mat(i,j);
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
        for (unsigned long i = 0; i < height(); ++i)
          (*this)(i,j) = rhs;

      return (*this);
    }

    /**
     * Subscript operator
     * @param i is the line selected
     * @param j is the column selected
     * @return the element in the matrix at the emplacement [i,j]
     */
    Data_Type& operator()(unsigned long i, unsigned long j)
    {
      assert(i < height());
      assert(j < width());
      return (*matrix)(i, columns(j));
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
      return (*matrix)(i, columns(j));
    }
    ///@}

  protected:
    /// The matrix used by the submatrix
    MatrixType* matrix;
    /// The lines used by this sub matrix
    ColumnVectorForColumn columns;
  };


  /// A view of a matrix with some columns and lines
  template<class Matrix_Type, class ColumnVectorForLine, class ColumnVectorForColumn>
  class SubColumnLineMatrix
  {
  public:
    typedef Matrix_Type MatrixType;
    /// Static dimensions
    enum{staticHeight = ColumnVectorForLine::staticHeight, staticWidth = ColumnVectorForColumn::staticHeight} ;

    typedef ColumnVectorForLine IterateLineVector;
    typedef ColumnVectorForColumn IterateColumnVector;
    /// Only to have a link on the type used in the matrix
    typedef typename MatrixType::Data_Type Data_Type;
    /// Typedef of self
    typedef SubColumnLineMatrix Self;
    /// Typedef of operation result
    typedef Matrix<Data_Type, staticHeight, staticWidth> Result;
    /// Typedef of comparison result
    typedef Matrix<bool, staticHeight, staticWidth> ComparisonResult;
    /// Typedef of transposition
    typedef Matrix<Data_Type, staticWidth, staticHeight> Transpose;
    /// Typedef of column selection
    typedef typename ColumnVector<Data_Type, staticHeight>::Type Column;
    /// Typedef of line selection
    typedef typename LineVector<Data_Type, staticWidth>::Type Line;
    /// Typedef for <<
    typedef typename std::ostream ostream;

    /// Typedef for data traits
    typedef DataTypeTraits<Data_Type> DataTraits;
    /// Typedef of iterator to an element
    typedef SubColumnLineMatrixIterator<Self> iterator;
    /// Typedef of const iterator to an element
    typedef SubColumnLineMatrixIterator<const Self> const_iterator;

    /**
     * Constructor of a submatrix with column view
     * @param mat is the matrix from which a submatrix is used
     * @param columns is a column matrix of the used columns
     * @param lines is a column matrix of the used lines
     */
    SubColumnLineMatrix(MatrixType& mat, const ColumnVectorForLine& columns, const ColumnVectorForColumn& lines)
      :matrix(&mat), columns(columns), lines(lines)
    {}

    /** @name Other operations (submatrix)*/
    /// @{
    /**
     * Returns the height of the matrix
     * @return the height of the matrix
     */
    unsigned long realHeight() const
    {
      return matrix->height();
    }
    /**
     * Returns the width of the matrix
     * @return the width of the matrix
     */
    unsigned long realWidth() const
    {
      return matrix->width();
    }
    /**
     * Returns the height of the matrix
     * @return the height of the matrix
     */
    unsigned long height() const
    {
      return lines.height();
    }
    /**
     * Returns the width of the matrix
     * @return the width of the matrix
     */
    unsigned long width() const
    {
      return columns.height();
    }
    /**
     * Returns the number of columns of the submatrix
     * @return the number of columns of the submatrix
     */
    unsigned long nbColumns() const
    {
      return columns.height();
    }
    /**
     * Returns the number of lines of the submatrix
     * @return the number of lines of the submatrix
     */
    unsigned long nbLines() const
    {
      return lines.height();
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
        for (unsigned long i = 0; i < height(); ++i)
          (*this)(i,j) = mat(i,j);
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
        for (unsigned long i = 0; i < height(); ++i)
          (*this)(i,j) = rhs;

      return (*this);
    }

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
     * Return an iterator to the beginning
     * @return an iterator to the beginning
     */
     iterator begin()
     {
       typename MatrixType::iterator index = matrix->begin();
       std::advance(index, lines(0) + columns(0) * matrix->height());
       return iterator(this, index, columns.begin(), columns.begin(), columns.end(), lines.begin(), lines.begin(), matrix);
     }

    /**
     * Return an iterator to the beginning
     * @return an iterator to the beginning
     */
     const_iterator begin() const
     {
       typename MatrixType::const_iterator index = matrix->cbegin();
       std::advance(index, lines(0) + columns(0) * matrix->height());
       return const_iterator(this, index, columns.begin(), columns.begin(), columns.end(), lines.begin(), lines.begin(), matrix);
     }

    /**
     * Return an iterator to the end
     * @return an iterator to the end
     */
     const_iterator end() const
     {
       typename MatrixType::const_iterator index = matrix->cbegin();
       std::advance(index, lines(0) + columns(width()) * matrix->height());
       return const_iterator(this, index, columns.end(), columns.begin(), columns.end(), lines.begin(), lines.begin(), matrix);
     }

    /**
     * Return an iterator to the end
     * @return an iterator to the end
     */
     iterator end()
     {
       typename MatrixType::iterator index = matrix->begin();
       std::advance(index, lines(0) + columns(width()) * matrix->height());
       return iterator(this, index, columns.end(), columns.begin(), columns.end(), lines.begin(), lines.begin(), matrix);
     }

    /**
     * Subscript operator
     * @param i is the line selected
     * @param j is the column selected
     * @return the element in the matrix at the emplacement [i,j]
     */
    Data_Type& operator()(unsigned long i, unsigned long j)
    {
      assert(i < height());
      assert(j < width());
      return (*matrix)(lines(i), columns(j));
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
      return (*matrix)(lines(i), columns(j));
    }
    ///@}

  protected:
    /// The matrix used by the submatrix
    MatrixType* matrix;
    /// The columns used by this sub matrix
    ColumnVectorForLine columns;
    /// The lines used by this sub matrix
    ColumnVectorForColumn lines;
  };
}

#endif
