/**
 * \file sub_matrix_lib.h
 * This header defines the sub matrix
 * - strassen multiplications
 * - operations on matrix
 */

#ifndef SUB_MATRIX_LIB
#define SUB_MATRIX_LIB

#include <boost/scoped_array.hpp>
#include <matrix/matrix_lib.h>
#include <matrix/matrix_traits.h>
#include <matrix/sub_matrix_iterator.h>

namespace Matrix
{
  /// The sub matrix class
  template<class Matrix_Type>
  struct SubMatrix
  {
    typedef Matrix_Type MatrixType;
    /// Only to have a link on the type used in the matrix
    typedef typename MatrixType::Data_Type Data_Type;
    /// Typedef of self
    typedef SubMatrix Self;
    /// Typedef of result
    typedef Matrix<Data_Type, 0U, 0U> Result;
    /// Typedef of transposition
    typedef Matrix<Data_Type, 0U, 0U> Transpose;
    /// Typedef of comparison result
    typedef Matrix<bool, 0U, 0U> ComparisonResult;
    /// Typedef of column selection
    typedef typename ColumnVector<typename MatrixType::Data_Type, 0U>::Type Column;
    /// Typedef of line selection
    typedef typename LineVector<typename MatrixType::Data_Type, 0U>::Type Line;
    /// Typedef for <<
    typedef typename std::ostream ostream;
    /// Typedef for >>
    typedef typename std::istream istream;
    /// Typedef for data traits
    typedef DataTypeTraits<Data_Type> DataTraits;
    /// Typedef of iterator to an element
    typedef SubMatrixIterator<Self> iterator;
    /// Typedef of const iterator to an element
    typedef SubMatrixIterator<const Self> const_iterator;

    /// Static dimensions
    enum{staticHeight = 0U, staticWidth = 0U} ;

    /** @name Constructors and destructors (submatrix)*/
    /// @{
    /**
     * Base constructor : whole matrix
     * @param mat is the matrix used for creating the submatrix
     */
    SubMatrix(MatrixType& mat)
      :matrix(&mat), m_lineOffset(0), m_columnOffset(0), m_lineSize(mat.width()), m_columnSize(mat.height())
    {
    }

    /**
     * Base constructor : part of a matrix
     * @param mat is the matrix used for creating the submatrix
     * @param lineOffset is the line offset in the matrix
     * @param columnOffset is the column offset in the matrix
     * @param lineSize is the line size in the matrix
     * @param columnSize is the column size in the matrix
     */
    SubMatrix(MatrixType& mat, unsigned long lineOffset, unsigned long columnOffset, unsigned long lineSize, unsigned long columnSize)
      :matrix(&mat), m_lineOffset(lineOffset), m_columnOffset(columnOffset), m_lineSize(lineSize), m_columnSize(columnSize)
    {
      assert(lineSize <= mat.width());
      assert(columnSize <= mat.height());
    }

    /**
     * Base constructor : part of a submatrix
     * @param mat is the submatrix used for creating the submatrix
     * @param lineOffset is the line offset in the submatrix
     * @param columnOffset is the column offset in the submatrix
     * @param lineSize is the line size in the submatrix
     * @param columnSize is the column size in the submatrix
     */
    SubMatrix(const Self mat, unsigned long lineOffset, unsigned long columnOffset, unsigned long lineSize, unsigned long columnSize)
      :matrix(mat.matrix), m_lineOffset(lineOffset + mat.lineOffset()), m_columnOffset(columnOffset + mat.columnOffset()), m_lineSize(lineSize), m_columnSize(columnSize)
    {
      assert(lineSize <= mat.width());
      assert(columnSize <= mat.height());
    }

    /**
     * Bogus constructor, should not be used
     */
    SubMatrix()
      :matrix(NULL), m_lineOffset(0), m_columnOffset(0), m_lineSize(0), m_columnSize(0)
    {
    }
    ///@}

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
      return m_columnSize;
    }
    /**
     * Returns the width of the matrix
     * @return the width of the matrix
     */
    unsigned long width() const
    {
      return m_lineSize;
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
     * Sets the height of the matrix
     * Does not modify the inner structure of the matrix !
     * @param height is the new height
     */
    void setRealHeight(unsigned long height)
    {
      matrix->setHeight(height);
    }

    /**
     * Sets the width of the matrix
     * Does not modify the inner structure of the matrix !
     * @param width is the new width
     */
    void setRealWidth(unsigned long width)
    {
      matrix->setWidth(width);
    }
    /**
     * Sets the height of the matrix
     * Does not modify the inner structure of the matrix !
     * @param columnSize is the new height
     */
    void setHeight(unsigned long columnSize)
    {
      m_columnSize = columnSize;
    }

    /**
     * Sets the width of the matrix
     * Does not modify the inner structure of the matrix !
     * @param lineSize is the new width
     */
    void setWidth(unsigned long lineSize)
    {
      m_lineSize = lineSize;
    }

    /**
     * Sets the number of lines of the matrix
     * @param height is the new height
     */
    void setNbLines(unsigned long height)
    {
      setHeight(height);
    }

    /**
     * Sets the number of columns of the matrix
     * @param width is the new height
     */
    void setNbColumns(unsigned long width)
    {
      setWidth(width);
    }

    /**
     * Gets the line size
     * @return the line size
     */
    unsigned long lineSize() const
    {
      return width();
    }

    /**
     * Gets the column size
     * @return the column size
     */
    unsigned long columnSize() const
    {
      return height();
    }

    /**
     * Gets the line offset
     * @return the line offset
     */
    unsigned long lineOffset() const
    {
      return m_lineOffset;
    }

    /**
     * Gets the column offset
     * @return the column offset
     */
    unsigned long columnOffset() const
    {
      return m_columnOffset;
    }

    /**
     * Sets the line size
     * @param lineSize is the new line size
     */
    void setLineSize(unsigned long lineSize)
    {
      setWidth(lineSize);
    }

    /**
     * Sets the column offsize
     * @param columnSize is the new column size
     */
    void setColumnSize(unsigned long columnSize)
    {
      setHeight(columnSize);
    }

    /**
     * Sets the line offset
     * @param lineOffset is the new line offset
     */
    void setLineOffset(unsigned long lineOffset)
    {
      m_lineOffset = lineOffset;
    }

    /**
     * Sets the column offset
     * @param columnOffset is the new column offset
     */
    void setColumnOffset(unsigned long columnOffset)
    {
      m_columnOffset = columnOffset;
    }

    /**
     * Assignement operator
     * @param mat is the submatrix that will be copied
     * @return the modified submatrix
     */
    template<class OtherMatrixType>
    Self& assign(const OtherMatrixType mat)
    {
      setWidth(mat.width());
      setHeight(mat.height());
      setLineOffset(mat.lineOffset());
      setColumnOffset(mat.columnOffset());
      matrix = mat.matrix;
      return *this;
    }

    /**
     * Return an iterator to the beginning
     * @return an iterator to the beginning
     */
    iterator begin()
    {
      typename MatrixType::iterator index = matrix->begin();
      std::advance(index, m_columnOffset + m_lineOffset * matrix->height());
      return iterator(this, index);
    }

    /**
     * Return an iterator to the beginning
     * @return an iterator to the beginning
     */
    const_iterator begin() const
    {
      typename MatrixType::const_iterator index = matrix->cbegin();
      std::advance(index, m_columnOffset + m_lineOffset * matrix->height());
      return const_iterator(this, index);
    }

    /**
     * Return an iterator to the end
     * @return an iterator to the end
     */
    const_iterator end() const
    {
      typename MatrixType::const_iterator index = matrix->cbegin();
      std::advance(index, m_columnOffset + (m_lineOffset + m_lineSize) * matrix->height());
      return const_iterator(this, index);
    }

    /**
     * Return an iterator to the end
     * @return an iterator to the end
     */
    iterator end()
    {
      typename MatrixType::iterator index = matrix->begin();
      std::advance(index, m_columnOffset + (m_lineOffset + m_lineSize) * matrix->height());
      return iterator(this, index);
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
      assert((i + columnOffset()) < realHeight());
      assert((j + lineOffset()) < realWidth());
      return (*matrix)((columnOffset() + i), (lineOffset() + j));
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
      if(((i + columnOffset()) >= realHeight())||((j + lineOffset()) >= realWidth()))
      {
        return DataTraits::zero((*matrix)(0, 0));
      }
      else
      {
        return (*matrix)((columnOffset() + i), (lineOffset() + j));
      }
    }

    /**
     * This function makes zero from elements that are small enough
     */
    void chopSmallElements()
    {

      for(unsigned long j = 0; j < width(); ++j)
      {
        for(unsigned long i = 0; i < height(); ++i)
        {
          Data_Type epsilon;
          if(i != j)
          {
            epsilon = (DataTraits::epsilon((*this)(i, j)) * (DataTraits::absolute((*this)(i, i)) + DataTraits::absolute((*this)(j, j)))) * width() * height();
            if(epsilon == DataTraits::zero((*this)(i, j)))
            {
              epsilon = DataTraits::epsilon((*this)(i, j)) * width() * height();
            }
            if(((*this)(i, j) <= epsilon) && (-epsilon <= (*this)(i, j)))
              (*this)(i,j) = DataTraits::zero((*this)(i, j));
          }
        }
      }
      for(unsigned long i = 0; i < width(); ++i)
      {
        Data_Type epsilon;
        epsilon = DataTraits::epsilon((*this)(i, i)) * (DataTraits::absolute((*this)(0, 0)) + DataTraits::absolute((*this)(height() - 1, width() - 1)));
        if(((*this)(i, i) <= epsilon) && (-epsilon <= (*this)(i, i)))
          (*this)(i,i) = DataTraits::zero((*this)(i, i));
      }
    }
    ///@}

    /** @name Shift (submatrix)*/
    /// @{
    /**
     * <<= operator : left shift by a constant
     * @param rhs is the shift asked
     * @return the modified submatrix
     */
    Self& operator<<=(const Data_Type rhs)
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
     * @return the modified submatrix
     */
    template<class OtherMatrixType>
    Self& operator<<=(const OtherMatrixType& rhs)
    {
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
     * @return the modified matsubmatrix
     */
    Self& operator>>=(const Data_Type rhs)
    {
      for(unsigned long j = 0; j < (width() - lineOffset()); ++j)
      {
        for(unsigned long i = 0; i < (height() - columnOffset()); ++i)
        {
          (*this)(i,j) >>= rhs;
        }
      }
      return (*this);
    }

    /**
     * >>= operator : right shift by a matrix
     * @param rhs is the matrix shift asked
     * @return the modified submatrix
     */
    template<class OtherMatrixType>
    Self& operator>>=(const OtherMatrixType& rhs)
    {
      assert(height() == rhs.height());
      assert(width() == rhs.width());
      for(unsigned long j = 0; j < (rhs.width() - rhs.lineOffset()); ++j)
      {
        for(unsigned long i = 0; i < (height() - columnOffset()); ++i)
        {
          (*this)(i,j) >>= rhs(i, j);
        }
      }
      return (*this);
    }
    ///@}

    /** @name Usual matrix operations (submatrix)*/
    /// @{
    /**
     * Tests if the submatrix is square or not
     * @return true if the submatrix is square
     */
    bool isSquare() const
    {
      return (columnSize() == lineSize());
    }
    ///@}
    /** @}
     */

  protected:
    MatrixType* matrix;
    /// Offset in the line of the parent Matrix
    unsigned long m_lineOffset;
    /// Offset in the column of the parent Matrix
    unsigned long m_columnOffset;
    /// Line size of the sub_matrix
    unsigned long m_lineSize;
    /// Column size of the sub_matrix
    unsigned long m_columnSize;
  };

  typedef SubMatrix<int> SubMatrixi;
  typedef SubMatrix<short> SubMatrixs;
  typedef SubMatrix<long> SubMatrixl;
  typedef SubMatrix<float> SubMatrixf;
  typedef SubMatrix<double> SubMatrixd;
}

#endif 
