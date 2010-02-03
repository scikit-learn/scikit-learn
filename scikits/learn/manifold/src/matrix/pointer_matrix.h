/**
 * \file pointer_matrix.h
 * Defines a matrix with the correct interface that takes a pointer as parameter
 */

#ifndef POINTERMATRIX
#define POINTERMATRIX

#include <matrix/matrix_traits.h>
#include <matrix/matrix_functions.h>
#include <matrix/container_iterator.h>

namespace Matrix
{
  /// A class that wrapps a pointer to an array
  template<class DataType>
  struct PointerMatrix
  {
    /// Only to have a link on the type used in the matrix
    typedef DataType Data_Type;
    /// Typedef of self
    typedef PointerMatrix Self;
    /// Typedef of operation result
    typedef Matrix<DataType, 0, 0> Result;
    /// Typedef of comparison result
    typedef Matrix<bool, 0, 0> ComparisonResult;
    /// Typedef of transposition
    typedef Matrix<DataType, 0, 0> Transpose;
    /// Typedef of pointer to an element
    typedef DataType* pointer;
    /// Typedef of iterator to an element
    typedef LinearIterator<Self> iterator;
    /// Typedef of const iterator to an element
    typedef LinearIterator<const Self> const_iterator;
    /// Typedef of column selection
    typedef typename ColumnVector<DataType, 0>::Type Column;
    /// Typedef of line selection
    typedef typename LineVector<DataType, 0>::Type Line;
    /// Typedef for <<
    typedef typename std::ostream ostream;
    /// Typedef for >>
    typedef typename std::istream istream;

    /// Typedef for data traits
    typedef DataTypeTraits<DataType> DataTraits;
    /// Typedef for own traits
    typedef DataTypeTraits<Self> OwnTraits;

    /// Static dimensions
    enum{staticHeight = 0, staticWidth = 0} ;

    /**
     * Constructor. Takes an array and uses it (warning, does not delete elements when destroyed)
     * @param height is the height of the new matrix
     * @param width is the width of the new matrix
     * @param elements is a pointer to the table of DataType elements
     */
    PointerMatrix(unsigned long height, unsigned long width, DataType* elements)
      :m_height(height), m_width(width), elements(elements)
    {
    }

    /**
     * Copy constructor
     * @param other is the matrix to be copied
     */
    PointerMatrix(const PointerMatrix& other)
      :m_height(other.m_height), m_width(other.m_width), elements(other.elements)
    {
    }

    /** @name Other operations */
    /// @{
    /**
     * Returns the height of the matrix
     * @return the height of the matrix
     */
    unsigned long height() const
    {
      return m_height;
    }
    /**
     * Returns the width of the matrix
     * @return the width of the matrix
     */
    unsigned long width() const
    {
      return m_width;
    }

    /**
     * Subscript operator
     * @param i is the line selected
     * @param j is the column selected
     * @return the element in the matrix at the emplacement [i,j]
     */
    DataType& operator()(unsigned long i, unsigned long j)
    {
      assert(i < height());
      assert(j < width());
      return elements[i + j * height()];
    }

    /**
     * Const subscript operator
     * @param i is the line selected
     * @param j is the column selected
     * @return the element in the matrix at the emplacement [i,j]
     */
    DataType operator()(unsigned long i, unsigned long j) const
    {
      assert(i < height());
      assert(j < width());
      return elements[i + j * height()];
    }

    /**
     * Return an iterator to the beginning
     * @return an iterator to the beginning
     */
    iterator begin()
    {
      return iterator(this, elements);
    }

    /**
     * Return an iterator to the beginning
     * @return an iterator to the beginning
     */
    const_iterator begin() const
    {
      return const_iterator(this, elements);
    }

    /**
     * Return an iterator to the beginning
     * @return an iterator to the beginning
     */
    const_iterator cbegin() const
    {
      return const_iterator(this, elements);
    }

    /**
     * Return an iterator to the end
     * @return an iterator to the end
     */
    const_iterator end() const
    {
      return const_iterator(this, elements + height() * width());
    }

    /**
     * Return an iterator to the end
     * @return an iterator to the end
     */
    iterator end()
    {
      return iterator(this, elements + height() * width());
    }

    /**
     * Assignmenent operator : assigning a new value
     * @param rhs is the new value in the matrix
     * @return the modified matrix
     */
    Self& operator=(const DataType rhs)
    {
      std::fill(begin(), end(), rhs);
      return (*this);
    }

// Putain de connard de GCC 4.0 de merde
#ifndef NOITERATORASSIGNEMENT
    /**
     * Assignmenent constructor
     * @param rhs is the copied matrix
     */
    template<class Assignable>
    typename boost::disable_if<has_const_iterator<Assignable>, Self&>::type operator=(const Assignable& rhs)
    {
      assert(rhs.height() == height());
      assert(rhs.width() == width());
      for(unsigned long i = 0; i < rhs.width(); ++i)
      {
        for(unsigned long j = 0; j < rhs.height(); ++j)
        {
          (*this)(j, i) = rhs(j, i);
        }
      }
      return (*this);
    }

    /**
     * Assignmenent constructor
     * @param rhs is the copied matrix
     */
    template<class Assignable>
    typename boost::enable_if<has_const_iterator<Assignable>, Self&>::type operator=(const Assignable& rhs)
    {
      assert(rhs.height() == height());
      assert(rhs.width() == width());
      std::copy(rhs.begin(), rhs.end(), begin());
      return (*this);
    }
#else
    /**
     * Assignmenent constructor
     * @param rhs is the copied matrix
     */
    template<class Assignable>
    Self& operator=(const Assignable& rhs)
    {
      assert(rhs.height() == height());
      assert(rhs.width() == width());
      for(unsigned long i = 0; i < rhs.width(); ++i)
      {
        for(unsigned long j = 0; j < rhs.height(); ++j)
        {
          (*this)(j, i) = rhs(j, i);
        }
      }
      return (*this);
    }
#endif
    ///@}

  private:
    /// Height of the matrix
    unsigned long m_height;
    /// Width of the matrix
    unsigned long m_width;
    /// Array of elements in the matrix
    DataType* elements;
  };
}

#endif
