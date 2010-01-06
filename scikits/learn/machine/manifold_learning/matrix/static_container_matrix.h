/**
 * \file static_container_matrix.h>
 * Contains the general static container with dimensions fixed at compile-time
 */

#ifndef STATICCONTAINERMATRIX
#define STATICCONTAINERMATRIX

#include <boost/static_assert.hpp>
#include <boost/utility/enable_if.hpp>
#include <cassert>
#include <matrix/matrix_functions.h>
#include <matrix/container_iterator.h>

namespace Matrix
{
#if (__GNUC__ == 4)
# if (__GNUC_MINOR__ == 0)
#  define NOITERATORASSIGNEMENT
# endif
#endif
#if (__GNUC__ < 4)
# define NOITERATORASSIGNEMENT
#endif

  /// The general matrix container
  template<class DataType, unsigned long Height, unsigned long Width>
  struct MatrixContainer
  {
  public:
    /// Typedef of self
    typedef MatrixContainer<DataType, Height, Width> Self;
    /// Typedef of internal data
    typedef DataType Data_Type;
    /// Typedef of internal data
    typedef DataType* MatrixContainerType;
    /// Typedef for data traits
    typedef DataTypeTraits<DataType> DataTraits;

    /// Typedef of pointer to an element
    typedef DataType* pointer;
    /// Typedef of iterator to an element
    typedef LinearIterator<Self> iterator;
    /// Typedef of const iterator to an element
    typedef LinearIterator<const Self> const_iterator;

    /** @name Constructors and destructors */
    /// @{
    /**
     * The constructor of the matrix
     * @param height is the number of lines
     * @param width is the number of columns
     * @param init_val is the initial value that has to be put in the matrix
     */
    MatrixContainer(unsigned long height, unsigned long width, DataType init_val)
    {
      assert(height == Height);
      assert(width == Width);
      std::fill(begin(), end(), init_val);
    }

    /**
     * The constructor of the matrix with functor
     * @param height is the number of lines
     * @param width is the number of columns
     * @param functor is the functor class used to fill the matrix
     */
    template<class Functor>
    MatrixContainer(Functor& functor, unsigned long height, unsigned long width)
    {
      assert(height == Height);
      assert(width == Width);
      std::generate<iterator, Functor&>(begin(), end(), functor);
    }

    /**
     * Copy constructor
     * @param other is the copied matrix
     */
    MatrixContainer(const Self& other)
    {
      std::copy(other.begin(), other.end(), begin());
    }

    /**
     * Copy constructor
     * @param other is the copied matrix
     */
    template<unsigned long otherHeight, unsigned long otherWidth>
    MatrixContainer(const MatrixContainer<DataType, otherHeight, otherWidth>& other)
    {
      assert(other.height() == Height);
      assert(other.width() == Width);
      std::copy(other.begin(), other.end(), begin());
    }

    /**
     * Constructor from a table of elements
     * @param height is the height of the new matrix
     * @param width is the width of the new matrix
     * @param elements is a pointer to the table of DataType elements
     */
    MatrixContainer(unsigned long height, unsigned long width, const DataType* elements)
    {
      assert(height == Height);
      assert(width == Width);
      std::copy(elements, elements + height * width, begin());
    }

    /**
     * Constructor from a matrix
     * @param other is the matrix
     * @param flag should not be used
     */
    template<class MatrixType>
    explicit MatrixContainer(const MatrixType& other, typename boost::disable_if<has_const_iterator<MatrixType>, bool>::type flag = true)
    {
      assert(other.height() == Height);
      assert(other.width() == Width);
      for(unsigned long j = 0; j < width(); ++j)
      {
        for(unsigned long i = 0; i < height(); ++i)
        {
          (*this)(i,j) = static_cast<typename MatrixType::Data_Type>(other(i,j));
        }
      }
    }

    /**
     * Constructor from a matrix
     * @param other is the matrix
     * @param flag should not be used
     */
    template<class MatrixType>
    MatrixContainer(const MatrixType& other, typename boost::enable_if<has_const_iterator<MatrixType>, bool>::type flag = true)
    {
      assert(other.height() == Height);
      assert(other.width() == Width);
      std::copy(other.begin(), other.end(), begin());
    }

    /**
     * Bogus constructor
     */
    MatrixContainer(void)
    {
    }
    ///@}

    /** @name Other operations */
    /// @{
    /**
     * Returns the height of the matrix
     * @return the height of the matrix
     */
    unsigned long height() const
    {
      return Height;
    }
    /**
     * Returns the width of the matrix
     * @return the width of the matrix
     */
    unsigned long width() const
    {
      return Width;
    }
    /**
     * Returns the number of columns of the matrix
     * @return the number of columns of the matrix
     */
    unsigned long nbColumns() const
    {
      return width();
    }
    /**
     * Returns the number of rows of the matrix
     * @return the number of rows of the matrix
     */
    unsigned long nbRows() const
    {
      return height();
    }
    /**
     * Returns the number of lines of the matrix
     * @return the number of lines of the matrix
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
    void setHeight(unsigned long height)
    {
      assert(false);
    }

    /**
     * Sets the width of the matrix
     * Does not modify the inner structure of the matrix !
     * @param width is the new height
     */
    void setWidth(unsigned long width)
    {
      assert(false);
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
     * Subscript operator
     * @param i is the line selected
     * @param j is the column selected
     * @return the element in the matrix at the emplacement [i,j]
     */
    DataType& operator()(unsigned long i, unsigned long j)
    {
      assert(i < height());
      assert(j < width());
      return mMatrix[i + j * height()];
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
      return mMatrix[i + j * height()];
    }

    /**
     * Return an iterator to the beginning
     * @return an iterator to the beginning
     */
    iterator begin()
    {
      return iterator(this, mMatrix);
    }

    /**
     * Return an iterator to the beginning
     * @return an iterator to the beginning
     */
    const_iterator begin() const
    {
      return const_iterator(this, mMatrix);
    }

    /**
     * Return an iterator to the beginning
     * @return an iterator to the beginning
     */
    const_iterator cbegin() const
    {
      return const_iterator(this, mMatrix);
    }

    /**
     * Return an iterator to the end
     * @return an iterator to the end
     */
    const_iterator end() const
    {
      return const_iterator(this, mMatrix + height() * width());
    }

    /**
     * Return an iterator to the end
     * @return an iterator to the end
     */
    iterator end()
    {
      return iterator(this, mMatrix + height() * width());
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

    /**
     * Assignement operator : assigning a matrix
     * @param rhs is the copied matrix
     * @return the new matrix
     */
    Self& operator=(const Self& rhs)
    {
      std::copy(rhs.begin(), rhs.end(), begin());
      return (*this);
    }

    /**
     * Assignement operator : assigning a matrix
     * @param rhs is the copied matrix
     * @return the new matrix
     */
    template<unsigned long otherHeight, unsigned long otherWidth>
    Self& operator=(const MatrixContainer<DataType, otherHeight, otherWidth>& rhs)
    {
      assert(rhs.height() == Height);
      assert(rhs.width() == Width);
      std::copy(rhs.begin(), rhs.end(), begin());
      return (*this);
    }

    /**
     * Swaps the content of two matrix
     * @param rhs is the other matrix
     * @return itself
     */
    template<unsigned long otherHeight, unsigned otherWidth>
    Self& swap(MatrixContainer<DataType, otherHeight, otherWidth>& rhs)
    {
      assert(rhs.height() == Height);
      assert(rhs.width() == Width);
      std::swap_ranges(begin(), end(), rhs.begin());
      return *this;
    }

// Putain de GCC 4.0 de merde
#ifndef NOITERATORASSIGNEMENT
    /**
     * Assignmenent constructor
     * @param rhs is the copied matrix
     */
    template<class Assignable>
    typename boost::disable_if<has_const_iterator<Assignable>, Self&>::type operator=(const Assignable& rhs)
    {
      assert(rhs.height() == Height);
      assert(rhs.width() == Width);
      for(unsigned long i = 0; i < Width; ++i)
      {
        for(unsigned long j = 0; j < Height; ++j)
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
      assert(rhs.height() == Height);
      assert(rhs.width() == Width);
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
      assert(rhs.height() == Height);
      assert(rhs.width() == Width);
      for(unsigned long i = 0; i < Width; ++i)
      {
        for(unsigned long j = 0; j < Height; ++j)
        {
          (*this)(j, i) = rhs(j, i);
        }
      }
      return (*this);
    }
#endif
    ///@}

    /**
     * Returns the dataholder
     * @return the smart pointer containing the data
     */
    const DataType* get() const
    {
      return &mMatrix[0];
    }

    /**
     * Returns the dataholder
     * @return the smart pointer containing the data
     */
    DataType* get()
    {
      return &mMatrix[0];
    }

  protected:
    /// Data storage
    DataType mMatrix[Width*Height];
  };
}

#endif
