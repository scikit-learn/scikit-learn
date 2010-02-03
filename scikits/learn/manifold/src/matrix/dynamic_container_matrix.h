/**
 * \file dynamic_container_matrix.h>
 * Contains the general dynamic container with dimensions not fixed at compile-time
 */

#ifndef DYNAMICCONTAINERMATRIX
#define DYNAMICCONTAINERMATRIX

#include <boost/scoped_array.hpp>
#include <boost/static_assert.hpp>
#include <cassert>
#include <matrix/container_iterator.h>

namespace Matrix
{
  /// The specialized matrix container for variable width
  template<class DataType, unsigned long Height>
  struct MatrixContainer<DataType, Height, 0U>
  {
  public:
    /// Typedef of self
    typedef MatrixContainer<DataType, Height, 0U> Self;
    /// Typedef of internal data
    typedef DataType Data_Type;
    /// Typedef of internal data
    typedef boost::scoped_array<DataType> MatrixContainerType;
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
      :m_width(width), mMatrix(new DataType[Height * m_width])
    {
      assert(height == Height);
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
      :m_width(width), mMatrix(new DataType[Height * m_width])
    {
      assert(height == Height);
      std::generate<iterator, Functor&>(begin(), end(), functor);
    }

    /**
     * Copy constructor
     * @param other is the copied matrix
     */
    MatrixContainer(const Self& other)
      :m_width(other.width()), mMatrix(new DataType[Height * width()])
    {
      std::copy(other.begin(), other.end(), begin());
    }

    /**
     * Copy constructor
     * @param other is the copied matrix
     */
    template<unsigned long Width>
    MatrixContainer(const MatrixContainer<DataType, 0U, Width>& other)
      :m_width(other.width()), mMatrix(new DataType[Height * width()])
    {
      std::copy(other.begin(), other.end(), begin());
    }

    /**
     * Copy constructor
     * @param other is the copied matrix
     */
    template<unsigned long Width>
    MatrixContainer(const MatrixContainer<DataType, Height, Width>& other)
      :m_width(other.width()), mMatrix(new DataType[Height * width()])
    {
      std::copy(other.begin(), other.end(), begin());
    }

    /**
     * Constructor from a table of elements
     * @param height is the height of the new matrix
     * @param width is the width of the new matrix
     * @param elements is a pointer to the table of DataType elements
     */
    MatrixContainer(unsigned long height, unsigned long width, const DataType* elements)
      :m_width(width), mMatrix(new DataType[Height * m_width])
    {
      std::copy(elements, elements + this->height() * this->width(), begin());
    }

    /**
     * Constructor from a matrix
     * @param other is the matrix
     * @param flag should not be used
     */
    template<class MatrixType>
    explicit MatrixContainer(const MatrixType& other, typename boost::disable_if<has_const_iterator<MatrixType>, bool>::type flag = true)
      :m_width(other.width()), mMatrix(new DataType[Height * m_width])
    {
      assert(other.height() == Height);
      for(unsigned long j = 0; j < width(); ++j)
      {
        for(unsigned long i = 0; i < height(); ++i)
        {
          (*this)(i, j) = static_cast<typename MatrixType::Data_Type>(other(i,j));
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
      :m_width(other.width()), mMatrix(new DataType[Height * m_width])
    {
      assert(other.height() == Height);
      std::copy(other.begin(), other.end(), begin());
    }

    /**
     * Bogus constructor
     */
    MatrixContainer(void)
      :m_width(0)
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
      return m_width;
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
     * @param width is the new width
     * \post width() = width
     */
    void setWidth(unsigned long width)
    {
      m_width = width;
      assert(width == this->width());
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
      return iterator(this, mMatrix.get());
    }

    /**
     * Return an iterator to the beginning
     * @return an iterator to the beginning
     */
    const_iterator begin() const
    {
      return const_iterator(this, mMatrix.get());
    }

    /**
     * Return an iterator to the beginning
     * @return an iterator to the beginning
     */
    const_iterator cbegin() const
    {
      return const_iterator(this, mMatrix.get());
    }

    /**
     * Return an iterator to the end
     * @return an iterator to the end
     */
    const_iterator end() const
    {
      return const_iterator(this, mMatrix.get() + height() * width());
    }

    /**
     * Return an iterator to the end
     * @return an iterator to the end
     */
    iterator end()
    {
      return iterator(this, mMatrix.get() + height() * width());
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
     * Assignmenent constructor
     * @param rhs is the copied matrix
     */
    template<unsigned long Width>
    Self& operator=(const MatrixContainer<DataType, 0U, Width>& rhs)
    {
      assert(rhs.height() == Height);
      setWidth(Width);
      resize();
      std::copy(rhs.begin(), rhs.end(), begin());
      return (*this);
    }

    /**
     * Assignmenent constructor
     * @param rhs is the copied matrix
     */
    template<unsigned long Width>
    Self& operator=(const MatrixContainer<DataType, Height, Width>& rhs)
    {
      assert(rhs.height() == Height);
      setWidth(Width);
      resize();
      std::copy(rhs.begin(), rhs.end(), begin());
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
      assert(rhs.height() == Height);
      setWidth(rhs.width());
      resize();
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
      assert(rhs.height() == Height);
      setWidth(rhs.width());
      resize();
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
      setWidth(rhs.width());
      resize();
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

    /**
     * Assignement operator : assigning a matrix
     * @param rhs is the copied matrix
     * @return the new matrix
     */
    Self& operator=(const Self& rhs)
    {
      assert(rhs.height() == Height);
      setWidth(rhs.width());
      resize();
      return (*this);
    }

    /**
     * Swaps the content of two matrix
     * @param rhs is the other matrix
     * @return itself
     */
    template<unsigned long otherHeight, unsigned long otherWidth>
    Self& swap(MatrixContainer<DataType, otherHeight, otherWidth>& rhs)
    {
      assert(rhs.height() == height());
      assert(rhs.width() == width());
      std::swap_ranges(begin(), end(), rhs.begin());
      return *this;
    }
    ///@}

    /**
     * Returns the dataholder
     * @return the smart pointer containing the data
     */
    const DataType* get() const
    {
      return mMatrix.get();
    }

    /**
     * Returns the dataholder
     * @return the smart pointer containing the data
     */
    DataType* get()
    {
      return mMatrix.get();
    }

  private:
    /// Initializes the indirection vector
    void resize()
    {
      mMatrix.reset(new DataType[height() * width()]);
    }

    /// Number of columns of the matrix
    unsigned long m_width;
  protected:
    /// Data storage
    boost::scoped_array<DataType> mMatrix;
  };

  /// The specialized matrix container for variable height
  template<class DataType, unsigned long Width>
  struct MatrixContainer<DataType, 0U, Width>
  {
  public:
    /// Typedef of self
    typedef MatrixContainer<DataType, 0U, Width> Self;
    /// Typedef of internal data
    typedef DataType Data_Type;
    /// Typedef of internal data
    typedef boost::scoped_array<DataType> MatrixContainerType;
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
      :m_height(height), mMatrix(new DataType[m_height * Width])
    {
      assert(width == Width);
      resize();
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
      :m_height(height), mMatrix(new DataType[m_height * Width])
    {
      assert(width == Width);
      resize();
      std::generate<iterator, Functor&>(begin(), end(), functor);
    }

    /**
     * Copy constructor
     * @param other is the copied matrix
     */
    MatrixContainer(const Self& other)
      :m_height(other.height()), mMatrix(new DataType[m_height * Width])
    {
      resize();
      std::copy(other.begin(), other.end(), begin());
    }

    /**
     * Copy constructor
     * @param other is the copied matrix
     */
    template<unsigned long Height>
    MatrixContainer(const MatrixContainer<DataType, Height, 0U>& other)
      :m_height(other.height()), mMatrix(new DataType[m_height * Width])
    {
      resize();
      std::copy(other.begin(), other.end(), begin());
    }

    /**
     * Copy constructor
     * @param other is the copied matrix
     */
    template<unsigned long Height>
    MatrixContainer(const MatrixContainer<DataType, Height, Width>& other)
      :m_height(other.height()), mMatrix(new DataType[height() * Width])
    {
      resize();
      std::copy(other.begin(), other.end(), begin());
    }

    /**
     * Constructor from a table of elements
     * @param height is the height of the new matrix
     * @param width is the width of the new matrix
     * @param elements is a pointer to the table of DataType elements
     */
    MatrixContainer(unsigned long height, unsigned long width, const DataType* elements)
      :m_height(height), mMatrix(new DataType[m_height * Width])
    {
      assert(width == Width);
      resize();
      std::copy(elements, elements + height * width, begin());
    }

    /**
     * Constructor from a matrix
     * @param other is the matrix
     * @param flag should not be used
     */
    template<class MatrixType>
    explicit MatrixContainer(const MatrixType& other, typename boost::disable_if<has_const_iterator<MatrixType>, bool>::type flag = true)
      :m_height(other.height()), mMatrix(new DataType[m_height * Width])
    {
      assert(other.width() == Width);
      resize();
      for(unsigned long j = 0; j < width(); ++j)
      {
        for(unsigned long i = 0; i < height(); ++i)
        {
          (*this)(i, j) = static_cast<typename MatrixType::Data_Type>(other(i,j));
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
      :m_height(other.height()), mMatrix(new DataType[m_height * Width])
    {
      assert(other.width() == Width);
      resize();
      std::copy(other.begin(), other.end(), begin());
    }

    /**
     * Bogus constructor
     */
    MatrixContainer(void)
      :m_height(0)
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
      return m_height;
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
     * \post height() = height
     */
    void setHeight(unsigned long height)
    {
      m_height = height;
      assert(height == this->height());
    }

    /**
     * Sets the width of the matrix
     * Does not modify the inner structure of the matrix !
     * @param width is the new width
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
      return iterator(this, mMatrix.get());
    }

    /**
     * Return an iterator to the beginning
     * @return an iterator to the beginning
     */
    const_iterator begin() const
    {
      return const_iterator(this, mMatrix.get());
    }

    /**
     * Return an iterator to the beginning
     * @return an iterator to the beginning
     */
    const_iterator cbegin() const
    {
      return const_iterator(this, mMatrix.get());
    }

    /**
     * Return an iterator to the end
     * @return an iterator to the end
     */
    const_iterator end() const
    {
      return const_iterator(this, mMatrix.get() + height() * width());
    }

    /**
     * Return an iterator to the end
     * @return an iterator to the end
     */
    iterator end()
    {
      return iterator(this, mMatrix.get() + height() * width());
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
     * Assignmenent constructor
     * @param rhs is the copied matrix
     */
    template<unsigned long Height>
    Self& operator=(const MatrixContainer<DataType, Height, 0U>& rhs)
    {
      assert(rhs.width() == Width);
      setHeight(Height);
      resize();
      std::copy(rhs.begin(), rhs.end(), begin());
      return (*this);
    }

    /**
     * Assignmenent constructor
     * @param rhs is the copied matrix
     */
    template<unsigned long Height>
    Self& operator=(const MatrixContainer<DataType, Height, Width>& rhs)
    {
      assert(rhs.width() == Width);
      setHeight(Height);
      resize();
      std::copy(rhs.begin(), rhs.end(), begin());
      return (*this);
    }

    /**
     * Assignement operator : assigning a matrix
     * @param rhs is the copied matrix
     * @return the new matrix
     */
    Self& operator=(const Self& rhs)
    {
      assert(rhs.width() == Width);
      setHeight(rhs.height());
      resize();
      std::copy(rhs.begin(), rhs.end(), begin());
      return (*this);
    }

    /**
     * Swaps the content of two matrix
     * @param rhs is the other matrix
     * @return itself
     */
    template<unsigned long otherHeight, unsigned long otherWidth>
    Self& swap(MatrixContainer<DataType, otherHeight, otherWidth>& rhs)
    {
      assert(rhs.height() == height());
      assert(rhs.width() == width());
      std::swap_ranges(begin(), end(), rhs.begin());
      return *this;
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
      assert(rhs.width() == Width);
      setHeight(rhs.height());
      resize();
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
      assert(rhs.width() == Width);
      setHeight(rhs.height());
      resize();
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
      assert(rhs.width() == Width);
      setHeight(rhs.height());
      resize();
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

    /**
     * Returns the dataholder
     * @return the smart pointer containing the data
     */
    const DataType* get() const
    {
      return mMatrix.get();
    }

    /**
     * Returns the dataholder
     * @return the smart pointer containing the data
     */
    DataType* get()
    {
      return mMatrix.get();
    }

  private:
    /// Initializes the indirection vector
    void resize()
    {
      mMatrix.reset(new DataType[height() * width()]);
    }

    /// Number of lines of the matrix
    unsigned long m_height;
  protected:
    /// Data storage
    boost::scoped_array<DataType> mMatrix;
  };

  /// The specialized variable matrix container
  template<class DataType>
  struct MatrixContainer<DataType, 0U, 0U>
  {
  public:
    /// Typedef of self
    typedef MatrixContainer<DataType, 0U, 0U> Self;
    /// Typedef of internal data
    typedef DataType Data_Type;
    /// Typedef of internal data
    typedef boost::scoped_array<DataType> MatrixContainerType;
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
      :m_height(height), m_width(width), mMatrix(new DataType[m_height * m_width])
    {
      resize();
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
      :m_height(height), m_width(width), mMatrix(new DataType[m_height * m_width])
    {
      resize();
      std::generate<iterator, Functor&>(begin(), end(), functor);
    }

    /**
     * Copy constructor
     * @param other is the copied matrix
     */
    MatrixContainer(const Self& other)
      :m_height(other.height()), m_width(other.width()), mMatrix(new DataType[m_height * m_width])
    {
      resize();
      std::copy(other.begin(), other.end(), begin());
    }
    /**
     * Copy constructor
     * @param other is the copied matrix
     */
    template<unsigned long Height, unsigned long Width>
    MatrixContainer(const MatrixContainer<DataType, Height, Width>& other)
      :m_height(other.height()), m_width(other.width()), mMatrix(new DataType[m_height * m_width])
    {
      resize();
      std::copy(other.begin(), other.end(), begin());
    }

    /**
     * Constructor from a table of elements
     * @param height is the height of the new matrix
     * @param width is the width of the new matrix
     * @param elements is a pointer to the table of DataType elements
     */
    MatrixContainer(unsigned long height, unsigned long width, const DataType* elements)
      :m_height(height), m_width(width), mMatrix(new DataType[m_height * m_width])
    {
      resize();
      std::copy(elements, elements + this->height() * this->width(), begin());
    }

    /**
     * Constructor from a matrix
     * @param other is the matrix
     * @param flag should not be used
     */
    template<class MatrixType>
    explicit MatrixContainer(const MatrixType& other, typename boost::disable_if<has_const_iterator<MatrixType>, bool>::type flag = true)
      :m_height(other.height()), m_width(other.width()), mMatrix(new DataType[m_height * m_width])
    {
      resize();
      for(unsigned long j = 0; j < width(); ++j)
      {
        for(unsigned long i = 0; i < height(); ++i)
        {
          (*this)(i, j) = static_cast<typename MatrixType::Data_Type>(other(i,j));
        }
      }
    }

    /**
     * Constructor from a matrix
     * @param other is the matrix
     * @param flag should not be used
     */
    template<class MatrixType>
    explicit MatrixContainer(const MatrixType& other, typename boost::enable_if<has_const_iterator<MatrixType>, bool>::type flag = false)
      :m_height(other.height()), m_width(other.width()), mMatrix(new DataType[m_height * m_width])
    {
      resize();
      std::copy(other.begin(), other.end(), begin());
    }

    /**
     * Bogus constructor
     */
    MatrixContainer(void)
      :m_height(0), m_width(0)
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
     * \post height() = height
     */
    void setHeight(unsigned long height)
    {
      m_height = height;
      assert(height == this->height());
    }

    /**
     * Sets the width of the matrix
     * Does not modify the inner structure of the matrix !
     * @param width is the new width
     * \post width() = width
     */
    void setWidth(unsigned long width)
    {
      m_width = width;
      assert(width == this->width());
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
      return iterator(this, mMatrix.get());
    }

    /**
     * Return an iterator to the beginning
     * @return an iterator to the beginning
     */
    const_iterator begin() const
    {
      return const_iterator(this, mMatrix.get());
    }

    /**
     * Return an iterator to the beginning
     * @return an iterator to the beginning
     */
    const_iterator cbegin() const
    {
      return const_iterator(this, mMatrix.get());
    }

    /**
     * Return an iterator to the end
     * @return an iterator to the end
     */
    const_iterator end() const
    {
      return const_iterator(this, mMatrix.get() + height() * width());
    }

    /**
     * Return an iterator to the end
     * @return an iterator to the end
     */
    iterator end()
    {
      return iterator(this, mMatrix.get() + height() * width());
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
     * Assignement constructor
     * @param rhs is the copied matrix
     */
     Self& operator=(const Self& rhs)
    {
      setHeight(rhs.height());
      setWidth(rhs.width());
      resize();
      std::copy(rhs.begin(), rhs.end(), begin());
      return (*this);
    }

    /**
     * Assignement constructor
     * @param rhs is the copied matrix
     */
    template<unsigned long Height, unsigned long Width>
     Self& operator=(const MatrixContainer<DataType, Height, Width>& rhs)
    {
      setHeight(rhs.height());
      setWidth(rhs.width());
      resize();
      std::copy(rhs.begin(), rhs.end(), begin());
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
      setHeight(rhs.height());
      setWidth(rhs.width());
      resize();
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
      setHeight(rhs.height());
      setWidth(rhs.width());
      resize();
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
      setHeight(rhs.height());
      setWidth(rhs.width());
      resize();
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

    /**
     * Swaps the content of two matrix
     * @param rhs is the other matrix
     * @return itself
     */
    template<unsigned long otherHeight, unsigned long otherWidth>
    Self& swap(MatrixContainer<DataType, otherHeight, otherWidth>& rhs)
    {
      assert(rhs.height() == height());
      assert(rhs.width() == width());
      std::swap_ranges(begin(), end(), rhs.begin());
      return *this;
    }
    ///@}

    /**
     * Returns the dataholder
     * @return the smart pointer containing the data
     */
    const DataType* get() const
    {
      return mMatrix.get();
    }

    /**
     * Returns the dataholder
     * @return the smart pointer containing the data
     */
    DataType* get()
    {
      return mMatrix.get();
    }

  private:
    /// Initializes the indirection vector
    void resize()
    {
      mMatrix.reset(new DataType[height() * width()]);
    }

    /// Number of lines of the matrix
    unsigned long m_height;
    /// Number of columns of the matrix
    unsigned long m_width;
  protected:
    /// Data storage
    boost::scoped_array<DataType> mMatrix;
  };
}

#endif
