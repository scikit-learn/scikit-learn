/**
 * \file sub_vector_matrix_iterator.h
 * Contains the iterator for a sub vector matrix
 */

#ifndef SUBVECTORMATRIXITERATOR
#define SUBVECTORMATRIXITERATOR

#include <iterator>
#include <boost/mpl/if.hpp>
#include <boost/type_traits.hpp>

namespace Matrix
{
  /// The standard iterator for sub line matrix
  template<class Container>
  class SubLineMatrixIterator : public boost::mpl::if_<typename boost::is_const<typename boost::remove_reference<Container>::type>::type, std::iterator<std::forward_iterator_tag, typename boost::add_const<typename Container::Data_Type>::type, ptrdiff_t, typename boost::add_const<typename Container::Data_Type>::type*, typename boost::add_const<typename Container::Data_Type>::type&>, std::iterator<std::forward_iterator_tag, typename Container::Data_Type, ptrdiff_t, typename Container::Data_Type*, typename Container::Data_Type&> >::type
  {
    typedef SubLineMatrixIterator Self;
    /// Type of the parent
    typedef typename boost::mpl::if_<typename boost::is_const<typename boost::remove_reference<Container>::type>::type, std::iterator<std::forward_iterator_tag, typename boost::add_const<typename Container::Data_Type>::type, ptrdiff_t, typename boost::add_const<typename Container::Data_Type>::type*, typename boost::add_const<typename Container::Data_Type>::type&>, std::iterator<std::forward_iterator_tag, typename Container::Data_Type, ptrdiff_t, typename Container::Data_Type*, typename Container::Data_Type&> >::type Parent;

    typedef SubLineMatrixIterator<typename boost::add_const<Container>::type> ConstSelf;
    typedef typename boost::mpl::if_<typename boost::is_const<typename boost::remove_reference<Container>::type>::type, typename Container::MatrixType::const_iterator, typename Container::MatrixType::iterator>::type ContainerIterator;

    typedef typename boost::mpl::if_<typename boost::is_const<typename boost::remove_reference<Container>::type>::type, typename Container::IterateLineVector::const_iterator, typename Container::IterateLineVector::iterator>::type LineIterator;

    typedef typename boost::add_const<typename Container::MatrixType>::type ConstMatrixType;
    /// Pointer to the container
    Container* containerPointer;
    /// Iterator to the parent
    ContainerIterator containerIterator;
    /// Iterator pointing to a specific line
    LineIterator line;
    /// Iterator to the first line
    LineIterator lineBegin_;
    /// Pointer to the inner container
    ConstMatrixType* matrix;
    /// Position in a column
    unsigned long position;
    /// Stride for operations
    unsigned long stride;
  public:

    /// Returns the position inside the container - not enough to check equality !
    ContainerIterator getPosition() const
    {
      return containerIterator;
    }

    /**
     * Constructor from a container and a position in the container
     * @param container is the container on which to iterate
     * @param containerIterator is a iterator in the container
     * @param line is the line that is currently used
     * @param lineBegin is the start of the line
     * @param position is the position in the line
     * @param stride is the step used in (in/de)crementation
     */
    SubLineMatrixIterator(Container* container, const ContainerIterator& containerIterator, const LineIterator& line, const LineIterator& lineBegin_, ConstMatrixType* matrix, unsigned long position = 0, unsigned long stride = 1U)
      :containerPointer(container), containerIterator(containerIterator), line(line), lineBegin_(lineBegin_), matrix(matrix), position(position), stride(stride)
    {
    }

    /// Preincrementation
    Self& operator++()
    {
      position += stride;
      if(position >= containerPointer->height())
      {
        position = position % containerPointer->height();
        LineIterator lastLine = line;
        line = lineBegin_;
        std::advance(line, position);
        std::advance(containerIterator, *line - *lastLine + matrix->height());
      }
      else
      {
        LineIterator lastLine = line;
        std::advance(line, stride);
        std::advance(containerIterator, *line - *lastLine);
      }
      return *this;
    }

    /// Equality operator
    template<class Iterator>
    bool operator!=(const Iterator& iterator)
    {
      return iterator.getPosition() != getPosition();
    }

    /**
     * Dereferencement operator
     * @return a reference to the value pointed by the iterator
     */
    typename Parent::reference operator*()
    {
      return *containerIterator;
    }

    /**
     * Dereferencement operator
     * @return the value pointed by the iterator
     */
    typename Parent::value_type operator*() const
    {
      return *containerIterator;
    }

    /**
     * Returns an line iterator pointing at the beginning of the line at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a line iterator pointing at the beginning of the current line
     */
    Self lineBegin(bool useDefault = true)
    {
      unsigned useStride = (useDefault?containerPointer->height():1U);
      return Self(containerPointer, containerIterator.lineBegin(false), line, lineBegin_, matrix, position, useStride);
    }

    /**
     * Returns an line iterator pointing at the beginning of the line at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a line iterator pointing at the beginning of the current line
     */
    ConstSelf lineBegin(bool useDefault = true) const
    {
      unsigned useStride = (useDefault?containerPointer->height():1U);
      return ConstSelf(containerPointer, containerIterator.lineBegin(false), line, lineBegin_, matrix, position, useStride);
    }

    /**
     * Returns an line iterator to the end of the line at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a line iterator pointing to the end of the current line
     */
    ConstSelf lineEnd(bool useDefault = true) const
    {
      unsigned useStride = (useDefault?containerPointer->height():1U);
      ContainerIterator currentIt = containerIterator.lineBegin(false);
      std::advance(currentIt, containerPointer->realHeight() * containerPointer->realWidth());
      return Self(containerPointer, currentIt, line, lineBegin_, matrix, position, useStride);
    }

    /**
     * Returns an column iterator pointing at the beginning of the column at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a column iterator pointing at the beginning of the current column
     */
    Self columnBegin(bool useDefault = true)
    {
      unsigned useStride = 1U;
      return Self(containerPointer, containerIterator.columnBegin(false), lineBegin_, lineBegin_, matrix, 0, useStride);
    }

    /**
     * Returns an column iterator pointing at the beginning of the column at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a column iterator pointing at the beginning of the current column
     */
    ConstSelf columnBegin(bool useDefault = true) const
    {
      unsigned useStride = 1U;
      return ConstSelf(containerPointer, containerIterator.columnBegin(false), lineBegin_, lineBegin_, matrix, 0, useStride);
    }

    /**
     * Returns an column iterator to the end of the line at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a column iterator pointing to the end of the current line
     */
    ConstSelf columnEnd(bool useDefault = true) const
    {
      unsigned useStride = 1U;
      ContainerIterator currentIt = containerIterator.lineBegin(false);
      std::advance(currentIt, containerPointer->height());
      return ConstSelf(containerPointer, currentIt, lineBegin_, lineBegin_, matrix, 0, useStride);
    }
  };

  /// The standard iterator for sub column matrix
  template<class Container>
  class SubColumnMatrixIterator : public boost::mpl::if_<typename boost::is_const<typename boost::remove_reference<Container>::type>::type, std::iterator<std::forward_iterator_tag, typename boost::add_const<typename Container::Data_Type>::type, ptrdiff_t, typename boost::add_const<typename Container::Data_Type>::type*, typename boost::add_const<typename Container::Data_Type>::type&>, std::iterator<std::forward_iterator_tag, typename Container::Data_Type, ptrdiff_t, typename Container::Data_Type*, typename Container::Data_Type&> >::type
  {
    typedef SubColumnMatrixIterator Self;
    /// Type of the parent
    typedef typename boost::mpl::if_<typename boost::is_const<typename boost::remove_reference<Container>::type>::type, std::iterator<std::forward_iterator_tag, typename boost::add_const<typename Container::Data_Type>::type, ptrdiff_t, typename boost::add_const<typename Container::Data_Type>::type*, typename boost::add_const<typename Container::Data_Type>::type&>, std::iterator<std::forward_iterator_tag, typename Container::Data_Type, ptrdiff_t, typename Container::Data_Type*, typename Container::Data_Type&> >::type Parent;

    typedef SubColumnMatrixIterator<typename boost::add_const<Container>::type> ConstSelf;
    typedef typename boost::mpl::if_<typename boost::is_const<typename boost::remove_reference<Container>::type>::type, typename Container::MatrixType::const_iterator, typename Container::MatrixType::iterator>::type ContainerIterator;

    typedef typename boost::mpl::if_<typename boost::is_const<typename boost::remove_reference<Container>::type>::type, typename Container::IterateColumnVector::const_iterator, typename Container::IterateColumnVector::iterator>::type ColumnIterator;

    /// Pointer to the container
    Container* containerPointer;
    /// Iterator to the parent
    ContainerIterator containerIterator;
    /// Iterator pointing to a specific column
    ColumnIterator column;
    /// Iterator pointing to a specific column
    ColumnIterator columnBegin_;
    /// Iterator pointing to a specific column
    ColumnIterator columnEnd_;
    /// Position in a column
    unsigned long position;
    /// Stride for operations
    unsigned long stride;
  public:

    /// Returns the position inside the container - not enough to check equality !
    std::pair<ColumnIterator, unsigned long> getPosition() const
    {
      return std::pair<ColumnIterator, unsigned long>(column, position);
    }

    /**
     * Constructor from a container and a position in the container
     * @param container is the container on which to iterate
     * @param containerIterator is a iterator in the container
     * @param column is the column that is currently used
     * @param position is the position in the column
     * @param stride is the step used in (in/de)crementation
     */
    SubColumnMatrixIterator(Container* container, const ContainerIterator& containerIterator, const ColumnIterator& column, const ColumnIterator& columnBegin_, const ColumnIterator& columnEnd_, unsigned long position = 0, unsigned long stride = 1U)
      :containerPointer(container), containerIterator(containerIterator), column(column), columnBegin_(columnBegin_), columnEnd_(columnEnd_), position(position), stride(stride)
    {
    }

    /// Preincrementation
    Self& operator++()
    {
      position += stride;
      std::advance(containerIterator, stride);
      if(position >= containerPointer->height())
      {
        position = position % containerPointer->height();
        ColumnIterator lastColumn = column;
        ++column;
        std::advance(containerIterator, (*column - *lastColumn - 1) * containerPointer->height());
      }
      return *this;
    }

    /// Equality operator
    template<class Iterator>
    bool operator!=(const Iterator& iterator)
    {
      return iterator.getPosition() != getPosition();
    }

    /**
     * Dereferencement operator
     * @return a reference to the value pointed by the iterator
     */
    typename Parent::reference operator*()
    {
      return *containerIterator;
    }

    /**
     * Dereferencement operator
     * @return the value pointed by the iterator
     */
    typename Parent::value_type operator*() const
    {
      return *containerIterator;
    }

    /**
     * Returns an line iterator pointing at the beginning of the line at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a line iterator pointing at the beginning of the current line
     */
    Self lineBegin(bool useDefault = true)
    {
      unsigned useStride = (useDefault?containerPointer->height():1U);
      ContainerIterator currentIt = containerIterator.lineBegin(false);
      std::advance(currentIt, *columnBegin_ * containerPointer->height());
      return Self(containerPointer, currentIt, columnBegin_, columnBegin_, columnEnd_, position, useStride);
    }

    /**
     * Returns an line iterator pointing at the beginning of the line at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a line iterator pointing at the beginning of the current line
     */
    ConstSelf lineBegin(bool useDefault = true) const
    {
      unsigned useStride = (useDefault?containerPointer->height():1U);
      ContainerIterator currentIt = containerIterator.lineBegin(false);
      std::advance(currentIt, *columnBegin_ * containerPointer->height());
      return ConstSelf(containerPointer, currentIt, columnBegin_, columnBegin_, columnEnd_, position, useStride);
    }

    /**
     * Returns an line iterator to the end of the line at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a line iterator pointing to the end of the current line
     */
    ConstSelf lineEnd(bool useDefault = true) const
    {
      unsigned useStride = (useDefault?containerPointer->height():1U);
      ContainerIterator currentIt = containerIterator.lineBegin(false);
      std::advance(currentIt, *columnEnd_ * containerPointer->height());
      return ConstSelf(containerPointer, currentIt, columnEnd_, columnBegin_, columnEnd_, position, useStride);
    }

    /**
     * Returns an column iterator pointing at the beginning of the column at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a column iterator pointing at the beginning of the current column
     */
    Self columnBegin(bool useDefault = true)
    {
      unsigned useStride = 1U;
      return Self(containerPointer, containerIterator.columnBegin(false), column, columnBegin_, columnEnd_, 0, useStride);
    }

    /**
     * Returns an column iterator pointing at the beginning of the column at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a column iterator pointing at the beginning of the current column
     */
    ConstSelf columnBegin(bool useDefault = true) const
    {
      unsigned useStride = 1U;
      return ConstSelf(containerPointer, containerIterator.columnBegin(false), column, columnBegin_, columnEnd_, 0, useStride);
    }

    /**
     * Returns an column iterator to the end of the line at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a column iterator pointing to the end of the current line
     */
    ConstSelf columnEnd(bool useDefault = true) const
    {
      unsigned useStride = 1U;
      ContainerIterator currentIt = containerIterator.columnBegin(false);
      std::advance(currentIt, containerPointer->height());
      ColumnIterator columnCopy = column;
      ++columnCopy;
      return ConstSelf(containerPointer, currentIt, columnCopy, columnBegin_, columnEnd_, 0, useStride);
    }
  };

  /// The standard iterator for sub column matrix
  template<class Container>
  class SubColumnLineMatrixIterator : public boost::mpl::if_<typename boost::is_const<typename boost::remove_reference<Container>::type>::type, std::iterator<std::forward_iterator_tag, typename boost::add_const<typename Container::Data_Type>::type, ptrdiff_t, typename boost::add_const<typename Container::Data_Type>::type*, typename boost::add_const<typename Container::Data_Type>::type&>, std::iterator<std::forward_iterator_tag, typename Container::Data_Type, ptrdiff_t, typename Container::Data_Type*, typename Container::Data_Type&> >::type
  {
    typedef SubColumnLineMatrixIterator Self;
    /// Type of the parent
    typedef typename boost::mpl::if_<typename boost::is_const<typename boost::remove_reference<Container>::type>::type, std::iterator<std::forward_iterator_tag, typename boost::add_const<typename Container::Data_Type>::type, ptrdiff_t, typename boost::add_const<typename Container::Data_Type>::type*, typename boost::add_const<typename Container::Data_Type>::type&>, std::iterator<std::forward_iterator_tag, typename Container::Data_Type, ptrdiff_t, typename Container::Data_Type*, typename Container::Data_Type&> >::type Parent;

    typedef SubColumnLineMatrixIterator<typename boost::add_const<Container>::type> ConstSelf;
    typedef typename boost::mpl::if_<typename boost::is_const<typename boost::remove_reference<Container>::type>::type, typename Container::MatrixType::const_iterator, typename Container::MatrixType::iterator>::type ContainerIterator;

    typedef typename boost::mpl::if_<typename boost::is_const<typename boost::remove_reference<Container>::type>::type, typename Container::IterateColumnVector::const_iterator, typename Container::IterateColumnVector::iterator>::type ColumnIterator;
    typedef typename boost::mpl::if_<typename boost::is_const<typename boost::remove_reference<Container>::type>::type, typename Container::IterateLineVector::const_iterator, typename Container::IterateLineVector::iterator>::type LineIterator;

    typedef typename boost::add_const<typename Container::MatrixType>::type ConstMatrixType;
    /// Pointer to the container
    Container* containerPointer;
    /// Iterator to the parent
    ContainerIterator containerIterator;
    /// Iterator pointing to a specific column
    ColumnIterator column;
    /// Iterator pointing to a specific column
    ColumnIterator columnBegin_;
    /// Iterator pointing to a specific column
    ColumnIterator columnEnd_;
    /// Iterator pointing to a specific line
    LineIterator line;
    /// Iterator to the first line
    LineIterator lineBegin_;
    /// Pointer to the inner container
    ConstMatrixType* matrix;
    /// Position in a column
    unsigned long position;
    /// Stride for operations
    unsigned long stride;
  public:

    /// Returns the position inside the container - not enough to check equality !
    std::pair<ColumnIterator, LineIterator> getPosition() const
    {
      return std::pair<ColumnIterator, LineIterator>(column, line);
    }

    /**
     * Constructor from a container and a position in the container
     * @param container is the container on which to iterate
     * @param containerIterator is a iterator in the container
     * @param column is the column that is currently used
     * @param line is the line that is currently used
     * @param lineBegin is the start of the line
     * @param position is the position in the column
     * @param stride is the step used in (in/de)crementation
     */
    SubColumnLineMatrixIterator(Container* container, const ContainerIterator& containerIterator, const ColumnIterator& column, const ColumnIterator& columnBegin_, const ColumnIterator& columnEnd_, const LineIterator& line, const LineIterator& lineBegin_, ConstMatrixType* matrix, unsigned long position = 0, unsigned long stride = 1U)
      :containerPointer(container), containerIterator(containerIterator), column(column), columnBegin_(columnBegin_), columnEnd_(columnEnd_), line(line), lineBegin_(lineBegin_), matrix(matrix), position(position), stride(stride)
    {
    }

    /// Preincrementation
    Self& operator++()
    {
      position += stride;
      if(position >= containerPointer->height())
      {
        position = position % containerPointer->height();
        ColumnIterator lastColumn = column;
        ++column;
        LineIterator lastLine = line;
        line = lineBegin_;
        std::advance(line, position);
        std::advance(containerIterator, *line - *lastLine + (*column - *lastColumn) * matrix->height());
      }
      else
      {
        LineIterator lastLine = line;
        std::advance(line, stride);
        std::advance(containerIterator, *line - *lastLine);
      }
      return *this;
    }

    /// Equality operator
    template<class Iterator>
    bool operator!=(const Iterator& iterator)
    {
      return iterator.getPosition() != getPosition();
    }

    /**
     * Dereferencement operator
     * @return a reference to the value pointed by the iterator
     */
    typename Parent::reference operator*()
    {
      return *containerIterator;
    }

    /**
     * Dereferencement operator
     * @return the value pointed by the iterator
     */
    typename Parent::value_type operator*() const
    {
      return *containerIterator;
    }

    /**
     * Returns an line iterator pointing at the beginning of the line at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a line iterator pointing at the beginning of the current line
     */
    Self lineBegin(bool useDefault = true)
    {
      unsigned useStride = (useDefault?containerPointer->height():1U);
      ContainerIterator currentIt = containerIterator.lineBegin(false);
      std::advance(currentIt, *columnBegin_ * containerPointer->realHeight());
      return Self(containerPointer, currentIt, columnBegin_, columnBegin_, columnEnd_, line, lineBegin_, matrix, position, useStride);
    }

    /**
     * Returns an line iterator pointing at the beginning of the line at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a line iterator pointing at the beginning of the current line
     */
    ConstSelf lineBegin(bool useDefault = true) const
    {
      unsigned useStride = (useDefault?containerPointer->height():1U);
      ContainerIterator currentIt = containerIterator.lineBegin(false);
      std::advance(currentIt, *columnBegin_ * containerPointer->realHeight());
      return ConstSelf(containerPointer, currentIt, columnBegin_, columnBegin_, columnEnd_, line, lineBegin_, matrix, position, useStride);
    }

    /**
     * Returns an line iterator to the end of the line at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a line iterator pointing to the end of the current line
     */
    ConstSelf lineEnd(bool useDefault = true) const
    {
      unsigned useStride = (useDefault?containerPointer->height():1U);
      ContainerIterator currentIt = containerIterator.lineBegin(false);
      std::advance(currentIt, *columnEnd_ * containerPointer->realHeight());
      return ConstSelf(containerPointer, currentIt, columnEnd_, columnBegin_, columnEnd_, line, lineBegin_, matrix, position, useStride);
    }

    /**
     * Returns an column iterator pointing at the beginning of the column at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a column iterator pointing at the beginning of the current column
     */
    Self columnBegin(bool useDefault = true)
    {
      unsigned useStride = 1U;
      return Self(containerPointer, containerIterator.columnBegin(false), column, columnBegin_, columnEnd_, lineBegin_, lineBegin_, matrix, 0, useStride);
    }

    /**
     * Returns an column iterator pointing at the beginning of the column at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a column iterator pointing at the beginning of the current column
     */
    ConstSelf columnBegin(bool useDefault = true) const
    {
      unsigned useStride = 1U;
      return ConstSelf(containerPointer, containerIterator.columnBegin(false), column, lineBegin_, lineBegin_, matrix, 0, useStride);
    }

    /**
     * Returns an column iterator to the end of the line at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a column iterator pointing to the end of the current line
     */
    ConstSelf columnEnd(bool useDefault = true) const
    {
      unsigned useStride = 1U;
      ContainerIterator currentIt = containerIterator.columnBegin(false);
      std::advance(currentIt, matrix->height());
      ColumnIterator nextColumn = column;
      ++nextColumn;
      return ConstSelf(containerPointer, currentIt, nextColumn, columnBegin_, columnEnd_, lineBegin_, lineBegin_, matrix, 0, useStride);
    }
  };
}

#endif
