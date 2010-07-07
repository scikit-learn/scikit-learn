/**
 * \file sub_matrix_iterator.h
 * Contains the iterator for a submatrix
 */

#ifndef SUBMATRIXITERATOR
#define SUBMATRIXITERATOR

#include <iterator>
#include <boost/mpl/if.hpp>
#include <boost/type_traits.hpp>

namespace Matrix
{
  /// The standard iterator for a submatrix
  template<class Container>
  class SubMatrixIterator : public boost::mpl::if_<typename boost::is_const<typename boost::remove_reference<Container>::type>::type, std::iterator<std::bidirectional_iterator_tag, typename boost::add_const<typename Container::Data_Type>::type, ptrdiff_t, typename boost::add_const<typename Container::Data_Type>::type*, typename boost::add_const<typename Container::Data_Type>::type&>, std::iterator<std::bidirectional_iterator_tag, typename Container::Data_Type, ptrdiff_t, typename Container::Data_Type*, typename Container::Data_Type&> >::type
  {
    typedef SubMatrixIterator Self;
    /// Type of the parent
    typedef typename boost::mpl::if_<typename boost::is_const<typename boost::remove_reference<Container>::type>::type, std::iterator<std::bidirectional_iterator_tag, typename boost::add_const<typename Container::Data_Type>::type, ptrdiff_t, typename boost::add_const<typename Container::Data_Type>::type*, typename boost::add_const<typename Container::Data_Type>::type&>, std::iterator<std::bidirectional_iterator_tag, typename Container::Data_Type, ptrdiff_t, typename Container::Data_Type*, typename Container::Data_Type&> >::type Parent;

    typedef SubMatrixIterator<typename boost::add_const<Container>::type> ConstSelf;
    typedef typename boost::mpl::if_<typename boost::is_const<typename boost::remove_reference<Container>::type>::type, typename Container::MatrixType::const_iterator, typename Container::MatrixType::iterator>::type ContainerIterator;
    typedef typename Container::MatrixType::const_iterator ConstContainerIterator;

    /// Iterator to the parent
    ContainerIterator containerIterator;
    /// Position in a column
    long position;
    /// Pointer to the container
    Container* containerPointer;
    /// Stride for operations
    unsigned long stride;
  public:
    /// Type of the data
    typedef typename Parent::value_type Data_Type;

    /// Returns the position inside the container
    ContainerIterator getPosition() const
    {
      return containerIterator;
    }

    /**
     * Constructor from a container and a position in the container
     * @param container is the container on which to iterate
     * @param containerIterator is a iterator in the container
     * @param stride is the step used in (in/de)crementation
     */
    SubMatrixIterator(Container* container, const ContainerIterator& containerIterator, long position = 0, unsigned long stride = 1U)
      :containerIterator(containerIterator), position(position), containerPointer(container), stride(stride)
    {
    }

    /// Preincrementation
    Self& operator++()
    {
      position += stride;
      std::advance(containerIterator, stride);
      if(position >= static_cast<long>(containerPointer->height()))
      {
        position = position % containerPointer->height();
        std::advance(containerIterator, containerPointer->realHeight() - containerPointer->height());
      }
      return *this;
    }

    /// Predecrementation
    Self& operator--()
    {
      position -= stride;
      std::advance(containerIterator, - stride);
      if(position < 0)
      {
        position = position % containerPointer->height();
        std::advance(containerIterator, - containerPointer->realHeight() + containerPointer->height());
      }
      return *this;
    }

    /// Equality operator
    template<class Iterator>
    bool operator!=(const Iterator& iterator) const
    {
      return iterator.getPosition() != getPosition();
    }

    /// Equality operator
    template<class Iterator>
    bool operator==(const Iterator& iterator) const
    {
      return iterator.getPosition() == getPosition();
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
      std::advance(currentIt, containerPointer->lineOffset() * containerPointer->realHeight());
      return Self(containerPointer, currentIt, position, useStride);
    }

    /**
     * Returns an line iterator pointing at the beginning of the line at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a line iterator pointing at the beginning of the current line
     */
    ConstSelf lineBegin(bool useDefault = true) const
    {
      unsigned useStride = (useDefault?containerPointer->height():1U);
      ConstContainerIterator currentIt = containerIterator.lineBegin(false);
      std::advance(currentIt, containerPointer->lineOffset() * containerPointer->realHeight());
      return ConstSelf(containerPointer, currentIt, position, useStride);
    }

    /**
     * Returns an line iterator to the end of the line at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a line iterator pointing to the end of the current line
     */
    ConstSelf lineEnd(bool useDefault = true) const
    {
      unsigned useStride = (useDefault?containerPointer->height():1U);
      ConstContainerIterator currentIt = containerIterator.lineBegin(false);
      std::advance(currentIt, (containerPointer->lineOffset() + containerPointer->width()) * containerPointer->realHeight());
      return ConstSelf(containerPointer, currentIt, position, useStride);
    }

    /**
     * Returns an column iterator pointing at the beginning of the column at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a column iterator pointing at the beginning of the current column
     */
    Self columnBegin(bool useDefault = true)
    {
      unsigned useStride = 1U;
      ContainerIterator currentIt = containerIterator.columnBegin(false);
      std::advance(currentIt, containerPointer->columnOffset());
      return Self(containerPointer, currentIt, 0, useStride);
    }

    /**
     * Returns an column iterator pointing at the beginning of the column at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a column iterator pointing at the beginning of the current column
     */
    ConstSelf columnBegin(bool useDefault = true) const
    {
      unsigned useStride = 1U;
      ConstContainerIterator currentIt = containerIterator.columnBegin(false);
      std::advance(currentIt, containerPointer->columnOffset());
      return Self(containerPointer, currentIt, 0, useStride);
    }

    /**
     * Returns an column iterator to the end of the line at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a column iterator pointing to the end of the current line
     */
    ConstSelf columnEnd(bool useDefault = true) const
    {
      unsigned useStride = 1U;
      ConstContainerIterator currentIt = containerIterator.columnBegin(false);
      std::advance(currentIt, containerPointer->columnOffset() + containerPointer->realHeight());
      return Self(containerPointer, currentIt, 0, useStride);
    }
  };
}
#endif
