/**
 * \file diagonal_matrix_iterator.h
 * Contains the iterator for diagonal only vision
 */

#ifndef DIAGONALITERATOR
#define DIAGONALITERATOR

#include <iterator>
#include <boost/mpl/if.hpp>
#include <boost/type_traits.hpp>

namespace Matrix
{
  /// The standard iterator for diagonal matrix
  template<class Container>
  class DiagonalIterator : public boost::mpl::if_<typename boost::is_const<typename boost::remove_reference<Container>::type>::type, std::iterator<std::random_access_iterator_tag, typename boost::add_const<typename Container::Data_Type>::type, ptrdiff_t, typename boost::add_const<typename Container::Data_Type>::type*, typename boost::add_const<typename Container::Data_Type>::type&>, std::iterator<std::random_access_iterator_tag, typename Container::Data_Type, ptrdiff_t, typename Container::Data_Type*, typename Container::Data_Type&> >::type
  {
    /// Type of the parent
    typedef typename boost::mpl::if_<typename boost::is_const<typename boost::remove_reference<Container>::type>::type, std::iterator<std::random_access_iterator_tag, typename boost::add_const<typename Container::Data_Type>::type, ptrdiff_t, typename boost::add_const<typename Container::Data_Type>::type*, typename boost::add_const<typename Container::Data_Type>::type&>, std::iterator<std::random_access_iterator_tag, typename Container::Data_Type, ptrdiff_t, typename Container::Data_Type*, typename Container::Data_Type&> >::type Parent;

    typedef DiagonalIterator<typename boost::add_const<Container>::type> ConstSelf;

    /// Position in the matrix
    unsigned long position;
    /// Pointer to the container
    Container* containerPointer;
    /// Stride for operations
    unsigned long stride;
  public:
    typedef DiagonalIterator Self;
    /// Type of the data
    typedef typename Parent::value_type Data_Type;

    /// Returns the position inside the container
    unsigned long getPosition() const
    {
      return position;
    }

    /**
     * Constructor from a container and a position in the container
     * @param container is the container on which to iterate
     * @param position is a pointer in the data
     * @param stride is the step used in (in/de)crementation
     */
    DiagonalIterator(Container* container, unsigned long position, unsigned long stride = 1U)
      :position(position), containerPointer(container), stride(stride)
    {
    }

    /// Preincrementation
    Self& operator++()
    {
      position += stride;
      return *this;
    }

    /// Predecrementation
    Self& operator--()
    {
      position -= stride;
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

    /// Less operator
    template<class Iterator>
    bool operator<(const Iterator& iterator) const
    {
      return getPosition() < iterator.getPosition();
    }

    /// Difference type
    template<class Iterator>
    typename Parent::difference_type operator-(const Iterator& iterator)
    {
      return getPosition() - iterator.getPosition();
    }

    /// + operator
    Self operator+(typename Parent::difference_type distance)
    {
      Self newIt(*this);
      newIt.position += distance * stride;
      return newIt;
    }

    /// += operator
    Self& operator+=(typename Parent::difference_type distance)
    {
      position += distance * stride;
      return *this;
    }

    /**
     * Dereferencement operator
     * @return the value pointed by the iterator
     */
    typename Parent::value_type operator*() const
    {
      unsigned long subposition = position % containerPointer->height();
      if(subposition == position / containerPointer->height())
        return (*containerPointer)(subposition, subposition);
      else
        return Container::DataTraits::zero((*containerPointer)(subposition, subposition));
    }

    /**
     * Returns an line iterator pointing at the beginning of the line at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a line iterator pointing at the beginning of the current line
     */
    Self lineBegin(bool useDefault = true)
    {
      unsigned useStride = (useDefault?containerPointer->height():1U);
      return Self(containerPointer, position % containerPointer->height(), useStride);
    }

    /**
     * Returns an line iterator pointing at the beginning of the line at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a line iterator pointing at the beginning of the current line
     */
    ConstSelf lineBegin(bool useDefault = true) const
    {
      unsigned useStride = (useDefault?containerPointer->height():1U);
      return ConstSelf(containerPointer, position % containerPointer->height(), useStride);
    }

    /**
     * Returns an line iterator to the end of the line at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a line iterator pointing to the end of the current line
     */
    ConstSelf lineEnd(bool useDefault = true) const
    {
      unsigned useStride = (useDefault?containerPointer->height():1U);
      return ConstSelf(containerPointer, position % containerPointer->height() + containerPointer->height() * containerPointer->width(), useStride);
    }

    /**
     * Returns an column iterator pointing at the beginning of the column at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a column iterator pointing at the beginning of the current column
     */
    Self columnBegin(bool useDefault = true)
    {
      unsigned useStride = 1U;
      return Self(containerPointer, position - position % containerPointer->height(), useStride);
    }

    /**
     * Returns an column iterator pointing at the beginning of the column at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a column iterator pointing at the beginning of the current column
     */
    ConstSelf columnBegin(bool useDefault = true) const
    {
      unsigned useStride = 1U;
      return ConstSelf(containerPointer, position - position % containerPointer->height(), useStride);
    }

    /**
     * Returns an column iterator to the end of the line at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a column iterator pointing to the end of the current line
     */
    ConstSelf columnEnd(bool useDefault = true) const
    {
      unsigned useStride = 1U;
      return ConstSelf(containerPointer, position - position % containerPointer->height() + containerPointer->height(), useStride);
    }
  };
}
#endif
