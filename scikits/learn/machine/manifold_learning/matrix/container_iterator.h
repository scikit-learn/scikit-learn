/**
 * \file container_iterator.h
 * Contains the standard container iterator
 */

#ifndef CONTAINERITERATOR
#define CONTAINERITERATOR

#include <iterator>
#include <boost/mpl/if.hpp>
#include <boost/type_traits.hpp>

namespace Matrix
{
  /// The standard iterator for linear container
  template<class Container>
  class LinearIterator : public boost::mpl::if_<typename boost::is_const<typename boost::remove_reference<Container>::type>::type, std::iterator<std::random_access_iterator_tag, typename boost::add_const<typename Container::Data_Type>::type, ptrdiff_t, typename boost::add_const<typename Container::Data_Type>::type*, typename boost::add_const<typename Container::Data_Type>::type&>, std::iterator<std::random_access_iterator_tag, typename Container::Data_Type, ptrdiff_t, typename Container::Data_Type*, typename Container::Data_Type&> >::type
  {
    /// Type of the parent
    typedef typename boost::mpl::if_<typename boost::is_const<typename boost::remove_reference<Container>::type>::type, std::iterator<std::random_access_iterator_tag, typename boost::add_const<typename Container::Data_Type>::type, ptrdiff_t, typename boost::add_const<typename Container::Data_Type>::type*, typename boost::add_const<typename Container::Data_Type>::type&>, std::iterator<std::random_access_iterator_tag, typename Container::Data_Type, ptrdiff_t, typename Container::Data_Type*, typename Container::Data_Type&> >::type Parent;

    typedef LinearIterator<typename boost::add_const<Container>::type> ConstSelf;

    /// Pointer to the data
    typename Parent::pointer dataPointer;
    /// Pointer to the container
    Container* containerPointer;
    /// Stride for operations
    unsigned long stride;
  public:
    typedef LinearIterator Self;
    /**
     * Return the inner pointer
     * @return the pointer in the iterator
     */
    typename Parent::pointer getPointer() const
    {
      return dataPointer;
    }
    /**
     * Return the inner pointer
     * @return the pointer in the iterator
     */
    Container* getContainer() const
    {
      return containerPointer;
    }

    /**
     * Constructor from a container and a position in the container
     * @param container is the container on which to iterate
     * @param position is a pointer in the data
     * @param stride is the step used in (in/de)crementation
     */
    LinearIterator(Container* container, typename Parent::pointer pointer, unsigned long stride = 1U)
      :dataPointer(pointer), containerPointer(container), stride(stride)
    {
    }

    /// Preincrementation
    Self& operator++()
    {
      dataPointer += stride;
      return *this;
    }

    /// Predecrementation
    Self& operator--()
    {
      dataPointer -= stride;
      return *this;
    }

    /// Equality operator
    template<class Iterator>
    bool operator!=(const Iterator& iterator) const
    {
      return iterator.getContainer() != getContainer() || iterator.getPointer() != getPointer();
    }
    /// Equality operator
    template<class Iterator>
    bool operator==(const Iterator& iterator) const
    {
      return iterator.getContainer() == getContainer() && iterator.getPointer() == getPointer();
    }

    /// Less operator
    template<class Iterator>
    bool operator<(const Iterator& iterator) const
    {
      return iterator.getContainer() == getContainer() && getPointer() < iterator.getPointer();
    }

    /// Difference type
    template<class Iterator>
    typename Parent::difference_type operator-(const Iterator& iterator)
    {
      return getPointer() - iterator.getPointer();
    }

    /// + operator
    Self operator+(typename Parent::difference_type distance)
    {
      Self newIt(*this);
      newIt += distance;
      return newIt;
    }

    /// += operator
    Self& operator+=(typename Parent::difference_type distance)
    {
      dataPointer += distance * stride;
      return (*this);
    }

    /**
     * Dereferencement operator
     * @return a reference to the value pointed by the iterator
     */
    typename Parent::reference operator*()
    {
//      assert(dataPointer - containerPointer->get() >= 0);
//      assert(dataPointer - containerPointer->get() < static_cast<typename Parent::difference_type>(containerPointer->height() * containerPointer->width()));
      return *dataPointer;
    }

    /**
     * Dereferencement operator
     * @return the value pointed by the iterator
     */
    typename Parent::value_type operator*() const
    {
      return *dataPointer;
    }

    /**
     * Returns an line iterator pointing at the beginning of the line at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a line iterator pointing at the beginning of the current line
     */
    Self lineBegin(bool useDefault = true)
    {
      unsigned useStride = (useDefault?containerPointer->height():1U);
      return Self(containerPointer, containerPointer->get() + (dataPointer - containerPointer->get()) % containerPointer->height(), useStride);
    }

    /**
     * Returns an line iterator pointing at the beginning of the line at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a line iterator pointing at the beginning of the current line
     */
    ConstSelf lineBegin(bool useDefault = true) const
    {
      unsigned useStride = (useDefault?containerPointer->height():1U);
      return ConstSelf(containerPointer, containerPointer->get() + (dataPointer - containerPointer->get()) % containerPointer->height(), useStride);
    }

    /**
     * Returns an line iterator to the end of the line at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a line iterator pointing to the end of the current line
     */
    ConstSelf lineEnd(bool useDefault = true) const
    {
      unsigned useStride = (useDefault?containerPointer->height():1U);
      return ConstSelf(containerPointer, containerPointer->get() + (dataPointer - containerPointer->get()) % containerPointer->height() + containerPointer->height() * containerPointer->width(), useStride);
    }

    /**
     * Returns an column iterator pointing at the beginning of the column at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a column iterator pointing at the beginning of the current column
     */
    Self columnBegin(bool useDefault = true)
    {
      unsigned useStride = 1U;
      return Self(containerPointer, dataPointer - (dataPointer - containerPointer->get()) % containerPointer->height(), useStride);
    }

    /**
     * Returns an column iterator pointing at the beginning of the column at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a column iterator pointing at the beginning of the current column
     */
    ConstSelf columnBegin(bool useDefault = true) const
    {
      unsigned useStride = 1U;
      return ConstSelf(containerPointer, dataPointer - (dataPointer - containerPointer->get()) % containerPointer->height(), useStride);
    }

    /**
     * Returns an column iterator to the end of the line at which the current iterator is
     * @param useDefault indicates if the the default stride is to be used
     * @return a column iterator pointing to the end of the current line
     */
    ConstSelf columnEnd(bool useDefault = true) const
    {
      unsigned useStride = 1U;
      return ConstSelf(containerPointer, dataPointer - (dataPointer - containerPointer->get()) % containerPointer->height() + containerPointer->height(), useStride);
    }
  };
}
#endif
