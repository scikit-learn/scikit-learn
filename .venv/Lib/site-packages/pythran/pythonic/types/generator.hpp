#ifndef PYTHONIC_TYPES_GENERATOR_HPP
#define PYTHONIC_TYPES_GENERATOR_HPP

#include "pythonic/include/types/generator.hpp"

#include "pythonic/builtins/StopIteration.hpp"
#include <iterator>

PYTHONIC_NS_BEGIN

namespace types
{
  template <class T>
  generator_iterator<T>::generator_iterator() : the_generator()
  {
    the_generator.__generator_state = -1;
  } // this represents the end

  template <class T>
  generator_iterator<T>::generator_iterator(T const &a_generator) : the_generator(a_generator)
  {
  }

  template <class T>
  generator_iterator<T> &generator_iterator<T>::operator++()
  {
    try {
      the_generator.next();
    } catch (types::StopIteration const &) {
      the_generator.__generator_state = -1;
    }
    return *this;
  }

  template <class T>
  typename T::result_type generator_iterator<T>::operator*() const
  {
    return *the_generator;
  }

  template <class T>
  bool generator_iterator<T>::operator!=(generator_iterator<T> const &other) const
  {
    assert(other.the_generator.__generator_state == -1 || the_generator.__generator_state == -1);
    return the_generator.__generator_state != other.the_generator.__generator_state;
  }

  template <class T>
  bool generator_iterator<T>::operator==(generator_iterator<T> const &other) const
  {
    assert(other.the_generator.__generator_state == -1 || the_generator.__generator_state == -1);
    return the_generator.__generator_state == other.the_generator.__generator_state;
  }

  template <class T>
  bool generator_iterator<T>::operator<(generator_iterator<T> const &other) const
  {
    assert(other.the_generator.__generator_state == -1 || the_generator.__generator_state == -1);
    return the_generator.__generator_state != other.the_generator.__generator_state;
  }
} // namespace types
PYTHONIC_NS_END

#endif
