#ifndef PYTHONIC_INCLUDE_TYPES_GENERATOR_HPP
#define PYTHONIC_INCLUDE_TYPES_GENERATOR_HPP

#include <cstddef>
#include <iterator>

PYTHONIC_NS_BEGIN

namespace types
{
  template <class T>
  struct generator_iterator
      : std::iterator<std::forward_iterator_tag, typename T::result_type, ptrdiff_t,
                      typename T::result_type *, typename T::result_type /* no ref */> {

    T the_generator;
    generator_iterator();
    generator_iterator(T const &a_generator);
    generator_iterator &operator++();
    typename T::result_type operator*() const;
    bool operator!=(generator_iterator<T> const &other) const;
    bool operator==(generator_iterator<T> const &other) const;
    bool operator<(generator_iterator<T> const &other) const;
  };
} // namespace types
PYTHONIC_NS_END

#endif
