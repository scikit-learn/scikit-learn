#ifndef PYTHONIC_NUMPY_COUNT_NONZERO_HPP
#define PYTHONIC_NUMPY_COUNT_NONZERO_HPP

#include "pythonic/include/numpy/count_nonzero.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class dtype, class E>
  auto _count_nonzero(E begin, E end, long &count, utils::int_<1>)
      -> std::enable_if_t<std::is_same<dtype, bool>::value>
  {
    for (; begin != end; ++begin)
      // Behaviour defined in the standard
      count += *begin;
  }

  template <class dtype, class E>
  auto _count_nonzero(E begin, E end, long &count, utils::int_<1>)
      -> std::enable_if_t<!std::is_same<dtype, bool>::value>
  {
    for (; begin != end; ++begin)
      if (*begin != static_cast<dtype>(0))
        ++count;
  }

  template <class dtype, class E, size_t N>
  void _count_nonzero(E begin, E end, long &count, utils::int_<N>)
  {
    for (; begin != end; ++begin)
      _count_nonzero<dtype>((*begin).begin(), (*begin).end(), count, utils::int_<N - 1>());
  }

  template <class E>
  long count_nonzero(E const &array)
  {
    long count(0);
    _count_nonzero<typename E::dtype>(array.begin(), array.end(), count, utils::int_<E::value>());
    return count;
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
