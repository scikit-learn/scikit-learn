#ifndef PYTHONIC_INCLUDE_NUMPY_COUNT_NONZERO_HPP
#define PYTHONIC_INCLUDE_NUMPY_COUNT_NONZERO_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class dtype, class E>
  auto _count_nonzero(E begin, E end, long &count, utils::int_<1>)
      -> std::enable_if_t<std::is_same<dtype, bool>::value>;

  template <class dtype, class E>
  auto _count_nonzero(E begin, E end, long &count, utils::int_<1>)
      -> std::enable_if_t<!std::is_same<dtype, bool>::value>;

  template <class dtype, class E, size_t N>
  void _count_nonzero(E begin, E end, long &count, utils::int_<N>);

  template <class E>
  long count_nonzero(E const &array);

  DEFINE_FUNCTOR(pythonic::numpy, count_nonzero);
} // namespace numpy
PYTHONIC_NS_END

#endif
