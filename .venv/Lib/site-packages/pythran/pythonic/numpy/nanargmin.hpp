#ifndef PYTHONIC_NUMPY_NANARGMIN_HPP
#define PYTHONIC_NUMPY_NANARGMIN_HPP

#include "pythonic/include/numpy/nanargmin.hpp"

#include "pythonic/builtins/ValueError.hpp"
#include "pythonic/numpy/isnan.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace
  {
    template <class E, class F>
    void _nanargmin(E begin, E end, F &min, long &index, long &where, utils::int_<1>)
    {
      for (; begin != end; ++begin, ++index) {
        auto curr = *begin;
        if (!functor::isnan()(curr) && curr < min) {
          min = curr;
          where = index;
        }
      }
    }

    template <class E, class F, size_t N>
    void _nanargmin(E begin, E end, F &min, long &index, long &where, utils::int_<N>)
    {
      for (; begin != end; ++begin)
        _nanargmin((*begin).begin(), (*begin).end(), min, index, where, utils::int_<N - 1>());
    }
  } // namespace

  template <class E>
  long nanargmin(E const &expr)
  {
    typename E::dtype min = std::numeric_limits<typename E::dtype>::infinity();
    long where = -1;
    long index = 0;
    _nanargmin(expr.begin(), expr.end(), min, index, where, utils::int_<E::value>());
    if (where >= 0)
      return where;
    else
      throw types::ValueError("empty sequence");
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
