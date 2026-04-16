#ifndef PYTHONIC_NUMPY_FLATNONZERO_HPP
#define PYTHONIC_NUMPY_FLATNONZERO_HPP

#include "pythonic/include/numpy/flatnonzero.hpp"

#include "pythonic/numpy/asarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace
  {
    template <class I, class O>
    void _flatnonzero(I begin, I end, O &out, long &i, utils::int_<1>)
    {
      for (; begin != end; ++begin, ++i)
        if (*begin)
          *out++ = i;
    }

    template <class I, class O, size_t N>
    void _flatnonzero(I begin, I end, O &out, long &i, utils::int_<N>)
    {
      for (; begin != end; ++begin)
        _flatnonzero((*begin).begin(), (*begin).end(), out, i, utils::int_<N - 1>());
    }
  } // namespace
  template <class E>
  types::ndarray<long, types::pshape<long>> flatnonzero(E const &expr)
  {
    long n = expr.flat_size();
    utils::shared_ref<types::raw_array<long>> buffer(n);
    long *iter = buffer->data;
    long i = 0;
    _flatnonzero(expr.begin(), expr.end(), iter, i, utils::int_<E::value>());
    types::pshape<long> shape = iter - buffer->data;
    return types::ndarray<long, types::pshape<long>>(std::move(buffer), shape);
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
