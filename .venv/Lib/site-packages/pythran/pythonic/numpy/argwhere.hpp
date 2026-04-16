#ifndef PYTHONIC_NUMPY_ARGWHERE_HPP
#define PYTHONIC_NUMPY_ARGWHERE_HPP

#include "pythonic/include/numpy/argwhere.hpp"

#include "pythonic/numpy/asarray.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  typename types::ndarray<long, types::array_tuple<long, 2>> argwhere(E const &expr)
  {
    constexpr long N = E::value;
    auto arr = asarray(expr);
    long sz = arr.flat_size();
    auto eshape = sutils::getshape(arr);

    utils::shared_ref<types::raw_array<long>> buffer(sz * N); // too much memory used
    long *buffer_iter = buffer->data;

    long real_sz = 0;
    auto iter = arr.fbegin();
    for (long i = 0; i < sz; ++i, ++iter) {
      if (*iter) {
        ++real_sz;
        long mult = 1;
        for (long j = N - 1; j > 0; j--) {
          buffer_iter[j] = (i / mult) % eshape[j];
          mult *= eshape[j];
        }
        buffer_iter[0] = i / mult;
        buffer_iter += N;
      }
    }
    types::array_tuple<long, 2> shape = {real_sz, N};
    return {buffer, shape};
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
