#ifndef PYTHONIC_NUMPY_DIFF_HPP
#define PYTHONIC_NUMPY_DIFF_HPP

#include "pythonic/include/numpy/diff.hpp"

#include "pythonic/numpy/asarray.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace details
  {
    template <class E>
    types::ndarray<typename E::dtype, types::array_tuple<long, E::value>> diff(E const &arr, long n,
                                                                               long axis)
    {
      auto shape = sutils::getshape(arr);
      auto stride = (axis == E::value - 1) ? arr.template shape<E::value - 1>()
                                           : std::accumulate(shape.begin() + axis + 1, shape.end(),
                                                             1L, std::multiplies<long>());
      --shape[axis];

      // this does not leak, but uses slightly too much memory
      auto out = arr.reshape(shape);

      auto iter = arr.fbegin();
      auto out_iter = out.fbegin();

      if (axis == E::value - 1) {
        for (long i = 0, sz = arr.flat_size(); i < sz; i += stride) {
          auto prev = *(iter + i);
          for (long k = 1; k < stride; ++k, ++out_iter) {
            auto nprev = *(iter + i + k);
            *(out_iter) = nprev - prev;
            prev = nprev;
          }
        }
      } else {
        iter += stride;
        for (auto out_end = out.fend(); out_iter != out_end; ++out_iter) {
          *out_iter = *iter++ - *out_iter;
        }
      }
      if (n == 1)
        return out;
      else
        return diff(out, n - 1, axis);
    }
  } // namespace details
  template <class E>
  types::ndarray<typename E::dtype, types::array_tuple<long, E::value>> diff(E const &expr, long n,
                                                                             long axis)
  {
    if (axis < 0)
      axis += E::value;
    // that's the only allocation that should happen
    return details::diff(array(expr), n, axis);
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
