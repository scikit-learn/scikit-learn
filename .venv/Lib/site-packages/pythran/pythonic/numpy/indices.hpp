#ifndef PYTHONIC_NUMPY_INDICES_HPP
#define PYTHONIC_NUMPY_INDICES_HPP

#include "pythonic/include/numpy/indices.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class pS, class dtype>
  types::ndarray<typename dtype::type,
                 sutils::push_front_t<pS, std::integral_constant<long, std::tuple_size<pS>::value>>>
  indices(pS const &shape, dtype)
  {
    auto constexpr N = std::tuple_size<pS>::value;
    sutils::push_front_t<pS, std::integral_constant<long, N>> oshape;
    sutils::scopy_shape<1, -1>(oshape, shape, std::make_index_sequence<N>());
    types::ndarray<typename dtype::type, sutils::push_front_t<pS, std::integral_constant<long, N>>>
        out(oshape, builtins::None);
    typename dtype::type *iters[N];
    for (size_t n = 0; n < N; ++n)
      iters[n] = out[n].buffer;
    size_t lens[N];
    lens[0] = out.flat_size() / std::get<0>(shape);

    auto ashape = sutils::array(shape);

    for (size_t n = 1; n < N; ++n)
      lens[n] = lens[n - 1] / ashape[n];
    for (long i = 0, n = out.flat_size() / N; i < n; ++i) {
      long mult = 1;
      for (long n = N - 1; n > 0; n--) {
        *(iters[n]++) = (i / mult) % ashape[n];
        mult *= ashape[n];
      }
      *(iters[0]++) = i / mult;
    }
    return out;
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
