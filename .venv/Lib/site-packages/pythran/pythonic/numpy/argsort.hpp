#ifndef PYTHONIC_NUMPY_ARGSORT_HPP
#define PYTHONIC_NUMPY_ARGSORT_HPP

#include "pythonic/include/numpy/argsort.hpp"
#include "pythonic/numpy/ndarray/sort.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  types::ndarray<long, types::array_tuple<long, 1>> argsort(E const &expr, types::none_type,
                                                            types::none_type)
  {
    auto out = functor::array{}(expr).flat();
    return argsort(out);
  }

  template <class T, class pS, class Sorter>
  types::ndarray<long, pS> _argsort(types::ndarray<T, pS> const &a, long axis, Sorter sorter)
  {
    constexpr auto N = std::tuple_size<pS>::value;
    if (axis < 0)
      axis += N;

    long const flat_size = a.flat_size();
    types::ndarray<long, pS> indices(a._shape, builtins::None);
    if (axis == N - 1) {
      size_t step = a.template shape<N - 1>();

      auto a_base = a.fbegin();
      for (long *iter_indices = indices.buffer, *end_indices = indices.buffer + flat_size;
           iter_indices != end_indices; iter_indices += step, a_base += step) {
        // fill with the original indices
        std::iota(iter_indices, iter_indices + step, 0L);
        // sort the index using the value from a
        sorter(iter_indices, iter_indices + step,
               [a_base](long i1, long i2) { return a_base[i1] < a_base[i2]; });
      }
    } else {
      auto out_shape = sutils::getshape(a);
      const long step =
          std::accumulate(out_shape.begin() + axis, out_shape.end(), 1L, std::multiplies<long>());
      long const buffer_size = out_shape[axis];
      const long stepper = step / out_shape[axis];
      const long n = flat_size / out_shape[axis];
      long ith = 0, nth = 0;

      long *buffer = utils::allocate<long>(buffer_size);
      long *buffer_start = buffer, *buffer_end = buffer + buffer_size;
      std::iota(buffer_start, buffer_end, 0L);
      for (long i = 0; i < n; i++) {
        auto a_base = a.fbegin() + ith;
        sorter(buffer, buffer + buffer_size, [a_base, stepper](long i1, long i2) {
          return a_base[i1 * stepper] < a_base[i2 * stepper];
        });

        for (long j = 0; j < buffer_size; ++j)
          indices.buffer[ith + j * stepper] = buffer[j];

        ith = step;
        if (ith >= flat_size) {
          ith = ++nth;
        }
      }
      utils::deallocate(buffer);
    }
    return indices;
  }

  template <class T, class pS>
  types::ndarray<long, pS> argsort(types::ndarray<T, pS> const &a, long axis, types::none_type)
  {
    return _argsort(a, axis, ndarray::quicksorter());
  }

  template <class T, class pS>
  types::ndarray<long, pS> argsort(types::ndarray<T, pS> const &a, long axis,
                                   types::str const &kind)
  {
    if (kind == "mergesort")
      return _argsort(a, axis, ndarray::mergesorter());
    else if (kind == "heapsort")
      return _argsort(a, axis, ndarray::heapsorter());
    else if (kind == "stable")
      return _argsort(a, axis, ndarray::stablesorter());
    return _argsort(a, axis, ndarray::quicksorter());
  }

  NUMPY_EXPR_TO_NDARRAY0_IMPL(argsort);
} // namespace numpy
PYTHONIC_NS_END

#endif
