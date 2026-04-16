#ifndef PYTHONIC_NUMPY_MEDIAN_HPP
#define PYTHONIC_NUMPY_MEDIAN_HPP

#include "pythonic/include/numpy/median.hpp"

#include "pythonic/numpy/asarray.hpp"
#include "pythonic/numpy/sort.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/allocate.hpp"
#include "pythonic/utils/functor.hpp"
#include <algorithm>

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace
  {

    template <class T_out, class T, class pS>
    std::enable_if_t<std::tuple_size<pS>::value != 1, void>
    _median(T_out *out, types::ndarray<T, pS> const &tmp, long axis)
    {
      auto tmp_shape = sutils::getshape(tmp);
      const long step =
          std::accumulate(tmp_shape.begin() + axis, tmp_shape.end(), 1L, std::multiplies<long>());
      long const buffer_size = tmp_shape[axis];
      T *buffer = utils::allocate<T>(buffer_size);
      const long stepper = step / tmp_shape[axis];
      const long n = tmp.flat_size() / tmp_shape[axis] * step;
      long ith = 0, nth = 0;
      for (long i = 0; i < n; i += step) {
        T *buffer_iter = buffer;
        T const *iter = tmp.buffer + ith;
        T const *iend = iter + step;
        while (iter != iend) {
          *buffer_iter++ = *iter;
          iter += stepper;
        }
        if (buffer_size % 2 == 1) {
          std::nth_element(buffer, buffer + buffer_size / 2, buffer_iter, ndarray::comparator<T>{});
          *out++ = buffer[buffer_size / 2];
        } else {
          std::nth_element(buffer, buffer + buffer_size / 2, buffer_iter, ndarray::comparator<T>{});
          auto t0 = buffer[buffer_size / 2];
          std::nth_element(buffer, buffer + buffer_size / 2 - 1, buffer + buffer_size / 2,
                           ndarray::comparator<T>{});
          auto t1 = buffer[buffer_size / 2 - 1];
          *out++ = (t0 + t1) / double(2);
        }
        ith += step;
        if (ith - nth == tmp.flat_size()) {
          ++nth;
          ith = nth;
        }
      }
      utils::deallocate(buffer);
    }
  } // namespace

  template <class T, class pS>
  decltype(std::declval<T>() + 1.) median(types::ndarray<T, pS> const &arr, types::none_type)
  {
    size_t n = arr.flat_size();
    T *tmp = utils::allocate<T>(n);
    std::copy(arr.buffer, arr.buffer + n, tmp);
    std::nth_element(tmp, tmp + n / 2, tmp + n, ndarray::comparator<T>{});
    T t0 = tmp[n / 2];
    if (n % 2 == 1) {
      utils::deallocate(tmp);
      return t0;
    } else {
      std::nth_element(tmp, tmp + n / 2 - 1, tmp + n / 2, ndarray::comparator<T>{});
      T t1 = tmp[n / 2 - 1];
      utils::deallocate(tmp);
      return (t0 + t1) / 2.;
    }
  }

  template <class T, class pS>
  std::enable_if_t<std::tuple_size<pS>::value != 1,
                   types::ndarray<decltype(std::declval<T>() + 1.),
                                  types::array_tuple<long, std::tuple_size<pS>::value - 1>>>
  median(types::ndarray<T, pS> const &arr, long axis)
  {
    constexpr auto N = std::tuple_size<pS>::value;
    if (axis < 0)
      axis += N;

    types::array_tuple<long, std::tuple_size<pS>::value - 1> shp;
    auto stmp = sutils::getshape(arr);
    auto next = std::copy(stmp.begin(), stmp.begin() + axis, shp.begin());
    std::copy(stmp.begin() + axis + 1, stmp.end(), next);

    types::ndarray<decltype(std::declval<T>() + 1.),
                   types::array_tuple<long, std::tuple_size<pS>::value - 1>>
        out(shp, types::none_type{});
    _median(out.buffer, arr, axis);
    return out;
  }

  template <class T, class pS>
  std::enable_if_t<std::tuple_size<pS>::value == 1, decltype(std::declval<T>() + 1.)>
  median(types::ndarray<T, pS> const &arr, long axis)
  {
    if (axis != 0)
      throw types::ValueError("axis out of bounds");
    return median(arr);
  }

  NUMPY_EXPR_TO_NDARRAY0_IMPL(median);
} // namespace numpy
PYTHONIC_NS_END

#endif
