#ifndef PYTHONIC_NUMPY_ROLL_HPP
#define PYTHONIC_NUMPY_ROLL_HPP

#include "pythonic/include/numpy/roll.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_conversion.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T, class pS>
  types::ndarray<T, pS> roll(types::ndarray<T, pS> const &expr, long shift)
  {
    long expr_fsize = expr.flat_size();
    if (expr_fsize == 0)
      return expr.copy();
    else if (shift < 0)
      shift += expr_fsize;
    else if (shift >= expr_fsize)
      shift %= expr_fsize;

    types::ndarray<T, pS> out(expr._shape, builtins::None);
    std::copy(expr.fbegin(), expr.fend() - shift,
              std::copy(expr.fend() - shift, expr.fend(), out.fbegin()));
    return out;
  }

  namespace
  {
    template <class To, class From, size_t N>
    To _roll(To to, From from, long shift, long axis, types::array_tuple<long, N> const &shape,
             utils::int_<N - 1>)
    {
      long dim = shape[N - 1];
      if (axis == N - 1) {
        const From split = from + (dim - shift);
        to = std::copy(split, from + dim, to);
        return std::copy(from, split, to);
      } else {
        return std::copy(from, from + dim, to);
      }
    }

    template <class To, class From, size_t N, size_t M>
    std::enable_if_t<M != N - 1, To> _roll(To to, From from, long shift, long axis,
                                           types::array_tuple<long, N> const &shape, utils::int_<M>)
    {
      long dim = shape[M];
      long offset =
          std::accumulate(shape.begin() + M + 1, shape.end(), 1L, std::multiplies<long>());
      if (axis == M) {
        const From split = from + (dim - shift) * offset;
        for (From iter = split, end = from + dim * offset; iter != end; iter += offset)
          to = _roll(to, iter, shift, axis, shape, utils::int_<M + 1>());
        for (From iter = from, end = split; iter != end; iter += offset)
          to = _roll(to, iter, shift, axis, shape, utils::int_<M + 1>());
      } else {
        for (From iter = from, end = from + dim * offset; iter != end; iter += offset)
          to = _roll(to, iter, shift, axis, shape, utils::int_<M + 1>());
      }
      return to;
    }
  } // namespace

  template <class T, class pS>
  types::ndarray<T, pS> roll(types::ndarray<T, pS> const &expr, long shift, long axis)
  {
    auto expr_shape = sutils::array(expr._shape);
    if (expr_shape[axis] == 0)
      return expr.copy();
    if (shift < 0)
      shift += expr_shape[axis];
    else if (shift >= expr_shape[axis])
      shift %= expr_shape[axis];

    types::ndarray<T, pS> out(expr._shape, builtins::None);
    _roll(out.fbegin(), expr.fbegin(), shift, axis, expr_shape, utils::int_<0>());
    return out;
  }

  namespace
  {
    template <class To, class From, size_t N>
    To _rolls(To to, From from, long shifts[N], types::array_tuple<long, N> const &shape,
              utils::int_<N - 1>)
    {
      long dim = shape[N - 1];
      if (long shift = shifts[N - 1]) {
        From from_split = from + (dim - shift);
        To next = std::copy(from_split, from + dim, to);
        return std::copy(from, from_split, next);
      } else {
        return std::copy(from, from + dim, to);
      }
    }

    template <class To, class From, size_t N, size_t M>
    std::enable_if_t<M != N - 1, To> _rolls(To to, From from, long shifts[N],
                                            types::array_tuple<long, N> const &shape,
                                            utils::int_<M>)
    {
      long dim = shape[M];
      long offset =
          std::accumulate(shape.begin() + M + 1, shape.end(), 1L, std::multiplies<long>());
      if (long shift = shifts[M]) {
        const From split = from + (dim - shift) * offset;
        for (From iter = split, end = from + dim * offset; iter != end; iter += offset)
          to = _rolls(to, iter, shifts, shape, utils::int_<M + 1>());
        for (From iter = from; iter != split; iter += offset)
          to = _rolls(to, iter, shifts, shape, utils::int_<M + 1>());
      } else {
        for (From iter = from, end = from + dim * offset; iter != end; iter += offset)
          to = _rolls(to, iter, shifts, shape, utils::int_<M + 1>());
      }
      return to;
    }
  } // namespace

  template <class T, class pS, size_t N>
  types::ndarray<T, pS> roll(types::ndarray<T, pS> const &expr, types::array_tuple<long, N> shifts,
                             types::array_tuple<long, N> axes)
  {
    constexpr long ndim = types::ndarray<T, pS>::value;
    auto expr_shape = sutils::array(expr._shape);
    long axes_shifts[ndim] = {0};
    for (size_t i = 0; i < N; ++i)
      axes_shifts[axes[i]] += shifts[i];

    for (size_t i = 0; i < N; ++i) {
      if (axes_shifts[i] < 0)
        axes_shifts[i] += expr_shape[i];
      else if (axes_shifts[i] >= expr_shape[i])
        axes_shifts[i] %= expr_shape[i];
    }

    types::ndarray<T, pS> out(expr._shape, builtins::None);
    _rolls(out.fbegin(), expr.fbegin(), axes_shifts, expr_shape, utils::int_<0>());
    return out;
  }

  NUMPY_EXPR_TO_NDARRAY0_IMPL(roll);
} // namespace numpy
PYTHONIC_NS_END

#endif
