#ifndef PYTHONIC_NUMPY_TRANSPOSE_HPP
#define PYTHONIC_NUMPY_TRANSPOSE_HPP

#include "pythonic/include/numpy/transpose.hpp"

#include "pythonic/builtins/ValueError.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/nested_container.hpp"
#include "pythonic/utils/numpy_conversion.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace
  {
    template <class T, class pS, class O, class Indices, class S, class Perm>
    O const *_transposer(types::ndarray<T, pS> &expr, O const *iter, Indices &indices,
                         S const &shape, Perm const &perm,
                         utils::int_<std::tuple_size<pS>::value - 1>)
    {
      for (long i = 0, n = shape[std::tuple_size<pS>::value - 1]; i < n; ++i) {
        indices[perm[std::tuple_size<pS>::value - 1]] = i;
        expr.fast(indices) = *iter++;
      }
      indices[perm[std::tuple_size<pS>::value - 1]] = 0;
      return iter;
    }

    template <class T, class pS, class O, class Indices, class S, class Perm, size_t I>
    std::enable_if_t<std::tuple_size<pS>::value - 1 != I, O const *>
    _transposer(types::ndarray<T, pS> &expr, O const *iter, Indices &indices, S const &shape,
                Perm const &perm, utils::int_<I>)
    {
      for (long i = 0, n = shape[I]; i < n; ++i) {
        indices[perm[I]] = i;
        iter = _transposer(expr, iter, indices, shape, perm, utils::int_<I + 1>());
      }
      indices[perm[I]] = 0;
      return iter;
    }
    template <class T, class pS>
    types::ndarray<T, types::array_tuple<long, std::tuple_size<pS>::value>>
    _transposer(types::ndarray<T, pS> const &a, long const l[std::tuple_size<pS>::value])
    {
      auto shape = sutils::getshape(a);
      types::array_tuple<long, std::tuple_size<pS>::value> shp;
      for (unsigned long i = 0; i < std::tuple_size<pS>::value; ++i)
        shp[i] = shape[l[i]];

      types::array_tuple<long, std::tuple_size<pS>::value> perm;
      for (std::size_t i = 0; i < std::tuple_size<pS>::value; ++i)
        perm[l[i]] = i;

      types::ndarray<T, types::array_tuple<long, std::tuple_size<pS>::value>> new_array(
          shp, builtins::None);

      auto const *iter = a.buffer;
      types::array_tuple<long, std::tuple_size<pS>::value> indices;
      _transposer(new_array, iter, indices, shape, perm, utils::int_<0>{});

      return new_array;
    }
  } // namespace

  template <class T, class pS>
  std::enable_if_t<(std::tuple_size<pS>::value > 2),
                   types::ndarray<T, types::array_tuple<long, std::tuple_size<pS>::value>>>
  transpose(types::ndarray<T, pS> const &a)
  {
    long t[std::tuple_size<pS>::value];
    for (unsigned long i = 0; i < std::tuple_size<pS>::value; ++i)
      t[std::tuple_size<pS>::value - 1 - i] = i;
    return _transposer(a, t);
  }

  template <class T, class pS, size_t M>
  types::ndarray<T, types::array_tuple<long, std::tuple_size<pS>::value>>
  transpose(types::ndarray<T, pS> const &a, types::array_tuple<long, M> const &t)
  {
    static_assert(std::tuple_size<pS>::value == M, "axes don't match array");

    long val = t[M - 1];
    if (val >= long(std::tuple_size<pS>::value))
      throw types::ValueError("invalid axis for this array");
    return _transposer(a, &t[0]);
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
