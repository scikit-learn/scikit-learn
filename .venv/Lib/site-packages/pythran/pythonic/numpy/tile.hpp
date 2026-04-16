#ifndef PYTHONIC_NUMPY_TILE_HPP
#define PYTHONIC_NUMPY_TILE_HPP

#include "pythonic/include/numpy/tile.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace
  {
    template <class I, class O>
    void _tile(I begin, I end, O &out, long rep, utils::int_<1>)
    {
      for (long i = 0; i < rep; ++i)
        out = std::copy(begin, end, out);
    }

    template <class I, class O, size_t N>
    void _tile(I begin, I end, O &out, long rep, utils::int_<N>)
    {
      for (; begin != end; ++begin)
        _tile((*begin).begin(), (*begin).end(), out, rep, utils::int_<N - 1>());
    }
  } // namespace

  template <class E>
  types::ndarray<typename E::dtype, types::array_tuple<long, E::value>> tile(E const &expr,
                                                                             long reps)
  {
    size_t n = expr.flat_size();
    types::ndarray<typename E::dtype, types::array_tuple<long, E::value>> out(
        types::array_tuple<long, 1>{{long(n * reps)}}, builtins::None);
    auto out_iter = out.fbegin();
    _tile(expr.begin(), expr.end(), out_iter, 1, utils::int_<E::value>());
    for (long i = 1; i < reps; ++i)
      out_iter = std::copy(out.fbegin(), out.fbegin() + n, out_iter);
    return out;
  }

  template <size_t Shift, class R, class S, size_t... Is>
  types::array_tuple<long, sizeof...(Is)> tile_init_shape(R const &reps, S const &expr_shape,
                                                          std::index_sequence<Is...>)
  {
    constexpr size_t M = S::value;
    return {{(reps[Is] * ((Is < Shift) ? 1 : expr_shape.template shape<(Is < M) ? Is : 0>()))...}};
  }

  template <class E, size_t N>
  types::ndarray<typename E::dtype, types::array_tuple<long, N>>
  tile(E const &expr, types::array_tuple<long, N> const &reps)
  {
    size_t n = expr.flat_size();
    types::array_tuple<long, N> shape =
        tile_init_shape<N - E::value>(reps, expr, std::make_index_sequence<N>());

    long last_rep = (E::value == N) ? std::get<N - 1>(reps) : 1;
    types::ndarray<typename E::dtype, types::array_tuple<long, N>> out(shape, builtins::None);
    auto out_iter = out.fbegin();
    _tile(expr.begin(), expr.end(), out_iter, last_rep, utils::int_<E::value>());

    size_t nreps = out.flat_size() / (n * last_rep);
    for (size_t i = 1; i < nreps; ++i)
      out_iter = std::copy(out.fbegin(), out.fbegin() + n, out_iter);
    return out;
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
