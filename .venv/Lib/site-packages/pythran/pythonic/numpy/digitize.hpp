#ifndef PYTHONIC_NUMPY_DIGITIZE_HPP
#define PYTHONIC_NUMPY_DIGITIZE_HPP

#include "pythonic/include/numpy/digitize.hpp"

#include "pythonic/builtins/None.hpp"
#include "pythonic/numpy/asarray.hpp"
#include "pythonic/operator_/gt.hpp"
#include "pythonic/operator_/lt.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace
  {
    template <class I, class O, class B, class Op>
    void _digitize(I begin, I end, O &out, B &bins, Op const &op, utils::int_<1>)
    {
      for (; begin != end; ++begin, ++out)
        *out = std::lower_bound(bins.begin(), bins.end(), *begin, op) - bins.begin();
    }

    template <class I, class O, class B, class Op, size_t N>
    void _digitize(I begin, I end, O &out, B &bins, Op const &op, utils::int_<N>)
    {
      for (; begin != end; ++begin)
        _digitize((*begin).begin(), (*begin).end(), out, bins, op, utils::int_<N - 1>());
    }
  } // namespace

  template <class E, class F>
  types::ndarray<long, types::pshape<long>> digitize(E const &expr, F const &b)
  {
    auto bins = asarray(b);
    bool is_increasing = bins.flat_size() > 1 && *bins.fbegin() < *(bins.fbegin() + 1);
    types::ndarray<long, types::pshape<long>> out(types::make_tuple(long(expr.flat_size())),
                                                  builtins::None);
    auto out_iter = out.fbegin();
    if (is_increasing)
      _digitize(expr.begin(), expr.end(), out_iter, bins, operator_::functor::lt(),
                utils::int_<E::value>());
    else
      _digitize(expr.begin(), expr.end(), out_iter, bins, operator_::functor::gt(),
                utils::int_<E::value>());
    return out;
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
