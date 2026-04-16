#ifndef PYTHONIC_NUMPY_PUT_HPP
#define PYTHONIC_NUMPY_PUT_HPP

#include "pythonic/include/numpy/put.hpp"

#include "pythonic/builtins/ValueError.hpp"
#include "pythonic/numpy/asarray.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_conversion.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class F, class T, class pS, class E>
  std::enable_if_t<types::is_numexpr_arg<F>::value, types::none_type>
  put(types::ndarray<T, pS> &expr, F const &ind, E const &v)
  {
    auto vind = asarray(ind);
    auto vv = asarray(v);
    for (long i = 0; i < ind.flat_size(); ++i) {
      auto val = *(vind.fbegin() + i);
      if (val >= expr.flat_size() || val < 0)
        throw types::ValueError("indice out of bound");
      *(expr.fbegin() + val) = *(vv.fbegin() + i % vv.flat_size());
    }
    return builtins::None;
  }

  template <class T, class pS>
  types::none_type put(types::ndarray<T, pS> &expr, long ind, T const &v)
  {
    if (ind >= expr.flat_size() || ind < 0)
      throw types::ValueError("indice out of bound");
    *(expr.fbegin() + ind) = v;
    return builtins::None;
  }

  template <class E, class M, class V>
  types::none_type put(E &, M const &, V const &)
  {
    throw std::runtime_error("put only partially implemented");
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
