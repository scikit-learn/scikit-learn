#ifndef PYTHONIC_NUMPY_FREXP_HPP
#define PYTHONIC_NUMPY_FREXP_HPP

#include "pythonic/include/numpy/frexp.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/types/traits.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_conversion.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T>
  std::enable_if_t<std::is_scalar<T>::value, std::tuple<T, int>> frexp(T val)
  {
    int exp;
    T significand = std::frexp(val, &exp);
    return std::make_tuple(significand, exp);
  }

  namespace
  {
    template <class E, class F, class G>
    void _frexp(E begin, E end, F significands_iter, G exps_iter, utils::int_<1>)
    {
      for (; begin != end; ++begin, ++significands_iter, ++exps_iter)
        *significands_iter = std::frexp(*begin, exps_iter);
    }

    template <class E, class F, class G, size_t N>
    void _frexp(E begin, E end, F significands_iter, G exps_iter, utils::int_<N>)
    {
      for (; begin != end; ++begin, ++significands_iter, ++exps_iter)
        _frexp((*begin).begin(), (*begin).end(), (*significands_iter).begin(), (*exps_iter).begin(),
               utils::int_<N - 1>());
    }
  } // namespace

  template <class E>
  std::enable_if_t<!types::is_dtype<E>::value,
                   std::tuple<types::ndarray<typename E::dtype, typename E::shape_t>,
                              types::ndarray<int, typename E::shape_t>>>
  frexp(E const &arr)
  {
    auto arr_shape = sutils::getshape(arr);
    types::ndarray<typename E::dtype, typename E::shape_t> significands(arr_shape, builtins::None);
    types::ndarray<int, typename E::shape_t> exps(arr_shape, builtins::None);
    _frexp(arr.begin(), arr.end(), significands.begin(), exps.begin(), utils::int_<E::value>());
    return std::make_tuple(significands, exps);
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
