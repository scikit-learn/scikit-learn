#ifndef PYTHONIC_INCLUDE_NUMPY_FFT_FFTN_HPP
#define PYTHONIC_INCLUDE_NUMPY_FFT_FFTN_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace fft
  {
    // without shape
    template <class T, class pS, class Axes = types::none_type, class Norm = types::none_type>
    types::ndarray<std::enable_if_t<std::is_integral<T>::value, std::complex<double>>,
                   types::array_tuple<long, std::tuple_size<pS>::value>>
    fftn(types::ndarray<T, pS> const &a, types::none_type s = {}, Axes const &axes = {},
         Norm const &norm = {});

    template <class T, class pS, class Axes = types::none_type, class Norm = types::none_type>
    types::ndarray<std::enable_if_t<std::is_floating_point<T>::value, std::complex<T>>,
                   types::array_tuple<long, std::tuple_size<pS>::value>>
    fftn(types::ndarray<T, pS> const &a, types::none_type s = {}, Axes const &axes = {},
         Norm const &norm = {});

    template <class T, class pS, class Axes = types::none_type, class Norm = types::none_type>
    types::ndarray<std::enable_if_t<types::is_complex<T>::value, T>,
                   types::array_tuple<long, std::tuple_size<pS>::value>>
    fftn(types::ndarray<T, pS> const &a, types::none_type s = {}, Axes const &axes = {},
         Norm const &norm = {});

    // with shape
    template <class T, class pS, class I, size_t N, class V, class Axes = types::none_type,
              class Norm = types::none_type>
    types::ndarray<std::enable_if_t<std::is_integral<T>::value, std::complex<double>>,
                   types::array_tuple<long, std::tuple_size<pS>::value>>
    fftn(types::ndarray<T, pS> const &a, types::array_base<I, N, V> const &s, Axes const &axes = {},
         Norm const &norm = {});

    template <class T, class pS, class I, size_t N, class V, class Axes = types::none_type,
              class Norm = types::none_type>
    types::ndarray<std::enable_if_t<std::is_floating_point<T>::value, std::complex<T>>,
                   types::array_tuple<long, std::tuple_size<pS>::value>>
    fftn(types::ndarray<T, pS> const &a, types::array_base<I, N, V> const &s, Axes const &axes = {},
         Norm const &norm = {});

    template <class T, class pS, class I, size_t N, class V, class Axes = types::none_type,
              class Norm = types::none_type>
    types::ndarray<std::enable_if_t<types::is_complex<T>::value, T>,
                   types::array_tuple<long, std::tuple_size<pS>::value>>
    fftn(types::ndarray<T, pS> const &a, types::array_base<I, N, V> const &s, Axes const &axes = {},
         Norm const &norm = {});

    NUMPY_EXPR_TO_NDARRAY0_DECL(fftn);
    DEFINE_FUNCTOR(pythonic::numpy::fft, fftn);
  } // namespace fft
} // namespace numpy
PYTHONIC_NS_END

#endif
