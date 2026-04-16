#ifndef PYTHONIC_NUMPY_FFT_IHFFT_HPP
#define PYTHONIC_NUMPY_FFT_IHFFT_HPP
#include "pythonic/builtins/None.hpp"
#include "pythonic/include/numpy/fft/ihfft.hpp"
#include "pythonic/include/utils/array_helper.hpp"
#include "pythonic/numpy/fft/c2c.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace fft
  {

    template <class T, class pS>
    types::ndarray<std::enable_if_t<std::is_floating_point<T>::value, std::complex<T>>,
                   types::array_tuple<long, std::tuple_size<pS>::value>>
    ihfft(types::ndarray<T, pS> const &in_array, types::none_type n, long axis,
          types::str const &norm)
    {
      return r2c(in_array, -1, axis, norm, false, false);
    }

    template <class T, class pS>
    types::ndarray<std::enable_if_t<std::is_floating_point<T>::value, std::complex<T>>,
                   types::array_tuple<long, std::tuple_size<pS>::value>>
    ihfft(types::ndarray<T, pS> const &in_array, types::none_type n, long axis,
          types::none_type norm)
    {
      return r2c(in_array, -1, axis, "", false, false);
    }

    template <class T, class pS>
    types::ndarray<std::enable_if_t<std::is_floating_point<T>::value, std::complex<T>>,
                   types::array_tuple<long, std::tuple_size<pS>::value>>
    ihfft(types::ndarray<T, pS> const &in_array, long n, long axis, types::none_type norm)
    {
      return r2c(in_array, n, axis, "", false, false);
    }

    template <class T, class pS>
    types::ndarray<std::enable_if_t<std::is_floating_point<T>::value, std::complex<T>>,
                   types::array_tuple<long, std::tuple_size<pS>::value>>
    ihfft(types::ndarray<T, pS> const &in_array, long n, long axis, types::str const &norm)
    {
      return r2c(in_array, n, axis, norm, false, false);
    }

    template <class T, class pS>
    types::ndarray<std::enable_if_t<std::is_integral<T>::value, std::complex<double>>,
                   types::array_tuple<long, std::tuple_size<pS>::value>>
    ihfft(types::ndarray<T, pS> const &in_array, types::none_type n, long axis,
          types::str const &norm)
    {
      auto tmp_array = _copy_to_double(in_array);
      return r2c(tmp_array, -1, axis, norm, false, false);
    }

    template <class T, class pS>
    types::ndarray<std::enable_if_t<std::is_integral<T>::value, std::complex<double>>,
                   types::array_tuple<long, std::tuple_size<pS>::value>>
    ihfft(types::ndarray<T, pS> const &in_array, types::none_type n, long axis,
          types::none_type norm)
    {
      auto tmp_array = _copy_to_double(in_array);
      return r2c(tmp_array, -1, axis, "", false, false);
    }

    template <class T, class pS>
    types::ndarray<std::enable_if_t<std::is_integral<T>::value, std::complex<double>>,
                   types::array_tuple<long, std::tuple_size<pS>::value>>
    ihfft(types::ndarray<T, pS> const &in_array, long n, long axis, types::none_type norm)
    {
      auto tmp_array = _copy_to_double(in_array);
      return r2c(tmp_array, n, axis, "", false, false);
    }

    template <class T, class pS>
    types::ndarray<std::enable_if_t<std::is_integral<T>::value, std::complex<double>>,
                   types::array_tuple<long, std::tuple_size<pS>::value>>
    ihfft(types::ndarray<T, pS> const &in_array, long n, long axis, types::str const &norm)
    {
      auto tmp_array = _copy_to_double(in_array);
      return r2c(tmp_array, n, axis, norm, false, false);
    }

    NUMPY_EXPR_TO_NDARRAY0_IMPL(ihfft);
  } // namespace fft
} // namespace numpy
PYTHONIC_NS_END

#endif
