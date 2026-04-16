#ifndef PYTHONIC_INCLUDE_NUMPY_FFT_HFFT_HPP
#define PYTHONIC_INCLUDE_NUMPY_FFT_HFFT_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

/**
 * **Noteable difference to numpy.fft.hfft:**
 * In contrast to numpy.fft.hfft this implementation preserves precision
 * of floating point and complex inputs, i.e. complex<float> input yields
 * complex<float> output. numpy.fft.fft always returns complex<double>, even for
 * long double input. This follows the same reasoning as given by numpy compiled
 * with intel_mkl (see here: https://github.com/IntelPython/mkl_fft/issues/10).
 * Conversion to double precision causes code to be slower and hurts use cases
 * where single precision preservation is desired, e.g. when interacting with
 *GPUs
 * or instruments. Moreover for the case of long double inputs, this avoids
 * loss of precision.
 **/

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace fft
  {

    template <class T, class pS>
    types::ndarray<T, types::array_tuple<long, std::tuple_size<pS>::value>>
    hfft(types::ndarray<std::complex<T>, pS> const &a, long n = -1, long axis = -1,
         types::str const &norm = {});

    template <class T, class pS>
    types::ndarray<T, types::array_tuple<long, std::tuple_size<pS>::value>>
    hfft(types::ndarray<std::complex<T>, pS> const &a, types::none_type n, long axis,
         types::str const &norm);

    template <class T, class pS>
    types::ndarray<T, types::array_tuple<long, std::tuple_size<pS>::value>>
    hfft(types::ndarray<std::complex<T>, pS> const &a, long n, long axis, types::none_type norm);

    template <class T, class pS>
    types::ndarray<T, types::array_tuple<long, std::tuple_size<pS>::value>>
    hfft(types::ndarray<std::complex<T>, pS> const &a, types::none_type n, long axis = -1,
         types::none_type norm = types::none_type{});

    template <class T, class pS>
    types::ndarray<std::enable_if_t<!types::is_complex<T>::value,
                                    std::conditional_t<std::is_integral<T>::value, double, T>>,
                   types::array_tuple<long, std::tuple_size<pS>::value>>
    hfft(types::ndarray<T, pS> const &a, long n = -1, long axis = -1, types::str const &norm = {});

    template <class T, class pS>
    types::ndarray<std::enable_if_t<!types::is_complex<T>::value,
                                    std::conditional_t<std::is_integral<T>::value, double, T>>,
                   types::array_tuple<long, std::tuple_size<pS>::value>>
    hfft(types::ndarray<T, pS> const &a, types::none_type n, long axis, types::str const &norm);

    template <class T, class pS>
    types::ndarray<std::enable_if_t<!types::is_complex<T>::value,
                                    std::conditional_t<std::is_integral<T>::value, double, T>>,
                   types::array_tuple<long, std::tuple_size<pS>::value>>
    hfft(types::ndarray<T, pS> const &a, long n, long axis, types::none_type norm);

    template <class T, class pS>
    types::ndarray<std::enable_if_t<!types::is_complex<T>::value,
                                    std::conditional_t<std::is_integral<T>::value, double, T>>,
                   types::array_tuple<long, std::tuple_size<pS>::value>>
    hfft(types::ndarray<T, pS> const &a, types::none_type n, long axis = -1,
         types::none_type norm = types::none_type{});

    NUMPY_EXPR_TO_NDARRAY0_DECL(hfft);
    DEFINE_FUNCTOR(pythonic::numpy::fft, hfft);
  } // namespace fft
} // namespace numpy
PYTHONIC_NS_END

#endif
