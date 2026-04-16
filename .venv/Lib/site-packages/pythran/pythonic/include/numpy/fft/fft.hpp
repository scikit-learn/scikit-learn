#ifndef PYTHONIC_INCLUDE_NUMPY_FFT_FFT_HPP
#define PYTHONIC_INCLUDE_NUMPY_FFT_FFT_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

/**
 * **Noteable difference to numpy.fft.fft:**
 * In contrast to numpy.fft.fft this implementation preserves precision
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

namespace numpy
{
  namespace fft
  {

    template <class T, class pS, class N = types::none_type, class Norm = types::none_type>
    types::ndarray<std::enable_if_t<types::is_complex<T>::value, T>,
                   types::array_tuple<long, std::tuple_size<pS>::value>>
    fft(types::ndarray<T, pS> const &a, N const &n = {}, long axis = -1, Norm const &norm = {});

    template <class T, class pS, class N = types::none_type, class Norm = types::none_type>
    types::ndarray<std::enable_if_t<std::is_floating_point<T>::value, std::complex<T>>,
                   types::array_tuple<long, std::tuple_size<pS>::value>>
    fft(types::ndarray<T, pS> const &a, N const &n = {}, long axis = -1, Norm const &norm = {});

    template <class T, class pS, class N = types::none_type, class Norm = types::none_type>
    types::ndarray<std::enable_if_t<std::is_integral<T>::value, std::complex<double>>,
                   types::array_tuple<long, std::tuple_size<pS>::value>>
    fft(types::ndarray<T, pS> const &a, N const &n = {}, long axis = -1, Norm const &norm = {});

    NUMPY_EXPR_TO_NDARRAY0_DECL(fft);
    DEFINE_FUNCTOR(pythonic::numpy::fft, fft);
  } // namespace fft
} // namespace numpy
PYTHONIC_NS_END

#endif
