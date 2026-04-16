#ifndef PYTHONIC_INCLUDE_NUMPY_FFT_C2C_HPP
#define PYTHONIC_INCLUDE_NUMPY_FFT_C2C_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace fft
  {

    template <class T, class pS>
    types::ndarray<std::complex<T>, types::array_tuple<long, std::tuple_size<pS>::value>>
    c2c(types::ndarray<std::complex<T>, pS> const &a, long n = -1, long axis = -1,
        types::str const &norm = {}, bool const forward = true);
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
