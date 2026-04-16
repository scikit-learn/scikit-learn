#ifndef PYTHONIC_NUMPY_VDOT_HPP
#define PYTHONIC_NUMPY_VDOT_HPP

#include "pythonic/include/numpy/vdot.hpp"

#include "pythonic/numpy/asarray.hpp"
#include "pythonic/numpy/conjugate.hpp"
#include "pythonic/numpy/dot.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class U, class V>
  auto vdot(U const &u, V const &v)
      -> decltype(functor::dot{}(functor::asarray{}(u).flat(), functor::asarray{}(v).flat()))
  {
    if (types::is_complex<typename U::dtype>::value &&
        types::is_complex<typename V::dtype>::value) {
      return functor::dot{}(functor::asarray{}(functor::conjugate{}(u)).flat(),
                            functor::asarray{}(v).flat());
    } else {
      return functor::dot{}(functor::asarray{}(u).flat(), functor::asarray{}(v).flat());
    }
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
