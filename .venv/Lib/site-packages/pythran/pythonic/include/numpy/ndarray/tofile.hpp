#ifndef PYTHONIC_INCLUDE_NUMPY_NDARRAY_TOFILE_HPP
#define PYTHONIC_INCLUDE_NUMPY_NDARRAY_TOFILE_HPP

#include "pythonic/types/ndarray.hpp"
#include "pythonic/types/str.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_conversion.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace ndarray
  {
    template <class T, class pS>
    void tofile(types::ndarray<T, pS> const &expr, types::str const &file_name,
                types::str const &sep = "", types::str const &format = "");

    NUMPY_EXPR_TO_NDARRAY0_DECL(tofile);
    DEFINE_FUNCTOR(pythonic::numpy::ndarray, tofile);
  } // namespace ndarray
} // namespace numpy
PYTHONIC_NS_END

#endif
