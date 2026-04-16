#ifndef PYTHONIC_NUMPY_APPEND_HPP
#define PYTHONIC_NUMPY_APPEND_HPP

#include "pythonic/include/numpy/append.hpp"

#include "pythonic/numpy/asarray.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T, class pS, class F>
  std::enable_if_t<!types::is_dtype<F>::value,
                   types::ndarray<typename __combined<T, typename types::dtype_of<F>::type>::type,
                                  types::pshape<long>>>
  append(types::ndarray<T, pS> const &nto, F const &data)
  {
    auto ndata = numpy::functor::asarray{}(data);
    long nsize = nto.flat_size() + ndata.flat_size();
    types::ndarray<typename __combined<T, typename types::dtype_of<F>::type>::type,
                   types::pshape<long>>
        out(types::pshape<long>(nsize), builtins::None);
    auto out_back = std::copy(nto.fbegin(), nto.fend(), out.fbegin());
    std::copy(ndata.fbegin(), ndata.fend(), out_back);
    return out;
  }
  template <class T, class pS, class F>
  std::enable_if_t<types::is_dtype<F>::value,
                   types::ndarray<typename __combined<T, typename types::dtype_of<F>::type>::type,
                                  types::pshape<long>>>
  append(types::ndarray<T, pS> const &nto, F const &data)
  {
    long nsize = nto.flat_size() + 1;
    types::ndarray<typename __combined<T, typename types::dtype_of<F>::type>::type,
                   types::pshape<long>>
        out(types::pshape<long>(nsize), builtins::None);
    auto out_back = std::copy(nto.fbegin(), nto.fend(), out.fbegin());
    *out_back = data;
    return out;
  }

  template <class T, class F>
  types::ndarray<typename __combined<typename types::dtype_of<T>::type,
                                     typename types::dtype_of<F>::type>::type,
                 types::pshape<long>>
  append(T const &to, F const &data)
  {
    return append(numpy::functor::asarray{}(to), data);
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
