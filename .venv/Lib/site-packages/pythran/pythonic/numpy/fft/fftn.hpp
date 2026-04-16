#ifndef PYTHONIC_NUMPY_FFT_FFTN_HPP
#define PYTHONIC_NUMPY_FFT_FFTN_HPP

#include "pythonic/builtins/None.hpp"
#include "pythonic/include/numpy/fft/fftn.hpp"
#include "pythonic/include/utils/array_helper.hpp"
#include "pythonic/numpy/fft/c2c.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace fft
  {
    namespace details
    {
      inline types::str normalize_norm(types::none_type const &)
      {
        return "backward";
      }
      template <class T>
      inline types::str normalize_norm(T const &norm)
      {
        return norm;
      }

      template <size_t K, size_t S>
          types::array_tuple < long,
          K<S ? K : S> normalize_axes(types::none_type const &)
      {
        if (S == 1)
          return {-1}; // FIXME: understand why this is needed
        types::array_tuple < long, K<S ? K : S> result;
        for (size_t i = 0; i < std::min(K, S); ++i)
          result[i] = (long)i;
        return result;
      }
      template <size_t K, size_t S, class T, size_t N, class V>
      types::array_base<T, N, V> const &normalize_axes(types::array_base<T, N, V> const &axes)
      {
        return axes;
      }
    } // namespace details

    // without shape

    template <class T, class pS, class Axes, class Norm>
    types::ndarray<std::enable_if_t<std::is_integral<T>::value, std::complex<double>>,
                   types::array_tuple<long, std::tuple_size<pS>::value>>
    fftn(types::ndarray<T, pS> const &in_array, types::none_type s, Axes const &axes,
         Norm const &norm)
    {
      auto tmp_array = _copy_to_double(in_array);
      return fftn(tmp_array, s, axes, norm);
    }

    template <class T, class pS, class Axes, class Norm>
    types::ndarray<std::enable_if_t<std::is_floating_point<T>::value, std::complex<T>>,
                   types::array_tuple<long, std::tuple_size<pS>::value>>
    fftn(types::ndarray<T, pS> const &in_array, types::none_type s, Axes const &axes,
         Norm const &norm)
    {
      auto naxes = details::normalize_axes<std::tuple_size<pS>::value, 1>(axes);
      auto nnorm = details::normalize_norm(norm);
      auto result = r2c(in_array, -1, naxes[0], nnorm.c_str(), true, true);
      for (size_t i = 1; i < naxes.size(); ++i)
        result = c2c(result, -1, naxes[i], nnorm.c_str(), true);
      return result;
    }

    template <class T, class pS, class Axes, class Norm>
    types::ndarray<std::enable_if_t<types::is_complex<T>::value, T>,
                   types::array_tuple<long, std::tuple_size<pS>::value>>
    fftn(types::ndarray<T, pS> const &in_array, types::none_type s, Axes const &axes,
         Norm const &norm)
    {
      auto naxes = details::normalize_axes<std::tuple_size<pS>::value, 1>(axes);
      auto nnorm = details::normalize_norm(norm);
      auto result = c2c(in_array, -1, naxes[0], nnorm.c_str(), true);
      for (size_t i = 1; i < naxes.size(); ++i)
        result = c2c(result, -1, naxes[i], nnorm.c_str(), true);
      return result;
    }

    // with shape
    template <class T, class pS, class I, size_t N, class V, class Axes, class Norm>
    types::ndarray<std::enable_if_t<std::is_integral<T>::value, std::complex<double>>,
                   types::array_tuple<long, std::tuple_size<pS>::value>>
    fftn(types::ndarray<T, pS> const &a, types::array_base<I, N, V> const &s, Axes const &axes,
         Norm const &norm)
    {
      auto tmp_array = _copy_to_double(a);
      return fftn(tmp_array, s, axes, norm);
    }

    template <class T, class pS, class I, size_t N, class V, class Axes, class Norm>
    types::ndarray<std::enable_if_t<std::is_floating_point<T>::value, std::complex<T>>,
                   types::array_tuple<long, std::tuple_size<pS>::value>>
    fftn(types::ndarray<T, pS> const &a, types::array_base<I, N, V> const &s, Axes const &axes,
         Norm const &norm)
    {
      auto nnorm = details::normalize_norm(norm);
      auto naxes = details::normalize_axes<std::tuple_size<pS>::value, N>(axes);
      size_t i = 0;
      auto out = r2c(a, s[i], naxes[i], nnorm.c_str(), true, true);
      for (++i; i < N; ++i) {
        out = c2c(out, s[i], naxes[i], nnorm.c_str(), true);
      }
      return out;
    }

    template <class T, class pS, class I, size_t N, class V, class Axes, class Norm>
    types::ndarray<std::enable_if_t<types::is_complex<T>::value, T>,
                   types::array_tuple<long, std::tuple_size<pS>::value>>
    fftn(types::ndarray<T, pS> const &a, types::array_base<I, N, V> const &s, Axes const &axes,
         Norm const &norm)
    {
      auto nnorm = details::normalize_norm(norm);
      auto naxes = details::normalize_axes<std::tuple_size<pS>::value, N>(axes);
      size_t i = 0;
      auto out = c2c(a, s[i], naxes[i], nnorm.c_str(), true);
      for (++i; i < N; ++i) {
        out = c2c(out, s[i], naxes[i], nnorm.c_str(), true);
      }
      return out;
    }

    NUMPY_EXPR_TO_NDARRAY0_IMPL(fftn);
  } // namespace fft
} // namespace numpy
PYTHONIC_NS_END

#endif
