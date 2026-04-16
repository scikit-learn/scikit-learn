#ifndef PYTHONIC_NUMPY_FFT_C2C_HPP
#define PYTHONIC_NUMPY_FFT_C2C_HPP

#include "pythonic/builtins/None.hpp"
#include "pythonic/include/numpy/fft/c2c.hpp"
#include "pythonic/include/utils/array_helper.hpp"
#include "pythonic/numpy/concatenate.hpp"
#include "pythonic/numpy/empty.hpp"
#include "pythonic/numpy/zeros.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

#include <array>
#include <cmath>
#include <cstring>

#include "pythonic/numpy/fft/pocketfft.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace fft
  {
    using pocketfft::shape_t;
    using pocketfft::stride_t;
    using ldbl_t = std::conditional_t<sizeof(long double) == sizeof(double), double, long double>;

    template <class T, class pS>
    types::ndarray<std::enable_if_t<std::is_integral<T>::value, double>,
                   types::array_tuple<long, std::tuple_size<pS>::value>>
    _copy_to_double(types::ndarray<T, pS> const &in_array)
    {
      auto out_shape = sutils::getshape(in_array);
      size_t l = in_array.flat_size();
      auto out_array = numpy::empty(out_shape, types::dtype_t<double>());
      std::copy(in_array.buffer, in_array.buffer + l, out_array.buffer);
      return out_array;
    }

    template <class T, class pS>
    types::ndarray<std::enable_if_t<std::is_floating_point<T>::value, std::complex<T>>,
                   types::array_tuple<long, std::tuple_size<pS>::value>>
    _copy_to_complex(types::ndarray<T, pS> const &in_array)
    {
      auto out_shape = sutils::getshape(in_array);
      size_t l = in_array.flat_size();
      auto out_array = numpy::empty(out_shape, types::dtype_t<typename std::complex<T>>());
      std::copy(in_array.buffer, in_array.buffer + l, out_array.buffer);
      return out_array;
    }

    template <class T, class pS>
    types::ndarray<std::enable_if_t<std::is_integral<T>::value, std::complex<double>>,
                   types::array_tuple<long, std::tuple_size<pS>::value>>
    _copy_to_complex(types::ndarray<T, pS> const &in_array)
    {
      auto out_shape = sutils::getshape(in_array);
      size_t l = in_array.flat_size();
      auto out_array = numpy::empty(out_shape, types::dtype_t<typename std::complex<double>>());
      std::copy(in_array.buffer, in_array.buffer + l, out_array.buffer);
      return out_array;
    }

    enum class Inorm : int {
      forward,
      explicit_forward,
      ortho,
      backward,
    };

    Inorm _get_inorm(types::str const &norm, bool forward)
    {
      if (norm == "ortho")
        return Inorm::ortho;
      if (norm == "forward")
        return Inorm::explicit_forward;
      return forward ? Inorm::forward : Inorm::backward;
    }

    template <typename T>
    T norm_fct(Inorm inorm, size_t N)
    {
      switch (inorm) {
      case Inorm::ortho:
        return T(1 / sqrt(ldbl_t(N)));
      case Inorm::backward:
        return T(1 / ldbl_t(N));
      case Inorm::explicit_forward:
        return T(1. / N);
      default:
        assert(false && "unreachable");
        return T(0);
      }
    }

    template <typename T>
    T norm_fct(Inorm inorm, const shape_t &shape, const shape_t &axes, size_t fct = 1,
               int delta = 0)
    {
      // Fast path
      if (inorm == Inorm::forward)
        return 1;

      size_t N(1);
      for (auto a : axes)
        N *= fct * size_t(int64_t(shape[a]) + delta);
      return norm_fct<T>(inorm, N);
    }

    template <class T, class pS>
    stride_t create_strides(types::ndarray<T, pS> const &in_array)
    {
      auto constexpr N = std::tuple_size<pS>::value;
      auto shape = sutils::getshape(in_array);
      stride_t strides = stride_t(N);
      strides[N - 1] = sizeof(T);
      std::transform(strides.rbegin(), strides.rend() - 1, shape.rbegin(), strides.rbegin() + 1,
                     std::multiplies<long>());
      return strides;
    }

    template <class T, class pS>
    types::ndarray<T, types::array_tuple<long, std::tuple_size<pS>::value>>
    _pad_in_array(types::ndarray<T, pS> const &in_array, long axis, long n)
    {
      types::ndarray<T, types::array_tuple<long, std::tuple_size<pS>::value>> extended_array;
      auto tmp_shape = sutils::getshape(in_array);
      tmp_shape[axis] = n;
      auto tmp_array = zeros(tmp_shape, types::dtype_t<T>());
      types::list<types::ndarray<T, types::array_tuple<long, std::tuple_size<pS>::value>>> bi(0);
      bi.push_back(in_array);
      bi.push_back(tmp_array);
      extended_array = concatenate(bi, axis);
      return extended_array;
    }

    template <class T, class pS>
    types::ndarray<T, types::array_tuple<long, std::tuple_size<pS>::value>>
    c2r(types::ndarray<std::complex<T>, pS> const &in_array, long n, long axis,
        types::str const &norm, bool forward)
    {
      auto constexpr N = std::tuple_size<pS>::value;
      Inorm inorm = _get_inorm(norm, forward);
      if (axis < 0)
        axis = N + axis;
      auto in_shape = sutils::getshape(in_array);
      long npts = in_shape[axis];
      if (n == -1)
        n = 2 * npts - 2;
      auto out_shape = sutils::getshape(in_array);
      out_shape[axis] = n;
      // Create output array.
      types::ndarray<T, types::array_tuple<long, std::tuple_size<pS>::value>> out_array(
          out_shape, builtins::None);
      std::complex<T> *d_in;
      types::ndarray<std::complex<T>, types::array_tuple<long, std::tuple_size<pS>::value>>
          extended_array;
      stride_t in_strides;
      auto out_strides = create_strides(out_array);
      if (n > 2 * npts - 2) {
        // extend array with zeros along axis direction
        extended_array = _pad_in_array(in_array, axis, n - (2 * npts - 2));
        in_strides = create_strides(extended_array);
        d_in = reinterpret_cast<std::complex<T> *>(extended_array.buffer);
      } else {
        in_strides =
            create_strides(in_array); // for cropped arrays we need to use different strides
        d_in = reinterpret_cast<std::complex<T> *>(in_array.buffer);
      }
      auto d_out = reinterpret_cast<T *>(out_array.buffer);
      // axes calculation is for 1D transform
      shape_t axes = shape_t(1);
      axes[0] = axis;
      shape_t shapes = shape_t(size_t(N));
      std::copy(out_shape.begin(), out_shape.begin() + N, shapes.begin());
      auto fct = norm_fct<T>(inorm, shapes, axes);
      pocketfft::c2r(shapes, in_strides, out_strides, axes, forward, d_in, d_out, fct, size_t(0));
      return out_array;
    }

    template <class T, class pS>
    types::ndarray<std::complex<T>, types::array_tuple<long, std::tuple_size<pS>::value>>
    c2c(types::ndarray<std::complex<T>, pS> const &in_array, long n, long axis,
        types::str const &norm, bool forward)
    {
      auto constexpr N = std::tuple_size<pS>::value;
      Inorm inorm = _get_inorm(norm, forward);
      if (axis < 0)
        axis = N + axis;
      auto in_shape = sutils::getshape(in_array);
      long npts = in_shape[axis];
      if (n == -1)
        n = npts;
      auto out_shape = sutils::getshape(in_array);
      out_shape[axis] = n;
      // Create output array.
      types::ndarray<std::complex<T>, types::array_tuple<long, std::tuple_size<pS>::value>>
          out_array(out_shape, builtins::None);
      std::complex<T> *d_in;
      types::ndarray<std::complex<T>, types::array_tuple<long, std::tuple_size<pS>::value>>
          extended_array;
      stride_t in_strides;
      if (n > npts) {
        // extend array with zeros along axis direction
        extended_array = _pad_in_array(in_array, axis, n - npts);
        d_in = reinterpret_cast<std::complex<T> *>(extended_array.buffer);
        in_strides = create_strides(extended_array); //
      } else {
        d_in = reinterpret_cast<std::complex<T> *>(in_array.buffer);
        in_strides =
            create_strides(in_array); // for cropped arrays we need to use different strides
      }
      auto d_out = reinterpret_cast<std::complex<T> *>(out_array.buffer);
      // axes calculation is for 1D transform
      shape_t axes = shape_t(1);
      axes[0] = axis;
      auto out_strides = create_strides(out_array);
      shape_t shapes = shape_t(size_t(N));
      for (size_t i = 0; i < N; ++i)
        shapes[i] = size_t(out_shape[i]);
      auto fct = norm_fct<T>(inorm, shapes, axes);
      pocketfft::c2c(shapes, in_strides, out_strides, axes, forward, d_in, d_out, fct, size_t(0));
      return out_array;
    }

    template <class T, class pS>
    types::ndarray<std::enable_if_t<std::is_floating_point<T>::value, std::complex<T>>,
                   types::array_tuple<long, std::tuple_size<pS>::value>>
    r2c(types::ndarray<T, pS> const &in_array, long n, long axis, types::str const &norm,
        bool forward, bool extend = true)
    {
      auto constexpr N = std::tuple_size<pS>::value;
      Inorm inorm = _get_inorm(norm, forward);
      if (axis < 0)
        axis = N + axis;
      auto in_shape = sutils::getshape(in_array);
      long npts = in_shape[axis];
      if (n == -1)
        n = npts;
      auto out_shape = sutils::getshape(in_array);
      if (extend) {
        out_shape[axis] = n;
      } else {
        out_shape[axis] = n / 2 + 1;
      }
      // Create output array.
      types::ndarray<std::complex<T>, types::array_tuple<long, std::tuple_size<pS>::value>>
          out_array(out_shape, builtins::None);
      T *d_in;
      types::ndarray<T, types::array_tuple<long, std::tuple_size<pS>::value>> extended_array;
      shape_t shapes = shape_t(size_t(N));
      stride_t in_strides;
      if (n > npts) {
        // extend array with zeros along axis direction
        extended_array = _pad_in_array(in_array, axis, n - npts);
        auto ext_shape = sutils::getshape(extended_array);
        std::copy(ext_shape.begin(), ext_shape.begin() + N, shapes.begin());
        d_in = reinterpret_cast<T *>(extended_array.buffer);
        in_strides = create_strides(extended_array);
      } else {
        d_in = reinterpret_cast<T *>(in_array.buffer);
        in_shape[axis] = n;
        std::copy(in_shape.begin(), in_shape.begin() + N, shapes.begin());
        in_strides =
            create_strides(in_array); // for cropped arrays we need to use different strides
      }
      auto d_out = reinterpret_cast<std::complex<T> *>(out_array.buffer);
      // axes calculation is for 1D transform
      shape_t axes = shape_t(1);
      axes[0] = axis;
      auto out_strides = create_strides(out_array);
      auto fct = norm_fct<T>(inorm, shapes, axes);
      pocketfft::r2c(shapes, in_strides, out_strides, axes, forward, d_in, d_out, fct, size_t(0));
      if (extend) {
        using namespace pocketfft::detail;
        ndarr<std::complex<T>> ares(out_array.buffer, shapes, out_strides);
        rev_iter iter(ares, axes);
        while (iter.remaining() > 0) {
          auto v = ares[iter.ofs()];
          ares[iter.rev_ofs()] = conj(v);
          iter.advance();
        }
      }
      return out_array;
    }
  } // namespace fft
} // namespace numpy
PYTHONIC_NS_END

#endif
