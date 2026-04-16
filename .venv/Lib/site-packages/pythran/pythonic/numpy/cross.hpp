#ifndef PYTHONIC_NUMPY_CROSS_HPP
#define PYTHONIC_NUMPY_CROSS_HPP

#include "pythonic/include/numpy/cross.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <size_t N, size_t En, size_t Fn>
  struct _cross {
    template <class Out, class E, class F>
    void operator()(Out obegin, Out oend, E ebegin, F fbegin)
    {
      while (obegin != oend) {
        _cross<N - 1, En, Fn>{}((*obegin).begin(), (*obegin).end(), (*ebegin).begin(),
                                (*fbegin).begin());
        ++obegin, ++ebegin, ++fbegin;
      }
    }
  };
  template <>
  struct _cross<1, 2, 2> {
    template <class Out, class E, class F>
    void operator()(Out obegin, Out oend, E ebegin, F fbegin)
    {
      auto e0 = *ebegin;
      ++ebegin;
      auto e1 = *ebegin;
      auto f0 = *fbegin;
      ++fbegin;
      auto f1 = *fbegin;
      *obegin = e0 * f1 - e1 * f0;
    }
  };
  template <>
  struct _cross<1, 2, 3> {
    template <class Out, class E, class F>
    void operator()(Out obegin, Out oend, E ebegin, F fbegin)
    {
      auto e0 = *ebegin;
      ++ebegin;
      auto e1 = *ebegin;
      decltype(e1) e2 = 0;
      auto f0 = *fbegin;
      ++fbegin;
      auto f1 = *fbegin;
      ++fbegin;
      auto f2 = *fbegin;
      *obegin = e1 * f2 - e2 * f1;
      ++obegin;
      *obegin = e2 * f0 - e0 * f2;
      ++obegin;
      *obegin = e0 * f1 - e1 * f0;
    }
  };
  template <>
  struct _cross<1, 3, 3> {
    template <class Out, class E, class F>
    void operator()(Out obegin, Out oend, E ebegin, F fbegin)
    {
      auto e0 = *ebegin;
      ++ebegin;
      auto e1 = *ebegin;
      ++ebegin;
      auto e2 = *ebegin;
      auto f0 = *fbegin;
      ++fbegin;
      auto f1 = *fbegin;
      ++fbegin;
      auto f2 = *fbegin;
      *obegin = e1 * f2 - e2 * f1;
      ++obegin;
      *obegin = e2 * f0 - e0 * f2;
      ++obegin;
      *obegin = e0 * f1 - e1 * f0;
    }
  };
  template <>
  struct _cross<1, 3, 2> {
    template <class Out, class E, class F>
    void operator()(Out obegin, Out oend, E ebegin, F fbegin)
    {
      auto e0 = *ebegin;
      ++ebegin;
      auto e1 = *ebegin;
      ++ebegin;
      auto e2 = *ebegin;
      auto f0 = *fbegin;
      ++fbegin;
      auto f1 = *fbegin;
      decltype(f1) f2 = 0;
      *obegin = e1 * f2 - e2 * f1;
      ++obegin;
      *obegin = e2 * f0 - e0 * f2;
      ++obegin;
      *obegin = e0 * f1 - e1 * f0;
    }
  };

  template <class E, class F>
  types::ndarray<typename __combined<typename E::dtype, typename F::dtype>::type,
                 types::array_tuple<long, E::value>>
  cross(E const &e, F const &f)
  {
    using dtype = typename __combined<typename E::dtype, typename F::dtype>::type;
    types::array_tuple<long, E::value> out_shape;
    sutils::copy_shape<0, 0>(out_shape, e, std::make_index_sequence<E::value - 1>());
    if (e.template shape<E::value - 1>() == 2) {
      if (f.template shape<F::value - 1>() == 2) {
        out_shape[E::value - 1] = 1;
        types::ndarray<dtype, types::array_tuple<long, E::value>> out{out_shape,
                                                                      types::none_type{}};
        _cross<E::value, 2, 2>{}(out.begin(), out.end(), e.begin(), f.begin());
        return out;
      } else {
        out_shape[E::value - 1] = 3;
        types::ndarray<dtype, types::array_tuple<long, E::value>> out{out_shape,
                                                                      types::none_type{}};
        _cross<E::value, 2, 3>{}(out.begin(), out.end(), e.begin(), f.begin());
        return out;
      }
    } else {
      if (f.template shape<F::value - 1>() == 2) {
        out_shape[E::value - 1] = 3;
        types::ndarray<dtype, types::array_tuple<long, E::value>> out{out_shape,
                                                                      types::none_type{}};
        _cross<E::value, 3, 2>{}(out.begin(), out.end(), e.begin(), f.begin());
        return out;
      } else {
        out_shape[E::value - 1] = 3;
        types::ndarray<dtype, types::array_tuple<long, E::value>> out{out_shape,
                                                                      types::none_type{}};
        _cross<E::value, 3, 3>{}(out.begin(), out.end(), e.begin(), f.begin());
        return out;
      }
    }
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
