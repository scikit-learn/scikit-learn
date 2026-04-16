#ifndef PYTHONIC_NUMPY_PARTIAL_SUM_HPP
#define PYTHONIC_NUMPY_PARTIAL_SUM_HPP

#include "pythonic/include/numpy/partial_sum.hpp"

#include "pythonic/builtins/ValueError.hpp"
#include "pythonic/types/ndarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  /**
   * The cast is perform to be numpy compliant
   *
   * a = numpy.array([1, 256])
   * In [10]: numpy.mod.accumulate(a, dtype=numpy.uint32)
   * Out[10]: array([1, 1], dtype=uint32)
   * In [11]: numpy.mod.accumulate(a, dtype=numpy.uint8)
   * Out[11]: array([1, 0], dtype=uint8)
   */
  namespace
  {

    template <class Op, size_t N, class A>
    struct _partial_sum {
      template <class E, class F>
      A operator()(E const &e, F &o)
      {
        auto it_begin = e.begin();
        A acc = _partial_sum<Op, N - 1, A>{}((*it_begin), o);
        ++it_begin;
        for (; it_begin < e.end(); ++it_begin)
          acc = _partial_sum<Op, N - 1, A>{}(*it_begin, o, acc);
        return acc;
      }
      template <class E, class F>
      A operator()(E const &e, F &o, A acc)
      {
        for (auto const &value : e)
          acc = _partial_sum<Op, N - 1, A>{}(value, o, acc);
        return acc;
      }
    };

    template <class Op, class A>
    struct _partial_sum<Op, 1, A> {
      template <class E, class F>
      A operator()(E const &e, F &o)
      {
        auto it_begin = e.begin();
        A acc = *it_begin;
        *o = acc;
        ++it_begin, ++o;
        for (; it_begin < e.end(); ++it_begin, ++o) {
          acc = Op{}(acc, (A)*it_begin);
          *o = acc;
        }
        return acc;
      }
      template <class E, class F>
      A operator()(E e, F &o, A acc)
      {
        for (auto const &value : e) {
          acc = Op{}(acc, (A)value);
          *o = acc;
          ++o;
        }
        return acc;
      }
    };
  } // namespace

  template <class Op, class E, class dtype>
  types::ndarray<typename dtype::type, types::pshape<long>> partial_sum(E const &expr, dtype d)
  {
    const long count = expr.flat_size();
    types::ndarray<typename dtype::type, types::pshape<long>> the_partial_sum{
        types::make_tuple(count), builtins::None};
    auto begin_it = the_partial_sum.begin();
    _partial_sum<Op, E::value, typename dtype::type>{}(expr, begin_it);
    return the_partial_sum;
  }

  template <class Op, class E, class dtype>
  auto partial_sum(E const &expr, long axis, dtype d)
      -> std::enable_if_t<E::value == 1, decltype(partial_sum<Op, E, dtype>(expr))>
  {
    if (axis != 0)
      throw types::ValueError("axis out of bounds");
    return partial_sum<Op, E, dtype>(expr);
  }

  template <class Op, class E, class dtype>
  std::enable_if_t<E::value != 1, partial_sum_type<Op, E, dtype>> partial_sum(E const &expr,
                                                                              long axis, dtype d)
  {
    if (axis < 0 || size_t(axis) >= E::value)
      throw types::ValueError("axis out of bounds");

    auto shape = sutils::getshape(expr);
    partial_sum_type<Op, E, dtype> the_partial_sum{shape, builtins::None};
    if (axis == 0) {
      auto it_begin = the_partial_sum.begin();
      _partial_sum<Op, 1, partial_sum_type2<Op, E, dtype>>{}(expr, it_begin);
    } else {
      std::transform(
          expr.begin(), expr.end(), the_partial_sum.begin(),
          [axis, d](typename std::iterator_traits<typename E::iterator>::value_type other) {
            return partial_sum<Op>(other, axis - 1, d);
          });
    }
    return the_partial_sum;
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
