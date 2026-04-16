#ifndef PYTHONIC_NUMPY_NDARRAY_SORT_HPP
#define PYTHONIC_NUMPY_NDARRAY_SORT_HPP

#include "pythonic/include/numpy/ndarray/sort.hpp"

#include <algorithm>

#include "pythonic/numpy/array.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/types/str.hpp"
#include "pythonic/utils/allocate.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/pdqsort.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace ndarray
  {
    namespace
    {
      struct quicksorter {
        template <class... Args>
        void operator()(Args &&...args)
        {
          pdqsort(std::forward<Args>(args)...);
        }
      };
      struct mergesorter {
        template <class It, class Cmp = std::less<typename It::value_type>>
        void operator()(It first, It last, Cmp cmp = Cmp())
        {
          if (last - first > 1) {
            It middle = first + (last - first) / 2;
            operator()(first, middle, cmp);
            operator()(middle, last, cmp);
            std::inplace_merge(first, middle, last, cmp);
          }
        }
      };
      struct heapsorter {
        template <class... Args>
        void operator()(Args &&...args)
        {
          return std::sort_heap(std::forward<Args>(args)...);
        }
      };
      struct stablesorter {
        template <class... Args>
        void operator()(Args &&...args)
        {
          return std::stable_sort(std::forward<Args>(args)...);
        }
      };

      template <class T>
      struct _comp : std::less<T> {
      };

      template <>
      struct _comp<float> {
        bool operator()(float x, float y)
        {
          if (__builtin_expect(std::isnan(y), false))
            return true;
          return x < y;
        }
      };

      template <>
      struct _comp<double> {
        bool operator()(double x, double y)
        {
          if (__builtin_expect(std::isnan(y), false))
            return true;
          return x < y;
        }
      };

      template <class T>
      struct _comp<std::complex<T>> {
        bool operator()(std::complex<T> const &i, std::complex<T> const &j) const
        {
          if (std::real(i) == std::real(j))
            return std::imag(i) < std::imag(j);
          else
            return std::real(i) < std::real(j);
        }
      };

      template <class T>
      using comparator = _comp<T>;

      template <class T, class pS, class Sorter>
      std::enable_if_t<std::tuple_size<pS>::value == 1, void> _sort(types::ndarray<T, pS> &out,
                                                                    long axis, Sorter sorter)
      {
        sorter(out.begin(), out.end(), comparator<T>{});
      }

      template <class T, class pS, class Sorter>
      std::enable_if_t<std::tuple_size<pS>::value != 1, void> _sort(types::ndarray<T, pS> &out,
                                                                    long axis, Sorter sorter)
      {
        constexpr auto N = std::tuple_size<pS>::value;
        if (axis < 0)
          axis += N;
        long const flat_size = out.flat_size();
        if (axis == N - 1) {
          const long step = out.template shape<N - 1>();
          for (T *out_iter = out.buffer, *end_iter = out.buffer + flat_size; out_iter != end_iter;
               out_iter += step)
            sorter(out_iter, out_iter + step, comparator<T>{});
        } else {
          auto out_shape = sutils::getshape(out);
          const long step = std::accumulate(out_shape.begin() + axis, out_shape.end(), 1L,
                                            std::multiplies<long>());
          long const buffer_size = out_shape[axis];
          const long stepper = step / out_shape[axis];
          const long n = flat_size / out_shape[axis];
          long ith = 0, nth = 0;
          T *buffer = utils::allocate<T>(buffer_size);
          for (long i = 0; i < n; i++) {
            for (long j = 0; j < buffer_size; ++j)
              buffer[j] = out.buffer[ith + j * stepper];
            sorter(buffer, buffer + buffer_size, comparator<T>{});
            for (long j = 0; j < buffer_size; ++j)
              out.buffer[ith + j * stepper] = buffer[j];

            ith += step;
            if (ith >= flat_size) {
              ith = ++nth;
            }
          }
          utils::deallocate(buffer);
        }
      }
    } // namespace

    template <class E>
    types::none_type sort(E &&expr, long axis, types::none_type)
    {
      _sort(expr, axis, quicksorter());
      return {};
    }

    template <class E>
    types::none_type sort(E &&expr, types::none_type)
    {
      _sort(expr, 0, quicksorter());
      return {};
    }

    template <class E>
    types::none_type sort(E &&expr, long axis, types::str const &kind)
    {
      if (kind == "quicksort")
        _sort(expr, axis, quicksorter());
      else if (kind == "mergesort")
        _sort(expr, axis, mergesorter());
      else if (kind == "heapsort")
        _sort(expr, axis, heapsorter());
      else if (kind == "stable")
        _sort(expr, axis, stablesorter());
      return {};
    }
  } // namespace ndarray
} // namespace numpy
PYTHONIC_NS_END

#endif
