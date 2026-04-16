#ifndef PYTHONIC_NUMPY_SETDIFF1D_HPP
#define PYTHONIC_NUMPY_SETDIFF1D_HPP

#include "pythonic/include/numpy/setdiff1d.hpp"

#include "pythonic/numpy/asarray.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/allocate.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/pdqsort.hpp"

#include <algorithm>
#include <set>

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace impl
  {

    template <typename InputIterator1, typename InputIterator2, typename OutputIterator>
    OutputIterator set_difference_unique(InputIterator1 first1, InputIterator1 last1,
                                         InputIterator2 first2, InputIterator2 last2,
                                         OutputIterator result)
    {
      while (first1 != last1 && first2 != last2) {
        auto const t1 = *first1;
        auto const t2 = *first2;
        if (t1 < t2) {
          *result = t1;
          while (*++first1 == t1)
            ;
          ++result;
        } else if (t2 < t1)
          while (*++first2 == t2)
            ;
        else {
          while (*++first1 == t1)
            ;
          while (*++first2 == t2)
            ;
        }
      }
      while (first1 != last1) {
        auto const t1 = *first1;
        *result = t1;
        while (*++first1 == t1)
          ;
        ++result;
      }
      return result;
    }
  } // namespace impl

  template <class T, class U>
  types::ndarray<typename __combined<typename types::dtype_of<T>::type,
                                     typename types::dtype_of<U>::type>::type,
                 types::pshape<long>>
  setdiff1d(T const &ar1, U const &ar2, bool assume_unique)
  {
    using dtype = typename __combined<typename types::dtype_of<T>::type,
                                      typename types::dtype_of<U>::type>::type;
    auto far1 = numpy::functor::array{}(ar1);
    auto far2 = numpy::functor::array{}(ar2);

    pdqsort(far1.fbegin(), far1.fend());
    pdqsort(far2.fbegin(), far2.fend());
    dtype *out = utils::allocate<dtype>(far1.flat_size() * far2.flat_size());

    dtype *out_last;
    if (assume_unique) {
      out_last = std::set_difference(far1.fbegin(), far1.fend(), far2.fbegin(), far2.fend(), out);
    } else {
      out_last =
          impl::set_difference_unique(far1.fbegin(), far1.fend(), far2.fbegin(), far2.fend(), out);
    }
    auto size = out_last - out;
    out = utils::reallocate(out, size);
    return {out, types::pshape<long>(size), types::ownership::owned};
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
