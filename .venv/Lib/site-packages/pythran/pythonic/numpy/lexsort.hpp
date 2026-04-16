#ifndef PYTHONIC_NUMPY_LEXSORT_HPP
#define PYTHONIC_NUMPY_LEXSORT_HPP

#include "pythonic/include/numpy/lexsort.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/pdqsort.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    template <size_t I>
    struct lexcmp_nth {
      template <class K>
      bool operator()(K const &keys, long i0, long i1) const
      {
        if (std::get<I - 1>(keys)[i0] < std::get<I - 1>(keys)[i1])
          return true;
        else if (std::get<I - 1>(keys)[i0] > std::get<I - 1>(keys)[i1])
          return false;
        else
          return lexcmp_nth<I - 1>{}(keys, i0, i1);
      }
    };
    template <>
    struct lexcmp_nth<0> {
      template <class K>
      bool operator()(K const &keys, long i0, long i1) const
      {
        return false;
      }
    };

    template <class K>
    struct lexcmp {
      K const &keys;
      lexcmp(K const &keys) : keys(keys)
      {
      }
      bool operator()(long i0, long i1)
      {
        return lexcmp_nth<std::tuple_size<K>::value>{}(keys, i0, i1);
      }
    };
  } // namespace details

  template <class pS>
  types::ndarray<long, types::pshape<long>> lexsort(pS const &keys)
  {
    long n = std::get<0>(keys).size();
    types::ndarray<long, types::pshape<long>> out(types::pshape<long>(n), builtins::None);
    // fill with the original indices
    std::iota(out.buffer, out.buffer + n, 0L);
    // then sort using keys as the comparator
    pdqsort(out.buffer, out.buffer + n, details::lexcmp<pS>(keys));
    return out;
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
