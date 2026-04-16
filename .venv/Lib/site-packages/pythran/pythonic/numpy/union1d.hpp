#ifndef PYTHONIC_NUMPY_UNION1D_HPP
#define PYTHONIC_NUMPY_UNION1D_HPP

#include "pythonic/include/numpy/union1d.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

#include <set>

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace
  {
    template <class I, class O>
    void _union1d(I begin, I end, O &out, utils::int_<1>)
    {
      for (; begin != end; ++begin)
        out.insert(*begin);
    }

    template <class I, class O, size_t N>
    void _union1d(I begin, I end, O &out, utils::int_<N>)
    {
      for (; begin != end; ++begin)
        _union1d((*begin).begin(), (*begin).end(), out, utils::int_<N - 1>());
    }
  } // namespace

  template <class E, class F>
  types::ndarray<typename __combined<typename E::dtype, typename F::dtype>::type,
                 types::pshape<long>>
  union1d(E const &e, F const &f)
  {
    std::set<typename __combined<typename E::dtype, typename F::dtype>::type> res;
    _union1d(e.begin(), e.end(), res, utils::int_<E::value>());
    _union1d(f.begin(), f.end(), res, utils::int_<F::value>());
    return {res};
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
