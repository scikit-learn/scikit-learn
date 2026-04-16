#ifndef PYTHONIC_NUMPY_SELECT_HPP
#define PYTHONIC_NUMPY_SELECT_HPP

#include "pythonic/include/numpy/select.hpp"

#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/int_.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace
  {
    // TODO It could certainly be represent as a numpy_***_expr as each
    // elements
    // is computed without information from neighbor.
    //
    template <class Ichoice, class Icond, class Iout, class Isel>
    long _select(Ichoice ibegin, Ichoice iend, Iout obegin, Isel sbegin, Icond cbegin, long size,
                 utils::int_<1>)
    {
      static_assert(std::is_same<Ichoice, int>::value, "");
      for (; ibegin != iend && size != 0; ++ibegin, ++obegin, ++sbegin, ++cbegin) {
        // If elements it not already selected && condition match, copy it!
        if (!*sbegin && *cbegin) {
          *obegin = *ibegin;
          *sbegin = true;
          size--;
        }
      }
      return size;
    }

    template <class Ichoice, class Icond, class Iout, class Isel, size_t N>
    long _select(Ichoice ibegin, Ichoice iend, Iout obegin, Isel sbegin, Icond cbegin, long size,
                 utils::int_<N>)
    {
      for (; ibegin != iend && size != 0; ++ibegin, ++obegin, ++sbegin, ++cbegin)
        size = _select((*ibegin).begin(), (*ibegin).end(), (*obegin).begin(), (*sbegin).begin(),
                       (*cbegin).begin(), size, utils::int_<N - 1>());
      return size;
    }
  } // namespace

  template <class C, class L>
  types::ndarray<typename L::dtype, types::array_tuple<long, L::value - 1>>
  select(C const &condlist, L const &choicelist, typename L::dtype _default)
  {
    constexpr size_t N = L::value - 1;
    auto &&choicelist0_shape = sutils::getshape(choicelist[0]);
    types::ndarray<typename L::dtype, types::array_tuple<long, N>> out(choicelist0_shape, _default);
    types::ndarray<typename L::dtype, types::array_tuple<long, N>> selected(choicelist0_shape,
                                                                            false);
    long size = selected.flat_size();
    for (long i = 0; i < condlist.size() && size != 0; i++)
      size = _select(choicelist[i].begin(), choicelist[i].end(), out.begin(), selected.begin(),
                     condlist.begin(), size, utils::int_<N>());
    return out;
  }

  template <class C, class L, class T>
  types::ndarray<typename L::dtype, sutils::pop_head_t<typename L::shape_t>>
  select_helper(C const &condlist, L const &choicelist, T _default)
  {
    types::ndarray<typename L::dtype, sutils::pop_head_t<typename L::shape_t>> out(
        sutils::getshape(choicelist[0]), _default);
    for (long i = 0; i < out.flat_size(); ++i)
      for (long j = 0; j < (long)condlist.size(); ++j)
        if (condlist[j].buffer[i]) {
          out.buffer[i] = choicelist[j].buffer[i];
          break;
        }
    return out;
  }

  template <class T, class TpS, class U, class UpS>
  std::enable_if_t<std::tuple_size<TpS>::value == std::tuple_size<UpS>::value,
                   types::ndarray<T, TpS>>
  select(types::list<types::ndarray<U, UpS>> const &condlist,
         types::list<types::ndarray<T, TpS>> const &choicelist, T _default)
  {
    return select_helper(condlist, choicelist, _default);
  }
  template <class T, class TpS, class U, class UpS, size_t M>
  std::enable_if_t<std::tuple_size<TpS>::value == std::tuple_size<UpS>::value,
                   types::ndarray<T, TpS>>
  select(types::static_list<types::ndarray<U, UpS>, M> const &condlist,
         types::static_list<types::ndarray<T, TpS>, M> const &choicelist, T _default)
  {
    return select_helper(condlist, choicelist, _default);
  }
  template <class T, class TpS, class U, class UpS, size_t M>
  std::enable_if_t<std::tuple_size<TpS>::value == std::tuple_size<UpS>::value,
                   types::ndarray<T, TpS>>
  select(types::static_list<types::ndarray<U, UpS>, M> const &condlist,
         types::list<types::ndarray<T, TpS>> const &choicelist, T _default)
  {
    return select_helper(condlist, choicelist, _default);
  }
  template <class T, class TpS, class U, class UpS, size_t M>
  std::enable_if_t<std::tuple_size<TpS>::value == std::tuple_size<UpS>::value,
                   types::ndarray<T, TpS>>
  select(types::list<types::ndarray<U, UpS>> const &condlist,
         types::static_list<types::ndarray<T, TpS>, M> const &choicelist, T _default)
  {
    return select_helper(condlist, choicelist, _default);
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
