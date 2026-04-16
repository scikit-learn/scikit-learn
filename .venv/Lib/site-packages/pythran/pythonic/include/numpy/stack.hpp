#ifndef PYTHONIC_INCLUDE_NUMPY_STACK_HPP
#define PYTHONIC_INCLUDE_NUMPY_STACK_HPP

#include <pythonic/include/types/ndarray.hpp>

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class ArraySequence>
  types::ndarray<typename ArraySequence::value_type::dtype,
                 types::array_tuple<long, ArraySequence::value_type::value + 1>>
  stack(ArraySequence const &args, long axis = 0);

  namespace details
  {
    template <class... Tys>
    using stack_helper_t = typename __combined<typename assignable<Tys>::type...>::type;
  }

  template <class... Tys>
  types::ndarray<typename details::stack_helper_t<Tys...>::dtype,
                 types::array_tuple<long, details::stack_helper_t<Tys...>::value + 1>>
  stack(std::tuple<Tys...> const &args, long axis = 0);

  DEFINE_FUNCTOR(pythonic::numpy, stack);
} // namespace numpy
PYTHONIC_NS_END

#endif
