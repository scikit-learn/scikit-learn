#ifndef PYTHONIC_NUMPY_STACK_HPP
#define PYTHONIC_NUMPY_STACK_HPP

#include "pythonic/builtins/ValueError.hpp"
#include "pythonic/builtins/len.hpp"
#include <pythonic/include/numpy/stack.hpp>
#include <pythonic/numpy/concatenate.hpp>

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class ArraySequence>
  types::ndarray<typename ArraySequence::value_type::dtype,
                 types::array_tuple<long, ArraySequence::value_type::value + 1>>
  stack(ArraySequence const &args, long axis)
  {
    if (builtins::len(args) == 0)
      throw pythonic::types::ValueError("need at least one array to stack");
    auto shape = sutils::getshape(args[0]);
    constexpr long N = std::tuple_size<decltype(shape)>::value; // The length of the shape
                                                                // array.
    auto values = sutils::array(shape); // You can't do shape[i] but you can do shape.array()[i]
    types::array_tuple<long, N + 1> new_shape; // A new array that's 1 element longer than shape.
    // Insert a "0" at the position indicated by axis.
    for (long i = 0; i < N + 1; i++) {
      if (i < axis)
        new_shape[i] = values[i];
      if (i == axis)
        new_shape[i] = 1;
      if (i > axis)
        new_shape[i] = values[i - 1];
    }

    // Create a new empty list.
    types::list<types::ndarray<typename ArraySequence::value_type::dtype,
                               types::array_tuple<long, ArraySequence::value_type::value + 1>>>
        bi(0);
    // Push the resized arrays into the list.
    for (auto &&arg : args) {
      bi.push_back(arg.reshape(new_shape));
    }
    // Call concatenate on this list.
    return concatenate(bi, axis);
  }
  template <size_t... Is, class... Tys>
  types::ndarray<typename details::stack_helper_t<Tys...>::dtype,
                 types::array_tuple<long, details::stack_helper_t<Tys...>::value + 1>>
  stack(std::tuple<Tys...> const &args, long axis, std::index_sequence<Is...>)
  {
    types::array_tuple<
        types::ndarray<typename details::stack_helper_t<Tys...>::dtype,
                       types::array_tuple<long, details::stack_helper_t<Tys...>::value>>,
        sizeof...(Tys)>
        vargs{{std::get<Is>(args)...}};
    return stack(vargs, axis);
  }

  template <class... Tys>
  types::ndarray<typename details::stack_helper_t<Tys...>::dtype,
                 types::array_tuple<long, details::stack_helper_t<Tys...>::value + 1>>
  stack(std::tuple<Tys...> const &args, long axis)
  {
    return stack(args, axis, std::make_index_sequence<sizeof...(Tys)>());
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
