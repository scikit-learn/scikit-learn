#ifndef PYTHONIC_INCLUDE_TYPES_NUMPY_UFUNC_HPP
#define PYTHONIC_INCLUDE_TYPES_NUMPY_UFUNC_HPP

#include "numpy/ufuncobject.h"

PYTHONIC_NS_BEGIN

namespace types
{
  namespace detail
  {
    template <class F, typename ResType, typename... ArgTypes, size_t... Is>
    void ufunc_wrapper(char *output, char **inputs, npy_intp n, npy_intp output_step,
                       const npy_intp *inputs_steps, std::index_sequence<Is...>)
    {
      for (npy_intp i = 0; i < n; ++i) {
        *(ResType *)output =
            F{}(*(std::tuple_element_t<Is, std::tuple<ArgTypes...>> *)(inputs[Is])...);
        output += output_step;
        (void)std::initializer_list<int>{((inputs[Is] += inputs_steps[Is]), 0)...};
      }
    }
  } // namespace detail
  template <class F, typename ResType, typename... ArgTypes>
  void ufunc_wrapper(char **args, npy_intp const *dimensions, npy_intp const *steps,
                     void * /*extra*/)
  {
    npy_intp output_step = steps[sizeof...(ArgTypes)];
    return detail::ufunc_wrapper<F, ResType, ArgTypes...>(
        args[sizeof...(ArgTypes)], args, dimensions[0], output_step, steps,
        std::make_index_sequence<sizeof...(ArgTypes)>());
  }
} // namespace types
PYTHONIC_NS_END

#endif
