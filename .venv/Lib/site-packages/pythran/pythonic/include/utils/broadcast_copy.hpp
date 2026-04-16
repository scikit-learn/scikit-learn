#ifndef PYTHONIC_INCLUDE_UTILS_BROADCAST_COPY_HPP
#define PYTHONIC_INCLUDE_UTILS_BROADCAST_COPY_HPP

#include "pythonic/include/types/tuple.hpp"

PYTHONIC_NS_BEGIN

namespace utils
{

  /* helper function to get the dimension of an array
   * yields 0 for scalar types
   */
  template <class T, typename EnableDefault = void>
  struct dim_of {
    static const size_t value = T::value;
  };

  template <class T, size_t N, class V>
  struct dim_of<types::array_base<T, N, V>, void> {
    static const size_t value = 1 + dim_of<T>::value;
  };

  template <class T>
  struct dim_of<T, std::enable_if_t<std::is_fundamental<T>::value>> {
    static const size_t value = 0;
  };

#define SPECIALIZE_DIM_OF(TYPE)                                                                    \
  template <>                                                                                      \
  struct dim_of<TYPE> {                                                                            \
    static const size_t value = 0;                                                                 \
  }

  SPECIALIZE_DIM_OF(std::complex<float>);
  SPECIALIZE_DIM_OF(std::complex<double>);

#undef SPECIALIZE_DIM_OF

  template <class E, class F, size_t N, int D, bool vector_form>
  E &broadcast_copy(E &self, F const &other);

  template <class Op, class E, class F, size_t N, int D, bool vector_form>
  E &broadcast_update(E &self, F const &other);
} // namespace utils
PYTHONIC_NS_END

#endif
