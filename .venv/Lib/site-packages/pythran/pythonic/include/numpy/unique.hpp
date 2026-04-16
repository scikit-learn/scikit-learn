#ifndef PYTHONIC_INCLUDE_NUMPY_UNIQUE_HPP
#define PYTHONIC_INCLUDE_NUMPY_UNIQUE_HPP

#include "pythonic/include/types/immediate.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/types/tuple.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  types::ndarray<typename E::dtype, types::pshape<long>> unique(E const &expr);

  template <class E>
  std::tuple<types::ndarray<typename E::dtype, types::pshape<long>>,
             types::ndarray<long, types::pshape<long>>>
  unique(E const &expr, types::true_immediate return_index);

  template <class E>
  types::ndarray<typename E::dtype, types::pshape<long>>
  unique(E const &expr, types::false_immediate return_index);

  template <class E>
  std::tuple<types::ndarray<typename E::dtype, types::pshape<long>>,
             types::ndarray<long, types::pshape<long>>>
  unique(E const &expr, types::false_immediate return_index, types::true_immediate return_inverse);

  template <class E>
  types::ndarray<typename E::dtype, types::pshape<long>>
  unique(E const &expr, types::false_immediate return_index, types::false_immediate return_inverse);

  template <class E>
  std::tuple<types::ndarray<typename E::dtype, types::pshape<long>>,
             types::ndarray<long, types::pshape<long>>>
  unique(E const &expr, types::true_immediate return_index, types::false_immediate return_inverse);

  template <class E>
  std::tuple<types::ndarray<typename E::dtype, types::pshape<long>>,
             types::ndarray<long, types::pshape<long>>, types::ndarray<long, types::pshape<long>>>
  unique(E const &expr, types::true_immediate return_index, types::true_immediate return_inverse);

  template <class E>
  std::tuple<types::ndarray<typename E::dtype, types::pshape<long>>,
             types::ndarray<long, types::pshape<long>>, types::ndarray<long, types::pshape<long>>,
             types::ndarray<long, types::pshape<long>>>
  unique(E const &expr, types::true_immediate return_index, types::true_immediate return_inverse,
         types::true_immediate return_counts);

  template <class E>
  std::tuple<types::ndarray<typename E::dtype, types::pshape<long>>,
             types::ndarray<long, types::pshape<long>>, types::ndarray<long, types::pshape<long>>>
  unique(E const &expr, types::true_immediate return_index, types::true_immediate return_inverse,
         types::false_immediate return_counts);

  template <class E>
  std::tuple<types::ndarray<typename E::dtype, types::pshape<long>>,
             types::ndarray<long, types::pshape<long>>>
  unique(E const &expr, types::true_immediate return_index, types::false_immediate return_inverse,
         types::false_immediate return_counts);

  template <class E>
  std::tuple<types::ndarray<typename E::dtype, types::pshape<long>>,
             types::ndarray<long, types::pshape<long>>, types::ndarray<long, types::pshape<long>>>
  unique(E const &expr, types::true_immediate return_index, types::false_immediate return_inverse,
         types::true_immediate return_counts);

  template <class E>
  std::tuple<types::ndarray<typename E::dtype, types::pshape<long>>,
             types::ndarray<long, types::pshape<long>>>
  unique(E const &expr, types::false_immediate return_index, types::true_immediate return_inverse,
         types::false_immediate return_counts);

  template <class E>
  std::tuple<types::ndarray<typename E::dtype, types::pshape<long>>,
             types::ndarray<long, types::pshape<long>>, types::ndarray<long, types::pshape<long>>>
  unique(E const &expr, types::false_immediate return_index, types::true_immediate return_inverse,
         types::true_immediate return_counts);

  template <class E>
  types::ndarray<typename E::dtype, types::pshape<long>>
  unique(E const &expr, types::false_immediate return_index, types::false_immediate return_inverse,
         types::false_immediate return_counts);

  template <class E>
  std::tuple<types::ndarray<typename E::dtype, types::pshape<long>>,
             types::ndarray<long, types::pshape<long>>>
  unique(E const &expr, types::false_immediate return_index, types::false_immediate return_inverse,
         types::true_immediate return_counts);

  DEFINE_FUNCTOR(pythonic::numpy, unique)
} // namespace numpy
PYTHONIC_NS_END

#endif
