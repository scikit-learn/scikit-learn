#ifndef PYTHONIC_INCLUDE_BUILTIN_FILTER_HPP
#define PYTHONIC_INCLUDE_BUILTIN_FILTER_HPP

#include "pythonic/include/itertools/common.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/iterator.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace details
  {

    template <typename Operator, typename List0>
    struct filter_iterator : std::iterator<std::forward_iterator_tag, typename List0::value_type> {
      using sequence_type = std::remove_cv_t<std::remove_reference_t<List0>>;

      Operator op;
      typename List0::iterator iter;
      // FIXME : iter_end should be const because filter should be evaluate
      // only once. Some tests doesn't work with it for now because of
      // uncorrect itertools.product implementation
      typename List0::iterator iter_end;

      bool test_filter(std::true_type);
      bool test_filter(std::false_type);

      filter_iterator() = default;
      filter_iterator(Operator _op, List0 &_seq);
      filter_iterator(itertools::npos, Operator _op, List0 &_seq);

      typename List0::value_type operator*() const;

      filter_iterator &operator++();
      void next_value();

      bool operator==(filter_iterator const &other) const;
      bool operator!=(filter_iterator const &other) const;
      bool operator<(filter_iterator const &other) const;
    };

    // Inherit from iterator_reminder to keep a reference on the iterator
    // and avoid a dangling reference
    // FIXME: It would be better to have a copy only if needed but Pythran
    // typing is not good enough for this as arguments have
    // remove_cv_t/remove_ref_t
    template <typename Operator, typename List0>
    struct filter : utils::iterator_reminder<false, List0>, filter_iterator<Operator, List0> {

      using value_type = typename List0::value_type;
      using iterator = filter_iterator<Operator, List0>;

      iterator end_iter;

      filter() = default;
      filter(Operator _op, List0 const &_seq);

      iterator &begin();
      iterator const &begin() const;
      iterator const &end() const;
    };
  } // namespace details

  template <typename Operator, typename List0>
  details::filter<std::remove_cv_t<std::remove_reference_t<Operator>>,
                  std::remove_cv_t<std::remove_reference_t<List0>>>
  filter(Operator &&_op, List0 &&_seq);

  DEFINE_FUNCTOR(pythonic::builtins, filter);
} // namespace builtins
PYTHONIC_NS_END

/* type inference stuff  {*/
#include "pythonic/include/types/combined.hpp"

template <class E, class Op, class T>
struct __combined<E, pythonic::builtins::details::filter<Op, T>> {
  using type = typename __combined<
      E, container<typename pythonic::builtins::details::filter<Op, T>::value_type>>::type;
};
/* } */

#endif
