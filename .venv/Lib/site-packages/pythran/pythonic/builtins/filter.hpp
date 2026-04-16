#ifndef PYTHONIC_BUILTIN_FILTER_HPP
#define PYTHONIC_BUILTIN_FILTER_HPP

#include "pythonic/include/builtins/filter.hpp"

#include "pythonic/itertools/common.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/iterator.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{
  namespace details
  {
    template <typename Operator, typename List0>
    bool filter_iterator<Operator, List0>::test_filter(std::false_type)
    {
      return op(*iter);
    }

    template <typename Operator, typename List0>
    bool filter_iterator<Operator, List0>::test_filter(std::true_type)
    {
      return *iter;
    }

    template <typename Operator, typename List0>
    filter_iterator<Operator, List0>::filter_iterator(Operator _op, List0 &_seq)
        : op(_op), iter(_seq.begin()), iter_end(_seq.end())
    {
      if (iter != iter_end && !test_filter(std::is_same<types::none_type, Operator>()))
        next_value();
    }

    template <typename Operator, typename List0>
    filter_iterator<Operator, List0>::filter_iterator(itertools::npos, Operator _op, List0 &_seq)
        : op(_op), iter(_seq.end()), iter_end(_seq.end())
    {
    }

    template <typename Operator, typename List0>
    typename List0::value_type filter_iterator<Operator, List0>::operator*() const
    {
      return *iter;
    }

    template <typename Operator, typename List0>
    filter_iterator<Operator, List0> &filter_iterator<Operator, List0>::operator++()
    {
      next_value();
      return *this;
    }

    template <typename Operator, typename List0>
    void filter_iterator<Operator, List0>::next_value()
    {
      while (++iter != iter_end) {
        if (test_filter(std::is_same<types::none_type, Operator>()))
          return;
      }
    }

    template <typename Operator, typename List0>
    bool filter_iterator<Operator, List0>::operator==(filter_iterator const &other) const
    {
      return !(iter != other.iter);
    }

    template <typename Operator, typename List0>
    bool filter_iterator<Operator, List0>::operator!=(filter_iterator const &other) const
    {
      return iter != other.iter;
    }

    template <typename Operator, typename List0>
    bool filter_iterator<Operator, List0>::operator<(filter_iterator const &other) const
    {
      return iter != other.iter;
    }

    template <typename Operator, typename List0>
    filter<Operator, List0>::filter(Operator _op, List0 const &_seq)
        : utils::iterator_reminder<false, List0>(_seq), iterator(_op, this->values),
          end_iter(itertools::npos(), _op, this->values)
    {
    }

    template <typename Operator, typename List0>
    typename filter<Operator, List0>::iterator &filter<Operator, List0>::begin()
    {
      return *this;
    }

    template <typename Operator, typename List0>
    typename filter<Operator, List0>::iterator const &filter<Operator, List0>::begin() const
    {
      return *this;
    }

    template <typename Operator, typename List0>
    typename filter<Operator, List0>::iterator const &filter<Operator, List0>::end() const
    {
      return end_iter;
    }
  } // namespace details

  template <typename Operator, typename List0>
  details::filter<std::remove_cv_t<std::remove_reference_t<Operator>>,
                  std::remove_cv_t<std::remove_reference_t<List0>>>
  filter(Operator &&_op, List0 &&_seq)
  {
    return {std::forward<Operator>(_op), std::forward<List0>(_seq)};
  }
} // namespace builtins
PYTHONIC_NS_END

#endif
