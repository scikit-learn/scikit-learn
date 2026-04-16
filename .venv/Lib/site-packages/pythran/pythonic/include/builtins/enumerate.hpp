#ifndef PYTHONIC_INCLUDE_BUILTIN_ENUMERATE_HPP
#define PYTHONIC_INCLUDE_BUILTIN_ENUMERATE_HPP

#include "pythonic/include/types/tuple.hpp"
#include "pythonic/include/utils/functor.hpp"

#include <iterator>

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace details
  {
    // FIXME return value may be a type::make_tuple
    template <class Iterator>
    using enumerate_iterator_base = std::iterator<
        typename std::iterator_traits<Iterator>::iterator_category,
        types::make_tuple_t<long, typename std::iterator_traits<Iterator>::value_type>>;

    template <class Iterator>
    struct enumerate_iterator : public enumerate_iterator_base<Iterator> {
      long value;
      Iterator iter;
      enumerate_iterator();
      enumerate_iterator(Iterator const &iter, long first);
      typename enumerate_iterator_base<Iterator>::value_type operator*() const
      {
        return types::make_tuple(value, *iter);
      }
      enumerate_iterator &operator++()
      {
        ++value, ++iter;
        return *this;
      }
      enumerate_iterator &operator+=(long n);
      bool operator!=(enumerate_iterator const &other) const;
      bool operator<(enumerate_iterator const &other) const;
      long operator-(enumerate_iterator const &other) const;
      bool operator==(enumerate_iterator const &it) const;
    };

    template <class Iterable>
    struct enumerate
        : private Iterable, /* to hold a reference on the iterable */
          public enumerate_iterator<typename Iterable::iterator> /* to be compatible with
                                                                    builtins.next*/
    {
      using iterator = enumerate_iterator<typename Iterable::iterator>;
      using iterator::operator*;
      using value_type = typename iterator::value_type;

      iterator end_iter;

      enumerate();
      enumerate(Iterable seq, long first);
      iterator &begin();
      iterator const &begin() const;
      iterator end() const;
    };
  } // namespace details

  template <class Iterable>
  details::enumerate<std::remove_cv_t<std::remove_reference_t<Iterable>>>
  enumerate(Iterable &&seq, long first = 0L);

  DEFINE_FUNCTOR(pythonic::builtins, enumerate);
} // namespace builtins
PYTHONIC_NS_END

#endif
