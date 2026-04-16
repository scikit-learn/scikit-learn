#ifndef PYTHONIC_INCLUDE_UTILS_ITERATOR_HPP
#define PYTHONIC_INCLUDE_UTILS_ITERATOR_HPP

PYTHONIC_NS_BEGIN

namespace utils
{

  template <class T>
  struct comparable_iterator : T {
    comparable_iterator();
    comparable_iterator(T const &t);
    bool operator<(comparable_iterator<T> other);
  };

  // Utility class to remind sequence we are iterating on to avoid dangling
  // reference
  template <bool as_tuple, class T, class... Others>
  struct iterator_reminder;

  template <class T, class... Others>
  struct iterator_reminder<true, T, Others...> {
    std::tuple<T, Others...> values;
    // FIXME : It works only because template arguments are ! references
    // so it trigger a copy.
    iterator_reminder() = default;
    iterator_reminder(T const &v, Others const &...o);
  };

  template <class T>
  struct iterator_reminder<false, T> {
    T values;
    iterator_reminder() = default;
    iterator_reminder(T const &v);
  };

  template <class T>
  struct iterator_reminder<true, T> {
    std::tuple<T> values;
    iterator_reminder() = default;
    iterator_reminder(T const &v);
  };

  /* Get the "minimum" of all iterators :
     - only random => random
     - at least one forward => forward
     */
  template <typename... Iters>
  struct iterator_min;

  template <typename T>
  struct iterator_min<T> {
    using type = typename std::iterator_traits<T>::iterator_category;
  };

  template <typename T, typename... Iters>
  struct iterator_min<T, Iters...> {
    using type =
        std::conditional_t<std::is_same<typename std::iterator_traits<T>::iterator_category,
                                        std::forward_iterator_tag>::value,
                           std::forward_iterator_tag, typename iterator_min<Iters...>::type>;
  };
} // namespace utils
PYTHONIC_NS_END

#endif
