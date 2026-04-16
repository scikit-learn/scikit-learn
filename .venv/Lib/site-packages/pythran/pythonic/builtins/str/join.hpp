#ifndef PYTHONIC_BUILTIN_STR_JOIN_HPP
#define PYTHONIC_BUILTIN_STR_JOIN_HPP

#include "pythonic/include/builtins/str/join.hpp"

#include "pythonic/builtins/len.hpp"
#include "pythonic/types/str.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace str
  {

    template <class S>
    types::str join(S const &s, types::str const &iterable)
    {
      long ssize =
          std::distance(std::begin(s), std::end(s)) - (std::is_same<S, types::str>::value ? 0 : 1);
      /* first iterate over iterable to gather sizes */
      size_t n = ssize * (iterable.size() - 1) + iterable.size();

      std::string out(n, 0);

      auto iter = iterable.chars().begin();
      auto oter = out.begin();
      if (iter != iterable.chars().end()) {
        *oter++ = *iter++;
        if (ssize)
          for (; iter != iterable.chars().end(); ++iter) {
            for (auto &&v : s)
              *oter++ = v.chars()[0];
            *oter++ = *iter;
          }
        else
          std::copy(iter, iterable.chars().end(), oter);
      }
      return {std::move(out)};
    }

    template <class S, class Iterable>
    std::enable_if_t<
        !std::is_same<std::remove_cv_t<std::remove_reference_t<Iterable>>, types::str>::value &&
            std::is_same<typename std::iterator_traits<typename std::remove_reference_t<
                             Iterable>::iterator>::iterator_category,
                         std::random_access_iterator_tag>::value,
        types::str>
    join(S const &s, Iterable &&iterable)
    {
      long ssize = builtins::functor::len{}(s);

      /* first iterate over iterable to gather sizes */
      long iterable_size = std::distance(iterable.begin(), iterable.end());
      if (iterable_size == 0)
        return "";
      size_t n = ssize * (iterable_size - 1);
      for (auto const &iter : iterable)
        n += builtins::len(iter);

      std::string out(n, 0);

      auto iter = iterable.begin();
      auto oter = out.begin();
      if (iter != iterable.end()) {
        auto tmp = *iter;
        auto const &stmp = tmp.chars();
        oter = std::copy(stmp.begin(), stmp.end(), oter);
        ++iter;
        if (ssize)
          for (; iter != iterable.end(); ++iter) {
            auto chars = s.chars();
            oter = std::copy(std::begin(chars), std::begin(chars) + ssize, oter);
            auto tmp = *iter;
            auto const &stmp = tmp.chars();
            oter = std::copy(stmp.begin(), stmp.end(), oter);
          }
        else
          for (; iter != iterable.end(); ++iter) {
            auto tmp = (*iter);
            auto const &stmp = tmp.chars();
            oter = std::copy(stmp.begin(), stmp.end(), oter);
          }
      }
      return {std::move(out)};
    }

    template <class S, class Iterable>
    std::enable_if_t<!std::is_same<typename std::iterator_traits<typename std::remove_reference_t<
                                       Iterable>::iterator>::iterator_category,
                                   std::random_access_iterator_tag>::value,
                     types::str>
    join(S const &s, Iterable &&iterable)
    {
      types::str out;
      auto iter = iterable.begin();
      if (iter != iterable.end()) {
        out += *iter;
        ++iter;
        for (; iter != iterable.end(); ++iter) {
          out += s;
          out += *iter;
        }
      }
      return out;
    }
  } // namespace str
} // namespace builtins
PYTHONIC_NS_END
#endif
