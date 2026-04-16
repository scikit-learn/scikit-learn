#ifndef PYTHONIC_INCLUDE_BUILTIN_PYTHRAN_STATIC_IF_HPP
#define PYTHONIC_INCLUDE_BUILTIN_PYTHRAN_STATIC_IF_HPP

#include "pythonic/include/builtins/pythran/is_none.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace pythran
  {
    namespace details
    {
      template <class IsSame>
      struct static_if;

      template <>
      struct static_if<types::true_type> {
        static_if(types::true_type)
        {
        }
        template <class F0, class F1>
        F0 operator()(F0 f0, F1 f1)
        {
          return f0;
        }
      };
      template <>
      struct static_if<types::false_type> {
        static_if(types::false_type)
        {
        }
        template <class F0, class F1>
        F1 operator()(F0 f0, F1 f1)
        {
          return f1;
        }
      };
      template <>
      struct static_if<bool> {
        bool state_;
        static_if(bool state) : state_(state)
        {
        }

        template <class F0, class F1>
        struct merged {
          bool state_;
          F0 f0;
          F1 f1;
          merged(bool state, F0 f0, F1 f1) : state_(state), f0(f0), f1(f1)
          {
          }
          template <class... Args>
          auto operator()(Args &&...args) const ->
              typename __combined<decltype(f0(std::forward<Args>(args)...)),
                                  decltype(f1(std::forward<Args>(args)...))>::type
          {
            if (state_)
              return f0(std::forward<Args>(args)...);
            else
              return f1(std::forward<Args>(args)...);
          }
        };

        template <class F0, class F1>
        merged<F0, F1> operator()(F0 f0, F1 f1)
        {
          return {state_, f0, f1};
        }
      };
    } // namespace details
    template <class T, class F0, class F1>
    auto static_if(T const &cond, F0 f0, F1 f1) -> decltype(details::static_if<T>{cond}(f0, f1));

    template <class F0, class F1>
    auto static_if(int const &cond, F0 f0, F1 f1)
        -> decltype(details::static_if<bool>{(bool)cond}(f0, f1))
    {
      return details::static_if<bool>{(bool)cond}(f0, f1);
    }

    DEFINE_FUNCTOR(pythonic::builtins::pythran, static_if);
  } // namespace pythran
} // namespace builtins
PYTHONIC_NS_END

#endif
