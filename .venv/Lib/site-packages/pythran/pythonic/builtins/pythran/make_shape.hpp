#ifndef PYTHONIC_BUILTIN_MAKE_SHAPE_HPP
#define PYTHONIC_BUILTIN_MAKE_SHAPE_HPP

PYTHONIC_NS_BEGIN

namespace builtins
{
  namespace pythran
  {
    template <class... Args>
    pythonic::types::pshape<Args...> make_shape(Args... args)
    {
      return {args...};
    }
  } // namespace pythran
} // namespace builtins
PYTHONIC_NS_END

#endif
