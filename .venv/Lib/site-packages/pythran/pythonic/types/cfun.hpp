#ifndef PYTHONIC_TYPES_CFUN_HPP
#define PYTHONIC_TYPES_CFUN_HPP

#include "pythonic/include/types/cfun.hpp"

PYTHONIC_NS_BEGIN

namespace types
{
  template <class ReturnType, class... ArgsType>
  cfun<ReturnType(ArgsType...)>::cfun(ReturnType (*fun)(ArgsType...)) : ptr(fun)
  {
  }

  template <class ReturnType, class... ArgsType>
  ReturnType cfun<ReturnType(ArgsType...)>::operator()(ArgsType... args) const
  {
    return (*ptr)(args...);
  }
} // namespace types
PYTHONIC_NS_END

#ifdef ENABLE_PYTHON_MODULE

PYTHONIC_NS_BEGIN

template <class R, class... Args>
PyObject *to_python<types::cfun<R(Args...)>>::convert(types::cfun<R(Args...)> const &v)
{
  return PyCapsule_New(v.ptr, nullptr, nullptr);
}

template <class R, class... Args>
bool from_python<types::cfun<R(Args...)>>::is_convertible(PyObject *obj)
{
  return PyCapsule_CheckExact(obj);
}

template <class R, class... Args>
types::cfun<R(Args...)> from_python<types::cfun<R(Args...)>>::convert(PyObject *obj)
{
  void *ptr = PyCapsule_GetPointer(obj, PyCapsule_GetName(obj) /* avoid the string check*/);
  return {reinterpret_cast<R (*)(Args...)>(ptr)};
}
PYTHONIC_NS_END

#endif

#endif
