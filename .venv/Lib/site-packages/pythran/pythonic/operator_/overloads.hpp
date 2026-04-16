#ifndef PYTHONIC_OPERATOR_OVERLOADS_HPP
#define PYTHONIC_OPERATOR_OVERLOADS_HPP

#include "pythonic/include/operator_/overloads.hpp"
#include <limits>

#define PYTHONIC_OPERATOR_OVERLOAD_IMPL(type, opname, op, overflow_check)                          \
  inline type opname(type a, type b)                                                               \
  {                                                                                                \
    assert((overflow_check) && "overflow check");                                                  \
    return a op b;                                                                                 \
  }

// workaround the fact that char and short computations are done using int in C,
// while they are done at their respective type in numpy
#define DEFINE_ALL_OPERATOR_OVERLOADS_NO_BOOL_IMPL(opname, op, overflow_check)                     \
  PYTHONIC_OPERATOR_OVERLOAD_IMPL(unsigned char, opname, op, true)                                 \
  PYTHONIC_OPERATOR_OVERLOAD_IMPL(char, opname, op, true)                                          \
  PYTHONIC_OPERATOR_OVERLOAD_IMPL(signed char, opname, op, overflow_check)                         \
  PYTHONIC_OPERATOR_OVERLOAD_IMPL(unsigned short, opname, op, true)                                \
  PYTHONIC_OPERATOR_OVERLOAD_IMPL(signed short, opname, op, overflow_check)                        \
  PYTHONIC_OPERATOR_OVERLOAD_IMPL(unsigned int, opname, op, true)                                  \
  PYTHONIC_OPERATOR_OVERLOAD_IMPL(signed int, opname, op, overflow_check)                          \
  PYTHONIC_OPERATOR_OVERLOAD_IMPL(unsigned long, opname, op, true)                                 \
  PYTHONIC_OPERATOR_OVERLOAD_IMPL(signed long, opname, op, overflow_check)                         \
  PYTHONIC_OPERATOR_OVERLOAD_IMPL(unsigned long long, opname, op, true)                            \
  PYTHONIC_OPERATOR_OVERLOAD_IMPL(signed long long, opname, op, overflow_check)

#define DEFINE_ALL_OPERATOR_OVERLOADS_IMPL(opname, op, overflow_check)                             \
  PYTHONIC_OPERATOR_OVERLOAD_IMPL(bool, opname, op, true)                                          \
  DEFINE_ALL_OPERATOR_OVERLOADS_NO_BOOL_IMPL(opname, op, overflow_check)

#endif
