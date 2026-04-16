#ifndef PYTHONIC_INCLUDE_OPERATOR_OVERLOADS_HPP
#define PYTHONIC_INCLUDE_OPERATOR_OVERLOADS_HPP

#define PYTHONIC_OPERATOR_OVERLOAD_DECL(type, opname, op) type opname(type a, type b);

// workaround the fact that char and short computations are done using int in C,
// while they are done at their respective type in numpy
#define DEFINE_ALL_OPERATOR_OVERLOADS_DECL(opname, op)                                             \
  PYTHONIC_OPERATOR_OVERLOAD_DECL(bool, opname, op)                                                \
  PYTHONIC_OPERATOR_OVERLOAD_DECL(unsigned char, opname, op)                                       \
  PYTHONIC_OPERATOR_OVERLOAD_DECL(char, opname, op)                                                \
  PYTHONIC_OPERATOR_OVERLOAD_DECL(signed char, opname, op)                                         \
  PYTHONIC_OPERATOR_OVERLOAD_DECL(unsigned short, opname, op)                                      \
  PYTHONIC_OPERATOR_OVERLOAD_DECL(signed short, opname, op)                                        \
  PYTHONIC_OPERATOR_OVERLOAD_DECL(unsigned int, opname, op)                                        \
  PYTHONIC_OPERATOR_OVERLOAD_DECL(signed int, opname, op)                                          \
  PYTHONIC_OPERATOR_OVERLOAD_DECL(unsigned long, opname, op)                                       \
  PYTHONIC_OPERATOR_OVERLOAD_DECL(signed long, opname, op)                                         \
  PYTHONIC_OPERATOR_OVERLOAD_DECL(unsigned long long, opname, op)                                  \
  PYTHONIC_OPERATOR_OVERLOAD_DECL(signed long long, opname, op)

#endif
