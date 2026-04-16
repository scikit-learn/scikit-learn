#ifndef PYTHONIC_TYPES_NUMPY_OPERATORS_HPP
#define PYTHONIC_TYPES_NUMPY_OPERATORS_HPP

#include "pythonic/include/types/numpy_operators.hpp"

#include "pythonic/numpy/bitwise_not.hpp"
#include "pythonic/numpy/mod.hpp"
#include "pythonic/operator_/add.hpp"
#include "pythonic/operator_/and_.hpp"
#include "pythonic/operator_/div.hpp"
#include "pythonic/operator_/eq.hpp"
#include "pythonic/operator_/ge.hpp"
#include "pythonic/operator_/gt.hpp"
#include "pythonic/operator_/le.hpp"
#include "pythonic/operator_/lshift.hpp"
#include "pythonic/operator_/lt.hpp"
#include "pythonic/operator_/mul.hpp"
#include "pythonic/operator_/ne.hpp"
#include "pythonic/operator_/neg.hpp"
#include "pythonic/operator_/not_.hpp"
#include "pythonic/operator_/or_.hpp"
#include "pythonic/operator_/pos.hpp"
#include "pythonic/operator_/rshift.hpp"
#include "pythonic/operator_/sub.hpp"
#include "pythonic/operator_/xor_.hpp"
#include "pythonic/types/numpy_broadcast.hpp"
#include "pythonic/types/numpy_op_helper.hpp"

PYTHONIC_NS_BEGIN
/* operators must live in the same namespace as the associated type */
namespace types
{
#define NUMPY_BINARY_FUNC_NAME operator+
#define NUMPY_BINARY_FUNC_SYM operator_::functor::add
#include "pythonic/types/numpy_binary_op.hpp"

#define NUMPY_BINARY_FUNC_NAME operator&
#define NUMPY_BINARY_FUNC_SYM operator_::functor::and_
#include "pythonic/types/numpy_binary_op.hpp"

#define NUMPY_UNARY_FUNC_NAME operator~
#define NUMPY_UNARY_FUNC_SYM numpy::functor::bitwise_not
#include "pythonic/types/numpy_unary_op.hpp"

#define NUMPY_BINARY_FUNC_NAME operator|
#define NUMPY_BINARY_FUNC_SYM operator_::functor::or_
#include "pythonic/types/numpy_binary_op.hpp"

#define NUMPY_BINARY_FUNC_NAME operator^
#define NUMPY_BINARY_FUNC_SYM operator_::functor::xor_
#include "pythonic/types/numpy_binary_op.hpp"

#define NUMPY_BINARY_FUNC_NAME operator/
#define NUMPY_BINARY_FUNC_SYM operator_::functor::div
#include "pythonic/types/numpy_binary_op.hpp"

#define NUMPY_BINARY_FUNC_NAME operator==
#define NUMPY_BINARY_FUNC_SYM operator_::functor::eq
#include "pythonic/types/numpy_binary_op.hpp"

#define NUMPY_BINARY_FUNC_NAME operator%
#define NUMPY_BINARY_FUNC_SYM numpy::functor::mod
#include "pythonic/types/numpy_binary_op.hpp"

#define NUMPY_BINARY_FUNC_NAME operator>
#define NUMPY_BINARY_FUNC_SYM operator_::functor::gt
#include "pythonic/types/numpy_binary_op.hpp"

#define NUMPY_BINARY_FUNC_NAME operator>=
#define NUMPY_BINARY_FUNC_SYM operator_::functor::ge
#include "pythonic/types/numpy_binary_op.hpp"

#define NUMPY_BINARY_FUNC_NAME operator<<
#define NUMPY_BINARY_FUNC_SYM operator_::functor::lshift
#include "pythonic/types/numpy_binary_op.hpp"

#define NUMPY_BINARY_FUNC_NAME operator<
#define NUMPY_BINARY_FUNC_SYM operator_::functor::lt
#include "pythonic/types/numpy_binary_op.hpp"

#define NUMPY_BINARY_FUNC_NAME operator<=
#define NUMPY_BINARY_FUNC_SYM operator_::functor::le
#include "pythonic/types/numpy_binary_op.hpp"

#define NUMPY_BINARY_FUNC_NAME operator*
#define NUMPY_BINARY_FUNC_SYM operator_::functor::mul
#include "pythonic/types/numpy_binary_op.hpp"

#define NUMPY_UNARY_FUNC_NAME operator-
#define NUMPY_UNARY_FUNC_SYM operator_::functor::neg
#include "pythonic/types/numpy_unary_op.hpp"

#define NUMPY_BINARY_FUNC_NAME operator!=
#define NUMPY_BINARY_FUNC_SYM operator_::functor::ne
#include "pythonic/types/numpy_binary_op.hpp"

#define NUMPY_UNARY_FUNC_NAME operator+
#define NUMPY_UNARY_FUNC_SYM operator_::functor::pos
#include "pythonic/types/numpy_unary_op.hpp"

#define NUMPY_UNARY_FUNC_NAME operator!
#define NUMPY_UNARY_FUNC_SYM operator_::functor::not_
#include "pythonic/types/numpy_unary_op.hpp"

#define NUMPY_BINARY_FUNC_NAME operator>>
#define NUMPY_BINARY_FUNC_SYM operator_::functor::rshift
#include "pythonic/types/numpy_binary_op.hpp"

#define NUMPY_BINARY_FUNC_NAME operator-
#define NUMPY_BINARY_FUNC_SYM operator_::functor::sub
#include "pythonic/types/numpy_binary_op.hpp"
} // namespace types
PYTHONIC_NS_END

#endif
