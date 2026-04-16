#ifndef PYTHONIC_INCLUDE_TYPES_NUMPY_OPERATORS_HPP
#define PYTHONIC_INCLUDE_TYPES_NUMPY_OPERATORS_HPP

#include "pythonic/include/numpy/bitwise_not.hpp"
#include "pythonic/include/numpy/mod.hpp"
#include "pythonic/include/operator_/add.hpp"
#include "pythonic/include/operator_/and_.hpp"
#include "pythonic/include/operator_/div.hpp"
#include "pythonic/include/operator_/eq.hpp"
#include "pythonic/include/operator_/ge.hpp"
#include "pythonic/include/operator_/gt.hpp"
#include "pythonic/include/operator_/le.hpp"
#include "pythonic/include/operator_/lshift.hpp"
#include "pythonic/include/operator_/lt.hpp"
#include "pythonic/include/operator_/mul.hpp"
#include "pythonic/include/operator_/ne.hpp"
#include "pythonic/include/operator_/neg.hpp"
#include "pythonic/include/operator_/not_.hpp"
#include "pythonic/include/operator_/or_.hpp"
#include "pythonic/include/operator_/pos.hpp"
#include "pythonic/include/operator_/rshift.hpp"
#include "pythonic/include/operator_/sub.hpp"
#include "pythonic/include/operator_/xor_.hpp"
#include "pythonic/include/types/numpy_broadcast.hpp"
#include "pythonic/include/types/numpy_op_helper.hpp"

PYTHONIC_NS_BEGIN
/* operators must live in the same namespace as the associated type */
namespace types
{
#define NUMPY_BINARY_FUNC_NAME operator+
#define NUMPY_BINARY_FUNC_SYM operator_::functor::add
#include "pythonic/include/types/numpy_binary_op.hpp"

#define NUMPY_BINARY_FUNC_NAME operator&
#define NUMPY_BINARY_FUNC_SYM operator_::functor::and_
#include "pythonic/include/types/numpy_binary_op.hpp"

#define NUMPY_UNARY_FUNC_NAME operator~
#define NUMPY_UNARY_FUNC_SYM numpy::functor::bitwise_not
#include "pythonic/include/types/numpy_unary_op.hpp"

#define NUMPY_BINARY_FUNC_NAME operator|
#define NUMPY_BINARY_FUNC_SYM operator_::functor::or_
#include "pythonic/include/types/numpy_binary_op.hpp"

#define NUMPY_BINARY_FUNC_NAME operator^
#define NUMPY_BINARY_FUNC_SYM operator_::functor::xor_
#include "pythonic/include/types/numpy_binary_op.hpp"

#define NUMPY_BINARY_FUNC_NAME operator/
#define NUMPY_BINARY_FUNC_SYM operator_::functor::div
#include "pythonic/include/types/numpy_binary_op.hpp"

#define NUMPY_BINARY_FUNC_NAME operator==
#define NUMPY_BINARY_FUNC_SYM operator_::functor::eq
#include "pythonic/include/types/numpy_binary_op.hpp"

#define NUMPY_BINARY_FUNC_NAME operator%
#define NUMPY_BINARY_FUNC_SYM numpy::functor::mod
#include "pythonic/include/types/numpy_binary_op.hpp"

#define NUMPY_BINARY_FUNC_NAME operator>
#define NUMPY_BINARY_FUNC_SYM operator_::functor::gt
#include "pythonic/include/types/numpy_binary_op.hpp"

#define NUMPY_BINARY_FUNC_NAME operator>=
#define NUMPY_BINARY_FUNC_SYM operator_::functor::ge
#include "pythonic/include/types/numpy_binary_op.hpp"

#define NUMPY_BINARY_FUNC_NAME operator<<
#define NUMPY_BINARY_FUNC_SYM operator_::functor::lshift
#include "pythonic/include/types/numpy_binary_op.hpp"

#define NUMPY_BINARY_FUNC_NAME operator<
#define NUMPY_BINARY_FUNC_SYM operator_::functor::lt
#include "pythonic/include/types/numpy_binary_op.hpp"

#define NUMPY_BINARY_FUNC_NAME operator<=
#define NUMPY_BINARY_FUNC_SYM operator_::functor::le
#include "pythonic/include/types/numpy_binary_op.hpp"

#define NUMPY_BINARY_FUNC_NAME operator*
#define NUMPY_BINARY_FUNC_SYM operator_::functor::mul
#include "pythonic/include/types/numpy_binary_op.hpp"

#define NUMPY_UNARY_FUNC_NAME operator-
#define NUMPY_UNARY_FUNC_SYM operator_::functor::neg
#include "pythonic/include/types/numpy_unary_op.hpp"

#define NUMPY_BINARY_FUNC_NAME operator!=
#define NUMPY_BINARY_FUNC_SYM operator_::functor::ne
#include "pythonic/include/types/numpy_binary_op.hpp"

#define NUMPY_UNARY_FUNC_NAME operator+
#define NUMPY_UNARY_FUNC_SYM operator_::functor::pos
#include "pythonic/include/types/numpy_unary_op.hpp"

#define NUMPY_UNARY_FUNC_NAME operator!
#define NUMPY_UNARY_FUNC_SYM operator_::functor::not_
#include "pythonic/include/types/numpy_unary_op.hpp"

#define NUMPY_BINARY_FUNC_NAME operator>>
#define NUMPY_BINARY_FUNC_SYM operator_::functor::rshift
#include "pythonic/include/types/numpy_binary_op.hpp"

#define NUMPY_BINARY_FUNC_NAME operator-
#define NUMPY_BINARY_FUNC_SYM operator_::functor::sub
#include "pythonic/include/types/numpy_binary_op.hpp"
} // namespace types
PYTHONIC_NS_END

#endif
