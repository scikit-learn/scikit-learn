# We declare these parser functions here to pass a C function type around.
# FIXME: do this in Python notation.

from .Scanning cimport PyrexScanner

ctypedef object (*p_sub_expr_func)(PyrexScanner obj)

cdef p_binop_expr(PyrexScanner s, ops, p_sub_expr_func p_sub_expr)
cdef p_rassoc_binop_expr(PyrexScanner s, unicode op, p_sub_expr_func p_subexpr)
