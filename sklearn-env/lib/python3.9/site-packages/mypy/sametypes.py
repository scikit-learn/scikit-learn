from typing import Sequence

from mypy.types import (
    Type, UnboundType, AnyType, NoneType, TupleType, TypedDictType,
    UnionType, CallableType, TypeVarType, Instance, TypeVisitor, ErasedType,
    Overloaded, PartialType, DeletedType, UninhabitedType, TypeType, LiteralType,
    ProperType, get_proper_type, TypeAliasType, ParamSpecType
)
from mypy.typeops import tuple_fallback, make_simplified_union


def is_same_type(left: Type, right: Type) -> bool:
    """Is 'left' the same type as 'right'?"""

    left = get_proper_type(left)
    right = get_proper_type(right)

    if isinstance(right, UnboundType):
        # Make unbound types same as anything else to reduce the number of
        # generated spurious error messages.
        return True
    else:
        # Simplify types to canonical forms.
        #
        # There are multiple possible union types that represent the same type,
        # such as Union[int, bool, str] and Union[int, str]. Also, some union
        # types can be simplified to non-union types such as Union[int, bool]
        # -> int. It would be nice if we always had simplified union types but
        # this is currently not the case, though it often is.
        left = simplify_union(left)
        right = simplify_union(right)

        return left.accept(SameTypeVisitor(right))


def simplify_union(t: Type) -> ProperType:
    t = get_proper_type(t)
    if isinstance(t, UnionType):
        return make_simplified_union(t.items)
    return t


def is_same_types(a1: Sequence[Type], a2: Sequence[Type]) -> bool:
    if len(a1) != len(a2):
        return False
    for i in range(len(a1)):
        if not is_same_type(a1[i], a2[i]):
            return False
    return True


class SameTypeVisitor(TypeVisitor[bool]):
    """Visitor for checking whether two types are the 'same' type."""

    def __init__(self, right: ProperType) -> None:
        self.right = right

    # visit_x(left) means: is left (which is an instance of X) the same type as
    # right?

    def visit_unbound_type(self, left: UnboundType) -> bool:
        return True

    def visit_any(self, left: AnyType) -> bool:
        return isinstance(self.right, AnyType)

    def visit_none_type(self, left: NoneType) -> bool:
        return isinstance(self.right, NoneType)

    def visit_uninhabited_type(self, t: UninhabitedType) -> bool:
        return isinstance(self.right, UninhabitedType)

    def visit_erased_type(self, left: ErasedType) -> bool:
        # We can get here when isinstance is used inside a lambda
        # whose type is being inferred. In any event, we have no reason
        # to think that an ErasedType will end up being the same as
        # any other type, except another ErasedType (for protocols).
        return isinstance(self.right, ErasedType)

    def visit_deleted_type(self, left: DeletedType) -> bool:
        return isinstance(self.right, DeletedType)

    def visit_instance(self, left: Instance) -> bool:
        return (isinstance(self.right, Instance) and
                left.type == self.right.type and
                is_same_types(left.args, self.right.args) and
                left.last_known_value == self.right.last_known_value)

    def visit_type_alias_type(self, left: TypeAliasType) -> bool:
        # Similar to protocols, two aliases with the same targets return False here,
        # but both is_subtype(t, s) and is_subtype(s, t) return True.
        return (isinstance(self.right, TypeAliasType) and
                left.alias == self.right.alias and
                is_same_types(left.args, self.right.args))

    def visit_type_var(self, left: TypeVarType) -> bool:
        return (isinstance(self.right, TypeVarType) and
                left.id == self.right.id)

    def visit_param_spec(self, left: ParamSpecType) -> bool:
        # Ignore upper bound since it's derived from flavor.
        return (isinstance(self.right, ParamSpecType) and
                left.id == self.right.id and left.flavor == self.right.flavor)

    def visit_callable_type(self, left: CallableType) -> bool:
        # FIX generics
        if isinstance(self.right, CallableType):
            cright = self.right
            return (is_same_type(left.ret_type, cright.ret_type) and
                    is_same_types(left.arg_types, cright.arg_types) and
                    left.arg_names == cright.arg_names and
                    left.arg_kinds == cright.arg_kinds and
                    left.is_type_obj() == cright.is_type_obj() and
                    left.is_ellipsis_args == cright.is_ellipsis_args)
        else:
            return False

    def visit_tuple_type(self, left: TupleType) -> bool:
        if isinstance(self.right, TupleType):
            return (is_same_type(tuple_fallback(left), tuple_fallback(self.right))
                    and is_same_types(left.items, self.right.items))
        else:
            return False

    def visit_typeddict_type(self, left: TypedDictType) -> bool:
        if isinstance(self.right, TypedDictType):
            if left.items.keys() != self.right.items.keys():
                return False
            for (_, left_item_type, right_item_type) in left.zip(self.right):
                if not is_same_type(left_item_type, right_item_type):
                    return False
            return True
        else:
            return False

    def visit_literal_type(self, left: LiteralType) -> bool:
        if isinstance(self.right, LiteralType):
            if left.value != self.right.value:
                return False
            return is_same_type(left.fallback, self.right.fallback)
        else:
            return False

    def visit_union_type(self, left: UnionType) -> bool:
        if isinstance(self.right, UnionType):
            # Check that everything in left is in right
            for left_item in left.items:
                if not any(is_same_type(left_item, right_item) for right_item in self.right.items):
                    return False

            # Check that everything in right is in left
            for right_item in self.right.items:
                if not any(is_same_type(right_item, left_item) for left_item in left.items):
                    return False

            return True
        else:
            return False

    def visit_overloaded(self, left: Overloaded) -> bool:
        if isinstance(self.right, Overloaded):
            return is_same_types(left.items, self.right.items)
        else:
            return False

    def visit_partial_type(self, left: PartialType) -> bool:
        # A partial type is not fully defined, so the result is indeterminate. We shouldn't
        # get here.
        raise RuntimeError

    def visit_type_type(self, left: TypeType) -> bool:
        if isinstance(self.right, TypeType):
            return is_same_type(left.item, self.right.item)
        else:
            return False
