from __future__ import annotations

from mypy.erasetype import erase_typevars
from mypy.nodes import TypeInfo
from mypy.types import (
    AnyType,
    Instance,
    ParamSpecType,
    ProperType,
    TupleType,
    Type,
    TypeOfAny,
    TypeVarLikeType,
    TypeVarTupleType,
    TypeVarType,
    UnpackType,
)


def fill_typevars(typ: TypeInfo) -> Instance | TupleType:
    """For a non-generic type, return instance type representing the type.

    For a generic G type with parameters T1, .., Tn, return G[T1, ..., Tn].
    """
    tvs: list[Type] = []
    # TODO: why do we need to keep both typ.type_vars and typ.defn.type_vars?
    for i in range(len(typ.defn.type_vars)):
        tv: TypeVarLikeType | UnpackType = typ.defn.type_vars[i]
        # Change the line number
        if isinstance(tv, TypeVarType):
            tv = tv.copy_modified(line=-1, column=-1)
        elif isinstance(tv, TypeVarTupleType):
            tv = UnpackType(
                TypeVarTupleType(
                    tv.name,
                    tv.fullname,
                    tv.id,
                    tv.upper_bound,
                    tv.tuple_fallback,
                    tv.default,
                    line=-1,
                    column=-1,
                )
            )
        else:
            assert isinstance(tv, ParamSpecType)
            tv = ParamSpecType(
                tv.name,
                tv.fullname,
                tv.id,
                tv.flavor,
                tv.upper_bound,
                tv.default,
                line=-1,
                column=-1,
            )
        tvs.append(tv)
    inst = Instance(typ, tvs)
    # TODO: do we need to also handle typeddict_type here and below?
    if typ.tuple_type is None:
        return inst
    return typ.tuple_type.copy_modified(fallback=inst)


def fill_typevars_with_any(typ: TypeInfo) -> Instance | TupleType:
    """Apply a correct number of Any's as type arguments to a type."""
    args: list[Type] = []
    for tv in typ.defn.type_vars:
        # Valid erasure for *Ts is *tuple[Any, ...], not just Any.
        if isinstance(tv, TypeVarTupleType):
            args.append(
                UnpackType(
                    tv.tuple_fallback.copy_modified(
                        args=[AnyType(TypeOfAny.special_form)]
                    )
                )
            )
        else:
            args.append(AnyType(TypeOfAny.special_form))
    inst = Instance(typ, args)
    if typ.tuple_type is None:
        return inst
    erased_tuple_type = erase_typevars(
        typ.tuple_type, {tv.id for tv in typ.defn.type_vars}
    )
    assert isinstance(erased_tuple_type, ProperType)
    if isinstance(erased_tuple_type, TupleType):
        return typ.tuple_type.copy_modified(fallback=inst)
    return inst


def has_no_typevars(typ: Type) -> bool:
    # We test if a type contains type variables by erasing all type variables
    # and comparing the result to the original type. We use comparison by equality that
    # in turn uses `__eq__` defined for types. Note: we can't use `is_same_type` because
    # it is not safe with unresolved forward references, while this function may be called
    # before forward references resolution patch pass. Note also that it is not safe to use
    # `is` comparison because `erase_typevars` doesn't preserve type identity.
    return typ == erase_typevars(typ)
