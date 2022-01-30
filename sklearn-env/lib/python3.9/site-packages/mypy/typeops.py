"""Miscellaneous type operations and helpers for use during type checking.

NOTE: These must not be accessed from mypy.nodes or mypy.types to avoid import
      cycles. These must not be called from the semantic analysis main pass
      since these may assume that MROs are ready.
"""

from typing import cast, Optional, List, Sequence, Set, Iterable, TypeVar, Dict, Tuple, Any
from typing_extensions import Type as TypingType
import itertools
import sys

from mypy.types import (
    TupleType, Instance, FunctionLike, Type, CallableType, TypeVarLikeType, Overloaded,
    TypeVarType, UninhabitedType, FormalArgument, UnionType, NoneType, TypedDictType,
    AnyType, TypeOfAny, TypeType, ProperType, LiteralType, get_proper_type, get_proper_types,
    copy_type, TypeAliasType, TypeQuery, ParamSpecType
)
from mypy.nodes import (
    FuncBase, FuncItem, FuncDef, OverloadedFuncDef, TypeInfo, ARG_STAR, ARG_STAR2, ARG_POS,
    Expression, StrExpr, Var, Decorator, SYMBOL_FUNCBASE_TYPES
)
from mypy.maptype import map_instance_to_supertype
from mypy.expandtype import expand_type_by_instance, expand_type

from mypy.typevars import fill_typevars

from mypy import state


def is_recursive_pair(s: Type, t: Type) -> bool:
    """Is this a pair of recursive type aliases?"""
    return (isinstance(s, TypeAliasType) and isinstance(t, TypeAliasType) and
            s.is_recursive and t.is_recursive)


def tuple_fallback(typ: TupleType) -> Instance:
    """Return fallback type for a tuple."""
    from mypy.join import join_type_list

    info = typ.partial_fallback.type
    if info.fullname != 'builtins.tuple':
        return typ.partial_fallback
    return Instance(info, [join_type_list(typ.items)])


def try_getting_instance_fallback(typ: ProperType) -> Optional[Instance]:
    """Returns the Instance fallback for this type if one exists.

    Otherwise, returns None.
    """
    if isinstance(typ, Instance):
        return typ
    elif isinstance(typ, TupleType):
        return tuple_fallback(typ)
    elif isinstance(typ, TypedDictType):
        return typ.fallback
    elif isinstance(typ, FunctionLike):
        return typ.fallback
    elif isinstance(typ, LiteralType):
        return typ.fallback
    else:
        return None


def type_object_type_from_function(signature: FunctionLike,
                                   info: TypeInfo,
                                   def_info: TypeInfo,
                                   fallback: Instance,
                                   is_new: bool) -> FunctionLike:
    # We first need to record all non-trivial (explicit) self types in __init__,
    # since they will not be available after we bind them. Note, we use explicit
    # self-types only in the defining class, similar to __new__ (but not exactly the same,
    # see comment in class_callable below). This is mostly useful for annotating library
    # classes such as subprocess.Popen.
    default_self = fill_typevars(info)
    if not is_new and not info.is_newtype:
        orig_self_types = [(it.arg_types[0] if it.arg_types and it.arg_types[0] != default_self
                            and it.arg_kinds[0] == ARG_POS else None) for it in signature.items]
    else:
        orig_self_types = [None] * len(signature.items)

    # The __init__ method might come from a generic superclass 'def_info'
    # with type variables that do not map identically to the type variables of
    # the class 'info' being constructed. For example:
    #
    #   class A(Generic[T]):
    #       def __init__(self, x: T) -> None: ...
    #   class B(A[List[T]]):
    #      ...
    #
    # We need to map B's __init__ to the type (List[T]) -> None.
    signature = bind_self(signature, original_type=default_self, is_classmethod=is_new)
    signature = cast(FunctionLike, map_type_from_supertype(signature, info, def_info))

    special_sig: Optional[str] = None
    if def_info.fullname == 'builtins.dict':
        # Special signature!
        special_sig = 'dict'

    if isinstance(signature, CallableType):
        return class_callable(signature, info, fallback, special_sig, is_new, orig_self_types[0])
    else:
        # Overloaded __init__/__new__.
        assert isinstance(signature, Overloaded)
        items: List[CallableType] = []
        for item, orig_self in zip(signature.items, orig_self_types):
            items.append(class_callable(item, info, fallback, special_sig, is_new, orig_self))
        return Overloaded(items)


def class_callable(init_type: CallableType, info: TypeInfo, type_type: Instance,
                   special_sig: Optional[str],
                   is_new: bool, orig_self_type: Optional[Type] = None) -> CallableType:
    """Create a type object type based on the signature of __init__."""
    variables: List[TypeVarLikeType] = []
    variables.extend(info.defn.type_vars)
    variables.extend(init_type.variables)

    from mypy.subtypes import is_subtype

    init_ret_type = get_proper_type(init_type.ret_type)
    orig_self_type = get_proper_type(orig_self_type)
    default_ret_type = fill_typevars(info)
    explicit_type = init_ret_type if is_new else orig_self_type
    if (
        isinstance(explicit_type, (Instance, TupleType))
        # Only use the declared return type from __new__ or declared self in __init__
        # if it is actually returning a subtype of what we would return otherwise.
        and is_subtype(explicit_type, default_ret_type, ignore_type_params=True)
    ):
        ret_type: Type = explicit_type
    else:
        ret_type = default_ret_type

    callable_type = init_type.copy_modified(
        ret_type=ret_type, fallback=type_type, name=None, variables=variables,
        special_sig=special_sig)
    c = callable_type.with_name(info.name)
    return c


def map_type_from_supertype(typ: Type,
                            sub_info: TypeInfo,
                            super_info: TypeInfo) -> Type:
    """Map type variables in a type defined in a supertype context to be valid
    in the subtype context. Assume that the result is unique; if more than
    one type is possible, return one of the alternatives.

    For example, assume

      class D(Generic[S]): ...
      class C(D[E[T]], Generic[T]): ...

    Now S in the context of D would be mapped to E[T] in the context of C.
    """
    # Create the type of self in subtype, of form t[a1, ...].
    inst_type = fill_typevars(sub_info)
    if isinstance(inst_type, TupleType):
        inst_type = tuple_fallback(inst_type)
    # Map the type of self to supertype. This gets us a description of the
    # supertype type variables in terms of subtype variables, i.e. t[t1, ...]
    # so that any type variables in tN are to be interpreted in subtype
    # context.
    inst_type = map_instance_to_supertype(inst_type, super_info)
    # Finally expand the type variables in type with those in the previously
    # constructed type. Note that both type and inst_type may have type
    # variables, but in type they are interpreted in supertype context while
    # in inst_type they are interpreted in subtype context. This works even if
    # the names of type variables in supertype and subtype overlap.
    return expand_type_by_instance(typ, inst_type)


def supported_self_type(typ: ProperType) -> bool:
    """Is this a supported kind of explicit self-types?

    Currently, this means a X or Type[X], where X is an instance or
    a type variable with an instance upper bound.
    """
    if isinstance(typ, TypeType):
        return supported_self_type(typ.item)
    return (isinstance(typ, TypeVarType) or
            (isinstance(typ, Instance) and typ != fill_typevars(typ.type)))


F = TypeVar('F', bound=FunctionLike)


def bind_self(method: F, original_type: Optional[Type] = None, is_classmethod: bool = False) -> F:
    """Return a copy of `method`, with the type of its first parameter (usually
    self or cls) bound to original_type.

    If the type of `self` is a generic type (T, or Type[T] for classmethods),
    instantiate every occurrence of type with original_type in the rest of the
    signature and in the return type.

    original_type is the type of E in the expression E.copy(). It is None in
    compatibility checks. In this case we treat it as the erasure of the
    declared type of self.

    This way we can express "the type of self". For example:

    T = TypeVar('T', bound='A')
    class A:
        def copy(self: T) -> T: ...

    class B(A): pass

    b = B().copy()  # type: B

    """
    if isinstance(method, Overloaded):
        return cast(F, Overloaded([bind_self(c, original_type, is_classmethod)
                                   for c in method.items]))
    assert isinstance(method, CallableType)
    func = method
    if not func.arg_types:
        # Invalid method, return something.
        return cast(F, func)
    if func.arg_kinds[0] == ARG_STAR:
        # The signature is of the form 'def foo(*args, ...)'.
        # In this case we shouldn't drop the first arg,
        # since func will be absorbed by the *args.

        # TODO: infer bounds on the type of *args?
        return cast(F, func)
    self_param_type = get_proper_type(func.arg_types[0])

    variables: Sequence[TypeVarLikeType] = []
    if func.variables and supported_self_type(self_param_type):
        from mypy.infer import infer_type_arguments

        if original_type is None:
            # TODO: type check method override (see #7861).
            original_type = erase_to_bound(self_param_type)
        original_type = get_proper_type(original_type)

        all_ids = func.type_var_ids()
        typeargs = infer_type_arguments(all_ids, self_param_type, original_type,
                                        is_supertype=True)
        if (is_classmethod
                # TODO: why do we need the extra guards here?
                and any(isinstance(get_proper_type(t), UninhabitedType) for t in typeargs)
                and isinstance(original_type, (Instance, TypeVarType, TupleType))):
            # In case we call a classmethod through an instance x, fallback to type(x)
            typeargs = infer_type_arguments(all_ids, self_param_type, TypeType(original_type),
                                            is_supertype=True)

        ids = [tid for tid in all_ids
               if any(tid == t.id for t in get_type_vars(self_param_type))]

        # Technically, some constrains might be unsolvable, make them <nothing>.
        to_apply = [t if t is not None else UninhabitedType() for t in typeargs]

        def expand(target: Type) -> Type:
            return expand_type(target, {id: to_apply[all_ids.index(id)] for id in ids})

        arg_types = [expand(x) for x in func.arg_types[1:]]
        ret_type = expand(func.ret_type)
        variables = [v for v in func.variables if v.id not in ids]
    else:
        arg_types = func.arg_types[1:]
        ret_type = func.ret_type
        variables = func.variables

    original_type = get_proper_type(original_type)
    if isinstance(original_type, CallableType) and original_type.is_type_obj():
        original_type = TypeType.make_normalized(original_type.ret_type)
    res = func.copy_modified(arg_types=arg_types,
                             arg_kinds=func.arg_kinds[1:],
                             arg_names=func.arg_names[1:],
                             variables=variables,
                             ret_type=ret_type,
                             bound_args=[original_type])
    return cast(F, res)


def erase_to_bound(t: Type) -> Type:
    # TODO: use value restrictions to produce a union?
    t = get_proper_type(t)
    if isinstance(t, TypeVarType):
        return t.upper_bound
    if isinstance(t, TypeType):
        if isinstance(t.item, TypeVarType):
            return TypeType.make_normalized(t.item.upper_bound)
    return t


def callable_corresponding_argument(typ: CallableType,
                                    model: FormalArgument) -> Optional[FormalArgument]:
    """Return the argument a function that corresponds to `model`"""

    by_name = typ.argument_by_name(model.name)
    by_pos = typ.argument_by_position(model.pos)
    if by_name is None and by_pos is None:
        return None
    if by_name is not None and by_pos is not None:
        if by_name == by_pos:
            return by_name
        # If we're dealing with an optional pos-only and an optional
        # name-only arg, merge them.  This is the case for all functions
        # taking both *args and **args, or a pair of functions like so:

        # def right(a: int = ...) -> None: ...
        # def left(__a: int = ..., *, a: int = ...) -> None: ...
        from mypy.subtypes import is_equivalent

        if (not (by_name.required or by_pos.required)
                and by_pos.name is None
                and by_name.pos is None
                and is_equivalent(by_name.typ, by_pos.typ)):
            return FormalArgument(by_name.name, by_pos.pos, by_name.typ, False)
    return by_name if by_name is not None else by_pos


def is_simple_literal(t: ProperType) -> bool:
    """
    Whether a type is a simple enough literal to allow for fast Union simplification

    For now this means enum or string
    """
    return isinstance(t, LiteralType) and (
            t.fallback.type.is_enum or t.fallback.type.fullname == 'builtins.str'
    )


def make_simplified_union(items: Sequence[Type],
                          line: int = -1, column: int = -1,
                          *, keep_erased: bool = False,
                          contract_literals: bool = True) -> ProperType:
    """Build union type with redundant union items removed.

    If only a single item remains, this may return a non-union type.

    Examples:

    * [int, str] -> Union[int, str]
    * [int, object] -> object
    * [int, int] -> int
    * [int, Any] -> Union[int, Any] (Any types are not simplified away!)
    * [Any, Any] -> Any

    Note: This must NOT be used during semantic analysis, since TypeInfos may not
          be fully initialized.
    The keep_erased flag is used for type inference against union types
    containing type variables. If set to True, keep all ErasedType items.
    """
    items = get_proper_types(items)
    while any(isinstance(typ, UnionType) for typ in items):
        all_items: List[ProperType] = []
        for typ in items:
            if isinstance(typ, UnionType):
                all_items.extend(get_proper_types(typ.items))
            else:
                all_items.append(typ)
        items = all_items

    from mypy.subtypes import is_proper_subtype

    removed: Set[int] = set()
    seen: Set[Tuple[str, str]] = set()

    # NB: having a separate fast path for Union of Literal and slow path for other things
    # would arguably be cleaner, however it breaks down when simplifying the Union of two
    # different enum types as try_expanding_enum_to_union works recursively and will
    # trigger intermediate simplifications that would render the fast path useless
    for i, item in enumerate(items):
        if i in removed:
            continue
        # Avoid slow nested for loop for Union of Literal of strings/enums (issue #9169)
        if is_simple_literal(item):
            assert isinstance(item, LiteralType)
            assert isinstance(item.value, str)
            k = (item.value, item.fallback.type.fullname)
            if k in seen:
                removed.add(i)
                continue

            # NB: one would naively expect that it would be safe to skip the slow path
            # always for literals. One would be sorely mistaken. Indeed, some simplifications
            # such as that of None/Optional when strict optional is false, do require that we
            # proceed with the slow path. Thankfully, all literals will have the same subtype
            # relationship to non-literal types, so we only need to do that walk for the first
            # literal, which keeps the fast path fast even in the presence of a mixture of
            # literals and other types.
            safe_skip = len(seen) > 0
            seen.add(k)
            if safe_skip:
                continue
        # Keep track of the truishness info for deleted subtypes which can be relevant
        cbt = cbf = False
        for j, tj in enumerate(items):
            # NB: we don't need to check literals as the fast path above takes care of that
            if (
                    i != j
                    and not is_simple_literal(tj)
                    and is_proper_subtype(tj, item, keep_erased_types=keep_erased)
                    and is_redundant_literal_instance(item, tj)  # XXX?
            ):
                # We found a redundant item in the union.
                removed.add(j)
                cbt = cbt or tj.can_be_true
                cbf = cbf or tj.can_be_false
        # if deleted subtypes had more general truthiness, use that
        if not item.can_be_true and cbt:
            items[i] = true_or_false(item)
        elif not item.can_be_false and cbf:
            items[i] = true_or_false(item)

    simplified_set = [items[i] for i in range(len(items)) if i not in removed]

    # If more than one literal exists in the union, try to simplify
    if (contract_literals and sum(isinstance(item, LiteralType) for item in simplified_set) > 1):
        simplified_set = try_contracting_literals_in_union(simplified_set)

    return UnionType.make_union(simplified_set, line, column)


def _get_type_special_method_bool_ret_type(t: Type) -> Optional[Type]:
    t = get_proper_type(t)

    if isinstance(t, Instance):
        bool_method = t.type.get("__bool__")
        if bool_method:
            callee = get_proper_type(bool_method.type)
            if isinstance(callee, CallableType):
                return callee.ret_type

    return None


def true_only(t: Type) -> ProperType:
    """
    Restricted version of t with only True-ish values
    """
    t = get_proper_type(t)

    if not t.can_be_true:
        # All values of t are False-ish, so there are no true values in it
        return UninhabitedType(line=t.line, column=t.column)
    elif not t.can_be_false:
        # All values of t are already True-ish, so true_only is idempotent in this case
        return t
    elif isinstance(t, UnionType):
        # The true version of a union type is the union of the true versions of its components
        new_items = [true_only(item) for item in t.items]
        can_be_true_items = [item for item in new_items if item.can_be_true]
        return make_simplified_union(can_be_true_items, line=t.line, column=t.column)
    else:
        ret_type = _get_type_special_method_bool_ret_type(t)

        if ret_type and ret_type.can_be_false and not ret_type.can_be_true:
            new_t = copy_type(t)
            new_t.can_be_true = False
            return new_t

        new_t = copy_type(t)
        new_t.can_be_false = False
        return new_t


def false_only(t: Type) -> ProperType:
    """
    Restricted version of t with only False-ish values
    """
    t = get_proper_type(t)

    if not t.can_be_false:
        if state.strict_optional:
            # All values of t are True-ish, so there are no false values in it
            return UninhabitedType(line=t.line)
        else:
            # When strict optional checking is disabled, everything can be
            # False-ish since anything can be None
            return NoneType(line=t.line)
    elif not t.can_be_true:
        # All values of t are already False-ish, so false_only is idempotent in this case
        return t
    elif isinstance(t, UnionType):
        # The false version of a union type is the union of the false versions of its components
        new_items = [false_only(item) for item in t.items]
        can_be_false_items = [item for item in new_items if item.can_be_false]
        return make_simplified_union(can_be_false_items, line=t.line, column=t.column)
    else:
        ret_type = _get_type_special_method_bool_ret_type(t)

        if ret_type and ret_type.can_be_true and not ret_type.can_be_false:
            new_t = copy_type(t)
            new_t.can_be_false = False
            return new_t

        new_t = copy_type(t)
        new_t.can_be_true = False
        return new_t


def true_or_false(t: Type) -> ProperType:
    """
    Unrestricted version of t with both True-ish and False-ish values
    """
    t = get_proper_type(t)

    if isinstance(t, UnionType):
        new_items = [true_or_false(item) for item in t.items]
        return make_simplified_union(new_items, line=t.line, column=t.column)

    new_t = copy_type(t)
    new_t.can_be_true = new_t.can_be_true_default()
    new_t.can_be_false = new_t.can_be_false_default()
    return new_t


def erase_def_to_union_or_bound(tdef: TypeVarLikeType) -> Type:
    # TODO(PEP612): fix for ParamSpecType
    if isinstance(tdef, ParamSpecType):
        return AnyType(TypeOfAny.from_error)
    assert isinstance(tdef, TypeVarType)
    if tdef.values:
        return make_simplified_union(tdef.values)
    else:
        return tdef.upper_bound


def erase_to_union_or_bound(typ: TypeVarType) -> ProperType:
    if typ.values:
        return make_simplified_union(typ.values)
    else:
        return get_proper_type(typ.upper_bound)


def function_type(func: FuncBase, fallback: Instance) -> FunctionLike:
    if func.type:
        assert isinstance(func.type, FunctionLike)
        return func.type
    else:
        # Implicit type signature with dynamic types.
        if isinstance(func, FuncItem):
            return callable_type(func, fallback)
        else:
            # Broken overloads can have self.type set to None.
            # TODO: should we instead always set the type in semantic analyzer?
            assert isinstance(func, OverloadedFuncDef)
            any_type = AnyType(TypeOfAny.from_error)
            dummy = CallableType([any_type, any_type],
                                 [ARG_STAR, ARG_STAR2],
                                 [None, None], any_type,
                                 fallback,
                                 line=func.line, is_ellipsis_args=True)
            # Return an Overloaded, because some callers may expect that
            # an OverloadedFuncDef has an Overloaded type.
            return Overloaded([dummy])


def callable_type(fdef: FuncItem, fallback: Instance,
                  ret_type: Optional[Type] = None) -> CallableType:
    # TODO: somewhat unfortunate duplication with prepare_method_signature in semanal
    if fdef.info and not fdef.is_static and fdef.arg_names:
        self_type: Type = fill_typevars(fdef.info)
        if fdef.is_class or fdef.name == '__new__':
            self_type = TypeType.make_normalized(self_type)
        args = [self_type] + [AnyType(TypeOfAny.unannotated)] * (len(fdef.arg_names)-1)
    else:
        args = [AnyType(TypeOfAny.unannotated)] * len(fdef.arg_names)

    return CallableType(
        args,
        fdef.arg_kinds,
        fdef.arg_names,
        ret_type or AnyType(TypeOfAny.unannotated),
        fallback,
        name=fdef.name,
        line=fdef.line,
        column=fdef.column,
        implicit=True,
        # We need this for better error messages, like missing `self` note:
        definition=fdef if isinstance(fdef, FuncDef) else None,
    )


def try_getting_str_literals(expr: Expression, typ: Type) -> Optional[List[str]]:
    """If the given expression or type corresponds to a string literal
    or a union of string literals, returns a list of the underlying strings.
    Otherwise, returns None.

    Specifically, this function is guaranteed to return a list with
    one or more strings if one one the following is true:

    1. 'expr' is a StrExpr
    2. 'typ' is a LiteralType containing a string
    3. 'typ' is a UnionType containing only LiteralType of strings
    """
    if isinstance(expr, StrExpr):
        return [expr.value]

    # TODO: See if we can eliminate this function and call the below one directly
    return try_getting_str_literals_from_type(typ)


def try_getting_str_literals_from_type(typ: Type) -> Optional[List[str]]:
    """If the given expression or type corresponds to a string Literal
    or a union of string Literals, returns a list of the underlying strings.
    Otherwise, returns None.

    For example, if we had the type 'Literal["foo", "bar"]' as input, this function
    would return a list of strings ["foo", "bar"].
    """
    return try_getting_literals_from_type(typ, str, "builtins.str")


def try_getting_int_literals_from_type(typ: Type) -> Optional[List[int]]:
    """If the given expression or type corresponds to an int Literal
    or a union of int Literals, returns a list of the underlying ints.
    Otherwise, returns None.

    For example, if we had the type 'Literal[1, 2, 3]' as input, this function
    would return a list of ints [1, 2, 3].
    """
    return try_getting_literals_from_type(typ, int, "builtins.int")


T = TypeVar('T')


def try_getting_literals_from_type(typ: Type,
                                   target_literal_type: TypingType[T],
                                   target_fullname: str) -> Optional[List[T]]:
    """If the given expression or type corresponds to a Literal or
    union of Literals where the underlying values corresponds to the given
    target type, returns a list of those underlying values. Otherwise,
    returns None.
    """
    typ = get_proper_type(typ)

    if isinstance(typ, Instance) and typ.last_known_value is not None:
        possible_literals: List[Type] = [typ.last_known_value]
    elif isinstance(typ, UnionType):
        possible_literals = list(typ.items)
    else:
        possible_literals = [typ]

    literals: List[T] = []
    for lit in get_proper_types(possible_literals):
        if isinstance(lit, LiteralType) and lit.fallback.type.fullname == target_fullname:
            val = lit.value
            if isinstance(val, target_literal_type):
                literals.append(val)
            else:
                return None
        else:
            return None
    return literals


def is_literal_type_like(t: Optional[Type]) -> bool:
    """Returns 'true' if the given type context is potentially either a LiteralType,
    a Union of LiteralType, or something similar.
    """
    t = get_proper_type(t)
    if t is None:
        return False
    elif isinstance(t, LiteralType):
        return True
    elif isinstance(t, UnionType):
        return any(is_literal_type_like(item) for item in t.items)
    elif isinstance(t, TypeVarType):
        return (is_literal_type_like(t.upper_bound)
                or any(is_literal_type_like(item) for item in t.values))
    else:
        return False


def get_enum_values(typ: Instance) -> List[str]:
    """Return the list of values for an Enum."""
    return [name for name, sym in typ.type.names.items() if isinstance(sym.node, Var)]


def is_singleton_type(typ: Type) -> bool:
    """Returns 'true' if this type is a "singleton type" -- if there exists
    exactly only one runtime value associated with this type.

    That is, given two values 'a' and 'b' that have the same type 't',
    'is_singleton_type(t)' returns True if and only if the expression 'a is b' is
    always true.

    Currently, this returns True when given NoneTypes, enum LiteralTypes and
    enum types with a single value.

    Note that other kinds of LiteralTypes cannot count as singleton types. For
    example, suppose we do 'a = 100000 + 1' and 'b = 100001'. It is not guaranteed
    that 'a is b' will always be true -- some implementations of Python will end up
    constructing two distinct instances of 100001.
    """
    typ = get_proper_type(typ)
    # TODO:
    # Also make this return True if the type corresponds to ... (ellipsis) or NotImplemented?
    return (
            isinstance(typ, NoneType)
            or (isinstance(typ, LiteralType)
                and (typ.is_enum_literal() or isinstance(typ.value, bool)))
            or (isinstance(typ, Instance) and typ.type.is_enum and len(get_enum_values(typ)) == 1)
    )


def try_expanding_sum_type_to_union(typ: Type, target_fullname: str) -> ProperType:
    """Attempts to recursively expand any enum Instances with the given target_fullname
    into a Union of all of its component LiteralTypes.

    For example, if we have:

        class Color(Enum):
            RED = 1
            BLUE = 2
            YELLOW = 3

        class Status(Enum):
            SUCCESS = 1
            FAILURE = 2
            UNKNOWN = 3

    ...and if we call `try_expanding_enum_to_union(Union[Color, Status], 'module.Color')`,
    this function will return Literal[Color.RED, Color.BLUE, Color.YELLOW, Status].
    """
    typ = get_proper_type(typ)

    if isinstance(typ, UnionType):
        items = [
            try_expanding_sum_type_to_union(item, target_fullname)
            for item in typ.relevant_items()
        ]
        return make_simplified_union(items, contract_literals=False)
    elif isinstance(typ, Instance) and typ.type.fullname == target_fullname:
        if typ.type.is_enum:
            new_items = []
            for name, symbol in typ.type.names.items():
                if not isinstance(symbol.node, Var):
                    continue
                # Skip "_order_" and "__order__", since Enum will remove it
                if name in ("_order_", "__order__"):
                    continue
                new_items.append(LiteralType(name, typ))
            # SymbolTables are really just dicts, and dicts are guaranteed to preserve
            # insertion order only starting with Python 3.7. So, we sort these for older
            # versions of Python to help make tests deterministic.
            #
            # We could probably skip the sort for Python 3.6 since people probably run mypy
            # only using CPython, but we might as well for the sake of full correctness.
            if sys.version_info < (3, 7):
                new_items.sort(key=lambda lit: lit.value)
            return make_simplified_union(new_items, contract_literals=False)
        elif typ.type.fullname == "builtins.bool":
            return make_simplified_union(
                [LiteralType(True, typ), LiteralType(False, typ)],
                contract_literals=False
            )

    return typ


def try_contracting_literals_in_union(types: Sequence[Type]) -> List[ProperType]:
    """Contracts any literal types back into a sum type if possible.

    Will replace the first instance of the literal with the sum type and
    remove all others.

    If we call `try_contracting_union(Literal[Color.RED, Color.BLUE, Color.YELLOW])`,
    this function will return Color.

    We also treat `Literal[True, False]` as `bool`.
    """
    proper_types = [get_proper_type(typ) for typ in types]
    sum_types: Dict[str, Tuple[Set[Any], List[int]]] = {}
    marked_for_deletion = set()
    for idx, typ in enumerate(proper_types):
        if isinstance(typ, LiteralType):
            fullname = typ.fallback.type.fullname
            if typ.fallback.type.is_enum or isinstance(typ.value, bool):
                if fullname not in sum_types:
                    sum_types[fullname] = (set(get_enum_values(typ.fallback))
                                           if typ.fallback.type.is_enum
                                           else set((True, False)),
                                           [])
                literals, indexes = sum_types[fullname]
                literals.discard(typ.value)
                indexes.append(idx)
                if not literals:
                    first, *rest = indexes
                    proper_types[first] = typ.fallback
                    marked_for_deletion |= set(rest)
    return list(itertools.compress(proper_types, [(i not in marked_for_deletion)
                                                  for i in range(len(proper_types))]))


def coerce_to_literal(typ: Type) -> Type:
    """Recursively converts any Instances that have a last_known_value or are
    instances of enum types with a single value into the corresponding LiteralType.
    """
    original_type = typ
    typ = get_proper_type(typ)
    if isinstance(typ, UnionType):
        new_items = [coerce_to_literal(item) for item in typ.items]
        return make_simplified_union(new_items)
    elif isinstance(typ, Instance):
        if typ.last_known_value:
            return typ.last_known_value
        elif typ.type.is_enum:
            enum_values = get_enum_values(typ)
            if len(enum_values) == 1:
                return LiteralType(value=enum_values[0], fallback=typ)
    return original_type


def get_type_vars(tp: Type) -> List[TypeVarType]:
    return tp.accept(TypeVarExtractor())


class TypeVarExtractor(TypeQuery[List[TypeVarType]]):
    def __init__(self) -> None:
        super().__init__(self._merge)

    def _merge(self, iter: Iterable[List[TypeVarType]]) -> List[TypeVarType]:
        out = []
        for item in iter:
            out.extend(item)
        return out

    def visit_type_var(self, t: TypeVarType) -> List[TypeVarType]:
        return [t]


def custom_special_method(typ: Type, name: str, check_all: bool = False) -> bool:
    """Does this type have a custom special method such as __format__() or __eq__()?

    If check_all is True ensure all items of a union have a custom method, not just some.
    """
    typ = get_proper_type(typ)
    if isinstance(typ, Instance):
        method = typ.type.get(name)
        if method and isinstance(method.node, (SYMBOL_FUNCBASE_TYPES, Decorator, Var)):
            if method.node.info:
                return not method.node.info.fullname.startswith('builtins.')
        return False
    if isinstance(typ, UnionType):
        if check_all:
            return all(custom_special_method(t, name, check_all) for t in typ.items)
        return any(custom_special_method(t, name) for t in typ.items)
    if isinstance(typ, TupleType):
        return custom_special_method(tuple_fallback(typ), name, check_all)
    if isinstance(typ, CallableType) and typ.is_type_obj():
        # Look up __method__ on the metaclass for class objects.
        return custom_special_method(typ.fallback, name, check_all)
    if isinstance(typ, AnyType):
        # Avoid false positives in uncertain cases.
        return True
    # TODO: support other types (see ExpressionChecker.has_member())?
    return False


def is_redundant_literal_instance(general: ProperType, specific: ProperType) -> bool:
    if not isinstance(general, Instance) or general.last_known_value is None:
        return True
    if isinstance(specific, Instance) and specific.last_known_value == general.last_known_value:
        return True
    if isinstance(specific, UninhabitedType):
        return True

    return False
