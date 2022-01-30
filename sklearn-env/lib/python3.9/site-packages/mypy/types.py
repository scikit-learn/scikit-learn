"""Classes for representing mypy types."""

import copy
import sys
from abc import abstractmethod

from typing import (
    Any, TypeVar, Dict, List, Tuple, cast, Set, Optional, Union, Iterable, NamedTuple,
    Sequence
)
from typing_extensions import ClassVar, Final, TYPE_CHECKING, overload, TypeAlias as _TypeAlias

from mypy.backports import OrderedDict
import mypy.nodes
from mypy import state
from mypy.nodes import (
    INVARIANT, SymbolNode, FuncDef,
    ArgKind, ARG_POS, ARG_STAR, ARG_STAR2,
)
from mypy.util import IdMapper
from mypy.bogus_type import Bogus


T = TypeVar('T')

JsonDict: _TypeAlias = Dict[str, Any]

# The set of all valid expressions that can currently be contained
# inside of a Literal[...].
#
# Literals can contain bytes and enum-values: we special-case both of these
# and store the value as a string. We rely on the fallback type that's also
# stored with the Literal to determine how a string is being used.
#
# TODO: confirm that we're happy with representing enums (and the
# other types) in the manner described above.
#
# Note: if we change the set of types included below, we must also
# make sure to audit the following methods:
#
# 1. types.LiteralType's serialize and deserialize methods: this method
#    needs to make sure it can convert the below types into JSON and back.
#
# 2. types.LiteralType's 'alue_repr` method: this method is ultimately used
#    by TypeStrVisitor's visit_literal_type to generate a reasonable
#    repr-able output.
#
# 3. server.astdiff.SnapshotTypeVisitor's visit_literal_type_method: this
#    method assumes that the following types supports equality checks and
#    hashability.
#
# Note: Although "Literal[None]" is a valid type, we internally always convert
# such a type directly into "None". So, "None" is not a valid parameter of
# LiteralType and is omitted from this list.
LiteralValue: _TypeAlias = Union[int, str, bool]


# If we only import type_visitor in the middle of the file, mypy
# breaks, and if we do it at the top, it breaks at runtime because of
# import cycle issues, so we do it at the top while typechecking and
# then again in the middle at runtime.
# We should be able to remove this once we are switched to the new
# semantic analyzer!
if TYPE_CHECKING:
    from mypy.type_visitor import (
        TypeVisitor as TypeVisitor,
        SyntheticTypeVisitor as SyntheticTypeVisitor,
    )

# Supported names of TypedDict type constructors.
TPDICT_NAMES: Final = (
    "typing.TypedDict",
    "typing_extensions.TypedDict",
    "mypy_extensions.TypedDict",
)

# Supported fallback instance type names for TypedDict types.
TPDICT_FB_NAMES: Final = (
    "typing._TypedDict",
    "typing_extensions._TypedDict",
    "mypy_extensions._TypedDict",
)

# A placeholder used for Bogus[...] parameters
_dummy: Final[Any] = object()


class TypeOfAny:
    """
    This class describes different types of Any. Each 'Any' can be of only one type at a time.
    """

    __slots__ = ()

    # Was this Any type inferred without a type annotation?
    unannotated: Final = 1
    # Does this Any come from an explicit type annotation?
    explicit: Final = 2
    # Does this come from an unfollowed import? See --disallow-any-unimported option
    from_unimported_type: Final = 3
    # Does this Any type come from omitted generics?
    from_omitted_generics: Final = 4
    # Does this Any come from an error?
    from_error: Final = 5
    # Is this a type that can't be represented in mypy's type system? For instance, type of
    # call to NewType...). Even though these types aren't real Anys, we treat them as such.
    # Also used for variables named '_'.
    special_form: Final = 6
    # Does this Any come from interaction with another Any?
    from_another_any: Final = 7
    # Does this Any come from an implementation limitation/bug?
    implementation_artifact: Final = 8
    # Does this Any come from use in the suggestion engine?  This is
    # used to ignore Anys inserted by the suggestion engine when
    # generating constraints.
    suggestion_engine: Final = 9


def deserialize_type(data: Union[JsonDict, str]) -> 'Type':
    if isinstance(data, str):
        return Instance.deserialize(data)
    classname = data['.class']
    method = deserialize_map.get(classname)
    if method is not None:
        return method(data)
    raise NotImplementedError('unexpected .class {}'.format(classname))


class Type(mypy.nodes.Context):
    """Abstract base class for all types."""

    __slots__ = ('can_be_true', 'can_be_false')
    # 'can_be_true' and 'can_be_false' mean whether the value of the
    # expression can be true or false in a boolean context. They are useful
    # when inferring the type of logic expressions like `x and y`.
    #
    # For example:
    #   * the literal `False` can't be true while `True` can.
    #   * a value with type `bool` can be true or false.
    #   * `None` can't be true
    #   * ...

    def __init__(self, line: int = -1, column: int = -1) -> None:
        super().__init__(line, column)
        self.can_be_true = self.can_be_true_default()
        self.can_be_false = self.can_be_false_default()

    def can_be_true_default(self) -> bool:
        return True

    def can_be_false_default(self) -> bool:
        return True

    def accept(self, visitor: 'TypeVisitor[T]') -> T:
        raise RuntimeError('Not implemented')

    def __repr__(self) -> str:
        return self.accept(TypeStrVisitor())

    def serialize(self) -> Union[JsonDict, str]:
        raise NotImplementedError('Cannot serialize {} instance'.format(self.__class__.__name__))

    @classmethod
    def deserialize(cls, data: JsonDict) -> 'Type':
        raise NotImplementedError('Cannot deserialize {} instance'.format(cls.__name__))


class TypeAliasType(Type):
    """A type alias to another type.

    NOTE: this is not being used yet, and the implementation is still incomplete.

    To support recursive type aliases we don't immediately expand a type alias
    during semantic analysis, but create an instance of this type that records the target alias
    definition node (mypy.nodes.TypeAlias) and type arguments (for generic aliases).

    This is very similar to how TypeInfo vs Instance interact, where a recursive class-based
    structure like
        class Node:
            value: int
            children: List[Node]
    can be represented in a tree-like manner.
    """

    __slots__ = ('alias', 'args', 'line', 'column', 'type_ref')

    def __init__(self, alias: Optional[mypy.nodes.TypeAlias], args: List[Type],
                 line: int = -1, column: int = -1) -> None:
        self.alias = alias
        self.args = args
        self.type_ref: Optional[str] = None
        super().__init__(line, column)

    def _expand_once(self) -> Type:
        """Expand to the target type exactly once.

        This doesn't do full expansion, i.e. the result can contain another
        (or even this same) type alias. Use this internal helper only when really needed,
        its public wrapper mypy.types.get_proper_type() is preferred.
        """
        assert self.alias is not None
        if self.alias.no_args:
            # We know that no_args=True aliases like L = List must have an instance
            # as their target.
            assert isinstance(self.alias.target, Instance)  # type: ignore[misc]
            return self.alias.target.copy_modified(args=self.args)
        return replace_alias_tvars(self.alias.target, self.alias.alias_tvars, self.args,
                                   self.line, self.column)

    def _partial_expansion(self) -> Tuple['ProperType', bool]:
        # Private method mostly for debugging and testing.
        unroller = UnrollAliasVisitor(set())
        unrolled = self.accept(unroller)
        assert isinstance(unrolled, ProperType)
        return unrolled, unroller.recursed

    def expand_all_if_possible(self) -> Optional['ProperType']:
        """Attempt a full expansion of the type alias (including nested aliases).

        If the expansion is not possible, i.e. the alias is (mutually-)recursive,
        return None.
        """
        unrolled, recursed = self._partial_expansion()
        if recursed:
            return None
        return unrolled

    @property
    def is_recursive(self) -> bool:
        assert self.alias is not None, 'Unfixed type alias'
        is_recursive = self.alias._is_recursive
        if is_recursive is None:
            is_recursive = self.expand_all_if_possible() is None
            # We cache the value on the underlying TypeAlias node as an optimization,
            # since the value is the same for all instances of the same alias.
            self.alias._is_recursive = is_recursive
        return is_recursive

    def can_be_true_default(self) -> bool:
        if self.alias is not None:
            return self.alias.target.can_be_true
        return super().can_be_true_default()

    def can_be_false_default(self) -> bool:
        if self.alias is not None:
            return self.alias.target.can_be_false
        return super().can_be_false_default()

    def accept(self, visitor: 'TypeVisitor[T]') -> T:
        return visitor.visit_type_alias_type(self)

    def __hash__(self) -> int:
        return hash((self.alias, tuple(self.args)))

    def __eq__(self, other: object) -> bool:
        # Note: never use this to determine subtype relationships, use is_subtype().
        if not isinstance(other, TypeAliasType):
            return NotImplemented
        return (self.alias == other.alias
                and self.args == other.args)

    def serialize(self) -> JsonDict:
        assert self.alias is not None
        data: JsonDict = {
            ".class": "TypeAliasType",
            "type_ref": self.alias.fullname,
            "args": [arg.serialize() for arg in self.args],
        }
        return data

    @classmethod
    def deserialize(cls, data: JsonDict) -> 'TypeAliasType':
        assert data['.class'] == 'TypeAliasType'
        args: List[Type] = []
        if 'args' in data:
            args_list = data['args']
            assert isinstance(args_list, list)
            args = [deserialize_type(arg) for arg in args_list]
        alias = TypeAliasType(None, args)
        alias.type_ref = data['type_ref']
        return alias

    def copy_modified(self, *,
                      args: Optional[List[Type]] = None) -> 'TypeAliasType':
        return TypeAliasType(
            self.alias,
            args if args is not None else self.args.copy(),
            self.line, self.column)


class TypeGuardedType(Type):
    """Only used by find_isinstance_check() etc."""

    __slots__ = ('type_guard',)

    def __init__(self, type_guard: Type):
        super().__init__(line=type_guard.line, column=type_guard.column)
        self.type_guard = type_guard

    def __repr__(self) -> str:
        return "TypeGuard({})".format(self.type_guard)


class RequiredType(Type):
    """Required[T] or NotRequired[T]. Only usable at top-level of a TypedDict definition."""

    def __init__(self, item: Type, *, required: bool) -> None:
        super().__init__(line=item.line, column=item.column)
        self.item = item
        self.required = required

    def __repr__(self) -> str:
        if self.required:
            return "Required[{}]".format(self.item)
        else:
            return "NotRequired[{}]".format(self.item)


class ProperType(Type):
    """Not a type alias.

    Every type except TypeAliasType must inherit from this type.
    """

    __slots__ = ()


class TypeVarId:
    # A type variable is uniquely identified by its raw id and meta level.

    # For plain variables (type parameters of generic classes and
    # functions) raw ids are allocated by semantic analysis, using
    # positive ids 1, 2, ... for generic class parameters and negative
    # ids -1, ... for generic function type arguments. This convention
    # is only used to keep type variable ids distinct when allocating
    # them; the type checker makes no distinction between class and
    # function type variables.

    # Metavariables are allocated unique ids starting from 1.
    raw_id: int = 0

    # Level of the variable in type inference. Currently either 0 for
    # declared types, or 1 for type inference metavariables.
    meta_level: int = 0

    # Class variable used for allocating fresh ids for metavariables.
    next_raw_id: ClassVar[int] = 1

    def __init__(self, raw_id: int, meta_level: int = 0) -> None:
        self.raw_id = raw_id
        self.meta_level = meta_level

    @staticmethod
    def new(meta_level: int) -> 'TypeVarId':
        raw_id = TypeVarId.next_raw_id
        TypeVarId.next_raw_id += 1
        return TypeVarId(raw_id, meta_level)

    def __repr__(self) -> str:
        return self.raw_id.__repr__()

    def __eq__(self, other: object) -> bool:
        if isinstance(other, TypeVarId):
            return (self.raw_id == other.raw_id and
                    self.meta_level == other.meta_level)
        else:
            return False

    def __ne__(self, other: object) -> bool:
        return not (self == other)

    def __hash__(self) -> int:
        return hash((self.raw_id, self.meta_level))

    def is_meta_var(self) -> bool:
        return self.meta_level > 0


class TypeVarLikeType(ProperType):

    __slots__ = ('name', 'fullname', 'id', 'upper_bound')

    name: str  # Name (may be qualified)
    fullname: str  # Fully qualified name
    id: TypeVarId
    upper_bound: Type

    def __init__(
        self, name: str, fullname: str, id: Union[TypeVarId, int], upper_bound: Type,
        line: int = -1, column: int = -1
    ) -> None:
        super().__init__(line, column)
        self.name = name
        self.fullname = fullname
        if isinstance(id, int):
            id = TypeVarId(id)
        self.id = id
        self.upper_bound = upper_bound

    def serialize(self) -> JsonDict:
        raise NotImplementedError

    @classmethod
    def deserialize(cls, data: JsonDict) -> 'TypeVarLikeType':
        raise NotImplementedError


class TypeVarType(TypeVarLikeType):
    """Type that refers to a type variable."""

    __slots__ = ('values', 'variance')

    values: List[Type]  # Value restriction, empty list if no restriction
    variance: int

    def __init__(self, name: str, fullname: str, id: Union[TypeVarId, int], values: List[Type],
                 upper_bound: Type, variance: int = INVARIANT, line: int = -1,
                 column: int = -1) -> None:
        super().__init__(name, fullname, id, upper_bound, line, column)
        assert values is not None, "No restrictions must be represented by empty list"
        self.values = values
        self.variance = variance

    @staticmethod
    def new_unification_variable(old: 'TypeVarType') -> 'TypeVarType':
        new_id = TypeVarId.new(meta_level=1)
        return TypeVarType(old.name, old.fullname, new_id, old.values,
                          old.upper_bound, old.variance, old.line, old.column)

    def accept(self, visitor: 'TypeVisitor[T]') -> T:
        return visitor.visit_type_var(self)

    def __hash__(self) -> int:
        return hash(self.id)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, TypeVarType):
            return NotImplemented
        return self.id == other.id

    def serialize(self) -> JsonDict:
        assert not self.id.is_meta_var()
        return {'.class': 'TypeVarType',
                'name': self.name,
                'fullname': self.fullname,
                'id': self.id.raw_id,
                'values': [v.serialize() for v in self.values],
                'upper_bound': self.upper_bound.serialize(),
                'variance': self.variance,
                }

    @classmethod
    def deserialize(cls, data: JsonDict) -> 'TypeVarType':
        assert data['.class'] == 'TypeVarType'
        return TypeVarType(
            data['name'],
            data['fullname'],
            data['id'],
            [deserialize_type(v) for v in data['values']],
            deserialize_type(data['upper_bound']),
            data['variance'],
        )


class ParamSpecFlavor:
    # Simple ParamSpec reference such as "P"
    BARE: Final = 0
    # P.args
    ARGS: Final = 1
    # P.kwargs
    KWARGS: Final = 2


class ParamSpecType(TypeVarLikeType):
    """Type that refers to a ParamSpec.

    A ParamSpec is a type variable that represents the parameter
    types, names and kinds of a callable (i.e., the signature without
    the return type).

    This can be one of these forms
     * P (ParamSpecFlavor.BARE)
     * P.args (ParamSpecFlavor.ARGS)
     * P.kwargs (ParamSpecFLavor.KWARGS)

    The upper_bound is really used as a fallback type -- it's shared
    with TypeVarType for simplicity. It can't be specified by the user
    and the value is directly derived from the flavor (currently
    always just 'object').
    """

    __slots__ = ('flavor',)

    flavor: int

    def __init__(
         self, name: str, fullname: str, id: Union[TypeVarId, int], flavor: int,
         upper_bound: Type, *, line: int = -1, column: int = -1
    ) -> None:
        super().__init__(name, fullname, id, upper_bound, line=line, column=column)
        self.flavor = flavor

    @staticmethod
    def new_unification_variable(old: 'ParamSpecType') -> 'ParamSpecType':
        new_id = TypeVarId.new(meta_level=1)
        return ParamSpecType(old.name, old.fullname, new_id, old.flavor, old.upper_bound,
                             line=old.line, column=old.column)

    def with_flavor(self, flavor: int) -> 'ParamSpecType':
        return ParamSpecType(self.name, self.fullname, self.id, flavor,
                             upper_bound=self.upper_bound)

    def accept(self, visitor: 'TypeVisitor[T]') -> T:
        return visitor.visit_param_spec(self)

    def name_with_suffix(self) -> str:
        n = self.name
        if self.flavor == ParamSpecFlavor.ARGS:
            return f'{n}.args'
        elif self.flavor == ParamSpecFlavor.KWARGS:
            return f'{n}.kwargs'
        return n

    def __hash__(self) -> int:
        return hash((self.id, self.flavor))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ParamSpecType):
            return NotImplemented
        # Upper bound can be ignored, since it's determined by flavor.
        return self.id == other.id and self.flavor == other.flavor

    def serialize(self) -> JsonDict:
        assert not self.id.is_meta_var()
        return {
            '.class': 'ParamSpecType',
            'name': self.name,
            'fullname': self.fullname,
            'id': self.id.raw_id,
            'flavor': self.flavor,
            'upper_bound': self.upper_bound.serialize(),
        }

    @classmethod
    def deserialize(cls, data: JsonDict) -> 'ParamSpecType':
        assert data['.class'] == 'ParamSpecType'
        return ParamSpecType(
            data['name'],
            data['fullname'],
            data['id'],
            data['flavor'],
            deserialize_type(data['upper_bound']),
        )


class UnboundType(ProperType):
    """Instance type that has not been bound during semantic analysis."""

    __slots__ = ('name', 'args', 'optional', 'empty_tuple_index',
                 'original_str_expr', 'original_str_fallback')

    def __init__(self,
                 name: Optional[str],
                 args: Optional[Sequence[Type]] = None,
                 line: int = -1,
                 column: int = -1,
                 optional: bool = False,
                 empty_tuple_index: bool = False,
                 original_str_expr: Optional[str] = None,
                 original_str_fallback: Optional[str] = None,
                 ) -> None:
        super().__init__(line, column)
        if not args:
            args = []
        assert name is not None
        self.name = name
        self.args = tuple(args)
        # Should this type be wrapped in an Optional?
        self.optional = optional
        # Special case for X[()]
        self.empty_tuple_index = empty_tuple_index
        # If this UnboundType was originally defined as a str or bytes, keep track of
        # the original contents of that string-like thing. This way, if this UnboundExpr
        # ever shows up inside of a LiteralType, we can determine whether that
        # Literal[...] is valid or not. E.g. Literal[foo] is most likely invalid
        # (unless 'foo' is an alias for another literal or something) and
        # Literal["foo"] most likely is.
        #
        # We keep track of the entire string instead of just using a boolean flag
        # so we can distinguish between things like Literal["foo"] vs
        # Literal["    foo   "].
        #
        # We also keep track of what the original base fallback type was supposed to be
        # so we don't have to try and recompute it later
        self.original_str_expr = original_str_expr
        self.original_str_fallback = original_str_fallback

    def copy_modified(self,
                      args: Bogus[Optional[Sequence[Type]]] = _dummy,
                      ) -> 'UnboundType':
        if args is _dummy:
            args = self.args
        return UnboundType(
            name=self.name,
            args=args,
            line=self.line,
            column=self.column,
            optional=self.optional,
            empty_tuple_index=self.empty_tuple_index,
            original_str_expr=self.original_str_expr,
            original_str_fallback=self.original_str_fallback,
        )

    def accept(self, visitor: 'TypeVisitor[T]') -> T:
        return visitor.visit_unbound_type(self)

    def __hash__(self) -> int:
        return hash((self.name, self.optional, tuple(self.args), self.original_str_expr))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, UnboundType):
            return NotImplemented
        return (self.name == other.name and self.optional == other.optional and
                self.args == other.args and self.original_str_expr == other.original_str_expr and
                self.original_str_fallback == other.original_str_fallback)

    def serialize(self) -> JsonDict:
        return {'.class': 'UnboundType',
                'name': self.name,
                'args': [a.serialize() for a in self.args],
                'expr': self.original_str_expr,
                'expr_fallback': self.original_str_fallback,
                }

    @classmethod
    def deserialize(cls, data: JsonDict) -> 'UnboundType':
        assert data['.class'] == 'UnboundType'
        return UnboundType(data['name'],
                           [deserialize_type(a) for a in data['args']],
                           original_str_expr=data['expr'],
                           original_str_fallback=data['expr_fallback'],
                           )


class CallableArgument(ProperType):
    """Represents a Arg(type, 'name') inside a Callable's type list.

    Note that this is a synthetic type for helping parse ASTs, not a real type.
    """

    __slots__ = ('typ', 'name', 'constructor')

    typ: Type
    name: Optional[str]
    constructor: Optional[str]

    def __init__(self, typ: Type, name: Optional[str], constructor: Optional[str],
                 line: int = -1, column: int = -1) -> None:
        super().__init__(line, column)
        self.typ = typ
        self.name = name
        self.constructor = constructor

    def accept(self, visitor: 'TypeVisitor[T]') -> T:
        assert isinstance(visitor, SyntheticTypeVisitor)
        return visitor.visit_callable_argument(self)

    def serialize(self) -> JsonDict:
        assert False, "Synthetic types don't serialize"


class TypeList(ProperType):
    """Information about argument types and names [...].

    This is used for the arguments of a Callable type, i.e. for
    [arg, ...] in Callable[[arg, ...], ret]. This is not a real type
    but a syntactic AST construct. UnboundTypes can also have TypeList
    types before they are processed into Callable types.
    """

    __slots__ = ('items',)

    items: List[Type]

    def __init__(self, items: List[Type], line: int = -1, column: int = -1) -> None:
        super().__init__(line, column)
        self.items = items

    def accept(self, visitor: 'TypeVisitor[T]') -> T:
        assert isinstance(visitor, SyntheticTypeVisitor)
        return visitor.visit_type_list(self)

    def serialize(self) -> JsonDict:
        assert False, "Synthetic types don't serialize"


class AnyType(ProperType):
    """The type 'Any'."""

    __slots__ = ('type_of_any', 'source_any', 'missing_import_name')

    def __init__(self,
                 type_of_any: int,
                 source_any: Optional['AnyType'] = None,
                 missing_import_name: Optional[str] = None,
                 line: int = -1,
                 column: int = -1) -> None:
        super().__init__(line, column)
        self.type_of_any = type_of_any
        # If this Any was created as a result of interacting with another 'Any', record the source
        # and use it in reports.
        self.source_any = source_any
        if source_any and source_any.source_any:
            self.source_any = source_any.source_any

        if source_any is None:
            self.missing_import_name = missing_import_name
        else:
            self.missing_import_name = source_any.missing_import_name

        # Only unimported type anys and anys from other anys should have an import name
        assert (missing_import_name is None or
                type_of_any in (TypeOfAny.from_unimported_type, TypeOfAny.from_another_any))
        # Only Anys that come from another Any can have source_any.
        assert type_of_any != TypeOfAny.from_another_any or source_any is not None
        # We should not have chains of Anys.
        assert not self.source_any or self.source_any.type_of_any != TypeOfAny.from_another_any

    @property
    def is_from_error(self) -> bool:
        return self.type_of_any == TypeOfAny.from_error

    def accept(self, visitor: 'TypeVisitor[T]') -> T:
        return visitor.visit_any(self)

    def copy_modified(self,
                      # Mark with Bogus because _dummy is just an object (with type Any)
                      type_of_any: Bogus[int] = _dummy,
                      original_any: Bogus[Optional['AnyType']] = _dummy,
                      ) -> 'AnyType':
        if type_of_any is _dummy:
            type_of_any = self.type_of_any
        if original_any is _dummy:
            original_any = self.source_any
        return AnyType(type_of_any=type_of_any, source_any=original_any,
                       missing_import_name=self.missing_import_name,
                       line=self.line, column=self.column)

    def __hash__(self) -> int:
        return hash(AnyType)

    def __eq__(self, other: object) -> bool:
        return isinstance(other, AnyType)

    def serialize(self) -> JsonDict:
        return {'.class': 'AnyType', 'type_of_any': self.type_of_any,
                'source_any': self.source_any.serialize() if self.source_any is not None else None,
                'missing_import_name': self.missing_import_name}

    @classmethod
    def deserialize(cls, data: JsonDict) -> 'AnyType':
        assert data['.class'] == 'AnyType'
        source = data['source_any']
        return AnyType(data['type_of_any'],
                       AnyType.deserialize(source) if source is not None else None,
                       data['missing_import_name'])


class UninhabitedType(ProperType):
    """This type has no members.

    This type is the bottom type.
    With strict Optional checking, it is the only common subtype between all
    other types, which allows `meet` to be well defined.  Without strict
    Optional checking, NoneType fills this role.

    In general, for any type T:
        join(UninhabitedType, T) = T
        meet(UninhabitedType, T) = UninhabitedType
        is_subtype(UninhabitedType, T) = True
    """

    __slots__ = ('ambiguous', 'is_noreturn',)

    is_noreturn: bool  # Does this come from a NoReturn?  Purely for error messages.
    # It is important to track whether this is an actual NoReturn type, or just a result
    # of ambiguous type inference, in the latter case we don't want to mark a branch as
    # unreachable in binder.
    ambiguous: bool  # Is this a result of inference for a variable without constraints?

    def __init__(self, is_noreturn: bool = False, line: int = -1, column: int = -1) -> None:
        super().__init__(line, column)
        self.is_noreturn = is_noreturn
        self.ambiguous = False

    def can_be_true_default(self) -> bool:
        return False

    def can_be_false_default(self) -> bool:
        return False

    def accept(self, visitor: 'TypeVisitor[T]') -> T:
        return visitor.visit_uninhabited_type(self)

    def __hash__(self) -> int:
        return hash(UninhabitedType)

    def __eq__(self, other: object) -> bool:
        return isinstance(other, UninhabitedType)

    def serialize(self) -> JsonDict:
        return {'.class': 'UninhabitedType',
                'is_noreturn': self.is_noreturn}

    @classmethod
    def deserialize(cls, data: JsonDict) -> 'UninhabitedType':
        assert data['.class'] == 'UninhabitedType'
        return UninhabitedType(is_noreturn=data['is_noreturn'])


class NoneType(ProperType):
    """The type of 'None'.

    This type can be written by users as 'None'.
    """

    __slots__ = ()

    def __init__(self, line: int = -1, column: int = -1) -> None:
        super().__init__(line, column)

    def can_be_true_default(self) -> bool:
        return False

    def __hash__(self) -> int:
        return hash(NoneType)

    def __eq__(self, other: object) -> bool:
        return isinstance(other, NoneType)

    def accept(self, visitor: 'TypeVisitor[T]') -> T:
        return visitor.visit_none_type(self)

    def serialize(self) -> JsonDict:
        return {'.class': 'NoneType'}

    @classmethod
    def deserialize(cls, data: JsonDict) -> 'NoneType':
        assert data['.class'] == 'NoneType'
        return NoneType()


# NoneType used to be called NoneTyp so to avoid needlessly breaking
# external plugins we keep that alias here.
NoneTyp = NoneType


class ErasedType(ProperType):
    """Placeholder for an erased type.

    This is used during type inference. This has the special property that
    it is ignored during type inference.
    """

    __slots__ = ()

    def accept(self, visitor: 'TypeVisitor[T]') -> T:
        return visitor.visit_erased_type(self)


class DeletedType(ProperType):
    """Type of deleted variables.

    These can be used as lvalues but not rvalues.
    """

    __slots__ = ('source',)

    source: Optional[str]  # May be None; name that generated this value

    def __init__(self, source: Optional[str] = None, line: int = -1, column: int = -1) -> None:
        super().__init__(line, column)
        self.source = source

    def accept(self, visitor: 'TypeVisitor[T]') -> T:
        return visitor.visit_deleted_type(self)

    def serialize(self) -> JsonDict:
        return {'.class': 'DeletedType',
                'source': self.source}

    @classmethod
    def deserialize(cls, data: JsonDict) -> 'DeletedType':
        assert data['.class'] == 'DeletedType'
        return DeletedType(data['source'])


# Fake TypeInfo to be used as a placeholder during Instance de-serialization.
NOT_READY: Final = mypy.nodes.FakeInfo("De-serialization failure: TypeInfo not fixed")


class Instance(ProperType):
    """An instance type of form C[T1, ..., Tn].

    The list of type variables may be empty.
    """

    __slots__ = ('type', 'args', 'erased', 'invalid', 'type_ref', 'last_known_value')

    def __init__(self, typ: mypy.nodes.TypeInfo, args: Sequence[Type],
                 line: int = -1, column: int = -1, erased: bool = False,
                 last_known_value: Optional['LiteralType'] = None) -> None:
        super().__init__(line, column)
        self.type = typ
        self.args = tuple(args)
        self.type_ref: Optional[str] = None

        # True if result of type variable substitution
        self.erased = erased

        # True if recovered after incorrect number of type arguments error
        self.invalid = False

        # This field keeps track of the underlying Literal[...] value associated with
        # this instance, if one is known.
        #
        # This field is set whenever possible within expressions, but is erased upon
        # variable assignment (see erasetype.remove_instance_last_known_values) unless
        # the variable is declared to be final.
        #
        # For example, consider the following program:
        #
        #     a = 1
        #     b: Final[int] = 2
        #     c: Final = 3
        #     print(a + b + c + 4)
        #
        # The 'Instance' objects associated with the expressions '1', '2', '3', and '4' will
        # have last_known_values of type Literal[1], Literal[2], Literal[3], and Literal[4]
        # respectively. However, the Instance object assigned to 'a' and 'b' will have their
        # last_known_value erased: variable 'a' is mutable; variable 'b' was declared to be
        # specifically an int.
        #
        # Or more broadly, this field lets this Instance "remember" its original declaration
        # when applicable. We want this behavior because we want implicit Final declarations
        # to act pretty much identically with constants: we should be able to replace any
        # places where we use some Final variable with the original value and get the same
        # type-checking behavior. For example, we want this program:
        #
        #    def expects_literal(x: Literal[3]) -> None: pass
        #    var: Final = 3
        #    expects_literal(var)
        #
        # ...to type-check in the exact same way as if we had written the program like this:
        #
        #    def expects_literal(x: Literal[3]) -> None: pass
        #    expects_literal(3)
        #
        # In order to make this work (especially with literal types), we need var's type
        # (an Instance) to remember the "original" value.
        #
        # Preserving this value within expressions is useful for similar reasons.
        #
        # Currently most of mypy will ignore this field and will continue to treat this type like
        # a regular Instance. We end up using this field only when we are explicitly within a
        # Literal context.
        self.last_known_value = last_known_value

    def accept(self, visitor: 'TypeVisitor[T]') -> T:
        return visitor.visit_instance(self)

    def __hash__(self) -> int:
        return hash((self.type, tuple(self.args), self.last_known_value))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Instance):
            return NotImplemented
        return (self.type == other.type
                and self.args == other.args
                and self.last_known_value == other.last_known_value)

    def serialize(self) -> Union[JsonDict, str]:
        assert self.type is not None
        type_ref = self.type.fullname
        if not self.args and not self.last_known_value:
            return type_ref
        data: JsonDict = {
            ".class": "Instance",
        }
        data["type_ref"] = type_ref
        data["args"] = [arg.serialize() for arg in self.args]
        if self.last_known_value is not None:
            data['last_known_value'] = self.last_known_value.serialize()
        return data

    @classmethod
    def deserialize(cls, data: Union[JsonDict, str]) -> 'Instance':
        if isinstance(data, str):
            inst = Instance(NOT_READY, [])
            inst.type_ref = data
            return inst
        assert data['.class'] == 'Instance'
        args: List[Type] = []
        if 'args' in data:
            args_list = data['args']
            assert isinstance(args_list, list)
            args = [deserialize_type(arg) for arg in args_list]
        inst = Instance(NOT_READY, args)
        inst.type_ref = data['type_ref']  # Will be fixed up by fixup.py later.
        if 'last_known_value' in data:
            inst.last_known_value = LiteralType.deserialize(data['last_known_value'])
        return inst

    def copy_modified(self, *,
                      args: Bogus[List[Type]] = _dummy,
                      erased: Bogus[bool] = _dummy,
                      last_known_value: Bogus[Optional['LiteralType']] = _dummy) -> 'Instance':
        return Instance(
            self.type,
            args if args is not _dummy else self.args,
            self.line,
            self.column,
            erased if erased is not _dummy else self.erased,
            last_known_value if last_known_value is not _dummy else self.last_known_value,
        )

    def has_readable_member(self, name: str) -> bool:
        return self.type.has_readable_member(name)


class FunctionLike(ProperType):
    """Abstract base class for function types."""

    __slots__ = ('fallback',)

    def __init__(self, line: int = -1, column: int = -1) -> None:
        super().__init__(line, column)
        self.can_be_false = False
        self.fallback: Instance

    @abstractmethod
    def is_type_obj(self) -> bool: pass

    @abstractmethod
    def type_object(self) -> mypy.nodes.TypeInfo: pass

    @property
    @abstractmethod
    def items(self) -> List['CallableType']: pass

    @abstractmethod
    def with_name(self, name: str) -> 'FunctionLike': pass

    @abstractmethod
    def get_name(self) -> Optional[str]: pass


FormalArgument = NamedTuple('FormalArgument', [
    ('name', Optional[str]),
    ('pos', Optional[int]),
    ('typ', Type),
    ('required', bool)])


class CallableType(FunctionLike):
    """Type of a non-overloaded callable object (such as function)."""

    __slots__ = ('arg_types',  # Types of function arguments
                 'arg_kinds',  # ARG_ constants
                 'arg_names',  # Argument names; None if not a keyword argument
                 'min_args',  # Minimum number of arguments; derived from arg_kinds
                 'ret_type',  # Return value type
                 'name',  # Name (may be None; for error messages and plugins)
                 'definition',  # For error messages.  May be None.
                 'variables',  # Type variables for a generic function
                 'is_ellipsis_args',  # Is this Callable[..., t] (with literal '...')?
                 'is_classmethod_class',  # Is this callable constructed for the benefit
                                          # of a classmethod's 'cls' argument?
                 'implicit',  # Was this type implicitly generated instead of explicitly
                              # specified by the user?
                 'special_sig',  # Non-None for signatures that require special handling
                                 # (currently only value is 'dict' for a signature similar to
                                 # 'dict')
                 'from_type_type',  # Was this callable generated by analyzing Type[...]
                                    # instantiation?
                 'bound_args',  # Bound type args, mostly unused but may be useful for
                                # tools that consume mypy ASTs
                 'def_extras',  # Information about original definition we want to serialize.
                                # This is used for more detailed error messages.
                 'type_guard',  # T, if -> TypeGuard[T] (ret_type is bool in this case).
                 )

    def __init__(self,
                 arg_types: Sequence[Type],
                 arg_kinds: List[ArgKind],
                 arg_names: Sequence[Optional[str]],
                 ret_type: Type,
                 fallback: Instance,
                 name: Optional[str] = None,
                 definition: Optional[SymbolNode] = None,
                 variables: Optional[Sequence[TypeVarLikeType]] = None,
                 line: int = -1,
                 column: int = -1,
                 is_ellipsis_args: bool = False,
                 implicit: bool = False,
                 special_sig: Optional[str] = None,
                 from_type_type: bool = False,
                 bound_args: Sequence[Optional[Type]] = (),
                 def_extras: Optional[Dict[str, Any]] = None,
                 type_guard: Optional[Type] = None,
                 ) -> None:
        super().__init__(line, column)
        assert len(arg_types) == len(arg_kinds) == len(arg_names)
        if variables is None:
            variables = []
        self.arg_types = list(arg_types)
        self.arg_kinds = arg_kinds
        self.arg_names = list(arg_names)
        self.min_args = arg_kinds.count(ARG_POS)
        self.ret_type = ret_type
        self.fallback = fallback
        assert not name or '<bound method' not in name
        self.name = name
        self.definition = definition
        self.variables = variables
        self.is_ellipsis_args = is_ellipsis_args
        self.implicit = implicit
        self.special_sig = special_sig
        self.from_type_type = from_type_type
        if not bound_args:
            bound_args = ()
        self.bound_args = bound_args
        if def_extras:
            self.def_extras = def_extras
        elif isinstance(definition, FuncDef):
            # This information would be lost if we don't have definition
            # after serialization, but it is useful in error messages.
            # TODO: decide how to add more info here (file, line, column)
            # without changing interface hash.
            self.def_extras = {
                'first_arg': (
                    definition.arguments[0].variable.name
                    if (getattr(definition, 'arguments', None)
                        and definition.arg_names
                        and definition.info
                        and not definition.is_static)
                    else None
                ),
            }
        else:
            self.def_extras = {}
        self.type_guard = type_guard

    def copy_modified(self,
                      arg_types: Bogus[Sequence[Type]] = _dummy,
                      arg_kinds: Bogus[List[ArgKind]] = _dummy,
                      arg_names: Bogus[List[Optional[str]]] = _dummy,
                      ret_type: Bogus[Type] = _dummy,
                      fallback: Bogus[Instance] = _dummy,
                      name: Bogus[Optional[str]] = _dummy,
                      definition: Bogus[SymbolNode] = _dummy,
                      variables: Bogus[Sequence[TypeVarLikeType]] = _dummy,
                      line: Bogus[int] = _dummy,
                      column: Bogus[int] = _dummy,
                      is_ellipsis_args: Bogus[bool] = _dummy,
                      implicit: Bogus[bool] = _dummy,
                      special_sig: Bogus[Optional[str]] = _dummy,
                      from_type_type: Bogus[bool] = _dummy,
                      bound_args: Bogus[List[Optional[Type]]] = _dummy,
                      def_extras: Bogus[Dict[str, Any]] = _dummy,
                      type_guard: Bogus[Optional[Type]] = _dummy,
                      ) -> 'CallableType':
        return CallableType(
            arg_types=arg_types if arg_types is not _dummy else self.arg_types,
            arg_kinds=arg_kinds if arg_kinds is not _dummy else self.arg_kinds,
            arg_names=arg_names if arg_names is not _dummy else self.arg_names,
            ret_type=ret_type if ret_type is not _dummy else self.ret_type,
            fallback=fallback if fallback is not _dummy else self.fallback,
            name=name if name is not _dummy else self.name,
            definition=definition if definition is not _dummy else self.definition,
            variables=variables if variables is not _dummy else self.variables,
            line=line if line is not _dummy else self.line,
            column=column if column is not _dummy else self.column,
            is_ellipsis_args=(
                is_ellipsis_args if is_ellipsis_args is not _dummy else self.is_ellipsis_args),
            implicit=implicit if implicit is not _dummy else self.implicit,
            special_sig=special_sig if special_sig is not _dummy else self.special_sig,
            from_type_type=from_type_type if from_type_type is not _dummy else self.from_type_type,
            bound_args=bound_args if bound_args is not _dummy else self.bound_args,
            def_extras=def_extras if def_extras is not _dummy else dict(self.def_extras),
            type_guard=type_guard if type_guard is not _dummy else self.type_guard,
        )

    def var_arg(self) -> Optional[FormalArgument]:
        """The formal argument for *args."""
        for position, (type, kind) in enumerate(zip(self.arg_types, self.arg_kinds)):
            if kind == ARG_STAR:
                return FormalArgument(None, position, type, False)
        return None

    def kw_arg(self) -> Optional[FormalArgument]:
        """The formal argument for **kwargs."""
        for position, (type, kind) in enumerate(zip(self.arg_types, self.arg_kinds)):
            if kind == ARG_STAR2:
                return FormalArgument(None, position, type, False)
        return None

    @property
    def is_var_arg(self) -> bool:
        """Does this callable have a *args argument?"""
        return ARG_STAR in self.arg_kinds

    @property
    def is_kw_arg(self) -> bool:
        """Does this callable have a **kwargs argument?"""
        return ARG_STAR2 in self.arg_kinds

    def is_type_obj(self) -> bool:
        return self.fallback.type.is_metaclass()

    def type_object(self) -> mypy.nodes.TypeInfo:
        assert self.is_type_obj()
        ret = get_proper_type(self.ret_type)
        if isinstance(ret, TypeVarType):
            ret = get_proper_type(ret.upper_bound)
        if isinstance(ret, TupleType):
            ret = ret.partial_fallback
        assert isinstance(ret, Instance)
        return ret.type

    def accept(self, visitor: 'TypeVisitor[T]') -> T:
        return visitor.visit_callable_type(self)

    def with_name(self, name: str) -> 'CallableType':
        """Return a copy of this type with the specified name."""
        return self.copy_modified(ret_type=self.ret_type, name=name)

    def get_name(self) -> Optional[str]:
        return self.name

    def max_possible_positional_args(self) -> int:
        """Returns maximum number of positional arguments this method could possibly accept.

        This takes into account *arg and **kwargs but excludes keyword-only args."""
        if self.is_var_arg or self.is_kw_arg:
            return sys.maxsize
        return sum([kind.is_positional() for kind in self.arg_kinds])

    def formal_arguments(self, include_star_args: bool = False) -> List[FormalArgument]:
        """Yields the formal arguments corresponding to this callable, ignoring *arg and **kwargs.

        To handle *args and **kwargs, use the 'callable.var_args' and 'callable.kw_args' fields,
        if they are not None.

        If you really want to include star args in the yielded output, set the
        'include_star_args' parameter to 'True'."""
        args = []
        done_with_positional = False
        for i in range(len(self.arg_types)):
            kind = self.arg_kinds[i]
            if kind.is_named() or kind.is_star():
                done_with_positional = True
            if not include_star_args and kind.is_star():
                continue

            required = kind.is_required()
            pos = None if done_with_positional else i
            arg = FormalArgument(
                self.arg_names[i],
                pos,
                self.arg_types[i],
                required
            )
            args.append(arg)
        return args

    def argument_by_name(self, name: Optional[str]) -> Optional[FormalArgument]:
        if name is None:
            return None
        seen_star = False
        for i, (arg_name, kind, typ) in enumerate(
                zip(self.arg_names, self.arg_kinds, self.arg_types)):
            # No more positional arguments after these.
            if kind.is_named() or kind.is_star():
                seen_star = True
            if kind.is_star():
                continue
            if arg_name == name:
                position = None if seen_star else i
                return FormalArgument(name, position, typ, kind.is_required())
        return self.try_synthesizing_arg_from_kwarg(name)

    def argument_by_position(self, position: Optional[int]) -> Optional[FormalArgument]:
        if position is None:
            return None
        if position >= len(self.arg_names):
            return self.try_synthesizing_arg_from_vararg(position)
        name, kind, typ = (
            self.arg_names[position],
            self.arg_kinds[position],
            self.arg_types[position],
        )
        if kind.is_positional():
            return FormalArgument(name, position, typ, kind == ARG_POS)
        else:
            return self.try_synthesizing_arg_from_vararg(position)

    def try_synthesizing_arg_from_kwarg(self,
                                        name: Optional[str]) -> Optional[FormalArgument]:
        kw_arg = self.kw_arg()
        if kw_arg is not None:
            return FormalArgument(name, None, kw_arg.typ, False)
        else:
            return None

    def try_synthesizing_arg_from_vararg(self,
                                         position: Optional[int]) -> Optional[FormalArgument]:
        var_arg = self.var_arg()
        if var_arg is not None:
            return FormalArgument(None, position, var_arg.typ, False)
        else:
            return None

    @property
    def items(self) -> List['CallableType']:
        return [self]

    def is_generic(self) -> bool:
        return bool(self.variables)

    def type_var_ids(self) -> List[TypeVarId]:
        a: List[TypeVarId] = []
        for tv in self.variables:
            a.append(tv.id)
        return a

    def param_spec(self) -> Optional[ParamSpecType]:
        """Return ParamSpec if callable can be called with one.

        A Callable accepting ParamSpec P args (*args, **kwargs) must have the
        two final parameters like this: *args: P.args, **kwargs: P.kwargs.
        """
        if len(self.arg_types) < 2:
            return None
        if self.arg_kinds[-2] != ARG_STAR or self.arg_kinds[-1] != ARG_STAR2:
            return None
        arg_type = self.arg_types[-2]
        if not isinstance(arg_type, ParamSpecType):
            return None
        return ParamSpecType(arg_type.name, arg_type.fullname, arg_type.id, ParamSpecFlavor.BARE,
                             arg_type.upper_bound)

    def expand_param_spec(self, c: 'CallableType') -> 'CallableType':
        return self.copy_modified(arg_types=self.arg_types[:-2] + c.arg_types,
                                  arg_kinds=self.arg_kinds[:-2] + c.arg_kinds,
                                  arg_names=self.arg_names[:-2] + c.arg_names)

    def __hash__(self) -> int:
        return hash((self.ret_type, self.is_type_obj(),
                     self.is_ellipsis_args, self.name,
                    tuple(self.arg_types), tuple(self.arg_names), tuple(self.arg_kinds)))

    def __eq__(self, other: object) -> bool:
        if isinstance(other, CallableType):
            return (self.ret_type == other.ret_type and
                    self.arg_types == other.arg_types and
                    self.arg_names == other.arg_names and
                    self.arg_kinds == other.arg_kinds and
                    self.name == other.name and
                    self.is_type_obj() == other.is_type_obj() and
                    self.is_ellipsis_args == other.is_ellipsis_args)
        else:
            return NotImplemented

    def serialize(self) -> JsonDict:
        # TODO: As an optimization, leave out everything related to
        # generic functions for non-generic functions.
        return {'.class': 'CallableType',
                'arg_types': [t.serialize() for t in self.arg_types],
                'arg_kinds': [int(x.value) for x in self.arg_kinds],
                'arg_names': self.arg_names,
                'ret_type': self.ret_type.serialize(),
                'fallback': self.fallback.serialize(),
                'name': self.name,
                # We don't serialize the definition (only used for error messages).
                'variables': [v.serialize() for v in self.variables],
                'is_ellipsis_args': self.is_ellipsis_args,
                'implicit': self.implicit,
                'bound_args': [(None if t is None else t.serialize())
                               for t in self.bound_args],
                'def_extras': dict(self.def_extras),
                'type_guard': self.type_guard.serialize() if self.type_guard is not None else None,
                }

    @classmethod
    def deserialize(cls, data: JsonDict) -> 'CallableType':
        assert data['.class'] == 'CallableType'
        # TODO: Set definition to the containing SymbolNode?
        return CallableType(
            [deserialize_type(t) for t in data['arg_types']],
            [ArgKind(x) for x in data['arg_kinds']],
            data['arg_names'],
            deserialize_type(data['ret_type']),
            Instance.deserialize(data['fallback']),
            name=data['name'],
            variables=[cast(TypeVarLikeType, deserialize_type(v)) for v in data['variables']],
            is_ellipsis_args=data['is_ellipsis_args'],
            implicit=data['implicit'],
            bound_args=[(None if t is None else deserialize_type(t)) for t in data['bound_args']],
            def_extras=data['def_extras'],
            type_guard=(deserialize_type(data['type_guard'])
                        if data['type_guard'] is not None else None),
        )


class Overloaded(FunctionLike):
    """Overloaded function type T1, ... Tn, where each Ti is CallableType.

    The variant to call is chosen based on static argument
    types. Overloaded function types can only be defined in stub
    files, and thus there is no explicit runtime dispatch
    implementation.
    """

    __slots__ = ('_items', 'fallback')

    _items: List[CallableType]  # Must not be empty

    def __init__(self, items: List[CallableType]) -> None:
        super().__init__(items[0].line, items[0].column)
        self._items = items
        self.fallback = items[0].fallback

    @property
    def items(self) -> List[CallableType]:
        return self._items

    def name(self) -> Optional[str]:
        return self.get_name()

    def is_type_obj(self) -> bool:
        # All the items must have the same type object status, so it's
        # sufficient to query only (any) one of them.
        return self._items[0].is_type_obj()

    def type_object(self) -> mypy.nodes.TypeInfo:
        # All the items must have the same type object, so it's sufficient to
        # query only (any) one of them.
        return self._items[0].type_object()

    def with_name(self, name: str) -> 'Overloaded':
        ni: List[CallableType] = []
        for it in self._items:
            ni.append(it.with_name(name))
        return Overloaded(ni)

    def get_name(self) -> Optional[str]:
        return self._items[0].name

    def accept(self, visitor: 'TypeVisitor[T]') -> T:
        return visitor.visit_overloaded(self)

    def __hash__(self) -> int:
        return hash(tuple(self.items))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Overloaded):
            return NotImplemented
        return self.items == other.items

    def serialize(self) -> JsonDict:
        return {'.class': 'Overloaded',
                'items': [t.serialize() for t in self.items],
                }

    @classmethod
    def deserialize(cls, data: JsonDict) -> 'Overloaded':
        assert data['.class'] == 'Overloaded'
        return Overloaded([CallableType.deserialize(t) for t in data['items']])


class TupleType(ProperType):
    """The tuple type Tuple[T1, ..., Tn] (at least one type argument).

    Instance variables:
        items: Tuple item types
        partial_fallback: The (imprecise) underlying instance type that is used
            for non-tuple methods. This is generally builtins.tuple[Any] for
            regular tuples, but it's different for named tuples and classes with
            a tuple base class. Use mypy.typeops.tuple_fallback to calculate the
            precise fallback type derived from item types.
        implicit: If True, derived from a tuple expression (t,....) instead of Tuple[t, ...]
    """

    __slots__ = ('items', 'partial_fallback', 'implicit')

    items: List[Type]
    partial_fallback: Instance
    implicit: bool

    def __init__(self, items: List[Type], fallback: Instance, line: int = -1,
                 column: int = -1, implicit: bool = False) -> None:
        super().__init__(line, column)
        self.items = items
        self.partial_fallback = fallback
        self.implicit = implicit
        self.can_be_true = len(self.items) > 0
        self.can_be_false = len(self.items) == 0

    def length(self) -> int:
        return len(self.items)

    def accept(self, visitor: 'TypeVisitor[T]') -> T:
        return visitor.visit_tuple_type(self)

    def __hash__(self) -> int:
        return hash((tuple(self.items), self.partial_fallback))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, TupleType):
            return NotImplemented
        return self.items == other.items and self.partial_fallback == other.partial_fallback

    def serialize(self) -> JsonDict:
        return {'.class': 'TupleType',
                'items': [t.serialize() for t in self.items],
                'partial_fallback': self.partial_fallback.serialize(),
                'implicit': self.implicit,
                }

    @classmethod
    def deserialize(cls, data: JsonDict) -> 'TupleType':
        assert data['.class'] == 'TupleType'
        return TupleType([deserialize_type(t) for t in data['items']],
                         Instance.deserialize(data['partial_fallback']),
                         implicit=data['implicit'])

    def copy_modified(self, *, fallback: Optional[Instance] = None,
                      items: Optional[List[Type]] = None) -> 'TupleType':
        if fallback is None:
            fallback = self.partial_fallback
        if items is None:
            items = self.items
        return TupleType(items, fallback, self.line, self.column)

    def slice(self, begin: Optional[int], end: Optional[int],
              stride: Optional[int]) -> 'TupleType':
        return TupleType(self.items[begin:end:stride], self.partial_fallback,
                         self.line, self.column, self.implicit)


class TypedDictType(ProperType):
    """Type of TypedDict object {'k1': v1, ..., 'kn': vn}.

    A TypedDict object is a dictionary with specific string (literal) keys. Each
    key has a value with a distinct type that depends on the key. TypedDict objects
    are normal dict objects at runtime.

    A TypedDictType can be either named or anonymous. If it's anonymous, its
    fallback will be typing_extensions._TypedDict (Instance). _TypedDict is a subclass
    of Mapping[str, object] and defines all non-mapping dict methods that TypedDict
    supports. Some dict methods are unsafe and not supported. _TypedDict isn't defined
    at runtime.

    If a TypedDict is named, its fallback will be an Instance of the named type
    (ex: "Point") whose TypeInfo has a typeddict_type that is anonymous. This
    is similar to how named tuples work.

    TODO: The fallback structure is perhaps overly complicated.
    """

    __slots__ = ('items', 'required_keys', 'fallback')

    items: "OrderedDict[str, Type]"  # item_name -> item_type
    required_keys: Set[str]
    fallback: Instance

    def __init__(self, items: 'OrderedDict[str, Type]', required_keys: Set[str],
                 fallback: Instance, line: int = -1, column: int = -1) -> None:
        super().__init__(line, column)
        self.items = items
        self.required_keys = required_keys
        self.fallback = fallback
        self.can_be_true = len(self.items) > 0
        self.can_be_false = len(self.required_keys) == 0

    def accept(self, visitor: 'TypeVisitor[T]') -> T:
        return visitor.visit_typeddict_type(self)

    def __hash__(self) -> int:
        return hash((frozenset(self.items.items()), self.fallback,
                     frozenset(self.required_keys)))

    def __eq__(self, other: object) -> bool:
        if isinstance(other, TypedDictType):
            if frozenset(self.items.keys()) != frozenset(other.items.keys()):
                return False
            for (_, left_item_type, right_item_type) in self.zip(other):
                if not left_item_type == right_item_type:
                    return False
            return self.fallback == other.fallback and self.required_keys == other.required_keys
        else:
            return NotImplemented

    def serialize(self) -> JsonDict:
        return {'.class': 'TypedDictType',
                'items': [[n, t.serialize()] for (n, t) in self.items.items()],
                'required_keys': sorted(self.required_keys),
                'fallback': self.fallback.serialize(),
                }

    @classmethod
    def deserialize(cls, data: JsonDict) -> 'TypedDictType':
        assert data['.class'] == 'TypedDictType'
        return TypedDictType(OrderedDict([(n, deserialize_type(t))
                                          for (n, t) in data['items']]),
                             set(data['required_keys']),
                             Instance.deserialize(data['fallback']))

    def is_anonymous(self) -> bool:
        return self.fallback.type.fullname in TPDICT_FB_NAMES

    def as_anonymous(self) -> 'TypedDictType':
        if self.is_anonymous():
            return self
        assert self.fallback.type.typeddict_type is not None
        return self.fallback.type.typeddict_type.as_anonymous()

    def copy_modified(self, *, fallback: Optional[Instance] = None,
                      item_types: Optional[List[Type]] = None,
                      required_keys: Optional[Set[str]] = None) -> 'TypedDictType':
        if fallback is None:
            fallback = self.fallback
        if item_types is None:
            items = self.items
        else:
            items = OrderedDict(zip(self.items, item_types))
        if required_keys is None:
            required_keys = self.required_keys
        return TypedDictType(items, required_keys, fallback, self.line, self.column)

    def create_anonymous_fallback(self, *, value_type: Type) -> Instance:
        anonymous = self.as_anonymous()
        return anonymous.fallback

    def names_are_wider_than(self, other: 'TypedDictType') -> bool:
        return len(other.items.keys() - self.items.keys()) == 0

    def zip(self, right: 'TypedDictType') -> Iterable[Tuple[str, Type, Type]]:
        left = self
        for (item_name, left_item_type) in left.items.items():
            right_item_type = right.items.get(item_name)
            if right_item_type is not None:
                yield (item_name, left_item_type, right_item_type)

    def zipall(self, right: 'TypedDictType') \
            -> Iterable[Tuple[str, Optional[Type], Optional[Type]]]:
        left = self
        for (item_name, left_item_type) in left.items.items():
            right_item_type = right.items.get(item_name)
            yield (item_name, left_item_type, right_item_type)
        for (item_name, right_item_type) in right.items.items():
            if item_name in left.items:
                continue
            yield (item_name, None, right_item_type)


class RawExpressionType(ProperType):
    """A synthetic type representing some arbitrary expression that does not cleanly
    translate into a type.

    This synthetic type is only used at the beginning stages of semantic analysis
    and should be completely removing during the process for mapping UnboundTypes to
    actual types: we either turn it into a LiteralType or an AnyType.

    For example, suppose `Foo[1]` is initially represented as the following:

        UnboundType(
            name='Foo',
            args=[
                RawExpressionType(value=1, base_type_name='builtins.int'),
            ],
        )

    As we perform semantic analysis, this type will transform into one of two
    possible forms.

    If 'Foo' was an alias for 'Literal' all along, this type is transformed into:

        LiteralType(value=1, fallback=int_instance_here)

    Alternatively, if 'Foo' is an unrelated class, we report an error and instead
    produce something like this:

        Instance(type=typeinfo_for_foo, args=[AnyType(TypeOfAny.from_error))

    If the "note" field is not None, the provided note will be reported alongside the
    error at this point.

    Note: if "literal_value" is None, that means this object is representing some
    expression that cannot possibly be a parameter of Literal[...]. For example,
    "Foo[3j]" would be represented as:

        UnboundType(
            name='Foo',
            args=[
                RawExpressionType(value=None, base_type_name='builtins.complex'),
            ],
        )
    """

    __slots__ = ('literal_value', 'base_type_name', 'note')

    def __init__(self,
                 literal_value: Optional[LiteralValue],
                 base_type_name: str,
                 line: int = -1,
                 column: int = -1,
                 note: Optional[str] = None,
                 ) -> None:
        super().__init__(line, column)
        self.literal_value = literal_value
        self.base_type_name = base_type_name
        self.note = note

    def simple_name(self) -> str:
        return self.base_type_name.replace("builtins.", "")

    def accept(self, visitor: 'TypeVisitor[T]') -> T:
        assert isinstance(visitor, SyntheticTypeVisitor)
        return visitor.visit_raw_expression_type(self)

    def serialize(self) -> JsonDict:
        assert False, "Synthetic types don't serialize"

    def __hash__(self) -> int:
        return hash((self.literal_value, self.base_type_name))

    def __eq__(self, other: object) -> bool:
        if isinstance(other, RawExpressionType):
            return (self.base_type_name == other.base_type_name
                    and self.literal_value == other.literal_value)
        else:
            return NotImplemented


class LiteralType(ProperType):
    """The type of a Literal instance. Literal[Value]

    A Literal always consists of:

    1. A native Python object corresponding to the contained inner value
    2. A fallback for this Literal. The fallback also corresponds to the
       parent type this Literal subtypes.

    For example, 'Literal[42]' is represented as
    'LiteralType(value=42, fallback=instance_of_int)'

    As another example, `Literal[Color.RED]` (where Color is an enum) is
    represented as `LiteralType(value="RED", fallback=instance_of_color)'.
    """
    __slots__ = ('value', 'fallback')

    def __init__(self, value: LiteralValue, fallback: Instance,
                 line: int = -1, column: int = -1) -> None:
        self.value = value
        super().__init__(line, column)
        self.fallback = fallback

    def can_be_false_default(self) -> bool:
        return not self.value

    def can_be_true_default(self) -> bool:
        return bool(self.value)

    def accept(self, visitor: 'TypeVisitor[T]') -> T:
        return visitor.visit_literal_type(self)

    def __hash__(self) -> int:
        return hash((self.value, self.fallback))

    def __eq__(self, other: object) -> bool:
        if isinstance(other, LiteralType):
            return self.fallback == other.fallback and self.value == other.value
        else:
            return NotImplemented

    def is_enum_literal(self) -> bool:
        return self.fallback.type.is_enum

    def value_repr(self) -> str:
        """Returns the string representation of the underlying type.

        This function is almost equivalent to running `repr(self.value)`,
        except it includes some additional logic to correctly handle cases
        where the value is a string, byte string, a unicode string, or an enum.
        """
        raw = repr(self.value)
        fallback_name = self.fallback.type.fullname

        # If this is backed by an enum,
        if self.is_enum_literal():
            return '{}.{}'.format(fallback_name, self.value)

        if fallback_name == 'builtins.bytes':
            # Note: 'builtins.bytes' only appears in Python 3, so we want to
            # explicitly prefix with a "b"
            return 'b' + raw
        elif fallback_name == 'builtins.unicode':
            # Similarly, 'builtins.unicode' only appears in Python 2, where we also
            # want to explicitly prefix
            return 'u' + raw
        else:
            # 'builtins.str' could mean either depending on context, but either way
            # we don't prefix: it's the "native" string. And of course, if value is
            # some other type, we just return that string repr directly.
            return raw

    def serialize(self) -> Union[JsonDict, str]:
        return {
            '.class': 'LiteralType',
            'value': self.value,
            'fallback': self.fallback.serialize(),
        }

    @classmethod
    def deserialize(cls, data: JsonDict) -> 'LiteralType':
        assert data['.class'] == 'LiteralType'
        return LiteralType(
            value=data['value'],
            fallback=Instance.deserialize(data['fallback']),
        )


class StarType(ProperType):
    """The star type *type_parameter.

    This is not a real type but a syntactic AST construct.
    """

    __slots__ = ('type',)

    type: Type

    def __init__(self, type: Type, line: int = -1, column: int = -1) -> None:
        super().__init__(line, column)
        self.type = type

    def accept(self, visitor: 'TypeVisitor[T]') -> T:
        assert isinstance(visitor, SyntheticTypeVisitor)
        return visitor.visit_star_type(self)

    def serialize(self) -> JsonDict:
        assert False, "Synthetic types don't serialize"


class UnionType(ProperType):
    """The union type Union[T1, ..., Tn] (at least one type argument)."""

    __slots__ = ('items', 'is_evaluated', 'uses_pep604_syntax')

    def __init__(self, items: Sequence[Type], line: int = -1, column: int = -1,
                 is_evaluated: bool = True, uses_pep604_syntax: bool = False) -> None:
        super().__init__(line, column)
        self.items = flatten_nested_unions(items)
        self.can_be_true = any(item.can_be_true for item in items)
        self.can_be_false = any(item.can_be_false for item in items)
        # is_evaluated should be set to false for type comments and string literals
        self.is_evaluated = is_evaluated
        # uses_pep604_syntax is True if Union uses OR syntax (X | Y)
        self.uses_pep604_syntax = uses_pep604_syntax

    def __hash__(self) -> int:
        return hash(frozenset(self.items))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, UnionType):
            return NotImplemented
        return frozenset(self.items) == frozenset(other.items)

    @overload
    @staticmethod
    def make_union(items: Sequence[ProperType],
                   line: int = -1, column: int = -1) -> ProperType: ...

    @overload
    @staticmethod
    def make_union(items: Sequence[Type], line: int = -1, column: int = -1) -> Type: ...

    @staticmethod
    def make_union(items: Sequence[Type], line: int = -1, column: int = -1) -> Type:
        if len(items) > 1:
            return UnionType(items, line, column)
        elif len(items) == 1:
            return items[0]
        else:
            return UninhabitedType()

    def length(self) -> int:
        return len(self.items)

    def accept(self, visitor: 'TypeVisitor[T]') -> T:
        return visitor.visit_union_type(self)

    def has_readable_member(self, name: str) -> bool:
        """For a tree of unions of instances, check whether all instances have a given member.

        TODO: Deal with attributes of TupleType etc.
        TODO: This should probably be refactored to go elsewhere.
        """
        return all((isinstance(x, UnionType) and x.has_readable_member(name)) or
                   (isinstance(x, Instance) and x.type.has_readable_member(name))
                   for x in get_proper_types(self.relevant_items()))

    def relevant_items(self) -> List[Type]:
        """Removes NoneTypes from Unions when strict Optional checking is off."""
        if state.strict_optional:
            return self.items
        else:
            return [i for i in get_proper_types(self.items) if not isinstance(i, NoneType)]

    def serialize(self) -> JsonDict:
        return {'.class': 'UnionType',
                'items': [t.serialize() for t in self.items],
                }

    @classmethod
    def deserialize(cls, data: JsonDict) -> 'UnionType':
        assert data['.class'] == 'UnionType'
        return UnionType([deserialize_type(t) for t in data['items']])


class PartialType(ProperType):
    """Type such as List[?] where type arguments are unknown, or partial None type.

    These are used for inferring types in multiphase initialization such as this:

      x = []       # x gets a partial type List[?], as item type is unknown
      x.append(1)  # partial type gets replaced with normal type List[int]

    Or with None:

      x = None  # x gets a partial type None
      if c:
          x = 1  # Infer actual type int for x
    """

    __slots__ = ('type', 'var', 'value_type')

    # None for the 'None' partial type; otherwise a generic class
    type: Optional[mypy.nodes.TypeInfo]
    var: mypy.nodes.Var
    # For partial defaultdict[K, V], the type V (K is unknown). If V is generic,
    # the type argument is Any and will be replaced later.
    value_type: Optional[Instance]

    def __init__(self,
                 type: 'Optional[mypy.nodes.TypeInfo]',
                 var: 'mypy.nodes.Var',
                 value_type: 'Optional[Instance]' = None) -> None:
        super().__init__()
        self.type = type
        self.var = var
        self.value_type = value_type

    def accept(self, visitor: 'TypeVisitor[T]') -> T:
        return visitor.visit_partial_type(self)


class EllipsisType(ProperType):
    """The type ... (ellipsis).

    This is not a real type but a syntactic AST construct, used in Callable[..., T], for example.

    A semantically analyzed type will never have ellipsis types.
    """

    __slots__ = ()

    def accept(self, visitor: 'TypeVisitor[T]') -> T:
        assert isinstance(visitor, SyntheticTypeVisitor)
        return visitor.visit_ellipsis_type(self)

    def serialize(self) -> JsonDict:
        assert False, "Synthetic types don't serialize"


class TypeType(ProperType):
    """For types like Type[User].

    This annotates variables that are class objects, constrained by
    the type argument.  See PEP 484 for more details.

    We may encounter expressions whose values are specific classes;
    those are represented as callables (possibly overloaded)
    corresponding to the class's constructor's signature and returning
    an instance of that class.  The difference with Type[C] is that
    those callables always represent the exact class given as the
    return type; Type[C] represents any class that's a subclass of C,
    and C may also be a type variable or a union (or Any).

    Many questions around subtype relationships between Type[C1] and
    def(...) -> C2 are answered by looking at the subtype
    relationships between C1 and C2, since Type[] is considered
    covariant.

    There's an unsolved problem with constructor signatures (also
    unsolved in PEP 484): calling a variable whose type is Type[C]
    assumes the constructor signature for C, even though a subclass of
    C might completely change the constructor signature.  For now we
    just assume that users of Type[C] are careful not to do that (in
    the future we might detect when they are violating that
    assumption).
    """

    __slots__ = ('item',)

    # This can't be everything, but it can be a class reference,
    # a generic class instance, a union, Any, a type variable...
    item: ProperType

    def __init__(self, item: Bogus[Union[Instance, AnyType, TypeVarType, TupleType, NoneType,
                                         CallableType]], *,
                 line: int = -1, column: int = -1) -> None:
        """To ensure Type[Union[A, B]] is always represented as Union[Type[A], Type[B]], item of
        type UnionType must be handled through make_normalized static method.
        """
        super().__init__(line, column)
        self.item = item

    @staticmethod
    def make_normalized(item: Type, *, line: int = -1, column: int = -1) -> ProperType:
        item = get_proper_type(item)
        if isinstance(item, UnionType):
            return UnionType.make_union(
                [TypeType.make_normalized(union_item) for union_item in item.items],
                line=line, column=column
            )
        return TypeType(item, line=line, column=column)  # type: ignore[arg-type]

    def accept(self, visitor: 'TypeVisitor[T]') -> T:
        return visitor.visit_type_type(self)

    def __hash__(self) -> int:
        return hash(self.item)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, TypeType):
            return NotImplemented
        return self.item == other.item

    def serialize(self) -> JsonDict:
        return {'.class': 'TypeType', 'item': self.item.serialize()}

    @classmethod
    def deserialize(cls, data: JsonDict) -> Type:
        assert data['.class'] == 'TypeType'
        return TypeType.make_normalized(deserialize_type(data['item']))


class PlaceholderType(ProperType):
    """Temporary, yet-unknown type during semantic analysis.

    This is needed when there's a reference to a type before the real symbol
    table entry of the target type is available (specifically, we use a
    temporary PlaceholderNode symbol node). Consider this example:

      class str(Sequence[str]): ...

    We use a PlaceholderType for the 'str' in 'Sequence[str]' since we can't create
    a TypeInfo for 'str' until all base classes have been resolved. We'll soon
    perform another analysis iteration which replaces the base class with a complete
    type without any placeholders. After semantic analysis, no placeholder types must
    exist.
    """

    __slots__ = ('fullname', 'args')

    def __init__(self, fullname: Optional[str], args: List[Type], line: int) -> None:
        super().__init__(line)
        self.fullname = fullname  # Must be a valid full name of an actual node (or None).
        self.args = args

    def accept(self, visitor: 'TypeVisitor[T]') -> T:
        assert isinstance(visitor, SyntheticTypeVisitor)
        return visitor.visit_placeholder_type(self)

    def serialize(self) -> str:
        # We should never get here since all placeholders should be replaced
        # during semantic analysis.
        assert False, "Internal error: unresolved placeholder type {}".format(self.fullname)


@overload
def get_proper_type(typ: None) -> None: ...
@overload
def get_proper_type(typ: Type) -> ProperType: ...


def get_proper_type(typ: Optional[Type]) -> Optional[ProperType]:
    """Get the expansion of a type alias type.

    If the type is already a proper type, this is a no-op. Use this function
    wherever a decision is made on a call like e.g. 'if isinstance(typ, UnionType): ...',
    because 'typ' in this case may be an alias to union. Note: if after making the decision
    on the isinstance() call you pass on the original type (and not one of its components)
    it is recommended to *always* pass on the unexpanded alias.
    """
    if typ is None:
        return None
    if isinstance(typ, TypeGuardedType):  # type: ignore[misc]
        typ = typ.type_guard
    while isinstance(typ, TypeAliasType):
        typ = typ._expand_once()
    assert isinstance(typ, ProperType), typ
    # TODO: store the name of original type alias on this type, so we can show it in errors.
    return typ


@overload
def get_proper_types(it: Iterable[Type]) -> List[ProperType]: ...  # type: ignore[misc]
@overload
def get_proper_types(it: Iterable[Optional[Type]]) -> List[Optional[ProperType]]: ...


def get_proper_types(it: Iterable[Optional[Type]]
                     ) -> Union[List[ProperType], List[Optional[ProperType]]]:
    return [get_proper_type(t) for t in it]


# We split off the type visitor base classes to another module
# to make it easier to gradually get modules working with mypyc.
# Import them here, after the types are defined.
# This is intended as a re-export also.
from mypy.type_visitor import (  # noqa
    TypeVisitor as TypeVisitor,
    SyntheticTypeVisitor as SyntheticTypeVisitor,
    TypeTranslator as TypeTranslator,
    TypeQuery as TypeQuery,
)


class TypeStrVisitor(SyntheticTypeVisitor[str]):
    """Visitor for pretty-printing types into strings.

    This is mostly for debugging/testing.

    Do not preserve original formatting.

    Notes:
     - Represent unbound types as Foo? or Foo?[...].
     - Represent the NoneType type as None.
    """

    def __init__(self, id_mapper: Optional[IdMapper] = None) -> None:
        self.id_mapper = id_mapper
        self.any_as_dots = False

    def visit_unbound_type(self, t: UnboundType) -> str:
        s = t.name + '?'
        if t.args:
            s += '[{}]'.format(self.list_str(t.args))
        return s

    def visit_type_list(self, t: TypeList) -> str:
        return '<TypeList {}>'.format(self.list_str(t.items))

    def visit_callable_argument(self, t: CallableArgument) -> str:
        typ = t.typ.accept(self)
        if t.name is None:
            return "{}({})".format(t.constructor, typ)
        else:
            return "{}({}, {})".format(t.constructor, typ, t.name)

    def visit_any(self, t: AnyType) -> str:
        if self.any_as_dots and t.type_of_any == TypeOfAny.special_form:
            return '...'
        return 'Any'

    def visit_none_type(self, t: NoneType) -> str:
        return "None"

    def visit_uninhabited_type(self, t: UninhabitedType) -> str:
        return "<nothing>"

    def visit_erased_type(self, t: ErasedType) -> str:
        return "<Erased>"

    def visit_deleted_type(self, t: DeletedType) -> str:
        if t.source is None:
            return "<Deleted>"
        else:
            return "<Deleted '{}'>".format(t.source)

    def visit_instance(self, t: Instance) -> str:
        if t.last_known_value and not t.args:
            # Instances with a literal fallback should never be generic. If they are,
            # something went wrong so we fall back to showing the full Instance repr.
            s = '{}?'.format(t.last_known_value)
        else:
            s = t.type.fullname or t.type.name or '<???>'

        if t.erased:
            s += '*'
        if t.args:
            s += '[{}]'.format(self.list_str(t.args))
        if self.id_mapper:
            s += '<{}>'.format(self.id_mapper.id(t.type))
        return s

    def visit_type_var(self, t: TypeVarType) -> str:
        if t.name is None:
            # Anonymous type variable type (only numeric id).
            s = '`{}'.format(t.id)
        else:
            # Named type variable type.
            s = '{}`{}'.format(t.name, t.id)
        if self.id_mapper and t.upper_bound:
            s += '(upper_bound={})'.format(t.upper_bound.accept(self))
        return s

    def visit_param_spec(self, t: ParamSpecType) -> str:
        if t.name is None:
            # Anonymous type variable type (only numeric id).
            s = f'`{t.id}'
        else:
            # Named type variable type.
            s = f'{t.name_with_suffix()}`{t.id}'
        return s

    def visit_callable_type(self, t: CallableType) -> str:
        param_spec = t.param_spec()
        if param_spec is not None:
            num_skip = 2
        else:
            num_skip = 0

        s = ''
        bare_asterisk = False
        for i in range(len(t.arg_types) - num_skip):
            if s != '':
                s += ', '
            if t.arg_kinds[i].is_named() and not bare_asterisk:
                s += '*, '
                bare_asterisk = True
            if t.arg_kinds[i] == ARG_STAR:
                s += '*'
            if t.arg_kinds[i] == ARG_STAR2:
                s += '**'
            name = t.arg_names[i]
            if name:
                s += name + ': '
            s += t.arg_types[i].accept(self)
            if t.arg_kinds[i].is_optional():
                s += ' ='

        if param_spec is not None:
            n = param_spec.name
            if s:
                s += ', '
            s += f'*{n}.args, **{n}.kwargs'

        s = '({})'.format(s)

        if not isinstance(get_proper_type(t.ret_type), NoneType):
            if t.type_guard is not None:
                s += ' -> TypeGuard[{}]'.format(t.type_guard.accept(self))
            else:
                s += ' -> {}'.format(t.ret_type.accept(self))

        if t.variables:
            vs = []
            for var in t.variables:
                if isinstance(var, TypeVarType):
                    # We reimplement TypeVarType.__repr__ here in order to support id_mapper.
                    if var.values:
                        vals = '({})'.format(', '.join(val.accept(self) for val in var.values))
                        vs.append('{} in {}'.format(var.name, vals))
                    elif not is_named_instance(var.upper_bound, 'builtins.object'):
                        vs.append('{} <: {}'.format(var.name, var.upper_bound.accept(self)))
                    else:
                        vs.append(var.name)
                else:
                    # For other TypeVarLikeTypes, just use the name
                    vs.append(var.name)
            s = '{} {}'.format('[{}]'.format(', '.join(vs)), s)

        return 'def {}'.format(s)

    def visit_overloaded(self, t: Overloaded) -> str:
        a = []
        for i in t.items:
            a.append(i.accept(self))
        return 'Overload({})'.format(', '.join(a))

    def visit_tuple_type(self, t: TupleType) -> str:
        s = self.list_str(t.items)
        if t.partial_fallback and t.partial_fallback.type:
            fallback_name = t.partial_fallback.type.fullname
            if fallback_name != 'builtins.tuple':
                return 'Tuple[{}, fallback={}]'.format(s, t.partial_fallback.accept(self))
        return 'Tuple[{}]'.format(s)

    def visit_typeddict_type(self, t: TypedDictType) -> str:
        def item_str(name: str, typ: str) -> str:
            if name in t.required_keys:
                return '{!r}: {}'.format(name, typ)
            else:
                return '{!r}?: {}'.format(name, typ)

        s = '{' + ', '.join(item_str(name, typ.accept(self))
                            for name, typ in t.items.items()) + '}'
        prefix = ''
        if t.fallback and t.fallback.type:
            if t.fallback.type.fullname not in TPDICT_FB_NAMES:
                prefix = repr(t.fallback.type.fullname) + ', '
        return 'TypedDict({}{})'.format(prefix, s)

    def visit_raw_expression_type(self, t: RawExpressionType) -> str:
        return repr(t.literal_value)

    def visit_literal_type(self, t: LiteralType) -> str:
        return 'Literal[{}]'.format(t.value_repr())

    def visit_star_type(self, t: StarType) -> str:
        s = t.type.accept(self)
        return '*{}'.format(s)

    def visit_union_type(self, t: UnionType) -> str:
        s = self.list_str(t.items)
        return 'Union[{}]'.format(s)

    def visit_partial_type(self, t: PartialType) -> str:
        if t.type is None:
            return '<partial None>'
        else:
            return '<partial {}[{}]>'.format(t.type.name,
                                             ', '.join(['?'] * len(t.type.type_vars)))

    def visit_ellipsis_type(self, t: EllipsisType) -> str:
        return '...'

    def visit_type_type(self, t: TypeType) -> str:
        return 'Type[{}]'.format(t.item.accept(self))

    def visit_placeholder_type(self, t: PlaceholderType) -> str:
        return '<placeholder {}>'.format(t.fullname)

    def visit_type_alias_type(self, t: TypeAliasType) -> str:
        if t.alias is not None:
            unrolled, recursed = t._partial_expansion()
            self.any_as_dots = recursed
            type_str = unrolled.accept(self)
            self.any_as_dots = False
            return type_str
        return '<alias (unfixed)>'

    def list_str(self, a: Iterable[Type]) -> str:
        """Convert items of an array to strings (pretty-print types)
        and join the results with commas.
        """
        res = []
        for t in a:
            res.append(t.accept(self))
        return ', '.join(res)


class UnrollAliasVisitor(TypeTranslator):
    def __init__(self, initial_aliases: Set[TypeAliasType]) -> None:
        self.recursed = False
        self.initial_aliases = initial_aliases

    def visit_type_alias_type(self, t: TypeAliasType) -> Type:
        if t in self.initial_aliases:
            self.recursed = True
            return AnyType(TypeOfAny.special_form)
        # Create a new visitor on encountering a new type alias, so that an alias like
        #     A = Tuple[B, B]
        #     B = int
        # will not be detected as recursive on the second encounter of B.
        subvisitor = UnrollAliasVisitor(self.initial_aliases | {t})
        result = get_proper_type(t).accept(subvisitor)
        if subvisitor.recursed:
            self.recursed = True
        return result


def strip_type(typ: Type) -> ProperType:
    """Make a copy of type without 'debugging info' (function name)."""
    typ = get_proper_type(typ)
    if isinstance(typ, CallableType):
        return typ.copy_modified(name=None)
    elif isinstance(typ, Overloaded):
        return Overloaded([cast(CallableType, strip_type(item))
                           for item in typ.items])
    else:
        return typ


def is_named_instance(t: Type, fullname: str) -> bool:
    t = get_proper_type(t)
    return isinstance(t, Instance) and t.type.fullname == fullname


TP = TypeVar('TP', bound=Type)


def copy_type(t: TP) -> TP:
    """
    Build a copy of the type; used to mutate the copy with truthiness information
    """
    return copy.copy(t)


class InstantiateAliasVisitor(TypeTranslator):
    def __init__(self, vars: List[str], subs: List[Type]) -> None:
        self.replacements = {v: s for (v, s) in zip(vars, subs)}

    def visit_type_alias_type(self, typ: TypeAliasType) -> Type:
        return typ.copy_modified(args=[t.accept(self) for t in typ.args])

    def visit_unbound_type(self, typ: UnboundType) -> Type:
        # TODO: stop using unbound type variables for type aliases.
        # Now that type aliases are very similar to TypeInfos we should
        # make type variable tracking similar as well. Maybe we can even support
        # upper bounds etc. for generic type aliases.
        if typ.name in self.replacements:
            return self.replacements[typ.name]
        return typ

    def visit_type_var(self, typ: TypeVarType) -> Type:
        if typ.name in self.replacements:
            return self.replacements[typ.name]
        return typ


def replace_alias_tvars(tp: Type, vars: List[str], subs: List[Type],
                        newline: int, newcolumn: int) -> Type:
    """Replace type variables in a generic type alias tp with substitutions subs
    resetting context. Length of subs should be already checked.
    """
    replacer = InstantiateAliasVisitor(vars, subs)
    new_tp = tp.accept(replacer)
    new_tp.line = newline
    new_tp.column = newcolumn
    return new_tp


class HasTypeVars(TypeQuery[bool]):
    def __init__(self) -> None:
        super().__init__(any)

    def visit_type_var(self, t: TypeVarType) -> bool:
        return True


def has_type_vars(typ: Type) -> bool:
    """Check if a type contains any type variables (recursively)."""
    return typ.accept(HasTypeVars())


def flatten_nested_unions(types: Iterable[Type],
                          handle_type_alias_type: bool = False) -> List[Type]:
    """Flatten nested unions in a type list."""
    # This and similar functions on unions can cause infinite recursion
    # if passed a "pathological" alias like A = Union[int, A] or similar.
    # TODO: ban such aliases in semantic analyzer.
    flat_items: List[Type] = []
    if handle_type_alias_type:
        types = get_proper_types(types)
    # TODO: avoid duplicate types in unions (e.g. using hash)
    for tp in types:
        if isinstance(tp, ProperType) and isinstance(tp, UnionType):
            flat_items.extend(flatten_nested_unions(tp.items,
                              handle_type_alias_type=handle_type_alias_type))
        else:
            flat_items.append(tp)
    return flat_items


def union_items(typ: Type) -> List[ProperType]:
    """Return the flattened items of a union type.

    For non-union types, return a list containing just the argument.
    """
    typ = get_proper_type(typ)
    if isinstance(typ, UnionType):
        items = []
        for item in typ.items:
            items.extend(union_items(item))
        return items
    else:
        return [typ]


def is_union_with_any(tp: Type) -> bool:
    """Is this a union with Any or a plain Any type?"""
    tp = get_proper_type(tp)
    if isinstance(tp, AnyType):
        return True
    if not isinstance(tp, UnionType):
        return False
    return any(is_union_with_any(t) for t in get_proper_types(tp.items))


def is_generic_instance(tp: Type) -> bool:
    tp = get_proper_type(tp)
    return isinstance(tp, Instance) and bool(tp.args)


def is_optional(t: Type) -> bool:
    t = get_proper_type(t)
    return isinstance(t, UnionType) and any(isinstance(get_proper_type(e), NoneType)
                                            for e in t.items)


def remove_optional(typ: Type) -> Type:
    typ = get_proper_type(typ)
    if isinstance(typ, UnionType):
        return UnionType.make_union([t for t in typ.items
                                     if not isinstance(get_proper_type(t), NoneType)])
    else:
        return typ


def is_literal_type(typ: ProperType, fallback_fullname: str, value: LiteralValue) -> bool:
    """Check if this type is a LiteralType with the given fallback type and value."""
    if isinstance(typ, Instance) and typ.last_known_value:
        typ = typ.last_known_value
    if not isinstance(typ, LiteralType):
        return False
    if typ.fallback.type.fullname != fallback_fullname:
        return False
    return typ.value == value


names: Final = globals().copy()
names.pop('NOT_READY', None)
deserialize_map: Final = {
    key: obj.deserialize
    for key, obj in names.items()
    if isinstance(obj, type) and issubclass(obj, Type) and obj is not Type
}


def callable_with_ellipsis(any_type: AnyType,
                           ret_type: Type,
                           fallback: Instance) -> CallableType:
    """Construct type Callable[..., ret_type]."""
    return CallableType([any_type, any_type],
                        [ARG_STAR, ARG_STAR2],
                        [None, None],
                        ret_type=ret_type,
                        fallback=fallback,
                        is_ellipsis_args=True)
