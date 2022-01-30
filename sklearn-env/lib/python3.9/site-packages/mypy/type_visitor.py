"""Type visitor classes.

This module defines the type visitors that are intended to be
subclassed by other code.  They have been separated out into their own
module to ease converting mypy to run under mypyc, since currently
mypyc-extension classes can extend interpreted classes but not the
other way around. Separating them out, then, allows us to compile
types before we can compile everything that uses a TypeVisitor.

The visitors are all re-exported from mypy.types and that is how
other modules refer to them.
"""

from abc import abstractmethod
from mypy.backports import OrderedDict
from typing import Generic, TypeVar, cast, Any, List, Callable, Iterable, Optional, Set, Sequence
from mypy_extensions import trait, mypyc_attr

T = TypeVar('T')

from mypy.types import (
    Type, AnyType, CallableType, Overloaded, TupleType, TypedDictType, LiteralType,
    RawExpressionType, Instance, NoneType, TypeType,
    UnionType, TypeVarType, PartialType, DeletedType, UninhabitedType, TypeVarLikeType,
    UnboundType, ErasedType, StarType, EllipsisType, TypeList, CallableArgument,
    PlaceholderType, TypeAliasType, ParamSpecType, get_proper_type
)


@trait
@mypyc_attr(allow_interpreted_subclasses=True)
class TypeVisitor(Generic[T]):
    """Visitor class for types (Type subclasses).

    The parameter T is the return type of the visit methods.
    """

    @abstractmethod
    def visit_unbound_type(self, t: UnboundType) -> T:
        pass

    @abstractmethod
    def visit_any(self, t: AnyType) -> T:
        pass

    @abstractmethod
    def visit_none_type(self, t: NoneType) -> T:
        pass

    @abstractmethod
    def visit_uninhabited_type(self, t: UninhabitedType) -> T:
        pass

    @abstractmethod
    def visit_erased_type(self, t: ErasedType) -> T:
        pass

    @abstractmethod
    def visit_deleted_type(self, t: DeletedType) -> T:
        pass

    @abstractmethod
    def visit_type_var(self, t: TypeVarType) -> T:
        pass

    @abstractmethod
    def visit_param_spec(self, t: ParamSpecType) -> T:
        pass

    @abstractmethod
    def visit_instance(self, t: Instance) -> T:
        pass

    @abstractmethod
    def visit_callable_type(self, t: CallableType) -> T:
        pass

    @abstractmethod
    def visit_overloaded(self, t: Overloaded) -> T:
        pass

    @abstractmethod
    def visit_tuple_type(self, t: TupleType) -> T:
        pass

    @abstractmethod
    def visit_typeddict_type(self, t: TypedDictType) -> T:
        pass

    @abstractmethod
    def visit_literal_type(self, t: LiteralType) -> T:
        pass

    @abstractmethod
    def visit_union_type(self, t: UnionType) -> T:
        pass

    @abstractmethod
    def visit_partial_type(self, t: PartialType) -> T:
        pass

    @abstractmethod
    def visit_type_type(self, t: TypeType) -> T:
        pass

    @abstractmethod
    def visit_type_alias_type(self, t: TypeAliasType) -> T:
        pass


@trait
@mypyc_attr(allow_interpreted_subclasses=True)
class SyntheticTypeVisitor(TypeVisitor[T]):
    """A TypeVisitor that also knows how to visit synthetic AST constructs.

       Not just real types."""

    @abstractmethod
    def visit_star_type(self, t: StarType) -> T:
        pass

    @abstractmethod
    def visit_type_list(self, t: TypeList) -> T:
        pass

    @abstractmethod
    def visit_callable_argument(self, t: CallableArgument) -> T:
        pass

    @abstractmethod
    def visit_ellipsis_type(self, t: EllipsisType) -> T:
        pass

    @abstractmethod
    def visit_raw_expression_type(self, t: RawExpressionType) -> T:
        pass

    @abstractmethod
    def visit_placeholder_type(self, t: PlaceholderType) -> T:
        pass


@mypyc_attr(allow_interpreted_subclasses=True)
class TypeTranslator(TypeVisitor[Type]):
    """Identity type transformation.

    Subclass this and override some methods to implement a non-trivial
    transformation.
    """

    def visit_unbound_type(self, t: UnboundType) -> Type:
        return t

    def visit_any(self, t: AnyType) -> Type:
        return t

    def visit_none_type(self, t: NoneType) -> Type:
        return t

    def visit_uninhabited_type(self, t: UninhabitedType) -> Type:
        return t

    def visit_erased_type(self, t: ErasedType) -> Type:
        return t

    def visit_deleted_type(self, t: DeletedType) -> Type:
        return t

    def visit_instance(self, t: Instance) -> Type:
        last_known_value: Optional[LiteralType] = None
        if t.last_known_value is not None:
            raw_last_known_value = t.last_known_value.accept(self)
            assert isinstance(raw_last_known_value, LiteralType)  # type: ignore
            last_known_value = raw_last_known_value
        return Instance(
            typ=t.type,
            args=self.translate_types(t.args),
            line=t.line,
            column=t.column,
            last_known_value=last_known_value,
        )

    def visit_type_var(self, t: TypeVarType) -> Type:
        return t

    def visit_param_spec(self, t: ParamSpecType) -> Type:
        return t

    def visit_partial_type(self, t: PartialType) -> Type:
        return t

    def visit_callable_type(self, t: CallableType) -> Type:
        return t.copy_modified(arg_types=self.translate_types(t.arg_types),
                               ret_type=t.ret_type.accept(self),
                               variables=self.translate_variables(t.variables))

    def visit_tuple_type(self, t: TupleType) -> Type:
        return TupleType(self.translate_types(t.items),
                         # TODO: This appears to be unsafe.
                         cast(Any, t.partial_fallback.accept(self)),
                         t.line, t.column)

    def visit_typeddict_type(self, t: TypedDictType) -> Type:
        items = OrderedDict([
            (item_name, item_type.accept(self))
            for (item_name, item_type) in t.items.items()
        ])
        return TypedDictType(items,
                             t.required_keys,
                             # TODO: This appears to be unsafe.
                             cast(Any, t.fallback.accept(self)),
                             t.line, t.column)

    def visit_literal_type(self, t: LiteralType) -> Type:
        fallback = t.fallback.accept(self)
        assert isinstance(fallback, Instance)  # type: ignore
        return LiteralType(
            value=t.value,
            fallback=fallback,
            line=t.line,
            column=t.column,
        )

    def visit_union_type(self, t: UnionType) -> Type:
        return UnionType(self.translate_types(t.items), t.line, t.column)

    def translate_types(self, types: Iterable[Type]) -> List[Type]:
        return [t.accept(self) for t in types]

    def translate_variables(self,
                            variables: Sequence[TypeVarLikeType]) -> Sequence[TypeVarLikeType]:
        return variables

    def visit_overloaded(self, t: Overloaded) -> Type:
        items: List[CallableType] = []
        for item in t.items:
            new = item.accept(self)
            assert isinstance(new, CallableType)  # type: ignore
            items.append(new)
        return Overloaded(items=items)

    def visit_type_type(self, t: TypeType) -> Type:
        return TypeType.make_normalized(t.item.accept(self), line=t.line, column=t.column)

    @abstractmethod
    def visit_type_alias_type(self, t: TypeAliasType) -> Type:
        # This method doesn't have a default implementation for type translators,
        # because type aliases are special: some information is contained in the
        # TypeAlias node, and we normally don't generate new nodes. Every subclass
        # must implement this depending on its semantics.
        pass


@mypyc_attr(allow_interpreted_subclasses=True)
class TypeQuery(SyntheticTypeVisitor[T]):
    """Visitor for performing queries of types.

    strategy is used to combine results for a series of types,
    common use cases involve a boolean query using `any` or `all`.

    Note: this visitor keeps an internal state (tracks type aliases to avoid
    recursion), so it should *never* be re-used for querying different types,
    create a new visitor instance instead.

    # TODO: check that we don't have existing violations of this rule.
    """

    def __init__(self, strategy: Callable[[Iterable[T]], T]) -> None:
        self.strategy = strategy
        # Keep track of the type aliases already visited. This is needed to avoid
        # infinite recursion on types like A = Union[int, List[A]].
        self.seen_aliases: Set[TypeAliasType] = set()

    def visit_unbound_type(self, t: UnboundType) -> T:
        return self.query_types(t.args)

    def visit_type_list(self, t: TypeList) -> T:
        return self.query_types(t.items)

    def visit_callable_argument(self, t: CallableArgument) -> T:
        return t.typ.accept(self)

    def visit_any(self, t: AnyType) -> T:
        return self.strategy([])

    def visit_uninhabited_type(self, t: UninhabitedType) -> T:
        return self.strategy([])

    def visit_none_type(self, t: NoneType) -> T:
        return self.strategy([])

    def visit_erased_type(self, t: ErasedType) -> T:
        return self.strategy([])

    def visit_deleted_type(self, t: DeletedType) -> T:
        return self.strategy([])

    def visit_type_var(self, t: TypeVarType) -> T:
        return self.query_types([t.upper_bound] + t.values)

    def visit_param_spec(self, t: ParamSpecType) -> T:
        return self.strategy([])

    def visit_partial_type(self, t: PartialType) -> T:
        return self.strategy([])

    def visit_instance(self, t: Instance) -> T:
        return self.query_types(t.args)

    def visit_callable_type(self, t: CallableType) -> T:
        # FIX generics
        return self.query_types(t.arg_types + [t.ret_type])

    def visit_tuple_type(self, t: TupleType) -> T:
        return self.query_types(t.items)

    def visit_typeddict_type(self, t: TypedDictType) -> T:
        return self.query_types(t.items.values())

    def visit_raw_expression_type(self, t: RawExpressionType) -> T:
        return self.strategy([])

    def visit_literal_type(self, t: LiteralType) -> T:
        return self.strategy([])

    def visit_star_type(self, t: StarType) -> T:
        return t.type.accept(self)

    def visit_union_type(self, t: UnionType) -> T:
        return self.query_types(t.items)

    def visit_overloaded(self, t: Overloaded) -> T:
        return self.query_types(t.items)

    def visit_type_type(self, t: TypeType) -> T:
        return t.item.accept(self)

    def visit_ellipsis_type(self, t: EllipsisType) -> T:
        return self.strategy([])

    def visit_placeholder_type(self, t: PlaceholderType) -> T:
        return self.query_types(t.args)

    def visit_type_alias_type(self, t: TypeAliasType) -> T:
        return get_proper_type(t).accept(self)

    def query_types(self, types: Iterable[Type]) -> T:
        """Perform a query for a list of types.

        Use the strategy to combine the results.
        Skip type aliases already visited types to avoid infinite recursion.
        """
        res: List[T] = []
        for t in types:
            if isinstance(t, TypeAliasType):
                # Avoid infinite recursion for recursive type aliases.
                # TODO: Ideally we should fire subvisitors here (or use caching) if we care
                #       about duplicates.
                if t in self.seen_aliases:
                    continue
                self.seen_aliases.add(t)
            res.append(t.accept(self))
        return self.strategy(res)
