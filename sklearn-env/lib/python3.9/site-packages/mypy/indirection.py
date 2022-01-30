from typing import Dict, Iterable, List, Optional, Set, Union

from mypy.types import TypeVisitor
import mypy.types as types
from mypy.util import split_module_names


def extract_module_names(type_name: Optional[str]) -> List[str]:
    """Returns the module names of a fully qualified type name."""
    if type_name is not None:
        # Discard the first one, which is just the qualified name of the type
        possible_module_names = split_module_names(type_name)
        return possible_module_names[1:]
    else:
        return []


class TypeIndirectionVisitor(TypeVisitor[Set[str]]):
    """Returns all module references within a particular type."""

    def __init__(self) -> None:
        self.cache: Dict[types.Type, Set[str]] = {}
        self.seen_aliases: Set[types.TypeAliasType] = set()

    def find_modules(self, typs: Iterable[types.Type]) -> Set[str]:
        self.seen_aliases.clear()
        return self._visit(typs)

    def _visit(self, typ_or_typs: Union[types.Type, Iterable[types.Type]]) -> Set[str]:
        typs = [typ_or_typs] if isinstance(typ_or_typs, types.Type) else typ_or_typs
        output: Set[str] = set()
        for typ in typs:
            if isinstance(typ, types.TypeAliasType):
                # Avoid infinite recursion for recursive type aliases.
                if typ in self.seen_aliases:
                    continue
                self.seen_aliases.add(typ)
            if typ in self.cache:
                modules = self.cache[typ]
            else:
                modules = typ.accept(self)
                self.cache[typ] = set(modules)
            output.update(modules)
        return output

    def visit_unbound_type(self, t: types.UnboundType) -> Set[str]:
        return self._visit(t.args)

    def visit_any(self, t: types.AnyType) -> Set[str]:
        return set()

    def visit_none_type(self, t: types.NoneType) -> Set[str]:
        return set()

    def visit_uninhabited_type(self, t: types.UninhabitedType) -> Set[str]:
        return set()

    def visit_erased_type(self, t: types.ErasedType) -> Set[str]:
        return set()

    def visit_deleted_type(self, t: types.DeletedType) -> Set[str]:
        return set()

    def visit_type_var(self, t: types.TypeVarType) -> Set[str]:
        return self._visit(t.values) | self._visit(t.upper_bound)

    def visit_param_spec(self, t: types.ParamSpecType) -> Set[str]:
        return set()

    def visit_instance(self, t: types.Instance) -> Set[str]:
        out = self._visit(t.args)
        if t.type:
            # Uses of a class depend on everything in the MRO,
            # as changes to classes in the MRO can add types to methods,
            # change property types, change the MRO itself, etc.
            for s in t.type.mro:
                out.update(split_module_names(s.module_name))
            if t.type.metaclass_type is not None:
                out.update(split_module_names(t.type.metaclass_type.type.module_name))
        return out

    def visit_callable_type(self, t: types.CallableType) -> Set[str]:
        out = self._visit(t.arg_types) | self._visit(t.ret_type)
        if t.definition is not None:
            out.update(extract_module_names(t.definition.fullname))
        return out

    def visit_overloaded(self, t: types.Overloaded) -> Set[str]:
        return self._visit(t.items) | self._visit(t.fallback)

    def visit_tuple_type(self, t: types.TupleType) -> Set[str]:
        return self._visit(t.items) | self._visit(t.partial_fallback)

    def visit_typeddict_type(self, t: types.TypedDictType) -> Set[str]:
        return self._visit(t.items.values()) | self._visit(t.fallback)

    def visit_literal_type(self, t: types.LiteralType) -> Set[str]:
        return self._visit(t.fallback)

    def visit_union_type(self, t: types.UnionType) -> Set[str]:
        return self._visit(t.items)

    def visit_partial_type(self, t: types.PartialType) -> Set[str]:
        return set()

    def visit_type_type(self, t: types.TypeType) -> Set[str]:
        return self._visit(t.item)

    def visit_type_alias_type(self, t: types.TypeAliasType) -> Set[str]:
        return self._visit(types.get_proper_type(t))
