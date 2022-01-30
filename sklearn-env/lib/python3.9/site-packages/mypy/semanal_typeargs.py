"""Verify properties of type arguments, like 'int' in C[int] being valid.

This must happen after semantic analysis since there can be placeholder
types until the end of semantic analysis, and these break various type
operations, including subtype checks.
"""

from typing import List, Optional, Set

from mypy.nodes import TypeInfo, Context, MypyFile, FuncItem, ClassDef, Block, FakeInfo
from mypy.types import (
    Type, Instance, TypeVarType, AnyType, get_proper_types, TypeAliasType, ParamSpecType,
    get_proper_type
)
from mypy.mixedtraverser import MixedTraverserVisitor
from mypy.subtypes import is_subtype
from mypy.sametypes import is_same_type
from mypy.errors import Errors
from mypy.scope import Scope
from mypy.options import Options
from mypy.errorcodes import ErrorCode
from mypy import message_registry, errorcodes as codes
from mypy.messages import format_type


class TypeArgumentAnalyzer(MixedTraverserVisitor):
    def __init__(self, errors: Errors, options: Options, is_typeshed_file: bool) -> None:
        self.errors = errors
        self.options = options
        self.is_typeshed_file = is_typeshed_file
        self.scope = Scope()
        # Should we also analyze function definitions, or only module top-levels?
        self.recurse_into_functions = True
        # Keep track of the type aliases already visited. This is needed to avoid
        # infinite recursion on types like A = Union[int, List[A]].
        self.seen_aliases: Set[TypeAliasType] = set()

    def visit_mypy_file(self, o: MypyFile) -> None:
        self.errors.set_file(o.path, o.fullname, scope=self.scope)
        with self.scope.module_scope(o.fullname):
            super().visit_mypy_file(o)

    def visit_func(self, defn: FuncItem) -> None:
        if not self.recurse_into_functions:
            return
        with self.scope.function_scope(defn):
            super().visit_func(defn)

    def visit_class_def(self, defn: ClassDef) -> None:
        with self.scope.class_scope(defn.info):
            super().visit_class_def(defn)

    def visit_block(self, o: Block) -> None:
        if not o.is_unreachable:
            super().visit_block(o)

    def visit_type_alias_type(self, t: TypeAliasType) -> None:
        super().visit_type_alias_type(t)
        if t in self.seen_aliases:
            # Avoid infinite recursion on recursive type aliases.
            # Note: it is fine to skip the aliases we have already seen in non-recursive types,
            # since errors there have already already reported.
            return
        self.seen_aliases.add(t)
        get_proper_type(t).accept(self)

    def visit_instance(self, t: Instance) -> None:
        # Type argument counts were checked in the main semantic analyzer pass. We assume
        # that the counts are correct here.
        info = t.type
        if isinstance(info, FakeInfo):
            return  # https://github.com/python/mypy/issues/11079
        for (i, arg), tvar in zip(enumerate(t.args), info.defn.type_vars):
            if isinstance(tvar, TypeVarType):
                if isinstance(arg, ParamSpecType):
                    # TODO: Better message
                    self.fail(f'Invalid location for ParamSpec "{arg.name}"', t)
                    continue
                if tvar.values:
                    if isinstance(arg, TypeVarType):
                        arg_values = arg.values
                        if not arg_values:
                            self.fail(
                                message_registry.INVALID_TYPEVAR_AS_TYPEARG.format(
                                    arg.name, info.name),
                                t, code=codes.TYPE_VAR)
                            continue
                    else:
                        arg_values = [arg]
                    self.check_type_var_values(info, arg_values, tvar.name, tvar.values, i + 1, t)
                if not is_subtype(arg, tvar.upper_bound):
                    self.fail(
                        message_registry.INVALID_TYPEVAR_ARG_BOUND.format(
                            format_type(arg), info.name, format_type(tvar.upper_bound)),
                        t, code=codes.TYPE_VAR)
        super().visit_instance(t)

    def check_type_var_values(self, type: TypeInfo, actuals: List[Type], arg_name: str,
                              valids: List[Type], arg_number: int, context: Context) -> None:
        for actual in get_proper_types(actuals):
            if (not isinstance(actual, AnyType) and
                    not any(is_same_type(actual, value)
                            for value in valids)):
                if len(actuals) > 1 or not isinstance(actual, Instance):
                    self.fail(
                        message_registry.INVALID_TYPEVAR_ARG_VALUE.format(type.name),
                        context, code=codes.TYPE_VAR)
                else:
                    class_name = '"{}"'.format(type.name)
                    actual_type_name = '"{}"'.format(actual.type.name)
                    self.fail(
                        message_registry.INCOMPATIBLE_TYPEVAR_VALUE.format(
                            arg_name, class_name, actual_type_name),
                        context,
                        code=codes.TYPE_VAR)

    def fail(self, msg: str, context: Context, *, code: Optional[ErrorCode] = None) -> None:
        self.errors.report(context.get_line(), context.get_column(), msg, code=code)
