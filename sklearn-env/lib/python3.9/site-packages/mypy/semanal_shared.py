"""Shared definitions used by different parts of semantic analysis."""

from abc import abstractmethod

from typing import Optional, List, Callable
from typing_extensions import Final
from mypy_extensions import trait

from mypy.nodes import (
    Context, SymbolTableNode, MypyFile, ImportedName, FuncDef, Node, TypeInfo, Expression, GDEF,
    SymbolNode, SymbolTable
)
from mypy.util import correct_relative_import
from mypy.types import (
    Type, FunctionLike, Instance, TupleType, TPDICT_FB_NAMES, ProperType, get_proper_type
)
from mypy.tvar_scope import TypeVarLikeScope
from mypy.errorcodes import ErrorCode
from mypy import join

# Priorities for ordering of patches within the "patch" phase of semantic analysis
# (after the main pass):

# Fix fallbacks (does joins)
PRIORITY_FALLBACKS: Final = 1


@trait
class SemanticAnalyzerCoreInterface:
    """A core abstract interface to generic semantic analyzer functionality.

    This is implemented by both semantic analyzer passes 2 and 3.
    """

    @abstractmethod
    def lookup_qualified(self, name: str, ctx: Context,
                         suppress_errors: bool = False) -> Optional[SymbolTableNode]:
        raise NotImplementedError

    @abstractmethod
    def lookup_fully_qualified(self, name: str) -> SymbolTableNode:
        raise NotImplementedError

    @abstractmethod
    def lookup_fully_qualified_or_none(self, name: str) -> Optional[SymbolTableNode]:
        raise NotImplementedError

    @abstractmethod
    def fail(self, msg: str, ctx: Context, serious: bool = False, *,
             blocker: bool = False, code: Optional[ErrorCode] = None) -> None:
        raise NotImplementedError

    @abstractmethod
    def note(self, msg: str, ctx: Context, *, code: Optional[ErrorCode] = None) -> None:
        raise NotImplementedError

    @abstractmethod
    def record_incomplete_ref(self) -> None:
        raise NotImplementedError

    @abstractmethod
    def defer(self) -> None:
        raise NotImplementedError

    @abstractmethod
    def is_incomplete_namespace(self, fullname: str) -> bool:
        """Is a module or class namespace potentially missing some definitions?"""
        raise NotImplementedError

    @property
    @abstractmethod
    def final_iteration(self) -> bool:
        """Is this the final iteration of semantic analysis?"""
        raise NotImplementedError

    @abstractmethod
    def is_future_flag_set(self, flag: str) -> bool:
        """Is the specific __future__ feature imported"""
        raise NotImplementedError

    @property
    @abstractmethod
    def is_stub_file(self) -> bool:
        raise NotImplementedError


@trait
class SemanticAnalyzerInterface(SemanticAnalyzerCoreInterface):
    """A limited abstract interface to some generic semantic analyzer pass 2 functionality.

    We use this interface for various reasons:

    * Looser coupling
    * Cleaner import graph
    * Less need to pass around callback functions
    """

    @abstractmethod
    def lookup(self, name: str, ctx: Context,
               suppress_errors: bool = False) -> Optional[SymbolTableNode]:
        raise NotImplementedError

    @abstractmethod
    def named_type(self, fullname: str,
                   args: Optional[List[Type]] = None) -> Instance:
        raise NotImplementedError

    @abstractmethod
    def named_type_or_none(self, fullname: str,
                           args: Optional[List[Type]] = None) -> Optional[Instance]:
        raise NotImplementedError

    @abstractmethod
    def accept(self, node: Node) -> None:
        raise NotImplementedError

    @abstractmethod
    def anal_type(self, t: Type, *,
                  tvar_scope: Optional[TypeVarLikeScope] = None,
                  allow_tuple_literal: bool = False,
                  allow_unbound_tvars: bool = False,
                  allow_required: bool = False,
                  report_invalid_types: bool = True) -> Optional[Type]:
        raise NotImplementedError

    @abstractmethod
    def basic_new_typeinfo(self, name: str, basetype_or_fallback: Instance, line: int) -> TypeInfo:
        raise NotImplementedError

    @abstractmethod
    def schedule_patch(self, priority: int, fn: Callable[[], None]) -> None:
        raise NotImplementedError

    @abstractmethod
    def add_symbol_table_node(self, name: str, stnode: SymbolTableNode) -> bool:
        """Add node to the current symbol table."""
        raise NotImplementedError

    @abstractmethod
    def current_symbol_table(self) -> SymbolTable:
        """Get currently active symbol table.

        May be module, class, or local namespace.
        """
        raise NotImplementedError

    @abstractmethod
    def add_symbol(self, name: str, node: SymbolNode, context: Context,
                   module_public: bool = True, module_hidden: bool = False,
                   can_defer: bool = True) -> bool:
        """Add symbol to the current symbol table."""
        raise NotImplementedError

    @abstractmethod
    def add_symbol_skip_local(self, name: str, node: SymbolNode) -> None:
        """Add symbol to the current symbol table, skipping locals.

        This is used to store symbol nodes in a symbol table that
        is going to be serialized (local namespaces are not serialized).
        See implementation docstring for more details.
        """
        raise NotImplementedError

    @abstractmethod
    def parse_bool(self, expr: Expression) -> Optional[bool]:
        raise NotImplementedError

    @abstractmethod
    def qualified_name(self, n: str) -> str:
        raise NotImplementedError

    @property
    @abstractmethod
    def is_typeshed_stub_file(self) -> bool:
        raise NotImplementedError

    @abstractmethod
    def is_func_scope(self) -> bool:
        raise NotImplementedError


def create_indirect_imported_name(file_node: MypyFile,
                                  module: str,
                                  relative: int,
                                  imported_name: str) -> Optional[SymbolTableNode]:
    """Create symbol table entry for a name imported from another module.

    These entries act as indirect references.
    """
    target_module, ok = correct_relative_import(
        file_node.fullname,
        relative,
        module,
        file_node.is_package_init_file())
    if not ok:
        return None
    target_name = '%s.%s' % (target_module, imported_name)
    link = ImportedName(target_name)
    # Use GDEF since this refers to a module-level definition.
    return SymbolTableNode(GDEF, link)


def set_callable_name(sig: Type, fdef: FuncDef) -> ProperType:
    sig = get_proper_type(sig)
    if isinstance(sig, FunctionLike):
        if fdef.info:
            if fdef.info.fullname in TPDICT_FB_NAMES:
                # Avoid exposing the internal _TypedDict name.
                class_name = 'TypedDict'
            else:
                class_name = fdef.info.name
            return sig.with_name(
                '{} of {}'.format(fdef.name, class_name))
        else:
            return sig.with_name(fdef.name)
    else:
        return sig


def calculate_tuple_fallback(typ: TupleType) -> None:
    """Calculate a precise item type for the fallback of a tuple type.

    This must be called only after the main semantic analysis pass, since joins
    aren't available before that.

    Note that there is an apparent chicken and egg problem with respect
    to verifying type arguments against bounds. Verifying bounds might
    require fallbacks, but we might use the bounds to calculate the
    fallbacks. In practice this is not a problem, since the worst that
    can happen is that we have invalid type argument values, and these
    can happen in later stages as well (they will generate errors, but
    we don't prevent their existence).
    """
    fallback = typ.partial_fallback
    assert fallback.type.fullname == 'builtins.tuple'
    fallback.args = (join.join_type_list(list(typ.items)),) + fallback.args[1:]
