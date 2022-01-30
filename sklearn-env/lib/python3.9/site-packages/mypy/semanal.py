"""The semantic analyzer.

Bind names to definitions and do various other simple consistency
checks.  Populate symbol tables.  The semantic analyzer also detects
special forms which reuse generic syntax such as NamedTuple and
cast().  Multiple analysis iterations may be needed to analyze forward
references and import cycles. Each iteration "fills in" additional
bindings and references until everything has been bound.

For example, consider this program:

  x = 1
  y = x

Here semantic analysis would detect that the assignment 'x = 1'
defines a new variable, the type of which is to be inferred (in a
later pass; type inference or type checking is not part of semantic
analysis).  Also, it would bind both references to 'x' to the same
module-level variable (Var) node.  The second assignment would also
be analyzed, and the type of 'y' marked as being inferred.

Semantic analysis of types is implemented in typeanal.py.

See semanal_main.py for the top-level logic.

Some important properties:

* After semantic analysis is complete, no PlaceholderNode and
  PlaceholderType instances should remain. During semantic analysis,
  if we encounter one of these, the current target should be deferred.

* A TypeInfo is only created once we know certain basic information about
  a type, such as the MRO, existence of a Tuple base class (e.g., for named
  tuples), and whether we have a TypedDict. We use a temporary
  PlaceholderNode node in the symbol table if some such information is
  missing.

* For assignments, we only add a non-placeholder symbol table entry once
  we know the sort of thing being defined (variable, NamedTuple, type alias,
  etc.).

* Every part of the analysis step must support multiple iterations over
  the same AST nodes, and each iteration must be able to fill in arbitrary
  things that were missing or incomplete in previous iterations.

* Changes performed by the analysis need to be reversible, since mypy
  daemon strips and reuses existing ASTs (to improve performance and/or
  reduce memory use).
"""

from contextlib import contextmanager

from typing import (
    List, Dict, Set, Tuple, cast, TypeVar, Union, Optional, Callable, Iterator, Iterable
)
from typing_extensions import Final, TypeAlias as _TypeAlias

from mypy.nodes import (
    MypyFile, TypeInfo, Node, AssignmentStmt, FuncDef, OverloadedFuncDef,
    ClassDef, Var, GDEF, FuncItem, Import, Expression, Lvalue,
    ImportFrom, ImportAll, Block, LDEF, NameExpr, MemberExpr,
    IndexExpr, TupleExpr, ListExpr, ExpressionStmt, ReturnStmt,
    RaiseStmt, AssertStmt, OperatorAssignmentStmt, WhileStmt,
    ForStmt, BreakStmt, ContinueStmt, IfStmt, TryStmt, WithStmt, DelStmt,
    GlobalDecl, SuperExpr, DictExpr, CallExpr, RefExpr, OpExpr, UnaryExpr,
    SliceExpr, CastExpr, RevealExpr, TypeApplication, Context, SymbolTable,
    SymbolTableNode, ListComprehension, GeneratorExpr,
    LambdaExpr, MDEF, Decorator, SetExpr, TypeVarExpr,
    StrExpr, BytesExpr, PrintStmt, ConditionalExpr, PromoteExpr,
    ComparisonExpr, StarExpr, ArgKind, ARG_POS, ARG_NAMED, type_aliases,
    YieldFromExpr, NamedTupleExpr, NonlocalDecl, SymbolNode,
    SetComprehension, DictionaryComprehension, TypeAlias, TypeAliasExpr,
    YieldExpr, ExecStmt, BackquoteExpr, ImportBase, AwaitExpr,
    IntExpr, FloatExpr, UnicodeExpr, TempNode, OverloadPart,
    PlaceholderNode, COVARIANT, CONTRAVARIANT, INVARIANT,
    get_nongen_builtins, get_member_expr_fullname, REVEAL_TYPE,
    REVEAL_LOCALS, is_final_node, TypedDictExpr, type_aliases_source_versions,
    typing_extensions_aliases,
    EnumCallExpr, RUNTIME_PROTOCOL_DECOS, FakeExpression, Statement, AssignmentExpr,
    ParamSpecExpr, EllipsisExpr, TypeVarLikeExpr, FuncBase, implicit_module_attrs,
)
from mypy.tvar_scope import TypeVarLikeScope
from mypy.typevars import fill_typevars
from mypy.visitor import NodeVisitor
from mypy.errors import Errors, report_internal_error
from mypy.messages import (
    best_matches, MessageBuilder, pretty_seq, SUGGESTED_TEST_FIXTURES, TYPES_FOR_UNIMPORTED_HINTS
)
from mypy.errorcodes import ErrorCode
from mypy import message_registry, errorcodes as codes
from mypy.types import (
    FunctionLike, UnboundType, TypeVarType, TupleType, UnionType, StarType,
    CallableType, Overloaded, Instance, Type, AnyType, LiteralType, LiteralValue,
    TypeTranslator, TypeOfAny, TypeType, NoneType, PlaceholderType, TPDICT_NAMES, ProperType,
    get_proper_type, get_proper_types, TypeAliasType, TypeVarLikeType
)
from mypy.typeops import function_type, get_type_vars
from mypy.type_visitor import TypeQuery
from mypy.typeanal import (
    TypeAnalyser, analyze_type_alias, no_subscript_builtin_alias,
    TypeVarLikeQuery, TypeVarLikeList, remove_dups, has_any_from_unimported_type,
    check_for_explicit_any, type_constructors, fix_instance_types
)
from mypy.exprtotype import expr_to_unanalyzed_type, TypeTranslationError
from mypy.options import Options
from mypy.plugin import (
    Plugin, ClassDefContext, SemanticAnalyzerPluginInterface,
    DynamicClassDefContext
)
from mypy.util import (
    correct_relative_import, unmangle, module_prefix, is_typeshed_file, unnamed_function,
)
from mypy.scope import Scope
from mypy.semanal_shared import (
    SemanticAnalyzerInterface, set_callable_name, calculate_tuple_fallback, PRIORITY_FALLBACKS
)
from mypy.semanal_namedtuple import NamedTupleAnalyzer
from mypy.semanal_typeddict import TypedDictAnalyzer
from mypy.semanal_enum import EnumCallAnalyzer, ENUM_BASES
from mypy.semanal_newtype import NewTypeAnalyzer
from mypy.reachability import (
    infer_reachability_of_if_statement, infer_condition_value, ALWAYS_FALSE, ALWAYS_TRUE,
    MYPY_TRUE, MYPY_FALSE
)
from mypy.mro import calculate_mro, MroError

T = TypeVar('T')


FUTURE_IMPORTS: Final = {
    '__future__.nested_scopes': 'nested_scopes',
    '__future__.generators': 'generators',
    '__future__.division': 'division',
    '__future__.absolute_import': 'absolute_import',
    '__future__.with_statement': 'with_statement',
    '__future__.print_function': 'print_function',
    '__future__.unicode_literals': 'unicode_literals',
    '__future__.barry_as_FLUFL': 'barry_as_FLUFL',
    '__future__.generator_stop': 'generator_stop',
    '__future__.annotations': 'annotations',
}


# Special cased built-in classes that are needed for basic functionality and need to be
# available very early on.
CORE_BUILTIN_CLASSES: Final = ["object", "bool", "function"]


# Used for tracking incomplete references
Tag: _TypeAlias = int


class SemanticAnalyzer(NodeVisitor[None],
                       SemanticAnalyzerInterface,
                       SemanticAnalyzerPluginInterface):
    """Semantically analyze parsed mypy files.

    The analyzer binds names and does various consistency checks for an
    AST. Note that type checking is performed as a separate pass.
    """

    __deletable__ = ['patches', 'options', 'cur_mod_node']

    # Module name space
    modules: Dict[str, MypyFile]
    # Global name space for current module
    globals: SymbolTable
    # Names declared using "global" (separate set for each scope)
    global_decls: List[Set[str]]
    # Names declared using "nonlocal" (separate set for each scope)
    nonlocal_decls: List[Set[str]]
    # Local names of function scopes; None for non-function scopes.
    locals: List[Optional[SymbolTable]]
    # Whether each scope is a comprehension scope.
    is_comprehension_stack: List[bool]
    # Nested block depths of scopes
    block_depth: List[int]
    # TypeInfo of directly enclosing class (or None)
    type: Optional[TypeInfo] = None
    # Stack of outer classes (the second tuple item contains tvars).
    type_stack: List[Optional[TypeInfo]]
    # Type variables bound by the current scope, be it class or function
    tvar_scope: TypeVarLikeScope
    # Per-module options
    options: Options

    # Stack of functions being analyzed
    function_stack: List[FuncItem]

    # Set to True if semantic analysis defines a name, or replaces a
    # placeholder definition. If some iteration makes no progress,
    # there can be at most one additional final iteration (see below).
    progress = False
    deferred = False  # Set to true if another analysis pass is needed
    incomplete = False  # Set to true if current module namespace is missing things
    # Is this the final iteration of semantic analysis (where we report
    # unbound names due to cyclic definitions and should not defer)?
    _final_iteration = False
    # These names couldn't be added to the symbol table due to incomplete deps.
    # Note that missing names are per module, _not_ per namespace. This means that e.g.
    # a missing name at global scope will block adding same name at a class scope.
    # This should not affect correctness and is purely a performance issue,
    # since it can cause unnecessary deferrals. These are represented as
    # PlaceholderNodes in the symbol table. We use this to ensure that the first
    # definition takes precedence even if it's incomplete.
    #
    # Note that a star import adds a special name '*' to the set, this blocks
    # adding _any_ names in the current file.
    missing_names: List[Set[str]]
    # Callbacks that will be called after semantic analysis to tweak things.
    patches: List[Tuple[int, Callable[[], None]]]
    loop_depth = 0         # Depth of breakable loops
    cur_mod_id = ''        # Current module id (or None) (phase 2)
    _is_stub_file = False   # Are we analyzing a stub file?
    _is_typeshed_stub_file = False  # Are we analyzing a typeshed stub file?
    imports: Set[str]  # Imported modules (during phase 2 analysis)
    # Note: some imports (and therefore dependencies) might
    # not be found in phase 1, for example due to * imports.
    errors: Errors  # Keeps track of generated errors
    plugin: Plugin  # Mypy plugin for special casing of library features
    statement: Optional[Statement] = None  # Statement/definition being analyzed
    future_import_flags: Set[str]

    # Mapping from 'async def' function definitions to their return type wrapped as a
    # 'Coroutine[Any, Any, T]'. Used to keep track of whether a function definition's
    # return type has already been wrapped, by checking if the function definition's
    # type is stored in this mapping and that it still matches.
    wrapped_coro_return_types: Dict[FuncDef, Type] = {}

    def __init__(self,
                 modules: Dict[str, MypyFile],
                 missing_modules: Set[str],
                 incomplete_namespaces: Set[str],
                 errors: Errors,
                 plugin: Plugin) -> None:
        """Construct semantic analyzer.

        We reuse the same semantic analyzer instance across multiple modules.

        Args:
            modules: Global modules dictionary
            missing_modules: Modules that could not be imported encountered so far
            incomplete_namespaces: Namespaces that are being populated during semantic analysis
                (can contain modules and classes within the current SCC; mutated by the caller)
            errors: Report analysis errors using this instance
        """
        self.locals = [None]
        self.is_comprehension_stack = [False]
        # Saved namespaces from previous iteration. Every top-level function/method body is
        # analyzed in several iterations until all names are resolved. We need to save
        # the local namespaces for the top level function and all nested functions between
        # these iterations. See also semanal_main.process_top_level_function().
        self.saved_locals: Dict[
            Union[FuncItem, GeneratorExpr, DictionaryComprehension], SymbolTable
        ] = {}
        self.imports = set()
        self.type = None
        self.type_stack = []
        # Are the namespaces of classes being processed complete?
        self.incomplete_type_stack: List[bool] = []
        self.tvar_scope = TypeVarLikeScope()
        self.function_stack = []
        self.block_depth = [0]
        self.loop_depth = 0
        self.errors = errors
        self.modules = modules
        self.msg = MessageBuilder(errors, modules)
        self.missing_modules = missing_modules
        self.missing_names = [set()]
        # These namespaces are still in process of being populated. If we encounter a
        # missing name in these namespaces, we need to defer the current analysis target,
        # since it's possible that the name will be there once the namespace is complete.
        self.incomplete_namespaces = incomplete_namespaces
        self.all_exports: List[str] = []
        # Map from module id to list of explicitly exported names (i.e. names in __all__).
        self.export_map: Dict[str, List[str]] = {}
        self.plugin = plugin
        # If True, process function definitions. If False, don't. This is used
        # for processing module top levels in fine-grained incremental mode.
        self.recurse_into_functions = True
        self.scope = Scope()

        # Trace line numbers for every file where deferral happened during analysis of
        # current SCC or top-level function.
        self.deferral_debug_context: List[Tuple[str, int]] = []

        self.future_import_flags: Set[str] = set()

    # mypyc doesn't properly handle implementing an abstractproperty
    # with a regular attribute so we make them properties
    @property
    def is_stub_file(self) -> bool:
        return self._is_stub_file

    @property
    def is_typeshed_stub_file(self) -> bool:
        return self._is_typeshed_stub_file

    @property
    def final_iteration(self) -> bool:
        return self._final_iteration

    #
    # Preparing module (performed before semantic analysis)
    #

    def prepare_file(self, file_node: MypyFile) -> None:
        """Prepare a freshly parsed file for semantic analysis."""
        if 'builtins' in self.modules:
            file_node.names['__builtins__'] = SymbolTableNode(GDEF,
                                                              self.modules['builtins'])
        if file_node.fullname == 'builtins':
            self.prepare_builtins_namespace(file_node)
        if file_node.fullname == 'typing':
            self.prepare_typing_namespace(file_node, type_aliases)
        if file_node.fullname == 'typing_extensions':
            self.prepare_typing_namespace(file_node, typing_extensions_aliases)

    def prepare_typing_namespace(self, file_node: MypyFile,
                                 aliases: Dict[str, str]) -> None:
        """Remove dummy alias definitions such as List = TypeAlias(object) from typing.

        They will be replaced with real aliases when corresponding targets are ready.
        """
        # This is all pretty unfortunate. typeshed now has a
        # sys.version_info check for OrderedDict, and we shouldn't
        # take it out, because it is correct and a typechecker should
        # use that as a source of truth. But instead we rummage
        # through IfStmts to remove the info first.  (I tried to
        # remove this whole machinery and ran into issues with the
        # builtins/typing import cycle.)
        def helper(defs: List[Statement]) -> None:
            for stmt in defs.copy():
                if isinstance(stmt, IfStmt):
                    for body in stmt.body:
                        helper(body.body)
                    if stmt.else_body:
                        helper(stmt.else_body.body)
                if (isinstance(stmt, AssignmentStmt) and len(stmt.lvalues) == 1 and
                        isinstance(stmt.lvalues[0], NameExpr)):
                    # Assignment to a simple name, remove it if it is a dummy alias.
                    if f'{file_node.fullname}.{stmt.lvalues[0].name}' in aliases:
                        defs.remove(stmt)

        helper(file_node.defs)

    def prepare_builtins_namespace(self, file_node: MypyFile) -> None:
        """Add certain special-cased definitions to the builtins module.

        Some definitions are too special or fundamental to be processed
        normally from the AST.
        """
        names = file_node.names

        # Add empty definition for core built-in classes, since they are required for basic
        # operation. These will be completed later on.
        for name in CORE_BUILTIN_CLASSES:
            cdef = ClassDef(name, Block([]))  # Dummy ClassDef, will be replaced later
            info = TypeInfo(SymbolTable(), cdef, 'builtins')
            info._fullname = 'builtins.%s' % name
            names[name] = SymbolTableNode(GDEF, info)

        bool_info = names['bool'].node
        assert isinstance(bool_info, TypeInfo)
        bool_type = Instance(bool_info, [])

        special_var_types: List[Tuple[str, Type]] = [
            ('None', NoneType()),
            # reveal_type is a mypy-only function that gives an error with
            # the type of its arg.
            ('reveal_type', AnyType(TypeOfAny.special_form)),
            # reveal_locals is a mypy-only function that gives an error with the types of
            # locals
            ('reveal_locals', AnyType(TypeOfAny.special_form)),
            ('True', bool_type),
            ('False', bool_type),
            ('__debug__', bool_type),
        ]

        for name, typ in special_var_types:
            v = Var(name, typ)
            v._fullname = 'builtins.%s' % name
            file_node.names[name] = SymbolTableNode(GDEF, v)

    #
    # Analyzing a target
    #

    def refresh_partial(self,
                        node: Union[MypyFile, FuncDef, OverloadedFuncDef],
                        patches: List[Tuple[int, Callable[[], None]]],
                        final_iteration: bool,
                        file_node: MypyFile,
                        options: Options,
                        active_type: Optional[TypeInfo] = None) -> None:
        """Refresh a stale target in fine-grained incremental mode."""
        self.patches = patches
        self.deferred = False
        self.incomplete = False
        self._final_iteration = final_iteration
        self.missing_names[-1] = set()

        with self.file_context(file_node, options, active_type):
            if isinstance(node, MypyFile):
                self.refresh_top_level(node)
            else:
                self.recurse_into_functions = True
                self.accept(node)
        del self.patches

    def refresh_top_level(self, file_node: MypyFile) -> None:
        """Reanalyze a stale module top-level in fine-grained incremental mode."""
        self.recurse_into_functions = False
        self.add_implicit_module_attrs(file_node)
        for d in file_node.defs:
            self.accept(d)
        if file_node.fullname == 'typing':
            self.add_builtin_aliases(file_node)
        if file_node.fullname == 'typing_extensions':
            self.add_typing_extension_aliases(file_node)
        self.adjust_public_exports()
        self.export_map[self.cur_mod_id] = self.all_exports
        self.all_exports = []

    def add_implicit_module_attrs(self, file_node: MypyFile) -> None:
        """Manually add implicit definitions of module '__name__' etc."""
        for name, t in implicit_module_attrs.items():
            # unicode docstrings should be accepted in Python 2
            if name == '__doc__':
                if self.options.python_version >= (3, 0):
                    typ: Type = UnboundType("__builtins__.str")
                else:
                    typ = UnionType([UnboundType('__builtins__.str'),
                                     UnboundType('__builtins__.unicode')])
            elif name == '__path__':
                if not file_node.is_package_init_file():
                    continue
                # Need to construct the type ourselves, to avoid issues with __builtins__.list
                # not being subscriptable or typing.List not getting bound
                sym = self.lookup_qualified("__builtins__.list", Context())
                if not sym:
                    continue
                node = sym.node
                assert isinstance(node, TypeInfo)
                typ = Instance(node, [self.str_type()])
            else:
                assert t is not None, 'type should be specified for {}'.format(name)
                typ = UnboundType(t)

            existing = file_node.names.get(name)
            if existing is not None and not isinstance(existing.node, PlaceholderNode):
                # Already exists.
                continue

            an_type = self.anal_type(typ)
            if an_type:
                var = Var(name, an_type)
                var._fullname = self.qualified_name(name)
                var.is_ready = True
                self.add_symbol(name, var, dummy_context())
            else:
                self.add_symbol(name,
                                PlaceholderNode(self.qualified_name(name), file_node, -1),
                                dummy_context())

    def add_builtin_aliases(self, tree: MypyFile) -> None:
        """Add builtin type aliases to typing module.

        For historical reasons, the aliases like `List = list` are not defined
        in typeshed stubs for typing module. Instead we need to manually add the
        corresponding nodes on the fly. We explicitly mark these aliases as normalized,
        so that a user can write `typing.List[int]`.
        """
        assert tree.fullname == 'typing'
        for alias, target_name in type_aliases.items():
            if type_aliases_source_versions[alias] > self.options.python_version:
                # This alias is not available on this Python version.
                continue
            name = alias.split('.')[-1]
            if name in tree.names and not isinstance(tree.names[name].node, PlaceholderNode):
                continue
            self.create_alias(tree, target_name, alias, name)

    def add_typing_extension_aliases(self, tree: MypyFile) -> None:
        """Typing extensions module does contain some type aliases.

        We need to analyze them as such, because in typeshed
        they are just defined as `_Alias()` call.
        Which is not supported natively.
        """
        assert tree.fullname == 'typing_extensions'

        for alias, target_name in typing_extensions_aliases.items():
            name = alias.split('.')[-1]
            if name in tree.names and isinstance(tree.names[name].node, TypeAlias):
                continue  # Do not reset TypeAliases on the second pass.

            # We need to remove any node that is there at the moment. It is invalid.
            tree.names.pop(name, None)

            # Now, create a new alias.
            self.create_alias(tree, target_name, alias, name)

    def create_alias(self, tree: MypyFile, target_name: str, alias: str, name: str) -> None:
        tag = self.track_incomplete_refs()
        n = self.lookup_fully_qualified_or_none(target_name)
        if n:
            if isinstance(n.node, PlaceholderNode):
                self.mark_incomplete(name, tree)
            else:
                # Found built-in class target. Create alias.
                target = self.named_type_or_none(target_name, [])
                assert target is not None
                # Transform List to List[Any], etc.
                fix_instance_types(target, self.fail, self.note, self.options.python_version)
                alias_node = TypeAlias(target, alias,
                                       line=-1, column=-1,  # there is no context
                                       no_args=True, normalized=True)
                self.add_symbol(name, alias_node, tree)
        elif self.found_incomplete_ref(tag):
            # Built-in class target may not ready yet -- defer.
            self.mark_incomplete(name, tree)
        else:
            # Test fixtures may be missing some builtin classes, which is okay.
            # Kill the placeholder if there is one.
            if name in tree.names:
                assert isinstance(tree.names[name].node, PlaceholderNode)
                del tree.names[name]

    def adjust_public_exports(self) -> None:
        """Adjust the module visibility of globals due to __all__."""
        if '__all__' in self.globals:
            for name, g in self.globals.items():
                # Being included in __all__ explicitly exports and makes public.
                if name in self.all_exports:
                    g.module_public = True
                    g.module_hidden = False
                # But when __all__ is defined, and a symbol is not included in it,
                # it cannot be public.
                else:
                    g.module_public = False

    @contextmanager
    def file_context(self,
                     file_node: MypyFile,
                     options: Options,
                     active_type: Optional[TypeInfo] = None) -> Iterator[None]:
        """Configure analyzer for analyzing targets within a file/class.

        Args:
            file_node: target file
            options: options specific to the file
            active_type: must be the surrounding class to analyze method targets
        """
        scope = self.scope
        self.options = options
        self.errors.set_file(file_node.path, file_node.fullname, scope=scope)
        self.cur_mod_node = file_node
        self.cur_mod_id = file_node.fullname
        with scope.module_scope(self.cur_mod_id):
            self._is_stub_file = file_node.path.lower().endswith('.pyi')
            self._is_typeshed_stub_file = is_typeshed_file(file_node.path)
            self.globals = file_node.names
            self.tvar_scope = TypeVarLikeScope()

            self.named_tuple_analyzer = NamedTupleAnalyzer(options, self)
            self.typed_dict_analyzer = TypedDictAnalyzer(options, self, self.msg)
            self.enum_call_analyzer = EnumCallAnalyzer(options, self)
            self.newtype_analyzer = NewTypeAnalyzer(options, self, self.msg)

            # Counter that keeps track of references to undefined things potentially caused by
            # incomplete namespaces.
            self.num_incomplete_refs = 0

            if active_type:
                self.incomplete_type_stack.append(False)
                scope.enter_class(active_type)
                self.enter_class(active_type.defn.info)
                for tvar in active_type.defn.type_vars:
                    self.tvar_scope.bind_existing(tvar)

            yield

            if active_type:
                scope.leave_class()
                self.leave_class()
                self.type = None
                self.incomplete_type_stack.pop()
        del self.options

    #
    # Functions
    #

    def visit_func_def(self, defn: FuncDef) -> None:
        self.statement = defn

        # Visit default values because they may contain assignment expressions.
        for arg in defn.arguments:
            if arg.initializer:
                arg.initializer.accept(self)

        defn.is_conditional = self.block_depth[-1] > 0

        # Set full names even for those definitions that aren't added
        # to a symbol table. For example, for overload items.
        defn._fullname = self.qualified_name(defn.name)

        # We don't add module top-level functions to symbol tables
        # when we analyze their bodies in the second phase on analysis,
        # since they were added in the first phase. Nested functions
        # get always added, since they aren't separate targets.
        if not self.recurse_into_functions or len(self.function_stack) > 0:
            if not defn.is_decorated and not defn.is_overload:
                self.add_function_to_symbol_table(defn)

        if not self.recurse_into_functions:
            return

        with self.scope.function_scope(defn):
            self.analyze_func_def(defn)

    def analyze_func_def(self, defn: FuncDef) -> None:
        self.function_stack.append(defn)

        if defn.type:
            assert isinstance(defn.type, CallableType)
            self.update_function_type_variables(defn.type, defn)
        self.function_stack.pop()

        if self.is_class_scope():
            # Method definition
            assert self.type is not None
            defn.info = self.type
            if defn.type is not None and defn.name in ('__init__', '__init_subclass__'):
                assert isinstance(defn.type, CallableType)
                if isinstance(get_proper_type(defn.type.ret_type), AnyType):
                    defn.type = defn.type.copy_modified(ret_type=NoneType())
            self.prepare_method_signature(defn, self.type)

        # Analyze function signature
        with self.tvar_scope_frame(self.tvar_scope.method_frame()):
            if defn.type:
                self.check_classvar_in_signature(defn.type)
                assert isinstance(defn.type, CallableType)
                # Signature must be analyzed in the surrounding scope so that
                # class-level imported names and type variables are in scope.
                analyzer = self.type_analyzer()
                tag = self.track_incomplete_refs()
                result = analyzer.visit_callable_type(defn.type, nested=False)
                # Don't store not ready types (including placeholders).
                if self.found_incomplete_ref(tag) or has_placeholder(result):
                    self.defer(defn)
                    return
                assert isinstance(result, ProperType)
                defn.type = result
                self.add_type_alias_deps(analyzer.aliases_used)
                self.check_function_signature(defn)
                if isinstance(defn, FuncDef):
                    assert isinstance(defn.type, CallableType)
                    defn.type = set_callable_name(defn.type, defn)

        self.analyze_arg_initializers(defn)
        self.analyze_function_body(defn)
        if (defn.is_coroutine and
                isinstance(defn.type, CallableType) and
                self.wrapped_coro_return_types.get(defn) != defn.type):
            if defn.is_async_generator:
                # Async generator types are handled elsewhere
                pass
            else:
                # A coroutine defined as `async def foo(...) -> T: ...`
                # has external return type `Coroutine[Any, Any, T]`.
                any_type = AnyType(TypeOfAny.special_form)
                ret_type = self.named_type_or_none('typing.Coroutine',
                                                   [any_type, any_type, defn.type.ret_type])
                assert ret_type is not None, "Internal error: typing.Coroutine not found"
                defn.type = defn.type.copy_modified(ret_type=ret_type)
                self.wrapped_coro_return_types[defn] = defn.type

    def prepare_method_signature(self, func: FuncDef, info: TypeInfo) -> None:
        """Check basic signature validity and tweak annotation of self/cls argument."""
        # Only non-static methods are special.
        functype = func.type
        if not func.is_static:
            if func.name in ['__init_subclass__', '__class_getitem__']:
                func.is_class = True
            if not func.arguments:
                self.fail('Method must have at least one argument', func)
            elif isinstance(functype, CallableType):
                self_type = get_proper_type(functype.arg_types[0])
                if isinstance(self_type, AnyType):
                    leading_type: Type = fill_typevars(info)
                    if func.is_class or func.name == '__new__':
                        leading_type = self.class_type(leading_type)
                    func.type = replace_implicit_first_type(functype, leading_type)

    def set_original_def(self, previous: Optional[Node], new: Union[FuncDef, Decorator]) -> bool:
        """If 'new' conditionally redefine 'previous', set 'previous' as original

        We reject straight redefinitions of functions, as they are usually
        a programming error. For example:

          def f(): ...
          def f(): ...  # Error: 'f' redefined
        """
        if isinstance(new, Decorator):
            new = new.func
        if (
            isinstance(previous, (FuncDef, Decorator))
            and unnamed_function(new.name)
            and unnamed_function(previous.name)
        ):
            return True
        if isinstance(previous, (FuncDef, Var, Decorator)) and new.is_conditional:
            new.original_def = previous
            return True
        else:
            return False

    def update_function_type_variables(self, fun_type: CallableType, defn: FuncItem) -> None:
        """Make any type variables in the signature of defn explicit.

        Update the signature of defn to contain type variable definitions
        if defn is generic.
        """
        with self.tvar_scope_frame(self.tvar_scope.method_frame()):
            a = self.type_analyzer()
            fun_type.variables = a.bind_function_type_variables(fun_type, defn)

    def visit_overloaded_func_def(self, defn: OverloadedFuncDef) -> None:
        self.statement = defn
        self.add_function_to_symbol_table(defn)

        if not self.recurse_into_functions:
            return

        # NB: Since _visit_overloaded_func_def will call accept on the
        # underlying FuncDefs, the function might get entered twice.
        # This is fine, though, because only the outermost function is
        # used to compute targets.
        with self.scope.function_scope(defn):
            self.analyze_overloaded_func_def(defn)

    def analyze_overloaded_func_def(self, defn: OverloadedFuncDef) -> None:
        # OverloadedFuncDef refers to any legitimate situation where you have
        # more than one declaration for the same function in a row.  This occurs
        # with a @property with a setter or a deleter, and for a classic
        # @overload.

        defn._fullname = self.qualified_name(defn.name)
        # TODO: avoid modifying items.
        defn.items = defn.unanalyzed_items.copy()

        first_item = defn.items[0]
        first_item.is_overload = True
        first_item.accept(self)

        if isinstance(first_item, Decorator) and first_item.func.is_property:
            # This is a property.
            first_item.func.is_overload = True
            self.analyze_property_with_multi_part_definition(defn)
            typ = function_type(first_item.func, self.named_type('builtins.function'))
            assert isinstance(typ, CallableType)
            types = [typ]
        else:
            # This is an a normal overload. Find the item signatures, the
            # implementation (if outside a stub), and any missing @overload
            # decorators.
            types, impl, non_overload_indexes = self.analyze_overload_sigs_and_impl(defn)
            defn.impl = impl
            if non_overload_indexes:
                self.handle_missing_overload_decorators(defn, non_overload_indexes,
                                                        some_overload_decorators=len(types) > 0)
            # If we found an implementation, remove it from the overload item list,
            # as it's special.
            if impl is not None:
                assert impl is defn.items[-1]
                defn.items = defn.items[:-1]
            elif not non_overload_indexes:
                self.handle_missing_overload_implementation(defn)

        if types:
            defn.type = Overloaded(types)
            defn.type.line = defn.line

        if not defn.items:
            # It was not a real overload after all, but function redefinition. We've
            # visited the redefinition(s) already.
            if not defn.impl:
                # For really broken overloads with no items and no implementation we need to keep
                # at least one item to hold basic information like function name.
                defn.impl = defn.unanalyzed_items[-1]
            return

        # We know this is an overload def. Infer properties and perform some checks.
        self.process_final_in_overload(defn)
        self.process_static_or_class_method_in_overload(defn)

    def analyze_overload_sigs_and_impl(
            self,
            defn: OverloadedFuncDef) -> Tuple[List[CallableType],
                                              Optional[OverloadPart],
                                              List[int]]:
        """Find overload signatures, the implementation, and items with missing @overload.

        Assume that the first was already analyzed. As a side effect:
        analyzes remaining items and updates 'is_overload' flags.
        """
        types = []
        non_overload_indexes = []
        impl: Optional[OverloadPart] = None
        for i, item in enumerate(defn.items):
            if i != 0:
                # Assume that the first item was already visited
                item.is_overload = True
                item.accept(self)
            # TODO: support decorated overloaded functions properly
            if isinstance(item, Decorator):
                callable = function_type(item.func, self.named_type('builtins.function'))
                assert isinstance(callable, CallableType)
                if not any(refers_to_fullname(dec, 'typing.overload')
                           for dec in item.decorators):
                    if i == len(defn.items) - 1 and not self.is_stub_file:
                        # Last item outside a stub is impl
                        impl = item
                    else:
                        # Oops it wasn't an overload after all. A clear error
                        # will vary based on where in the list it is, record
                        # that.
                        non_overload_indexes.append(i)
                else:
                    item.func.is_overload = True
                    types.append(callable)
            elif isinstance(item, FuncDef):
                if i == len(defn.items) - 1 and not self.is_stub_file:
                    impl = item
                else:
                    non_overload_indexes.append(i)
        return types, impl, non_overload_indexes

    def handle_missing_overload_decorators(self,
                                           defn: OverloadedFuncDef,
                                           non_overload_indexes: List[int],
                                           some_overload_decorators: bool) -> None:
        """Generate errors for overload items without @overload.

        Side effect: remote non-overload items.
        """
        if some_overload_decorators:
            # Some of them were overloads, but not all.
            for idx in non_overload_indexes:
                if self.is_stub_file:
                    self.fail("An implementation for an overloaded function "
                              "is not allowed in a stub file", defn.items[idx])
                else:
                    self.fail("The implementation for an overloaded function "
                              "must come last", defn.items[idx])
        else:
            for idx in non_overload_indexes[1:]:
                self.name_already_defined(defn.name, defn.items[idx], defn.items[0])
            if defn.impl:
                self.name_already_defined(defn.name, defn.impl, defn.items[0])
        # Remove the non-overloads
        for idx in reversed(non_overload_indexes):
            del defn.items[idx]

    def handle_missing_overload_implementation(self, defn: OverloadedFuncDef) -> None:
        """Generate error about missing overload implementation (only if needed)."""
        if not self.is_stub_file:
            if self.type and self.type.is_protocol and not self.is_func_scope():
                # An overloaded protocol method doesn't need an implementation.
                for item in defn.items:
                    if isinstance(item, Decorator):
                        item.func.is_abstract = True
                    else:
                        item.is_abstract = True
            else:
                self.fail(
                    "An overloaded function outside a stub file must have an implementation",
                    defn)

    def process_final_in_overload(self, defn: OverloadedFuncDef) -> None:
        """Detect the @final status of an overloaded function (and perform checks)."""
        # If the implementation is marked as @final (or the first overload in
        # stubs), then the whole overloaded definition if @final.
        if any(item.is_final for item in defn.items):
            # We anyway mark it as final because it was probably the intention.
            defn.is_final = True
            # Only show the error once per overload
            bad_final = next(ov for ov in defn.items if ov.is_final)
            if not self.is_stub_file:
                self.fail("@final should be applied only to overload implementation",
                          bad_final)
            elif any(item.is_final for item in defn.items[1:]):
                bad_final = next(ov for ov in defn.items[1:] if ov.is_final)
                self.fail("In a stub file @final must be applied only to the first overload",
                          bad_final)
        if defn.impl is not None and defn.impl.is_final:
            defn.is_final = True

    def process_static_or_class_method_in_overload(self, defn: OverloadedFuncDef) -> None:
        class_status = []
        static_status = []
        for item in defn.items:
            if isinstance(item, Decorator):
                inner = item.func
            elif isinstance(item, FuncDef):
                inner = item
            else:
                assert False, "The 'item' variable is an unexpected type: {}".format(type(item))
            class_status.append(inner.is_class)
            static_status.append(inner.is_static)

        if defn.impl is not None:
            if isinstance(defn.impl, Decorator):
                inner = defn.impl.func
            elif isinstance(defn.impl, FuncDef):
                inner = defn.impl
            else:
                assert False, "Unexpected impl type: {}".format(type(defn.impl))
            class_status.append(inner.is_class)
            static_status.append(inner.is_static)

        if len(set(class_status)) != 1:
            self.msg.overload_inconsistently_applies_decorator('classmethod', defn)
        elif len(set(static_status)) != 1:
            self.msg.overload_inconsistently_applies_decorator('staticmethod', defn)
        else:
            defn.is_class = class_status[0]
            defn.is_static = static_status[0]

    def analyze_property_with_multi_part_definition(self, defn: OverloadedFuncDef) -> None:
        """Analyze a property defined using multiple methods (e.g., using @x.setter).

        Assume that the first method (@property) has already been analyzed.
        """
        defn.is_property = True
        items = defn.items
        first_item = cast(Decorator, defn.items[0])
        deleted_items = []
        for i, item in enumerate(items[1:]):
            if isinstance(item, Decorator):
                if len(item.decorators) == 1:
                    node = item.decorators[0]
                    if isinstance(node, MemberExpr):
                        if node.name == 'setter':
                            # The first item represents the entire property.
                            first_item.var.is_settable_property = True
                            # Get abstractness from the original definition.
                            item.func.is_abstract = first_item.func.is_abstract
                else:
                    self.fail("Decorated property not supported", item)
                item.func.accept(self)
            else:
                self.fail('Unexpected definition for property "{}"'.format(first_item.func.name),
                          item)
                deleted_items.append(i + 1)
        for i in reversed(deleted_items):
            del items[i]

    def add_function_to_symbol_table(self, func: Union[FuncDef, OverloadedFuncDef]) -> None:
        if self.is_class_scope():
            assert self.type is not None
            func.info = self.type
        func._fullname = self.qualified_name(func.name)
        self.add_symbol(func.name, func, func)

    def analyze_arg_initializers(self, defn: FuncItem) -> None:
        with self.tvar_scope_frame(self.tvar_scope.method_frame()):
            # Analyze default arguments
            for arg in defn.arguments:
                if arg.initializer:
                    arg.initializer.accept(self)

    def analyze_function_body(self, defn: FuncItem) -> None:
        is_method = self.is_class_scope()
        with self.tvar_scope_frame(self.tvar_scope.method_frame()):
            # Bind the type variables again to visit the body.
            if defn.type:
                a = self.type_analyzer()
                a.bind_function_type_variables(cast(CallableType, defn.type), defn)
            self.function_stack.append(defn)
            with self.enter(defn):
                for arg in defn.arguments:
                    self.add_local(arg.variable, defn)

                # The first argument of a non-static, non-class method is like 'self'
                # (though the name could be different), having the enclosing class's
                # instance type.
                if is_method and not defn.is_static and not defn.is_class and defn.arguments:
                    defn.arguments[0].variable.is_self = True

                defn.body.accept(self)
            self.function_stack.pop()

    def check_classvar_in_signature(self, typ: ProperType) -> None:
        if isinstance(typ, Overloaded):
            for t in typ.items:  # type: ProperType
                self.check_classvar_in_signature(t)
            return
        if not isinstance(typ, CallableType):
            return
        for t in get_proper_types(typ.arg_types) + [get_proper_type(typ.ret_type)]:
            if self.is_classvar(t):
                self.fail_invalid_classvar(t)
                # Show only one error per signature
                break

    def check_function_signature(self, fdef: FuncItem) -> None:
        sig = fdef.type
        assert isinstance(sig, CallableType)
        if len(sig.arg_types) < len(fdef.arguments):
            self.fail('Type signature has too few arguments', fdef)
            # Add dummy Any arguments to prevent crashes later.
            num_extra_anys = len(fdef.arguments) - len(sig.arg_types)
            extra_anys = [AnyType(TypeOfAny.from_error)] * num_extra_anys
            sig.arg_types.extend(extra_anys)
        elif len(sig.arg_types) > len(fdef.arguments):
            self.fail('Type signature has too many arguments', fdef, blocker=True)

    def visit_decorator(self, dec: Decorator) -> None:
        self.statement = dec
        # TODO: better don't modify them at all.
        dec.decorators = dec.original_decorators.copy()
        dec.func.is_conditional = self.block_depth[-1] > 0
        if not dec.is_overload:
            self.add_symbol(dec.name, dec, dec)
        dec.func._fullname = self.qualified_name(dec.name)
        for d in dec.decorators:
            d.accept(self)
        removed: List[int] = []
        no_type_check = False
        for i, d in enumerate(dec.decorators):
            # A bunch of decorators are special cased here.
            if refers_to_fullname(d, 'abc.abstractmethod'):
                removed.append(i)
                dec.func.is_abstract = True
                self.check_decorated_function_is_method('abstractmethod', dec)
            elif (refers_to_fullname(d, 'asyncio.coroutines.coroutine') or
                  refers_to_fullname(d, 'types.coroutine')):
                removed.append(i)
                dec.func.is_awaitable_coroutine = True
            elif refers_to_fullname(d, 'builtins.staticmethod'):
                removed.append(i)
                dec.func.is_static = True
                dec.var.is_staticmethod = True
                self.check_decorated_function_is_method('staticmethod', dec)
            elif refers_to_fullname(d, 'builtins.classmethod'):
                removed.append(i)
                dec.func.is_class = True
                dec.var.is_classmethod = True
                self.check_decorated_function_is_method('classmethod', dec)
            elif (refers_to_fullname(d, 'builtins.property') or
                  refers_to_fullname(d, 'abc.abstractproperty') or
                  refers_to_fullname(d, 'functools.cached_property')):
                removed.append(i)
                dec.func.is_property = True
                dec.var.is_property = True
                if refers_to_fullname(d, 'abc.abstractproperty'):
                    dec.func.is_abstract = True
                elif refers_to_fullname(d, 'functools.cached_property'):
                    dec.var.is_settable_property = True
                self.check_decorated_function_is_method('property', dec)
                if len(dec.func.arguments) > 1:
                    self.fail('Too many arguments', dec.func)
            elif refers_to_fullname(d, 'typing.no_type_check'):
                dec.var.type = AnyType(TypeOfAny.special_form)
                no_type_check = True
            elif (refers_to_fullname(d, 'typing.final') or
                  refers_to_fullname(d, 'typing_extensions.final')):
                if self.is_class_scope():
                    assert self.type is not None, "No type set at class scope"
                    if self.type.is_protocol:
                        self.msg.protocol_members_cant_be_final(d)
                    else:
                        dec.func.is_final = True
                        dec.var.is_final = True
                    removed.append(i)
                else:
                    self.fail("@final cannot be used with non-method functions", d)
        for i in reversed(removed):
            del dec.decorators[i]
        if (not dec.is_overload or dec.var.is_property) and self.type:
            dec.var.info = self.type
            dec.var.is_initialized_in_class = True
        if not no_type_check and self.recurse_into_functions:
            dec.func.accept(self)
        if dec.decorators and dec.var.is_property:
            self.fail('Decorated property not supported', dec)

    def check_decorated_function_is_method(self, decorator: str,
                                           context: Context) -> None:
        if not self.type or self.is_func_scope():
            self.fail('"%s" used with a non-method' % decorator, context)

    #
    # Classes
    #

    def visit_class_def(self, defn: ClassDef) -> None:
        self.statement = defn
        self.incomplete_type_stack.append(not defn.info)
        with self.tvar_scope_frame(self.tvar_scope.class_frame()):
            self.analyze_class(defn)
        self.incomplete_type_stack.pop()

    def analyze_class(self, defn: ClassDef) -> None:
        fullname = self.qualified_name(defn.name)
        if not defn.info and not self.is_core_builtin_class(defn):
            # Add placeholder so that self-references in base classes can be
            # resolved.  We don't want this to cause a deferral, since if there
            # are no incomplete references, we'll replace this with a TypeInfo
            # before returning.
            placeholder = PlaceholderNode(fullname, defn, defn.line, becomes_typeinfo=True)
            self.add_symbol(defn.name, placeholder, defn, can_defer=False)

        tag = self.track_incomplete_refs()

        # Restore base classes after previous iteration (things like Generic[T] might be removed).
        defn.base_type_exprs.extend(defn.removed_base_type_exprs)
        defn.removed_base_type_exprs.clear()

        self.update_metaclass(defn)

        bases = defn.base_type_exprs
        bases, tvar_defs, is_protocol = self.clean_up_bases_and_infer_type_variables(
            defn, bases, context=defn)

        for tvd in tvar_defs:
            if (isinstance(tvd, TypeVarType)
                    and any(has_placeholder(t) for t in [tvd.upper_bound] + tvd.values)):
                # Some type variable bounds or values are not ready, we need
                # to re-analyze this class.
                self.defer()

        self.analyze_class_keywords(defn)
        result = self.analyze_base_classes(bases)

        if result is None or self.found_incomplete_ref(tag):
            # Something was incomplete. Defer current target.
            self.mark_incomplete(defn.name, defn)
            return

        base_types, base_error = result
        if any(isinstance(base, PlaceholderType) for base, _ in base_types):
            # We need to know the TypeInfo of each base to construct the MRO. Placeholder types
            # are okay in nested positions, since they can't affect the MRO.
            self.mark_incomplete(defn.name, defn)
            return

        is_typeddict, info = self.typed_dict_analyzer.analyze_typeddict_classdef(defn)
        if is_typeddict:
            for decorator in defn.decorators:
                decorator.accept(self)
                if isinstance(decorator, RefExpr):
                    if decorator.fullname in ('typing.final',
                                              'typing_extensions.final'):
                        self.fail("@final cannot be used with TypedDict", decorator)
            if info is None:
                self.mark_incomplete(defn.name, defn)
            else:
                self.prepare_class_def(defn, info)
            return

        if self.analyze_namedtuple_classdef(defn):
            return

        # Create TypeInfo for class now that base classes and the MRO can be calculated.
        self.prepare_class_def(defn)

        defn.type_vars = tvar_defs
        defn.info.type_vars = [tvar.name for tvar in tvar_defs]
        if base_error:
            defn.info.fallback_to_any = True

        with self.scope.class_scope(defn.info):
            self.configure_base_classes(defn, base_types)
            defn.info.is_protocol = is_protocol
            self.analyze_metaclass(defn)
            defn.info.runtime_protocol = False
            for decorator in defn.decorators:
                self.analyze_class_decorator(defn, decorator)
            self.analyze_class_body_common(defn)

    def is_core_builtin_class(self, defn: ClassDef) -> bool:
        return self.cur_mod_id == 'builtins' and defn.name in CORE_BUILTIN_CLASSES

    def analyze_class_body_common(self, defn: ClassDef) -> None:
        """Parts of class body analysis that are common to all kinds of class defs."""
        self.enter_class(defn.info)
        defn.defs.accept(self)
        self.apply_class_plugin_hooks(defn)
        self.leave_class()

    def analyze_namedtuple_classdef(self, defn: ClassDef) -> bool:
        """Check if this class can define a named tuple."""
        if defn.info and defn.info.is_named_tuple:
            # Don't reprocess everything. We just need to process methods defined
            # in the named tuple class body.
            is_named_tuple, info = True, defn.info  # type: bool, Optional[TypeInfo]
        else:
            is_named_tuple, info = self.named_tuple_analyzer.analyze_namedtuple_classdef(
                defn, self.is_stub_file)
        if is_named_tuple:
            if info is None:
                self.mark_incomplete(defn.name, defn)
            else:
                self.prepare_class_def(defn, info)
                with self.scope.class_scope(defn.info):
                    with self.named_tuple_analyzer.save_namedtuple_body(info):
                        self.analyze_class_body_common(defn)
            return True
        return False

    def apply_class_plugin_hooks(self, defn: ClassDef) -> None:
        """Apply a plugin hook that may infer a more precise definition for a class."""
        def get_fullname(expr: Expression) -> Optional[str]:
            if isinstance(expr, CallExpr):
                return get_fullname(expr.callee)
            elif isinstance(expr, IndexExpr):
                return get_fullname(expr.base)
            elif isinstance(expr, RefExpr):
                if expr.fullname:
                    return expr.fullname
                # If we don't have a fullname look it up. This happens because base classes are
                # analyzed in a different manner (see exprtotype.py) and therefore those AST
                # nodes will not have full names.
                sym = self.lookup_type_node(expr)
                if sym:
                    return sym.fullname
            return None

        for decorator in defn.decorators:
            decorator_name = get_fullname(decorator)
            if decorator_name:
                hook = self.plugin.get_class_decorator_hook(decorator_name)
                if hook:
                    hook(ClassDefContext(defn, decorator, self))

        if defn.metaclass:
            metaclass_name = get_fullname(defn.metaclass)
            if metaclass_name:
                hook = self.plugin.get_metaclass_hook(metaclass_name)
                if hook:
                    hook(ClassDefContext(defn, defn.metaclass, self))

        for base_expr in defn.base_type_exprs:
            base_name = get_fullname(base_expr)
            if base_name:
                hook = self.plugin.get_base_class_hook(base_name)
                if hook:
                    hook(ClassDefContext(defn, base_expr, self))

    def analyze_class_keywords(self, defn: ClassDef) -> None:
        for value in defn.keywords.values():
            value.accept(self)

    def enter_class(self, info: TypeInfo) -> None:
        # Remember previous active class
        self.type_stack.append(self.type)
        self.locals.append(None)  # Add class scope
        self.is_comprehension_stack.append(False)
        self.block_depth.append(-1)  # The class body increments this to 0
        self.type = info
        self.missing_names.append(set())

    def leave_class(self) -> None:
        """ Restore analyzer state. """
        self.block_depth.pop()
        self.locals.pop()
        self.is_comprehension_stack.pop()
        self.type = self.type_stack.pop()
        self.missing_names.pop()

    def analyze_class_decorator(self, defn: ClassDef, decorator: Expression) -> None:
        decorator.accept(self)
        if isinstance(decorator, RefExpr):
            if decorator.fullname in RUNTIME_PROTOCOL_DECOS:
                if defn.info.is_protocol:
                    defn.info.runtime_protocol = True
                else:
                    self.fail('@runtime_checkable can only be used with protocol classes',
                              defn)
            elif decorator.fullname in ('typing.final',
                                        'typing_extensions.final'):
                defn.info.is_final = True

    def clean_up_bases_and_infer_type_variables(
            self,
            defn: ClassDef,
            base_type_exprs: List[Expression],
            context: Context) -> Tuple[List[Expression],
                                       List[TypeVarLikeType],
                                       bool]:
        """Remove extra base classes such as Generic and infer type vars.

        For example, consider this class:

          class Foo(Bar, Generic[T]): ...

        Now we will remove Generic[T] from bases of Foo and infer that the
        type variable 'T' is a type argument of Foo.

        Note that this is performed *before* semantic analysis.

        Returns (remaining base expressions, inferred type variables, is protocol).
        """
        removed: List[int] = []
        declared_tvars: TypeVarLikeList = []
        is_protocol = False
        for i, base_expr in enumerate(base_type_exprs):
            self.analyze_type_expr(base_expr)

            try:
                base = self.expr_to_unanalyzed_type(base_expr)
            except TypeTranslationError:
                # This error will be caught later.
                continue
            result = self.analyze_class_typevar_declaration(base)
            if result is not None:
                if declared_tvars:
                    self.fail('Only single Generic[...] or Protocol[...] can be in bases', context)
                removed.append(i)
                tvars = result[0]
                is_protocol |= result[1]
                declared_tvars.extend(tvars)
            if isinstance(base, UnboundType):
                sym = self.lookup_qualified(base.name, base)
                if sym is not None and sym.node is not None:
                    if (sym.node.fullname in ('typing.Protocol', 'typing_extensions.Protocol') and
                            i not in removed):
                        # also remove bare 'Protocol' bases
                        removed.append(i)
                        is_protocol = True

        all_tvars = self.get_all_bases_tvars(base_type_exprs, removed)
        if declared_tvars:
            if len(remove_dups(declared_tvars)) < len(declared_tvars):
                self.fail("Duplicate type variables in Generic[...] or Protocol[...]", context)
            declared_tvars = remove_dups(declared_tvars)
            if not set(all_tvars).issubset(set(declared_tvars)):
                self.fail("If Generic[...] or Protocol[...] is present"
                          " it should list all type variables", context)
                # In case of error, Generic tvars will go first
                declared_tvars = remove_dups(declared_tvars + all_tvars)
        else:
            declared_tvars = all_tvars
        for i in reversed(removed):
            # We need to actually remove the base class expressions like Generic[T],
            # mostly because otherwise they will create spurious dependencies in fine
            # grained incremental mode.
            defn.removed_base_type_exprs.append(defn.base_type_exprs[i])
            del base_type_exprs[i]
        tvar_defs: List[TypeVarLikeType] = []
        for name, tvar_expr in declared_tvars:
            tvar_def = self.tvar_scope.bind_new(name, tvar_expr)
            tvar_defs.append(tvar_def)
        return base_type_exprs, tvar_defs, is_protocol

    def analyze_class_typevar_declaration(
        self,
        base: Type
    ) -> Optional[Tuple[TypeVarLikeList, bool]]:
        """Analyze type variables declared using Generic[...] or Protocol[...].

        Args:
            base: Non-analyzed base class

        Return None if the base class does not declare type variables. Otherwise,
        return the type variables.
        """
        if not isinstance(base, UnboundType):
            return None
        unbound = base
        sym = self.lookup_qualified(unbound.name, unbound)
        if sym is None or sym.node is None:
            return None
        if (sym.node.fullname == 'typing.Generic' or
                sym.node.fullname == 'typing.Protocol' and base.args or
                sym.node.fullname == 'typing_extensions.Protocol' and base.args):
            is_proto = sym.node.fullname != 'typing.Generic'
            tvars: TypeVarLikeList = []
            for arg in unbound.args:
                tag = self.track_incomplete_refs()
                tvar = self.analyze_unbound_tvar(arg)
                if tvar:
                    tvars.append(tvar)
                elif not self.found_incomplete_ref(tag):
                    self.fail('Free type variable expected in %s[...]' %
                              sym.node.name, base)
            return tvars, is_proto
        return None

    def analyze_unbound_tvar(self, t: Type) -> Optional[Tuple[str, TypeVarLikeExpr]]:
        if not isinstance(t, UnboundType):
            return None
        unbound = t
        sym = self.lookup_qualified(unbound.name, unbound)
        if sym and isinstance(sym.node, PlaceholderNode):
            self.record_incomplete_ref()
        if sym and isinstance(sym.node, ParamSpecExpr):
            if sym.fullname and not self.tvar_scope.allow_binding(sym.fullname):
                # It's bound by our type variable scope
                return None
            return unbound.name, sym.node
        if sym is None or not isinstance(sym.node, TypeVarExpr):
            return None
        elif sym.fullname and not self.tvar_scope.allow_binding(sym.fullname):
            # It's bound by our type variable scope
            return None
        else:
            assert isinstance(sym.node, TypeVarExpr)
            return unbound.name, sym.node

    def get_all_bases_tvars(self,
                            base_type_exprs: List[Expression],
                            removed: List[int]) -> TypeVarLikeList:
        """Return all type variable references in bases."""
        tvars: TypeVarLikeList = []
        for i, base_expr in enumerate(base_type_exprs):
            if i not in removed:
                try:
                    base = self.expr_to_unanalyzed_type(base_expr)
                except TypeTranslationError:
                    # This error will be caught later.
                    continue
                base_tvars = base.accept(TypeVarLikeQuery(self.lookup_qualified, self.tvar_scope))
                tvars.extend(base_tvars)
        return remove_dups(tvars)

    def prepare_class_def(self, defn: ClassDef, info: Optional[TypeInfo] = None) -> None:
        """Prepare for the analysis of a class definition.

        Create an empty TypeInfo and store it in a symbol table, or if the 'info'
        argument is provided, store it instead (used for magic type definitions).
        """
        if not defn.info:
            defn.fullname = self.qualified_name(defn.name)
            # TODO: Nested classes
            info = info or self.make_empty_type_info(defn)
            defn.info = info
            info.defn = defn
            if not self.is_func_scope():
                info._fullname = self.qualified_name(defn.name)
            else:
                info._fullname = info.name
        self.add_symbol(defn.name, defn.info, defn)
        if self.is_nested_within_func_scope():
            # We need to preserve local classes, let's store them
            # in globals under mangled unique names
            #
            # TODO: Putting local classes into globals breaks assumptions in fine-grained
            #       incremental mode and we should avoid it. In general, this logic is too
            #       ad-hoc and needs to be removed/refactored.
            if '@' not in defn.info._fullname:
                local_name = defn.info.name + '@' + str(defn.line)
                if defn.info.is_named_tuple:
                    # Module is already correctly set in _fullname for named tuples.
                    defn.info._fullname += '@' + str(defn.line)
                else:
                    defn.info._fullname = self.cur_mod_id + '.' + local_name
            else:
                # Preserve name from previous fine-grained incremental run.
                local_name = defn.info.name
            defn.fullname = defn.info._fullname
            self.globals[local_name] = SymbolTableNode(GDEF, defn.info)

    def make_empty_type_info(self, defn: ClassDef) -> TypeInfo:
        if (self.is_module_scope()
                and self.cur_mod_id == 'builtins'
                and defn.name in CORE_BUILTIN_CLASSES):
            # Special case core built-in classes. A TypeInfo was already
            # created for it before semantic analysis, but with a dummy
            # ClassDef. Patch the real ClassDef object.
            info = self.globals[defn.name].node
            assert isinstance(info, TypeInfo)
        else:
            info = TypeInfo(SymbolTable(), defn, self.cur_mod_id)
            info.set_line(defn)
        return info

    def get_name_repr_of_expr(self, expr: Expression) -> Optional[str]:
        """Try finding a short simplified textual representation of a base class expression."""
        if isinstance(expr, NameExpr):
            return expr.name
        if isinstance(expr, MemberExpr):
            return get_member_expr_fullname(expr)
        if isinstance(expr, IndexExpr):
            return self.get_name_repr_of_expr(expr.base)
        if isinstance(expr, CallExpr):
            return self.get_name_repr_of_expr(expr.callee)
        return None

    def analyze_base_classes(
            self,
            base_type_exprs: List[Expression]) -> Optional[Tuple[List[Tuple[ProperType,
                                                                            Expression]],
                                                                 bool]]:
        """Analyze base class types.

        Return None if some definition was incomplete. Otherwise, return a tuple
        with these items:

         * List of (analyzed type, original expression) tuples
         * Boolean indicating whether one of the bases had a semantic analysis error
        """
        is_error = False
        bases = []
        for base_expr in base_type_exprs:
            if (isinstance(base_expr, RefExpr) and
                    base_expr.fullname in ('typing.NamedTuple',) + TPDICT_NAMES):
                # Ignore magic bases for now.
                continue

            try:
                base = self.expr_to_analyzed_type(base_expr, allow_placeholder=True)
            except TypeTranslationError:
                name = self.get_name_repr_of_expr(base_expr)
                if isinstance(base_expr, CallExpr):
                    msg = 'Unsupported dynamic base class'
                else:
                    msg = 'Invalid base class'
                if name:
                    msg += ' "{}"'.format(name)
                self.fail(msg, base_expr)
                is_error = True
                continue
            if base is None:
                return None
            base = get_proper_type(base)
            bases.append((base, base_expr))
        return bases, is_error

    def configure_base_classes(self,
                               defn: ClassDef,
                               bases: List[Tuple[ProperType, Expression]]) -> None:
        """Set up base classes.

        This computes several attributes on the corresponding TypeInfo defn.info
        related to the base classes: defn.info.bases, defn.info.mro, and
        miscellaneous others (at least tuple_type, fallback_to_any, and is_enum.)
        """
        base_types: List[Instance] = []
        info = defn.info

        info.tuple_type = None
        for base, base_expr in bases:
            if isinstance(base, TupleType):
                actual_base = self.configure_tuple_base_class(defn, base, base_expr)
                base_types.append(actual_base)
            elif isinstance(base, Instance):
                if base.type.is_newtype:
                    self.fail('Cannot subclass "NewType"', defn)
                if self.enum_has_final_values(base):
                    # This means that are trying to subclass a non-default
                    # Enum class, with defined members. This is not possible.
                    # In runtime, it will raise. We need to mark this type as final.
                    # However, methods can be defined on a type: only values can't.
                    # We also don't count values with annotations only.
                    base.type.is_final = True
                base_types.append(base)
            elif isinstance(base, AnyType):
                if self.options.disallow_subclassing_any:
                    if isinstance(base_expr, (NameExpr, MemberExpr)):
                        msg = 'Class cannot subclass "{}" (has type "Any")'.format(base_expr.name)
                    else:
                        msg = 'Class cannot subclass value of type "Any"'
                    self.fail(msg, base_expr)
                info.fallback_to_any = True
            else:
                msg = 'Invalid base class'
                name = self.get_name_repr_of_expr(base_expr)
                if name:
                    msg += ' "{}"'.format(name)
                self.fail(msg, base_expr)
                info.fallback_to_any = True
            if self.options.disallow_any_unimported and has_any_from_unimported_type(base):
                if isinstance(base_expr, (NameExpr, MemberExpr)):
                    prefix = "Base type {}".format(base_expr.name)
                else:
                    prefix = "Base type"
                self.msg.unimported_type_becomes_any(prefix, base, base_expr)
            check_for_explicit_any(base, self.options, self.is_typeshed_stub_file, self.msg,
                                   context=base_expr)

        # Add 'object' as implicit base if there is no other base class.
        if not base_types and defn.fullname != 'builtins.object':
            base_types.append(self.object_type())

        info.bases = base_types

        # Calculate the MRO.
        if not self.verify_base_classes(defn):
            self.set_dummy_mro(defn.info)
            return
        self.calculate_class_mro(defn, self.object_type)

    def enum_has_final_values(self, base: Instance) -> bool:
        if (
            base.type.is_enum
            and base.type.fullname not in ENUM_BASES
            and base.type.names
            and base.type.defn
        ):
            for sym in base.type.names.values():
                if isinstance(sym.node, (FuncBase, Decorator)):
                    continue  # A method
                if not isinstance(sym.node, Var):
                    return True  # Can be a class
                if self.is_stub_file or sym.node.has_explicit_value:
                    # Corner case: assignments like `x: int` are fine in `.py` files.
                    # But, not is `.pyi` files, because we don't know
                    # if there's aactually a value or not.
                    return True
        return False

    def configure_tuple_base_class(self,
                                   defn: ClassDef,
                                   base: TupleType,
                                   base_expr: Expression) -> Instance:
        info = defn.info

        # There may be an existing valid tuple type from previous semanal iterations.
        # Use equality to check if it is the case.
        if info.tuple_type and info.tuple_type != base:
            self.fail("Class has two incompatible bases derived from tuple", defn)
            defn.has_incompatible_baseclass = True
        info.tuple_type = base
        if isinstance(base_expr, CallExpr):
            defn.analyzed = NamedTupleExpr(base.partial_fallback.type)
            defn.analyzed.line = defn.line
            defn.analyzed.column = defn.column

        if base.partial_fallback.type.fullname == 'builtins.tuple':
            # Fallback can only be safely calculated after semantic analysis, since base
            # classes may be incomplete. Postpone the calculation.
            self.schedule_patch(PRIORITY_FALLBACKS, lambda: calculate_tuple_fallback(base))

        return base.partial_fallback

    def set_dummy_mro(self, info: TypeInfo) -> None:
        # Give it an MRO consisting of just the class itself and object.
        info.mro = [info, self.object_type().type]
        info.bad_mro = True

    def calculate_class_mro(self, defn: ClassDef,
                            obj_type: Optional[Callable[[], Instance]] = None) -> None:
        """Calculate method resolution order for a class.

        `obj_type` may be omitted in the third pass when all classes are already analyzed.
        It exists just to fill in empty base class list during second pass in case of
        an import cycle.
        """
        try:
            calculate_mro(defn.info, obj_type)
        except MroError:
            self.fail('Cannot determine consistent method resolution '
                      'order (MRO) for "%s"' % defn.name, defn)
            self.set_dummy_mro(defn.info)
        # Allow plugins to alter the MRO to handle the fact that `def mro()`
        # on metaclasses permits MRO rewriting.
        if defn.fullname:
            hook = self.plugin.get_customize_class_mro_hook(defn.fullname)
            if hook:
                hook(ClassDefContext(defn, FakeExpression(), self))

    def update_metaclass(self, defn: ClassDef) -> None:
        """Lookup for special metaclass declarations, and update defn fields accordingly.

        * __metaclass__ attribute in Python 2
        * six.with_metaclass(M, B1, B2, ...)
        * @six.add_metaclass(M)
        * future.utils.with_metaclass(M, B1, B2, ...)
        * past.utils.with_metaclass(M, B1, B2, ...)
        """

        # Look for "__metaclass__ = <metaclass>" in Python 2
        python2_meta_expr: Optional[Expression] = None
        if self.options.python_version[0] == 2:
            for body_node in defn.defs.body:
                if isinstance(body_node, ClassDef) and body_node.name == "__metaclass__":
                    self.fail("Metaclasses defined as inner classes are not supported", body_node)
                    break
                elif isinstance(body_node, AssignmentStmt) and len(body_node.lvalues) == 1:
                    lvalue = body_node.lvalues[0]
                    if isinstance(lvalue, NameExpr) and lvalue.name == "__metaclass__":
                        python2_meta_expr = body_node.rvalue

        # Look for six.with_metaclass(M, B1, B2, ...)
        with_meta_expr: Optional[Expression] = None
        if len(defn.base_type_exprs) == 1:
            base_expr = defn.base_type_exprs[0]
            if isinstance(base_expr, CallExpr) and isinstance(base_expr.callee, RefExpr):
                base_expr.accept(self)
                if (base_expr.callee.fullname in {'six.with_metaclass',
                                                  'future.utils.with_metaclass',
                                                  'past.utils.with_metaclass'}
                        and len(base_expr.args) >= 1
                        and all(kind == ARG_POS for kind in base_expr.arg_kinds)):
                    with_meta_expr = base_expr.args[0]
                    defn.base_type_exprs = base_expr.args[1:]

        # Look for @six.add_metaclass(M)
        add_meta_expr: Optional[Expression] = None
        for dec_expr in defn.decorators:
            if isinstance(dec_expr, CallExpr) and isinstance(dec_expr.callee, RefExpr):
                dec_expr.callee.accept(self)
                if (dec_expr.callee.fullname == 'six.add_metaclass'
                    and len(dec_expr.args) == 1
                        and dec_expr.arg_kinds[0] == ARG_POS):
                    add_meta_expr = dec_expr.args[0]
                    break

        metas = {defn.metaclass, python2_meta_expr, with_meta_expr, add_meta_expr} - {None}
        if len(metas) == 0:
            return
        if len(metas) > 1:
            self.fail("Multiple metaclass definitions", defn)
            return
        defn.metaclass = metas.pop()

    def verify_base_classes(self, defn: ClassDef) -> bool:
        info = defn.info
        cycle = False
        for base in info.bases:
            baseinfo = base.type
            if self.is_base_class(info, baseinfo):
                self.fail('Cycle in inheritance hierarchy', defn)
                cycle = True
            if baseinfo.fullname == 'builtins.bool':
                self.fail('"%s" is not a valid base class' %
                          baseinfo.name, defn, blocker=True)
                return False
        dup = find_duplicate(info.direct_base_classes())
        if dup:
            self.fail('Duplicate base class "%s"' % dup.name, defn, blocker=True)
            return False
        return not cycle

    def is_base_class(self, t: TypeInfo, s: TypeInfo) -> bool:
        """Determine if t is a base class of s (but do not use mro)."""
        # Search the base class graph for t, starting from s.
        worklist = [s]
        visited = {s}
        while worklist:
            nxt = worklist.pop()
            if nxt == t:
                return True
            for base in nxt.bases:
                if base.type not in visited:
                    worklist.append(base.type)
                    visited.add(base.type)
        return False

    def analyze_metaclass(self, defn: ClassDef) -> None:
        if defn.metaclass:
            metaclass_name = None
            if isinstance(defn.metaclass, NameExpr):
                metaclass_name = defn.metaclass.name
            elif isinstance(defn.metaclass, MemberExpr):
                metaclass_name = get_member_expr_fullname(defn.metaclass)
            if metaclass_name is None:
                self.fail('Dynamic metaclass not supported for "%s"' % defn.name, defn.metaclass)
                return
            sym = self.lookup_qualified(metaclass_name, defn.metaclass)
            if sym is None:
                # Probably a name error - it is already handled elsewhere
                return
            if isinstance(sym.node, Var) and isinstance(get_proper_type(sym.node.type), AnyType):
                # 'Any' metaclass -- just ignore it.
                #
                # TODO: A better approach would be to record this information
                #       and assume that the type object supports arbitrary
                #       attributes, similar to an 'Any' base class.
                return
            if isinstance(sym.node, PlaceholderNode):
                self.defer(defn)
                return
            if not isinstance(sym.node, TypeInfo) or sym.node.tuple_type is not None:
                self.fail('Invalid metaclass "%s"' % metaclass_name, defn.metaclass)
                return
            if not sym.node.is_metaclass():
                self.fail('Metaclasses not inheriting from "type" are not supported',
                          defn.metaclass)
                return
            inst = fill_typevars(sym.node)
            assert isinstance(inst, Instance)
            defn.info.declared_metaclass = inst
        defn.info.metaclass_type = defn.info.calculate_metaclass_type()
        if any(info.is_protocol for info in defn.info.mro):
            if (not defn.info.metaclass_type or
                    defn.info.metaclass_type.type.fullname == 'builtins.type'):
                # All protocols and their subclasses have ABCMeta metaclass by default.
                # TODO: add a metaclass conflict check if there is another metaclass.
                abc_meta = self.named_type_or_none('abc.ABCMeta', [])
                if abc_meta is not None:  # May be None in tests with incomplete lib-stub.
                    defn.info.metaclass_type = abc_meta
        if defn.info.metaclass_type is None:
            # Inconsistency may happen due to multiple baseclasses even in classes that
            # do not declare explicit metaclass, but it's harder to catch at this stage
            if defn.metaclass is not None:
                self.fail('Inconsistent metaclass structure for "%s"' % defn.name, defn)
        else:
            if defn.info.metaclass_type.type.has_base('enum.EnumMeta'):
                defn.info.is_enum = True
                if defn.type_vars:
                    self.fail("Enum class cannot be generic", defn)

    #
    # Imports
    #

    def visit_import(self, i: Import) -> None:
        self.statement = i
        for id, as_id in i.ids:
            # Modules imported in a stub file without using 'import X as X' won't get exported
            # When implicit re-exporting is disabled, we have the same behavior as stubs.
            use_implicit_reexport = not self.is_stub_file and self.options.implicit_reexport
            if as_id is not None:
                base_id = id
                imported_id = as_id
                module_public = use_implicit_reexport or id.split(".")[-1] == as_id
            else:
                base_id = id.split('.')[0]
                imported_id = base_id
                module_public = use_implicit_reexport
            self.add_module_symbol(base_id, imported_id, context=i, module_public=module_public,
                                   module_hidden=not module_public)

    def visit_import_from(self, imp: ImportFrom) -> None:
        self.statement = imp
        module_id = self.correct_relative_import(imp)
        module = self.modules.get(module_id)
        for id, as_id in imp.names:
            fullname = module_id + '.' + id
            self.set_future_import_flags(fullname)
            if module is None:
                node = None
            elif module_id == self.cur_mod_id and fullname in self.modules:
                # Submodule takes precedence over definition in surround package, for
                # compatibility with runtime semantics in typical use cases. This
                # could more precisely model runtime semantics by taking into account
                # the line number beyond which the local definition should take
                # precedence, but doesn't seem to be important in most use cases.
                node = SymbolTableNode(GDEF, self.modules[fullname])
            else:
                if id == as_id == '__all__' and module_id in self.export_map:
                    self.all_exports[:] = self.export_map[module_id]
                node = module.names.get(id)

            missing_submodule = False
            imported_id = as_id or id

            # Modules imported in a stub file without using 'from Y import X as X' will
            # not get exported.
            # When implicit re-exporting is disabled, we have the same behavior as stubs.
            use_implicit_reexport = not self.is_stub_file and self.options.implicit_reexport
            module_public = use_implicit_reexport or (as_id is not None and id == as_id)

            # If the module does not contain a symbol with the name 'id',
            # try checking if it's a module instead.
            if not node:
                mod = self.modules.get(fullname)
                if mod is not None:
                    kind = self.current_symbol_kind()
                    node = SymbolTableNode(kind, mod)
                elif fullname in self.missing_modules:
                    missing_submodule = True
            # If it is still not resolved, check for a module level __getattr__
            if (module and not node and (module.is_stub or self.options.python_version >= (3, 7))
                    and '__getattr__' in module.names):
                # We store the fullname of the original definition so that we can
                # detect whether two imported names refer to the same thing.
                fullname = module_id + '.' + id
                gvar = self.create_getattr_var(module.names['__getattr__'], imported_id, fullname)
                if gvar:
                    self.add_symbol(
                        imported_id, gvar, imp, module_public=module_public,
                        module_hidden=not module_public
                    )
                    continue

            if node and not node.module_hidden:
                self.process_imported_symbol(
                    node, module_id, id, imported_id, fullname, module_public, context=imp
                )
            elif module and not missing_submodule:
                # Target module exists but the imported name is missing or hidden.
                self.report_missing_module_attribute(
                    module_id, id, imported_id, module_public=module_public,
                    module_hidden=not module_public, context=imp
                )
            else:
                # Import of a missing (sub)module.
                self.add_unknown_imported_symbol(
                    imported_id, imp, target_name=fullname, module_public=module_public,
                    module_hidden=not module_public
                )

    def process_imported_symbol(self,
                                node: SymbolTableNode,
                                module_id: str,
                                id: str,
                                imported_id: str,
                                fullname: str,
                                module_public: bool,
                                context: ImportBase) -> None:
        module_hidden = not module_public and fullname not in self.modules

        if isinstance(node.node, PlaceholderNode):
            if self.final_iteration:
                self.report_missing_module_attribute(
                    module_id, id, imported_id, module_public=module_public,
                    module_hidden=module_hidden, context=context
                )
                return
            else:
                # This might become a type.
                self.mark_incomplete(imported_id, node.node,
                                     module_public=module_public,
                                     module_hidden=module_hidden,
                                     becomes_typeinfo=True)
        existing_symbol = self.globals.get(imported_id)
        if (existing_symbol and not isinstance(existing_symbol.node, PlaceholderNode) and
                not isinstance(node.node, PlaceholderNode)):
            # Import can redefine a variable. They get special treatment.
            if self.process_import_over_existing_name(
                    imported_id, existing_symbol, node, context):
                return
        if existing_symbol and isinstance(node.node, PlaceholderNode):
            # Imports are special, some redefinitions are allowed, so wait until
            # we know what is the new symbol node.
            return
        # NOTE: we take the original node even for final `Var`s. This is to support
        # a common pattern when constants are re-exported (same applies to import *).
        self.add_imported_symbol(imported_id, node, context,
                                 module_public=module_public,
                                 module_hidden=module_hidden)

    def report_missing_module_attribute(
        self, import_id: str, source_id: str, imported_id: str, module_public: bool,
        module_hidden: bool, context: Node
    ) -> None:
        # Missing attribute.
        if self.is_incomplete_namespace(import_id):
            # We don't know whether the name will be there, since the namespace
            # is incomplete. Defer the current target.
            self.mark_incomplete(
                imported_id, context, module_public=module_public, module_hidden=module_hidden
            )
            return
        message = 'Module "{}" has no attribute "{}"'.format(import_id, source_id)
        # Suggest alternatives, if any match is found.
        module = self.modules.get(import_id)
        if module:
            if not self.options.implicit_reexport and source_id in module.names.keys():
                message = ('Module "{}" does not explicitly export attribute "{}"'
                           '; implicit reexport disabled'.format(import_id, source_id))
            else:
                alternatives = set(module.names.keys()).difference({source_id})
                matches = best_matches(source_id, alternatives)[:3]
                if matches:
                    suggestion = "; maybe {}?".format(pretty_seq(matches, "or"))
                    message += "{}".format(suggestion)
        self.fail(message, context, code=codes.ATTR_DEFINED)
        self.add_unknown_imported_symbol(
            imported_id, context, target_name=None, module_public=module_public,
            module_hidden=not module_public
        )

        if import_id == 'typing':
            # The user probably has a missing definition in a test fixture. Let's verify.
            fullname = 'builtins.{}'.format(source_id.lower())
            if (self.lookup_fully_qualified_or_none(fullname) is None and
                    fullname in SUGGESTED_TEST_FIXTURES):
                # Yes. Generate a helpful note.
                self.msg.add_fixture_note(fullname, context)

    def process_import_over_existing_name(self,
                                          imported_id: str, existing_symbol: SymbolTableNode,
                                          module_symbol: SymbolTableNode,
                                          import_node: ImportBase) -> bool:
        if existing_symbol.node is module_symbol.node:
            # We added this symbol on previous iteration.
            return False
        if (existing_symbol.kind in (LDEF, GDEF, MDEF) and
                isinstance(existing_symbol.node, (Var, FuncDef, TypeInfo, Decorator, TypeAlias))):
            # This is a valid import over an existing definition in the file. Construct a dummy
            # assignment that we'll use to type check the import.
            lvalue = NameExpr(imported_id)
            lvalue.kind = existing_symbol.kind
            lvalue.node = existing_symbol.node
            rvalue = NameExpr(imported_id)
            rvalue.kind = module_symbol.kind
            rvalue.node = module_symbol.node
            if isinstance(rvalue.node, TypeAlias):
                # Suppress bogus errors from the dummy assignment if rvalue is an alias.
                # Otherwise mypy may complain that alias is invalid in runtime context.
                rvalue.is_alias_rvalue = True
            assignment = AssignmentStmt([lvalue], rvalue)
            for node in assignment, lvalue, rvalue:
                node.set_line(import_node)
            import_node.assignments.append(assignment)
            return True
        return False

    def correct_relative_import(self, node: Union[ImportFrom, ImportAll]) -> str:
        import_id, ok = correct_relative_import(self.cur_mod_id, node.relative, node.id,
                                                self.cur_mod_node.is_package_init_file())
        if not ok:
            self.fail("Relative import climbs too many namespaces", node)
        return import_id

    def visit_import_all(self, i: ImportAll) -> None:
        i_id = self.correct_relative_import(i)
        if i_id in self.modules:
            m = self.modules[i_id]
            if self.is_incomplete_namespace(i_id):
                # Any names could be missing from the current namespace if the target module
                # namespace is incomplete.
                self.mark_incomplete('*', i)
            for name, node in m.names.items():
                fullname = i_id + '.' + name
                self.set_future_import_flags(fullname)
                if node is None:
                    continue
                # if '__all__' exists, all nodes not included have had module_public set to
                # False, and we can skip checking '_' because it's been explicitly included.
                if node.module_public and (not name.startswith('_') or '__all__' in m.names):
                    if isinstance(node.node, MypyFile):
                        # Star import of submodule from a package, add it as a dependency.
                        self.imports.add(node.node.fullname)
                    existing_symbol = self.lookup_current_scope(name)
                    if existing_symbol and not isinstance(node.node, PlaceholderNode):
                        # Import can redefine a variable. They get special treatment.
                        if self.process_import_over_existing_name(
                                name, existing_symbol, node, i):
                            continue
                    # In stub files, `from x import *` always reexports the symbols.
                    # In regular files, only if implicit reexports are enabled.
                    module_public = self.is_stub_file or self.options.implicit_reexport
                    self.add_imported_symbol(name, node, i,
                                             module_public=module_public,
                                             module_hidden=not module_public)

        else:
            # Don't add any dummy symbols for 'from x import *' if 'x' is unknown.
            pass

    #
    # Assignment
    #

    def visit_assignment_expr(self, s: AssignmentExpr) -> None:
        s.value.accept(self)
        self.analyze_lvalue(s.target, escape_comprehensions=True, has_explicit_value=True)

    def visit_assignment_stmt(self, s: AssignmentStmt) -> None:
        self.statement = s

        # Special case assignment like X = X.
        if self.analyze_identity_global_assignment(s):
            return

        tag = self.track_incomplete_refs()
        s.rvalue.accept(self)
        if self.found_incomplete_ref(tag) or self.should_wait_rhs(s.rvalue):
            # Initializer couldn't be fully analyzed. Defer the current node and give up.
            # Make sure that if we skip the definition of some local names, they can't be
            # added later in this scope, since an earlier definition should take precedence.
            for expr in names_modified_by_assignment(s):
                self.mark_incomplete(expr.name, expr)
            return

        # The r.h.s. is now ready to be classified, first check if it is a special form:
        special_form = False
        # * type alias
        if self.check_and_set_up_type_alias(s):
            s.is_alias_def = True
            special_form = True
        # * type variable definition
        elif self.process_typevar_declaration(s):
            special_form = True
        elif self.process_paramspec_declaration(s):
            special_form = True
        # * type constructors
        elif self.analyze_namedtuple_assign(s):
            special_form = True
        elif self.analyze_typeddict_assign(s):
            special_form = True
        elif self.newtype_analyzer.process_newtype_declaration(s):
            special_form = True
        elif self.analyze_enum_assign(s):
            special_form = True

        if special_form:
            self.record_special_form_lvalue(s)
            return
        # Clear the alias flag if assignment turns out not a special form after all. It
        # may be set to True while there were still placeholders due to forward refs.
        s.is_alias_def = False

        # OK, this is a regular assignment, perform the necessary analysis steps.
        s.is_final_def = self.unwrap_final(s)
        self.analyze_lvalues(s)
        self.check_final_implicit_def(s)
        self.check_classvar(s)
        self.process_type_annotation(s)
        self.apply_dynamic_class_hook(s)
        self.store_final_status(s)
        if not s.type:
            self.process_module_assignment(s.lvalues, s.rvalue, s)
        self.process__all__(s)
        self.process__deletable__(s)
        self.process__slots__(s)

    def analyze_identity_global_assignment(self, s: AssignmentStmt) -> bool:
        """Special case 'X = X' in global scope.

        This allows supporting some important use cases.

        Return true if special casing was applied.
        """
        if not isinstance(s.rvalue, NameExpr) or len(s.lvalues) != 1:
            # Not of form 'X = X'
            return False
        lvalue = s.lvalues[0]
        if not isinstance(lvalue, NameExpr) or s.rvalue.name != lvalue.name:
            # Not of form 'X = X'
            return False
        if self.type is not None or self.is_func_scope():
            # Not in global scope
            return False
        # It's an assignment like 'X = X' in the global scope.
        name = lvalue.name
        sym = self.lookup(name, s)
        if sym is None:
            if self.final_iteration:
                # Fall back to normal assignment analysis.
                return False
            else:
                self.defer()
                return True
        else:
            if sym.node is None:
                # Something special -- fall back to normal assignment analysis.
                return False
            if name not in self.globals:
                # The name is from builtins. Add an alias to the current module.
                self.add_symbol(name, sym.node, s)
            if not isinstance(sym.node, PlaceholderNode):
                for node in s.rvalue, lvalue:
                    node.node = sym.node
                    node.kind = GDEF
                    node.fullname = sym.node.fullname
            return True

    def should_wait_rhs(self, rv: Expression) -> bool:
        """Can we already classify this r.h.s. of an assignment or should we wait?

        This returns True if we don't have enough information to decide whether
        an assignment is just a normal variable definition or a special form.
        Always return False if this is a final iteration. This will typically cause
        the lvalue to be classified as a variable plus emit an error.
        """
        if self.final_iteration:
            # No chance, nothing has changed.
            return False
        if isinstance(rv, NameExpr):
            n = self.lookup(rv.name, rv)
            if n and isinstance(n.node, PlaceholderNode) and not n.node.becomes_typeinfo:
                return True
        elif isinstance(rv, MemberExpr):
            fname = get_member_expr_fullname(rv)
            if fname:
                n = self.lookup_qualified(fname, rv, suppress_errors=True)
                if n and isinstance(n.node, PlaceholderNode) and not n.node.becomes_typeinfo:
                    return True
        elif isinstance(rv, IndexExpr) and isinstance(rv.base, RefExpr):
            return self.should_wait_rhs(rv.base)
        elif isinstance(rv, CallExpr) and isinstance(rv.callee, RefExpr):
            # This is only relevant for builtin SCC where things like 'TypeVar'
            # may be not ready.
            return self.should_wait_rhs(rv.callee)
        return False

    def can_be_type_alias(self, rv: Expression, allow_none: bool = False) -> bool:
        """Is this a valid r.h.s. for an alias definition?

        Note: this function should be only called for expressions where self.should_wait_rhs()
        returns False.
        """
        if isinstance(rv, RefExpr) and self.is_type_ref(rv, bare=True):
            return True
        if isinstance(rv, IndexExpr) and self.is_type_ref(rv.base, bare=False):
            return True
        if self.is_none_alias(rv):
            return True
        if allow_none and isinstance(rv, NameExpr) and rv.fullname == 'builtins.None':
            return True
        if (isinstance(rv, OpExpr)
                and rv.op == '|'
                and self.can_be_type_alias(rv.left, allow_none=True)
                and self.can_be_type_alias(rv.right, allow_none=True)):
            return True
        return False

    def is_type_ref(self, rv: Expression, bare: bool = False) -> bool:
        """Does this expression refer to a type?

        This includes:
          * Special forms, like Any or Union
          * Classes (except subscripted enums)
          * Other type aliases
          * PlaceholderNodes with becomes_typeinfo=True (these can be not ready class
            definitions, and not ready aliases).

        If bare is True, this is not a base of an index expression, so some special
        forms are not valid (like a bare Union).

        Note: This method should be only used in context of a type alias definition.
        This method can only return True for RefExprs, to check if C[int] is a valid
        target for type alias call this method on expr.base (i.e. on C in C[int]).
        See also can_be_type_alias().
        """
        if not isinstance(rv, RefExpr):
            return False
        if isinstance(rv.node, TypeVarExpr):
            self.fail('Type variable "{}" is invalid as target for type alias'.format(
                rv.fullname), rv)
            return False

        if bare:
            # These three are valid even if bare, for example
            # A = Tuple is just equivalent to A = Tuple[Any, ...].
            valid_refs = {'typing.Any', 'typing.Tuple', 'typing.Callable'}
        else:
            valid_refs = type_constructors

        if isinstance(rv.node, TypeAlias) or rv.fullname in valid_refs:
            return True
        if isinstance(rv.node, TypeInfo):
            if bare:
                return True
            # Assignment color = Color['RED'] defines a variable, not an alias.
            return not rv.node.is_enum
        if isinstance(rv.node, Var):
            return rv.node.fullname in (
                'typing.NoReturn',
                'typing_extensions.NoReturn',
                'mypy_extensions.NoReturn',
            )

        if isinstance(rv, NameExpr):
            n = self.lookup(rv.name, rv)
            if n and isinstance(n.node, PlaceholderNode) and n.node.becomes_typeinfo:
                return True
        elif isinstance(rv, MemberExpr):
            fname = get_member_expr_fullname(rv)
            if fname:
                # The r.h.s. for variable definitions may not be a type reference but just
                # an instance attribute, so suppress the errors.
                n = self.lookup_qualified(fname, rv, suppress_errors=True)
                if n and isinstance(n.node, PlaceholderNode) and n.node.becomes_typeinfo:
                    return True
        return False

    def is_none_alias(self, node: Expression) -> bool:
        """Is this a r.h.s. for a None alias?

        We special case the assignments like Void = type(None), to allow using
        Void in type annotations.
        """
        if isinstance(node, CallExpr):
            if (isinstance(node.callee, NameExpr) and len(node.args) == 1 and
                    isinstance(node.args[0], NameExpr)):
                call = self.lookup_qualified(node.callee.name, node.callee)
                arg = self.lookup_qualified(node.args[0].name, node.args[0])
                if (call is not None and call.node and call.node.fullname == 'builtins.type' and
                        arg is not None and arg.node and arg.node.fullname == 'builtins.None'):
                    return True
        return False

    def record_special_form_lvalue(self, s: AssignmentStmt) -> None:
        """Record minimal necessary information about l.h.s. of a special form.

        This exists mostly for compatibility with the old semantic analyzer.
        """
        lvalue = s.lvalues[0]
        assert isinstance(lvalue, NameExpr)
        lvalue.is_special_form = True
        if self.current_symbol_kind() == GDEF:
            lvalue.fullname = self.qualified_name(lvalue.name)
        lvalue.kind = self.current_symbol_kind()

    def analyze_enum_assign(self, s: AssignmentStmt) -> bool:
        """Check if s defines an Enum."""
        if isinstance(s.rvalue, CallExpr) and isinstance(s.rvalue.analyzed, EnumCallExpr):
            # Already analyzed enum -- nothing to do here.
            return True
        return self.enum_call_analyzer.process_enum_call(s, self.is_func_scope())

    def analyze_namedtuple_assign(self, s: AssignmentStmt) -> bool:
        """Check if s defines a namedtuple."""
        if isinstance(s.rvalue, CallExpr) and isinstance(s.rvalue.analyzed, NamedTupleExpr):
            return True  # This is a valid and analyzed named tuple definition, nothing to do here.
        if len(s.lvalues) != 1 or not isinstance(s.lvalues[0], (NameExpr, MemberExpr)):
            return False
        lvalue = s.lvalues[0]
        name = lvalue.name
        internal_name, info = self.named_tuple_analyzer.check_namedtuple(s.rvalue, name,
                                                                         self.is_func_scope())
        if internal_name is None:
            return False
        if isinstance(lvalue, MemberExpr):
            self.fail("NamedTuple type as an attribute is not supported", lvalue)
            return False
        if internal_name != name:
            self.fail('First argument to namedtuple() should be "{}", not "{}"'.format(
                name, internal_name), s.rvalue, code=codes.NAME_MATCH)
            return True
        # Yes, it's a valid namedtuple, but defer if it is not ready.
        if not info:
            self.mark_incomplete(name, lvalue, becomes_typeinfo=True)
        return True

    def analyze_typeddict_assign(self, s: AssignmentStmt) -> bool:
        """Check if s defines a typed dict."""
        if isinstance(s.rvalue, CallExpr) and isinstance(s.rvalue.analyzed, TypedDictExpr):
            return True  # This is a valid and analyzed typed dict definition, nothing to do here.
        if len(s.lvalues) != 1 or not isinstance(s.lvalues[0], (NameExpr, MemberExpr)):
            return False
        lvalue = s.lvalues[0]
        name = lvalue.name
        is_typed_dict, info = self.typed_dict_analyzer.check_typeddict(s.rvalue, name,
                                                                       self.is_func_scope())
        if not is_typed_dict:
            return False
        if isinstance(lvalue, MemberExpr):
            self.fail("TypedDict type as attribute is not supported", lvalue)
            return False
        # Yes, it's a valid typed dict, but defer if it is not ready.
        if not info:
            self.mark_incomplete(name, lvalue, becomes_typeinfo=True)
        return True

    def analyze_lvalues(self, s: AssignmentStmt) -> None:
        # We cannot use s.type, because analyze_simple_literal_type() will set it.
        explicit = s.unanalyzed_type is not None
        if self.is_final_type(s.unanalyzed_type):
            # We need to exclude bare Final.
            assert isinstance(s.unanalyzed_type, UnboundType)
            if not s.unanalyzed_type.args:
                explicit = False

        if s.rvalue:
            if isinstance(s.rvalue, TempNode):
                has_explicit_value = not s.rvalue.no_rhs
            else:
                has_explicit_value = True
        else:
            has_explicit_value = False

        for lval in s.lvalues:
            self.analyze_lvalue(lval,
                                explicit_type=explicit,
                                is_final=s.is_final_def,
                                has_explicit_value=has_explicit_value)

    def apply_dynamic_class_hook(self, s: AssignmentStmt) -> None:
        if not isinstance(s.rvalue, CallExpr):
            return
        fname = None
        call = s.rvalue
        while True:
            if isinstance(call.callee, RefExpr):
                fname = call.callee.fullname
            # check if method call
            if fname is None and isinstance(call.callee, MemberExpr):
                callee_expr = call.callee.expr
                if isinstance(callee_expr, RefExpr) and callee_expr.fullname:
                    method_name = call.callee.name
                    fname = callee_expr.fullname + '.' + method_name
                elif isinstance(callee_expr, CallExpr):
                    # check if chain call
                    call = callee_expr
                    continue
            break
        if not fname:
            return
        hook = self.plugin.get_dynamic_class_hook(fname)
        if not hook:
            return
        for lval in s.lvalues:
            if not isinstance(lval, NameExpr):
                continue
            hook(DynamicClassDefContext(call, lval.name, self))

    def unwrap_final(self, s: AssignmentStmt) -> bool:
        """Strip Final[...] if present in an assignment.

        This is done to invoke type inference during type checking phase for this
        assignment. Also, Final[...] doesn't affect type in any way -- it is rather an
        access qualifier for given `Var`.

        Also perform various consistency checks.

        Returns True if Final[...] was present.
        """
        if not s.unanalyzed_type or not self.is_final_type(s.unanalyzed_type):
            return False
        assert isinstance(s.unanalyzed_type, UnboundType)
        if len(s.unanalyzed_type.args) > 1:
            self.fail("Final[...] takes at most one type argument", s.unanalyzed_type)
        invalid_bare_final = False
        if not s.unanalyzed_type.args:
            s.type = None
            if isinstance(s.rvalue, TempNode) and s.rvalue.no_rhs:
                invalid_bare_final = True
                self.fail("Type in Final[...] can only be omitted if there is an initializer", s)
        else:
            s.type = s.unanalyzed_type.args[0]

        if s.type is not None and self.is_classvar(s.type):
            self.fail("Variable should not be annotated with both ClassVar and Final", s)
            return False

        if len(s.lvalues) != 1 or not isinstance(s.lvalues[0], RefExpr):
            self.fail("Invalid final declaration", s)
            return False
        lval = s.lvalues[0]
        assert isinstance(lval, RefExpr)

        # Reset inferred status if it was set due to simple literal rvalue on previous iteration.
        # TODO: this is a best-effort quick fix, we should avoid the need to manually sync this,
        # see https://github.com/python/mypy/issues/6458.
        if lval.is_new_def:
            lval.is_inferred_def = s.type is None

        if self.loop_depth > 0:
            self.fail("Cannot use Final inside a loop", s)
        if self.type and self.type.is_protocol:
            self.msg.protocol_members_cant_be_final(s)
        if (isinstance(s.rvalue, TempNode) and s.rvalue.no_rhs and
                not self.is_stub_file and not self.is_class_scope()):
            if not invalid_bare_final:  # Skip extra error messages.
                self.msg.final_without_value(s)
        return True

    def check_final_implicit_def(self, s: AssignmentStmt) -> None:
        """Do basic checks for final declaration on self in __init__.

        Additional re-definition checks are performed by `analyze_lvalue`.
        """
        if not s.is_final_def:
            return
        lval = s.lvalues[0]
        assert isinstance(lval, RefExpr)
        if isinstance(lval, MemberExpr):
            if not self.is_self_member_ref(lval):
                self.fail("Final can be only applied to a name or an attribute on self", s)
                s.is_final_def = False
                return
            else:
                assert self.function_stack
                if self.function_stack[-1].name != '__init__':
                    self.fail("Can only declare a final attribute in class body or __init__", s)
                    s.is_final_def = False
                    return

    def store_final_status(self, s: AssignmentStmt) -> None:
        """If this is a locally valid final declaration, set the corresponding flag on `Var`."""
        if s.is_final_def:
            if len(s.lvalues) == 1 and isinstance(s.lvalues[0], RefExpr):
                node = s.lvalues[0].node
                if isinstance(node, Var):
                    node.is_final = True
                    node.final_value = self.unbox_literal(s.rvalue)
                    if (self.is_class_scope() and
                            (isinstance(s.rvalue, TempNode) and s.rvalue.no_rhs)):
                        node.final_unset_in_class = True
        else:
            for lval in self.flatten_lvalues(s.lvalues):
                # Special case: we are working with an `Enum`:
                #
                #   class MyEnum(Enum):
                #       key = 'some value'
                #
                # Here `key` is implicitly final. In runtime, code like
                #
                #     MyEnum.key = 'modified'
                #
                # will fail with `AttributeError: Cannot reassign members.`
                # That's why we need to replicate this.
                if (isinstance(lval, NameExpr) and
                        isinstance(self.type, TypeInfo) and
                        self.type.is_enum):
                    cur_node = self.type.names.get(lval.name, None)
                    if (cur_node and isinstance(cur_node.node, Var) and
                            not (isinstance(s.rvalue, TempNode) and s.rvalue.no_rhs)):
                        cur_node.node.is_final = True
                        s.is_final_def = True

                # Special case: deferred initialization of a final attribute in __init__.
                # In this case we just pretend this is a valid final definition to suppress
                # errors about assigning to final attribute.
                if isinstance(lval, MemberExpr) and self.is_self_member_ref(lval):
                    assert self.type, "Self member outside a class"
                    cur_node = self.type.names.get(lval.name, None)
                    if cur_node and isinstance(cur_node.node, Var) and cur_node.node.is_final:
                        assert self.function_stack
                        top_function = self.function_stack[-1]
                        if (top_function.name == '__init__' and
                                cur_node.node.final_unset_in_class and
                                not cur_node.node.final_set_in_init and
                                not (isinstance(s.rvalue, TempNode) and s.rvalue.no_rhs)):
                            cur_node.node.final_set_in_init = True
                            s.is_final_def = True

    def flatten_lvalues(self, lvalues: List[Expression]) -> List[Expression]:
        res: List[Expression] = []
        for lv in lvalues:
            if isinstance(lv, (TupleExpr, ListExpr)):
                res.extend(self.flatten_lvalues(lv.items))
            else:
                res.append(lv)
        return res

    def unbox_literal(self, e: Expression) -> Optional[Union[int, float, bool, str]]:
        if isinstance(e, (IntExpr, FloatExpr, StrExpr)):
            return e.value
        elif isinstance(e, NameExpr) and e.name in ('True', 'False'):
            return True if e.name == 'True' else False
        return None

    def process_type_annotation(self, s: AssignmentStmt) -> None:
        """Analyze type annotation or infer simple literal type."""
        if s.type:
            lvalue = s.lvalues[-1]
            allow_tuple_literal = isinstance(lvalue, TupleExpr)
            analyzed = self.anal_type(s.type, allow_tuple_literal=allow_tuple_literal)
            # Don't store not ready types (including placeholders).
            if analyzed is None or has_placeholder(analyzed):
                return
            s.type = analyzed
            if (self.type and self.type.is_protocol and isinstance(lvalue, NameExpr) and
                    isinstance(s.rvalue, TempNode) and s.rvalue.no_rhs):
                if isinstance(lvalue.node, Var):
                    lvalue.node.is_abstract_var = True
        else:
            if (self.type and self.type.is_protocol and
                    self.is_annotated_protocol_member(s) and not self.is_func_scope()):
                self.fail('All protocol members must have explicitly declared types', s)
            # Set the type if the rvalue is a simple literal (even if the above error occurred).
            if len(s.lvalues) == 1 and isinstance(s.lvalues[0], RefExpr):
                if s.lvalues[0].is_inferred_def:
                    s.type = self.analyze_simple_literal_type(s.rvalue, s.is_final_def)
        if s.type:
            # Store type into nodes.
            for lvalue in s.lvalues:
                self.store_declared_types(lvalue, s.type)

    def is_annotated_protocol_member(self, s: AssignmentStmt) -> bool:
        """Check whether a protocol member is annotated.

        There are some exceptions that can be left unannotated, like ``__slots__``."""
        return any(
            (
                isinstance(lv, NameExpr)
                and lv.name != '__slots__'
                and lv.is_inferred_def
            )
            for lv in s.lvalues
        )

    def analyze_simple_literal_type(self, rvalue: Expression, is_final: bool) -> Optional[Type]:
        """Return builtins.int if rvalue is an int literal, etc.

        If this is a 'Final' context, we return "Literal[...]" instead."""
        if self.options.semantic_analysis_only or self.function_stack:
            # Skip this if we're only doing the semantic analysis pass.
            # This is mostly to avoid breaking unit tests.
            # Also skip inside a function; this is to avoid confusing
            # the code that handles dead code due to isinstance()
            # inside type variables with value restrictions (like
            # AnyStr).
            return None
        if isinstance(rvalue, FloatExpr):
            return self.named_type_or_none('builtins.float')

        value: Optional[LiteralValue] = None
        type_name: Optional[str] = None
        if isinstance(rvalue, IntExpr):
            value, type_name = rvalue.value, 'builtins.int'
        if isinstance(rvalue, StrExpr):
            value, type_name = rvalue.value, 'builtins.str'
        if isinstance(rvalue, BytesExpr):
            value, type_name = rvalue.value, 'builtins.bytes'
        if isinstance(rvalue, UnicodeExpr):
            value, type_name = rvalue.value, 'builtins.unicode'

        if type_name is not None:
            assert value is not None
            typ = self.named_type_or_none(type_name)
            if typ and is_final:
                return typ.copy_modified(last_known_value=LiteralType(
                    value=value,
                    fallback=typ,
                    line=typ.line,
                    column=typ.column,
                ))
            return typ

        return None

    def analyze_alias(self, rvalue: Expression,
                      allow_placeholder: bool = False) -> Tuple[Optional[Type], List[str],
                                                                Set[str], List[str]]:
        """Check if 'rvalue' is a valid type allowed for aliasing (e.g. not a type variable).

        If yes, return the corresponding type, a list of
        qualified type variable names for generic aliases, a set of names the alias depends on,
        and a list of type variables if the alias is generic.
        An schematic example for the dependencies:
            A = int
            B = str
            analyze_alias(Dict[A, B])[2] == {'__main__.A', '__main__.B'}
        """
        dynamic = bool(self.function_stack and self.function_stack[-1].is_dynamic())
        global_scope = not self.type and not self.function_stack
        res = analyze_type_alias(rvalue,
                                 self,
                                 self.tvar_scope,
                                 self.plugin,
                                 self.options,
                                 self.is_typeshed_stub_file,
                                 allow_placeholder=allow_placeholder,
                                 in_dynamic_func=dynamic,
                                 global_scope=global_scope)
        typ: Optional[Type] = None
        if res:
            typ, depends_on = res
            found_type_vars = typ.accept(TypeVarLikeQuery(self.lookup_qualified, self.tvar_scope))
            alias_tvars = [name for (name, node) in found_type_vars]
            qualified_tvars = [node.fullname for (name, node) in found_type_vars]
        else:
            alias_tvars = []
            depends_on = set()
            qualified_tvars = []
        return typ, alias_tvars, depends_on, qualified_tvars

    def check_and_set_up_type_alias(self, s: AssignmentStmt) -> bool:
        """Check if assignment creates a type alias and set it up as needed.

        Return True if it is a type alias (even if the target is not ready),
        or False otherwise.

        Note: the resulting types for subscripted (including generic) aliases
        are also stored in rvalue.analyzed.
        """
        lvalue = s.lvalues[0]
        if len(s.lvalues) > 1 or not isinstance(lvalue, NameExpr):
            # First rule: Only simple assignments like Alias = ... create aliases.
            return False

        pep_613 = False
        if s.unanalyzed_type is not None and isinstance(s.unanalyzed_type, UnboundType):
            lookup = self.lookup(s.unanalyzed_type.name, s, suppress_errors=True)
            if lookup and lookup.fullname in ("typing.TypeAlias", "typing_extensions.TypeAlias"):
                pep_613 = True
        if s.unanalyzed_type is not None and not pep_613:
            # Second rule: Explicit type (cls: Type[A] = A) always creates variable, not alias.
            # unless using PEP 613 `cls: TypeAlias = A`
            return False

        existing = self.current_symbol_table().get(lvalue.name)
        # Third rule: type aliases can't be re-defined. For example:
        #     A: Type[float] = int
        #     A = float  # OK, but this doesn't define an alias
        #     B = int
        #     B = float  # Error!
        # Don't create an alias in these cases:
        if (existing
                and (isinstance(existing.node, Var)  # existing variable
                     or (isinstance(existing.node, TypeAlias)
                         and not s.is_alias_def)  # existing alias
                     or (isinstance(existing.node, PlaceholderNode)
                         and existing.node.node.line < s.line))):  # previous incomplete definition
            # TODO: find a more robust way to track the order of definitions.
            # Note: if is_alias_def=True, this is just a node from previous iteration.
            if isinstance(existing.node, TypeAlias) and not s.is_alias_def:
                self.fail('Cannot assign multiple types to name "{}"'
                          ' without an explicit "Type[...]" annotation'
                          .format(lvalue.name), lvalue)
            return False

        non_global_scope = self.type or self.is_func_scope()
        if isinstance(s.rvalue, RefExpr) and non_global_scope and not pep_613:
            # Fourth rule (special case): Non-subscripted right hand side creates a variable
            # at class and function scopes. For example:
            #
            #   class Model:
            #       ...
            #   class C:
            #       model = Model # this is automatically a variable with type 'Type[Model]'
            #
            # without this rule, this typical use case will require a lot of explicit
            # annotations (see the second rule).
            return False
        rvalue = s.rvalue
        if not self.can_be_type_alias(rvalue) and not pep_613:
            return False

        if existing and not isinstance(existing.node, (PlaceholderNode, TypeAlias)):
            # Cannot redefine existing node as type alias.
            return False

        res: Optional[Type] = None
        if self.is_none_alias(rvalue):
            res = NoneType()
            alias_tvars, depends_on, qualified_tvars = \
                [], set(), []  # type: List[str], Set[str], List[str]
        else:
            tag = self.track_incomplete_refs()
            res, alias_tvars, depends_on, qualified_tvars = \
                self.analyze_alias(rvalue, allow_placeholder=True)
            if not res:
                return False
            # TODO: Maybe we only need to reject top-level placeholders, similar
            #       to base classes.
            if self.found_incomplete_ref(tag) or has_placeholder(res):
                # Since we have got here, we know this must be a type alias (incomplete refs
                # may appear in nested positions), therefore use becomes_typeinfo=True.
                self.mark_incomplete(lvalue.name, rvalue, becomes_typeinfo=True)
                return True
        self.add_type_alias_deps(depends_on)
        # In addition to the aliases used, we add deps on unbound
        # type variables, since they are erased from target type.
        self.add_type_alias_deps(qualified_tvars)
        # The above are only direct deps on other aliases.
        # For subscripted aliases, type deps from expansion are added in deps.py
        # (because the type is stored).
        check_for_explicit_any(res, self.options, self.is_typeshed_stub_file, self.msg,
                               context=s)
        # When this type alias gets "inlined", the Any is not explicit anymore,
        # so we need to replace it with non-explicit Anys.
        if not has_placeholder(res):
            res = make_any_non_explicit(res)
        # Note: with the new (lazy) type alias representation we only need to set no_args to True
        # if the expected number of arguments is non-zero, so that aliases like A = List work.
        # However, eagerly expanding aliases like Text = str is a nice performance optimization.
        no_args = isinstance(res, Instance) and not res.args  # type: ignore[misc]
        fix_instance_types(res, self.fail, self.note, self.options.python_version)
        # Aliases defined within functions can't be accessed outside
        # the function, since the symbol table will no longer
        # exist. Work around by expanding them eagerly when used.
        eager = self.is_func_scope()
        alias_node = TypeAlias(res,
                               self.qualified_name(lvalue.name),
                               s.line,
                               s.column,
                               alias_tvars=alias_tvars,
                               no_args=no_args,
                               eager=eager)
        if isinstance(s.rvalue, (IndexExpr, CallExpr)):  # CallExpr is for `void = type(None)`
            s.rvalue.analyzed = TypeAliasExpr(alias_node)
            s.rvalue.analyzed.line = s.line
            # we use the column from resulting target, to get better location for errors
            s.rvalue.analyzed.column = res.column
        elif isinstance(s.rvalue, RefExpr):
            s.rvalue.is_alias_rvalue = True

        if existing:
            # An alias gets updated.
            updated = False
            if isinstance(existing.node, TypeAlias):
                if existing.node.target != res:
                    # Copy expansion to the existing alias, this matches how we update base classes
                    # for a TypeInfo _in place_ if there are nested placeholders.
                    existing.node.target = res
                    existing.node.alias_tvars = alias_tvars
                    existing.node.no_args = no_args
                    updated = True
            else:
                # Otherwise just replace existing placeholder with type alias.
                existing.node = alias_node
                updated = True
            if updated:
                if self.final_iteration:
                    self.cannot_resolve_name(lvalue.name, 'name', s)
                    return True
                else:
                    self.progress = True
                    # We need to defer so that this change can get propagated to base classes.
                    self.defer(s)
        else:
            self.add_symbol(lvalue.name, alias_node, s)
        if isinstance(rvalue, RefExpr) and isinstance(rvalue.node, TypeAlias):
            alias_node.normalized = rvalue.node.normalized
        return True

    def analyze_lvalue(self,
                       lval: Lvalue,
                       nested: bool = False,
                       explicit_type: bool = False,
                       is_final: bool = False,
                       escape_comprehensions: bool = False,
                       has_explicit_value: bool = False) -> None:
        """Analyze an lvalue or assignment target.

        Args:
            lval: The target lvalue
            nested: If true, the lvalue is within a tuple or list lvalue expression
            explicit_type: Assignment has type annotation
            escape_comprehensions: If we are inside a comprehension, set the variable
                in the enclosing scope instead. This implements
                https://www.python.org/dev/peps/pep-0572/#scope-of-the-target
        """
        if escape_comprehensions:
            assert isinstance(lval, NameExpr), "assignment expression target must be NameExpr"
        if isinstance(lval, NameExpr):
            self.analyze_name_lvalue(
                lval, explicit_type, is_final,
                escape_comprehensions,
                has_explicit_value=has_explicit_value,
            )
        elif isinstance(lval, MemberExpr):
            self.analyze_member_lvalue(lval, explicit_type, is_final)
            if explicit_type and not self.is_self_member_ref(lval):
                self.fail('Type cannot be declared in assignment to non-self '
                          'attribute', lval)
        elif isinstance(lval, IndexExpr):
            if explicit_type:
                self.fail('Unexpected type declaration', lval)
            lval.accept(self)
        elif isinstance(lval, TupleExpr):
            self.analyze_tuple_or_list_lvalue(lval, explicit_type)
        elif isinstance(lval, StarExpr):
            if nested:
                self.analyze_lvalue(lval.expr, nested, explicit_type)
            else:
                self.fail('Starred assignment target must be in a list or tuple', lval)
        else:
            self.fail('Invalid assignment target', lval)

    def analyze_name_lvalue(self,
                            lvalue: NameExpr,
                            explicit_type: bool,
                            is_final: bool,
                            escape_comprehensions: bool,
                            has_explicit_value: bool) -> None:
        """Analyze an lvalue that targets a name expression.

        Arguments are similar to "analyze_lvalue".
        """
        if lvalue.node:
            # This has been bound already in a previous iteration.
            return

        name = lvalue.name
        if self.is_alias_for_final_name(name):
            if is_final:
                self.fail("Cannot redefine an existing name as final", lvalue)
            else:
                self.msg.cant_assign_to_final(name, self.type is not None, lvalue)

        kind = self.current_symbol_kind()
        names = self.current_symbol_table(escape_comprehensions=escape_comprehensions)
        existing = names.get(name)

        outer = self.is_global_or_nonlocal(name)
        if kind == MDEF and isinstance(self.type, TypeInfo) and self.type.is_enum:
            # Special case: we need to be sure that `Enum` keys are unique.
            if existing:
                self.fail('Attempted to reuse member name "{}" in Enum definition "{}"'.format(
                    name, self.type.name,
                ), lvalue)

        if (not existing or isinstance(existing.node, PlaceholderNode)) and not outer:
            # Define new variable.
            var = self.make_name_lvalue_var(lvalue, kind, not explicit_type, has_explicit_value)
            added = self.add_symbol(name, var, lvalue, escape_comprehensions=escape_comprehensions)
            # Only bind expression if we successfully added name to symbol table.
            if added:
                lvalue.is_new_def = True
                lvalue.is_inferred_def = True
                lvalue.kind = kind
                lvalue.node = var
                if kind == GDEF:
                    lvalue.fullname = var._fullname
                else:
                    lvalue.fullname = lvalue.name
                if self.is_func_scope():
                    if unmangle(name) == '_':
                        # Special case for assignment to local named '_': always infer 'Any'.
                        typ = AnyType(TypeOfAny.special_form)
                        self.store_declared_types(lvalue, typ)
            if is_final and self.is_final_redefinition(kind, name):
                self.fail("Cannot redefine an existing name as final", lvalue)
        else:
            self.make_name_lvalue_point_to_existing_def(lvalue, explicit_type, is_final)

    def is_final_redefinition(self, kind: int, name: str) -> bool:
        if kind == GDEF:
            return self.is_mangled_global(name) and not self.is_initial_mangled_global(name)
        elif kind == MDEF and self.type:
            return unmangle(name) + "'" in self.type.names
        return False

    def is_alias_for_final_name(self, name: str) -> bool:
        if self.is_func_scope():
            if not name.endswith("'"):
                # Not a mangled name -- can't be an alias
                return False
            name = unmangle(name)
            assert self.locals[-1] is not None, "No locals at function scope"
            existing = self.locals[-1].get(name)
            return existing is not None and is_final_node(existing.node)
        elif self.type is not None:
            orig_name = unmangle(name) + "'"
            if name == orig_name:
                return False
            existing = self.type.names.get(orig_name)
            return existing is not None and is_final_node(existing.node)
        else:
            orig_name = unmangle(name) + "'"
            if name == orig_name:
                return False
            existing = self.globals.get(orig_name)
            return existing is not None and is_final_node(existing.node)

    def make_name_lvalue_var(
        self, lvalue: NameExpr, kind: int, inferred: bool, has_explicit_value: bool,
    ) -> Var:
        """Return a Var node for an lvalue that is a name expression."""
        v = Var(lvalue.name)
        v.set_line(lvalue)
        v.is_inferred = inferred
        if kind == MDEF:
            assert self.type is not None
            v.info = self.type
            v.is_initialized_in_class = True
        if kind != LDEF:
            v._fullname = self.qualified_name(lvalue.name)
        else:
            # fullanme should never stay None
            v._fullname = lvalue.name
        v.is_ready = False  # Type not inferred yet
        v.has_explicit_value = has_explicit_value
        return v

    def make_name_lvalue_point_to_existing_def(
            self,
            lval: NameExpr,
            explicit_type: bool,
            is_final: bool) -> None:
        """Update an lvalue to point to existing definition in the same scope.

        Arguments are similar to "analyze_lvalue".

        Assume that an existing name exists.
        """
        if is_final:
            # Redefining an existing name with final is always an error.
            self.fail("Cannot redefine an existing name as final", lval)
        original_def = self.lookup(lval.name, lval, suppress_errors=True)
        if original_def is None and self.type and not self.is_func_scope():
            # Workaround to allow "x, x = ..." in class body.
            original_def = self.type.get(lval.name)
        if explicit_type:
            # Don't re-bind if there is a type annotation.
            self.name_already_defined(lval.name, lval, original_def)
        else:
            # Bind to an existing name.
            if original_def:
                self.bind_name_expr(lval, original_def)
            else:
                self.name_not_defined(lval.name, lval)
            self.check_lvalue_validity(lval.node, lval)

    def analyze_tuple_or_list_lvalue(self, lval: TupleExpr,
                                     explicit_type: bool = False) -> None:
        """Analyze an lvalue or assignment target that is a list or tuple."""
        items = lval.items
        star_exprs = [item for item in items if isinstance(item, StarExpr)]

        if len(star_exprs) > 1:
            self.fail('Two starred expressions in assignment', lval)
        else:
            if len(star_exprs) == 1:
                star_exprs[0].valid = True
            for i in items:
                self.analyze_lvalue(
                    lval=i,
                    nested=True,
                    explicit_type=explicit_type,
                    # Lists and tuples always have explicit values defined:
                    # `a, b, c = value`
                    has_explicit_value=True,
                )

    def analyze_member_lvalue(self, lval: MemberExpr, explicit_type: bool, is_final: bool) -> None:
        """Analyze lvalue that is a member expression.

        Arguments:
            lval: The target lvalue
            explicit_type: Assignment has type annotation
            is_final: Is the target final
        """
        if lval.node:
            # This has been bound already in a previous iteration.
            return
        lval.accept(self)
        if self.is_self_member_ref(lval):
            assert self.type, "Self member outside a class"
            cur_node = self.type.names.get(lval.name)
            node = self.type.get(lval.name)
            if cur_node and is_final:
                # Overrides will be checked in type checker.
                self.fail("Cannot redefine an existing name as final", lval)
            # On first encounter with this definition, if this attribute was defined before
            # with an inferred type and it's marked with an explicit type now, give an error.
            if (not lval.node and cur_node and isinstance(cur_node.node, Var) and
                    cur_node.node.is_inferred and explicit_type):
                self.attribute_already_defined(lval.name, lval, cur_node)
            # If the attribute of self is not defined in superclasses, create a new Var, ...
            if (node is None
                    or (isinstance(node.node, Var) and node.node.is_abstract_var)
                    # ... also an explicit declaration on self also creates a new Var.
                    # Note that `explicit_type` might has been erased for bare `Final`,
                    # so we also check if `is_final` is passed.
                    or (cur_node is None and (explicit_type or is_final))):
                if self.type.is_protocol and node is None:
                    self.fail("Protocol members cannot be defined via assignment to self", lval)
                else:
                    # Implicit attribute definition in __init__.
                    lval.is_new_def = True
                    lval.is_inferred_def = True
                    v = Var(lval.name)
                    v.set_line(lval)
                    v._fullname = self.qualified_name(lval.name)
                    v.info = self.type
                    v.is_ready = False
                    v.explicit_self_type = explicit_type or is_final
                    lval.def_var = v
                    lval.node = v
                    # TODO: should we also set lval.kind = MDEF?
                    self.type.names[lval.name] = SymbolTableNode(MDEF, v, implicit=True)
        self.check_lvalue_validity(lval.node, lval)

    def is_self_member_ref(self, memberexpr: MemberExpr) -> bool:
        """Does memberexpr to refer to an attribute of self?"""
        if not isinstance(memberexpr.expr, NameExpr):
            return False
        node = memberexpr.expr.node
        return isinstance(node, Var) and node.is_self

    def check_lvalue_validity(self, node: Union[Expression, SymbolNode, None],
                              ctx: Context) -> None:
        if isinstance(node, TypeVarExpr):
            self.fail('Invalid assignment target', ctx)
        elif isinstance(node, TypeInfo):
            self.fail(message_registry.CANNOT_ASSIGN_TO_TYPE, ctx)

    def store_declared_types(self, lvalue: Lvalue, typ: Type) -> None:
        if isinstance(typ, StarType) and not isinstance(lvalue, StarExpr):
            self.fail('Star type only allowed for starred expressions', lvalue)
        if isinstance(lvalue, RefExpr):
            lvalue.is_inferred_def = False
            if isinstance(lvalue.node, Var):
                var = lvalue.node
                var.type = typ
                var.is_ready = True
            # If node is not a variable, we'll catch it elsewhere.
        elif isinstance(lvalue, TupleExpr):
            typ = get_proper_type(typ)
            if isinstance(typ, TupleType):
                if len(lvalue.items) != len(typ.items):
                    self.fail('Incompatible number of tuple items', lvalue)
                    return
                for item, itemtype in zip(lvalue.items, typ.items):
                    self.store_declared_types(item, itemtype)
            else:
                self.fail('Tuple type expected for multiple variables',
                          lvalue)
        elif isinstance(lvalue, StarExpr):
            # Historical behavior for the old parser
            if isinstance(typ, StarType):
                self.store_declared_types(lvalue.expr, typ.type)
            else:
                self.store_declared_types(lvalue.expr, typ)
        else:
            # This has been flagged elsewhere as an error, so just ignore here.
            pass

    def process_typevar_declaration(self, s: AssignmentStmt) -> bool:
        """Check if s declares a TypeVar; it yes, store it in symbol table.

        Return True if this looks like a type variable declaration (but maybe
        with errors), otherwise return False.
        """
        call = self.get_typevarlike_declaration(s, ("typing.TypeVar",))
        if not call:
            return False

        lvalue = s.lvalues[0]
        assert isinstance(lvalue, NameExpr)
        if s.type:
            self.fail("Cannot declare the type of a type variable", s)
            return False

        name = lvalue.name
        if not self.check_typevarlike_name(call, name, s):
            return False

        # Constraining types
        n_values = call.arg_kinds[1:].count(ARG_POS)
        values = self.analyze_value_types(call.args[1:1 + n_values])

        res = self.process_typevar_parameters(call.args[1 + n_values:],
                                              call.arg_names[1 + n_values:],
                                              call.arg_kinds[1 + n_values:],
                                              n_values,
                                              s)
        if res is None:
            return False
        variance, upper_bound = res

        existing = self.current_symbol_table().get(name)
        if existing and not (isinstance(existing.node, PlaceholderNode) or
                             # Also give error for another type variable with the same name.
                             (isinstance(existing.node, TypeVarExpr) and
                              existing.node is call.analyzed)):
            self.fail('Cannot redefine "%s" as a type variable' % name, s)
            return False

        if self.options.disallow_any_unimported:
            for idx, constraint in enumerate(values, start=1):
                if has_any_from_unimported_type(constraint):
                    prefix = "Constraint {}".format(idx)
                    self.msg.unimported_type_becomes_any(prefix, constraint, s)

            if has_any_from_unimported_type(upper_bound):
                prefix = "Upper bound of type variable"
                self.msg.unimported_type_becomes_any(prefix, upper_bound, s)

        for t in values + [upper_bound]:
            check_for_explicit_any(t, self.options, self.is_typeshed_stub_file, self.msg,
                                   context=s)

        # mypyc suppresses making copies of a function to check each
        # possible type, so set the upper bound to Any to prevent that
        # from causing errors.
        if values and self.options.mypyc:
            upper_bound = AnyType(TypeOfAny.implementation_artifact)

        # Yes, it's a valid type variable definition! Add it to the symbol table.
        if not call.analyzed:
            type_var = TypeVarExpr(name, self.qualified_name(name),
                                   values, upper_bound, variance)
            type_var.line = call.line
            call.analyzed = type_var
        else:
            assert isinstance(call.analyzed, TypeVarExpr)
            if call.analyzed.values != values or call.analyzed.upper_bound != upper_bound:
                self.progress = True
            call.analyzed.upper_bound = upper_bound
            call.analyzed.values = values

        self.add_symbol(name, call.analyzed, s)
        return True

    def check_typevarlike_name(self, call: CallExpr, name: str, context: Context) -> bool:
        """Checks that the name of a TypeVar or ParamSpec matches its variable."""
        name = unmangle(name)
        assert isinstance(call.callee, RefExpr)
        typevarlike_type = (
            call.callee.name if isinstance(call.callee, NameExpr) else call.callee.fullname
        )
        if len(call.args) < 1:
            self.fail("Too few arguments for {}()".format(typevarlike_type), context)
            return False
        if (not isinstance(call.args[0], (StrExpr, BytesExpr, UnicodeExpr))
                or not call.arg_kinds[0] == ARG_POS):
            self.fail("{}() expects a string literal as first argument".format(typevarlike_type),
                      context)
            return False
        elif call.args[0].value != name:
            msg = 'String argument 1 "{}" to {}(...) does not match variable name "{}"'
            self.fail(msg.format(call.args[0].value, typevarlike_type, name), context)
            return False
        return True

    def get_typevarlike_declaration(self, s: AssignmentStmt,
                                    typevarlike_types: Tuple[str, ...]) -> Optional[CallExpr]:
        """Returns the call expression if `s` is a declaration of `typevarlike_type`
        (TypeVar or ParamSpec), or None otherwise.
        """
        if len(s.lvalues) != 1 or not isinstance(s.lvalues[0], NameExpr):
            return None
        if not isinstance(s.rvalue, CallExpr):
            return None
        call = s.rvalue
        callee = call.callee
        if not isinstance(callee, RefExpr):
            return None
        if callee.fullname not in typevarlike_types:
            return None
        return call

    def process_typevar_parameters(self, args: List[Expression],
                                   names: List[Optional[str]],
                                   kinds: List[ArgKind],
                                   num_values: int,
                                   context: Context) -> Optional[Tuple[int, Type]]:
        has_values = (num_values > 0)
        covariant = False
        contravariant = False
        upper_bound: Type = self.object_type()
        for param_value, param_name, param_kind in zip(args, names, kinds):
            if not param_kind.is_named():
                self.fail(message_registry.TYPEVAR_UNEXPECTED_ARGUMENT, context)
                return None
            if param_name == 'covariant':
                if (isinstance(param_value, NameExpr)
                        and param_value.name in ('True', 'False')):
                    covariant = param_value.name == 'True'
                else:
                    self.fail(message_registry.TYPEVAR_VARIANCE_DEF.format(
                        'covariant'), context)
                    return None
            elif param_name == 'contravariant':
                if (isinstance(param_value, NameExpr)
                        and param_value.name in ('True', 'False')):
                    contravariant = param_value.name == 'True'
                else:
                    self.fail(message_registry.TYPEVAR_VARIANCE_DEF.format(
                        'contravariant'), context)
                    return None
            elif param_name == 'bound':
                if has_values:
                    self.fail("TypeVar cannot have both values and an upper bound", context)
                    return None
                try:
                    # We want to use our custom error message below, so we suppress
                    # the default error message for invalid types here.
                    analyzed = self.expr_to_analyzed_type(param_value,
                                                          allow_placeholder=True,
                                                          report_invalid_types=False)
                    if analyzed is None:
                        # Type variables are special: we need to place them in the symbol table
                        # soon, even if upper bound is not ready yet. Otherwise avoiding
                        # a "deadlock" in this common pattern would be tricky:
                        #     T = TypeVar('T', bound=Custom[Any])
                        #     class Custom(Generic[T]):
                        #         ...
                        analyzed = PlaceholderType(None, [], context.line)
                    upper_bound = get_proper_type(analyzed)
                    if isinstance(upper_bound, AnyType) and upper_bound.is_from_error:
                        self.fail(message_registry.TYPEVAR_BOUND_MUST_BE_TYPE, param_value)
                        # Note: we do not return 'None' here -- we want to continue
                        # using the AnyType as the upper bound.
                except TypeTranslationError:
                    self.fail(message_registry.TYPEVAR_BOUND_MUST_BE_TYPE, param_value)
                    return None
            elif param_name == 'values':
                # Probably using obsolete syntax with values=(...). Explain the current syntax.
                self.fail('TypeVar "values" argument not supported', context)
                self.fail("Use TypeVar('T', t, ...) instead of TypeVar('T', values=(t, ...))",
                          context)
                return None
            else:
                self.fail('{}: "{}"'.format(
                    message_registry.TYPEVAR_UNEXPECTED_ARGUMENT, param_name,
                ), context)
                return None

        if covariant and contravariant:
            self.fail("TypeVar cannot be both covariant and contravariant", context)
            return None
        elif num_values == 1:
            self.fail("TypeVar cannot have only a single constraint", context)
            return None
        elif covariant:
            variance = COVARIANT
        elif contravariant:
            variance = CONTRAVARIANT
        else:
            variance = INVARIANT
        return variance, upper_bound

    def process_paramspec_declaration(self, s: AssignmentStmt) -> bool:
        """Checks if s declares a ParamSpec; if yes, store it in symbol table.

        Return True if this looks like a ParamSpec (maybe with errors), otherwise return False.

        In the future, ParamSpec may accept bounds and variance arguments, in which
        case more aggressive sharing of code with process_typevar_declaration should be pursued.
        """
        call = self.get_typevarlike_declaration(
            s, ("typing_extensions.ParamSpec", "typing.ParamSpec")
        )
        if not call:
            return False

        lvalue = s.lvalues[0]
        assert isinstance(lvalue, NameExpr)
        if s.type:
            self.fail("Cannot declare the type of a parameter specification", s)
            return False

        name = lvalue.name
        if not self.check_typevarlike_name(call, name, s):
            return False

        # PEP 612 reserves the right to define bound, covariant and contravariant arguments to
        # ParamSpec in a later PEP. If and when that happens, we should do something
        # on the lines of process_typevar_parameters

        if not call.analyzed:
            paramspec_var = ParamSpecExpr(
                name, self.qualified_name(name), self.object_type(), INVARIANT
            )
            paramspec_var.line = call.line
            call.analyzed = paramspec_var
        else:
            assert isinstance(call.analyzed, ParamSpecExpr)
        self.add_symbol(name, call.analyzed, s)
        return True

    def basic_new_typeinfo(self, name: str,
                           basetype_or_fallback: Instance,
                           line: int) -> TypeInfo:
        if self.is_func_scope() and not self.type and '@' not in name:
            name += '@' + str(line)
        class_def = ClassDef(name, Block([]))
        if self.is_func_scope() and not self.type:
            # Full names of generated classes should always be prefixed with the module names
            # even if they are nested in a function, since these classes will be (de-)serialized.
            # (Note that the caller should append @line to the name to avoid collisions.)
            # TODO: clean this up, see #6422.
            class_def.fullname = self.cur_mod_id + '.' + self.qualified_name(name)
        else:
            class_def.fullname = self.qualified_name(name)

        info = TypeInfo(SymbolTable(), class_def, self.cur_mod_id)
        class_def.info = info
        mro = basetype_or_fallback.type.mro
        if not mro:
            # Forward reference, MRO should be recalculated in third pass.
            mro = [basetype_or_fallback.type, self.object_type().type]
        info.mro = [info] + mro
        info.bases = [basetype_or_fallback]
        return info

    def analyze_value_types(self, items: List[Expression]) -> List[Type]:
        """Analyze types from values expressions in type variable definition."""
        result: List[Type] = []
        for node in items:
            try:
                analyzed = self.anal_type(self.expr_to_unanalyzed_type(node),
                                          allow_placeholder=True)
                if analyzed is None:
                    # Type variables are special: we need to place them in the symbol table
                    # soon, even if some value is not ready yet, see process_typevar_parameters()
                    # for an example.
                    analyzed = PlaceholderType(None, [], node.line)
                result.append(analyzed)
            except TypeTranslationError:
                self.fail('Type expected', node)
                result.append(AnyType(TypeOfAny.from_error))
        return result

    def check_classvar(self, s: AssignmentStmt) -> None:
        """Check if assignment defines a class variable."""
        lvalue = s.lvalues[0]
        if len(s.lvalues) != 1 or not isinstance(lvalue, RefExpr):
            return
        if not s.type or not self.is_classvar(s.type):
            return
        if self.is_class_scope() and isinstance(lvalue, NameExpr):
            node = lvalue.node
            if isinstance(node, Var):
                node.is_classvar = True
            analyzed = self.anal_type(s.type)
            if analyzed is not None and get_type_vars(analyzed):
                # This means that we have a type var defined inside of a ClassVar.
                # This is not allowed by PEP526.
                # See https://github.com/python/mypy/issues/11538
                self.fail(message_registry.CLASS_VAR_WITH_TYPEVARS, s)
        elif not isinstance(lvalue, MemberExpr) or self.is_self_member_ref(lvalue):
            # In case of member access, report error only when assigning to self
            # Other kinds of member assignments should be already reported
            self.fail_invalid_classvar(lvalue)

    def is_classvar(self, typ: Type) -> bool:
        if not isinstance(typ, UnboundType):
            return False
        sym = self.lookup_qualified(typ.name, typ)
        if not sym or not sym.node:
            return False
        return sym.node.fullname == 'typing.ClassVar'

    def is_final_type(self, typ: Optional[Type]) -> bool:
        if not isinstance(typ, UnboundType):
            return False
        sym = self.lookup_qualified(typ.name, typ)
        if not sym or not sym.node:
            return False
        return sym.node.fullname in ('typing.Final', 'typing_extensions.Final')

    def fail_invalid_classvar(self, context: Context) -> None:
        self.fail(message_registry.CLASS_VAR_OUTSIDE_OF_CLASS, context)

    def process_module_assignment(self, lvals: List[Lvalue], rval: Expression,
                                  ctx: AssignmentStmt) -> None:
        """Propagate module references across assignments.

        Recursively handles the simple form of iterable unpacking; doesn't
        handle advanced unpacking with *rest, dictionary unpacking, etc.

        In an expression like x = y = z, z is the rval and lvals will be [x,
        y].

        """
        if (isinstance(rval, (TupleExpr, ListExpr))
                and all(isinstance(v, TupleExpr) for v in lvals)):
            # rval and all lvals are either list or tuple, so we are dealing
            # with unpacking assignment like `x, y = a, b`. Mypy didn't
            # understand our all(isinstance(...)), so cast them as TupleExpr
            # so mypy knows it is safe to access their .items attribute.
            seq_lvals = cast(List[TupleExpr], lvals)
            # given an assignment like:
            #     (x, y) = (m, n) = (a, b)
            # we now have:
            #     seq_lvals = [(x, y), (m, n)]
            #     seq_rval = (a, b)
            # We now zip this into:
            #     elementwise_assignments = [(a, x, m), (b, y, n)]
            # where each elementwise assignment includes one element of rval and the
            # corresponding element of each lval. Basically we unpack
            #     (x, y) = (m, n) = (a, b)
            # into elementwise assignments
            #     x = m = a
            #     y = n = b
            # and then we recursively call this method for each of those assignments.
            # If the rval and all lvals are not all of the same length, zip will just ignore
            # extra elements, so no error will be raised here; mypy will later complain
            # about the length mismatch in type-checking.
            elementwise_assignments = zip(rval.items, *[v.items for v in seq_lvals])
            for rv, *lvs in elementwise_assignments:
                self.process_module_assignment(lvs, rv, ctx)
        elif isinstance(rval, RefExpr):
            rnode = self.lookup_type_node(rval)
            if rnode and isinstance(rnode.node, MypyFile):
                for lval in lvals:
                    if not isinstance(lval, RefExpr):
                        continue
                    # respect explicitly annotated type
                    if (isinstance(lval.node, Var) and lval.node.type is not None):
                        continue

                    # We can handle these assignments to locals and to self
                    if isinstance(lval, NameExpr):
                        lnode = self.current_symbol_table().get(lval.name)
                    elif isinstance(lval, MemberExpr) and self.is_self_member_ref(lval):
                        assert self.type is not None
                        lnode = self.type.names.get(lval.name)
                    else:
                        continue

                    if lnode:
                        if isinstance(lnode.node, MypyFile) and lnode.node is not rnode.node:
                            assert isinstance(lval, (NameExpr, MemberExpr))
                            self.fail(
                                'Cannot assign multiple modules to name "{}" '
                                'without explicit "types.ModuleType" annotation'.format(lval.name),
                                ctx)
                        # never create module alias except on initial var definition
                        elif lval.is_inferred_def:
                            assert rnode.node is not None
                            lnode.node = rnode.node

    def process__all__(self, s: AssignmentStmt) -> None:
        """Export names if argument is a __all__ assignment."""
        if (len(s.lvalues) == 1 and isinstance(s.lvalues[0], NameExpr) and
                s.lvalues[0].name == '__all__' and s.lvalues[0].kind == GDEF and
                isinstance(s.rvalue, (ListExpr, TupleExpr))):
            self.add_exports(s.rvalue.items)

    def process__deletable__(self, s: AssignmentStmt) -> None:
        if not self.options.mypyc:
            return
        if (len(s.lvalues) == 1 and isinstance(s.lvalues[0], NameExpr) and
                s.lvalues[0].name == '__deletable__' and s.lvalues[0].kind == MDEF):
            rvalue = s.rvalue
            if not isinstance(rvalue, (ListExpr, TupleExpr)):
                self.fail('"__deletable__" must be initialized with a list or tuple expression', s)
                return
            items = rvalue.items
            attrs = []
            for item in items:
                if not isinstance(item, StrExpr):
                    self.fail('Invalid "__deletable__" item; string literal expected', item)
                else:
                    attrs.append(item.value)
            assert self.type
            self.type.deletable_attributes = attrs

    def process__slots__(self, s: AssignmentStmt) -> None:
        """
        Processing ``__slots__`` if defined in type.

        See: https://docs.python.org/3/reference/datamodel.html#slots
        """
        # Later we can support `__slots__` defined as `__slots__ = other = ('a', 'b')`
        if (isinstance(self.type, TypeInfo) and
                len(s.lvalues) == 1 and isinstance(s.lvalues[0], NameExpr) and
                s.lvalues[0].name == '__slots__' and s.lvalues[0].kind == MDEF):

            # We understand `__slots__` defined as string, tuple, list, set, and dict:
            if not isinstance(s.rvalue, (StrExpr, ListExpr, TupleExpr, SetExpr, DictExpr)):
                # For example, `__slots__` can be defined as a variable,
                # we don't support it for now.
                return

            if any(p.slots is None for p in self.type.mro[1:-1]):
                # At least one type in mro (excluding `self` and `object`)
                # does not have concrete `__slots__` defined. Ignoring.
                return

            concrete_slots = True
            rvalue: List[Expression] = []
            if isinstance(s.rvalue, StrExpr):
                rvalue.append(s.rvalue)
            elif isinstance(s.rvalue, (ListExpr, TupleExpr, SetExpr)):
                rvalue.extend(s.rvalue.items)
            else:
                # We have a special treatment of `dict` with possible `{**kwargs}` usage.
                # In this case we consider all `__slots__` to be non-concrete.
                for key, _ in s.rvalue.items:
                    if concrete_slots and key is not None:
                        rvalue.append(key)
                    else:
                        concrete_slots = False

            slots = []
            for item in rvalue:
                # Special case for `'__dict__'` value:
                # when specified it will still allow any attribute assignment.
                if isinstance(item, StrExpr) and item.value != '__dict__':
                    slots.append(item.value)
                else:
                    concrete_slots = False
            if not concrete_slots:
                # Some slot items are dynamic, we don't want any false positives,
                # so, we just pretend that this type does not have any slots at all.
                return

            # We need to copy all slots from super types:
            for super_type in self.type.mro[1:-1]:
                assert super_type.slots is not None
                slots.extend(super_type.slots)
            self.type.slots = set(slots)

    #
    # Misc statements
    #

    def visit_block(self, b: Block) -> None:
        if b.is_unreachable:
            return
        self.block_depth[-1] += 1
        for s in b.body:
            self.accept(s)
        self.block_depth[-1] -= 1

    def visit_block_maybe(self, b: Optional[Block]) -> None:
        if b:
            self.visit_block(b)

    def visit_expression_stmt(self, s: ExpressionStmt) -> None:
        self.statement = s
        s.expr.accept(self)

    def visit_return_stmt(self, s: ReturnStmt) -> None:
        self.statement = s
        if not self.is_func_scope():
            self.fail('"return" outside function', s)
        if s.expr:
            s.expr.accept(self)

    def visit_raise_stmt(self, s: RaiseStmt) -> None:
        self.statement = s
        if s.expr:
            s.expr.accept(self)
        if s.from_expr:
            s.from_expr.accept(self)

    def visit_assert_stmt(self, s: AssertStmt) -> None:
        self.statement = s
        if s.expr:
            s.expr.accept(self)
        if s.msg:
            s.msg.accept(self)

    def visit_operator_assignment_stmt(self,
                                       s: OperatorAssignmentStmt) -> None:
        self.statement = s
        s.lvalue.accept(self)
        s.rvalue.accept(self)
        if (isinstance(s.lvalue, NameExpr) and s.lvalue.name == '__all__' and
                s.lvalue.kind == GDEF and isinstance(s.rvalue, (ListExpr, TupleExpr))):
            self.add_exports(s.rvalue.items)

    def visit_while_stmt(self, s: WhileStmt) -> None:
        self.statement = s
        s.expr.accept(self)
        self.loop_depth += 1
        s.body.accept(self)
        self.loop_depth -= 1
        self.visit_block_maybe(s.else_body)

    def visit_for_stmt(self, s: ForStmt) -> None:
        self.statement = s
        s.expr.accept(self)

        # Bind index variables and check if they define new names.
        self.analyze_lvalue(s.index, explicit_type=s.index_type is not None)
        if s.index_type:
            if self.is_classvar(s.index_type):
                self.fail_invalid_classvar(s.index)
            allow_tuple_literal = isinstance(s.index, TupleExpr)
            analyzed = self.anal_type(s.index_type, allow_tuple_literal=allow_tuple_literal)
            if analyzed is not None:
                self.store_declared_types(s.index, analyzed)
                s.index_type = analyzed

        self.loop_depth += 1
        self.visit_block(s.body)
        self.loop_depth -= 1

        self.visit_block_maybe(s.else_body)

    def visit_break_stmt(self, s: BreakStmt) -> None:
        self.statement = s
        if self.loop_depth == 0:
            self.fail('"break" outside loop', s, serious=True, blocker=True)

    def visit_continue_stmt(self, s: ContinueStmt) -> None:
        self.statement = s
        if self.loop_depth == 0:
            self.fail('"continue" outside loop', s, serious=True, blocker=True)

    def visit_if_stmt(self, s: IfStmt) -> None:
        self.statement = s
        infer_reachability_of_if_statement(s, self.options)
        for i in range(len(s.expr)):
            s.expr[i].accept(self)
            self.visit_block(s.body[i])
        self.visit_block_maybe(s.else_body)

    def visit_try_stmt(self, s: TryStmt) -> None:
        self.statement = s
        self.analyze_try_stmt(s, self)

    def analyze_try_stmt(self, s: TryStmt, visitor: NodeVisitor[None]) -> None:
        s.body.accept(visitor)
        for type, var, handler in zip(s.types, s.vars, s.handlers):
            if type:
                type.accept(visitor)
            if var:
                self.analyze_lvalue(var)
            handler.accept(visitor)
        if s.else_body:
            s.else_body.accept(visitor)
        if s.finally_body:
            s.finally_body.accept(visitor)

    def visit_with_stmt(self, s: WithStmt) -> None:
        self.statement = s
        types: List[Type] = []

        if s.unanalyzed_type:
            assert isinstance(s.unanalyzed_type, ProperType)
            actual_targets = [t for t in s.target if t is not None]
            if len(actual_targets) == 0:
                # We have a type for no targets
                self.fail('Invalid type comment: "with" statement has no targets', s)
            elif len(actual_targets) == 1:
                # We have one target and one type
                types = [s.unanalyzed_type]
            elif isinstance(s.unanalyzed_type, TupleType):
                # We have multiple targets and multiple types
                if len(actual_targets) == len(s.unanalyzed_type.items):
                    types = s.unanalyzed_type.items.copy()
                else:
                    # But it's the wrong number of items
                    self.fail('Incompatible number of types for "with" targets', s)
            else:
                # We have multiple targets and one type
                self.fail('Multiple types expected for multiple "with" targets', s)

        new_types: List[Type] = []
        for e, n in zip(s.expr, s.target):
            e.accept(self)
            if n:
                self.analyze_lvalue(n, explicit_type=s.unanalyzed_type is not None)

                # Since we have a target, pop the next type from types
                if types:
                    t = types.pop(0)
                    if self.is_classvar(t):
                        self.fail_invalid_classvar(n)
                    allow_tuple_literal = isinstance(n, TupleExpr)
                    analyzed = self.anal_type(t, allow_tuple_literal=allow_tuple_literal)
                    if analyzed is not None:
                        # TODO: Deal with this better
                        new_types.append(analyzed)
                        self.store_declared_types(n, analyzed)

        s.analyzed_types = new_types

        self.visit_block(s.body)

    def visit_del_stmt(self, s: DelStmt) -> None:
        self.statement = s
        s.expr.accept(self)
        if not self.is_valid_del_target(s.expr):
            self.fail('Invalid delete target', s)

    def is_valid_del_target(self, s: Expression) -> bool:
        if isinstance(s, (IndexExpr, NameExpr, MemberExpr)):
            return True
        elif isinstance(s, (TupleExpr, ListExpr)):
            return all(self.is_valid_del_target(item) for item in s.items)
        else:
            return False

    def visit_global_decl(self, g: GlobalDecl) -> None:
        self.statement = g
        for name in g.names:
            if name in self.nonlocal_decls[-1]:
                self.fail('Name "{}" is nonlocal and global'.format(name), g)
            self.global_decls[-1].add(name)

    def visit_nonlocal_decl(self, d: NonlocalDecl) -> None:
        self.statement = d
        if not self.is_func_scope():
            self.fail("nonlocal declaration not allowed at module level", d)
        else:
            for name in d.names:
                for table in reversed(self.locals[:-1]):
                    if table is not None and name in table:
                        break
                else:
                    self.fail('No binding for nonlocal "{}" found'.format(name), d)

                if self.locals[-1] is not None and name in self.locals[-1]:
                    self.fail('Name "{}" is already defined in local '
                              'scope before nonlocal declaration'.format(name), d)

                if name in self.global_decls[-1]:
                    self.fail('Name "{}" is nonlocal and global'.format(name), d)
                self.nonlocal_decls[-1].add(name)

    def visit_print_stmt(self, s: PrintStmt) -> None:
        self.statement = s
        for arg in s.args:
            arg.accept(self)
        if s.target:
            s.target.accept(self)

    def visit_exec_stmt(self, s: ExecStmt) -> None:
        self.statement = s
        s.expr.accept(self)
        if s.globals:
            s.globals.accept(self)
        if s.locals:
            s.locals.accept(self)

    #
    # Expressions
    #

    def visit_name_expr(self, expr: NameExpr) -> None:
        n = self.lookup(expr.name, expr)
        if n:
            self.bind_name_expr(expr, n)

    def bind_name_expr(self, expr: NameExpr, sym: SymbolTableNode) -> None:
        """Bind name expression to a symbol table node."""
        if isinstance(sym.node, TypeVarExpr) and self.tvar_scope.get_binding(sym):
            self.fail('"{}" is a type variable and only valid in type '
                      'context'.format(expr.name), expr)
        elif isinstance(sym.node, PlaceholderNode):
            self.process_placeholder(expr.name, 'name', expr)
        else:
            expr.kind = sym.kind
            expr.node = sym.node
            expr.fullname = sym.fullname

    def visit_super_expr(self, expr: SuperExpr) -> None:
        if not self.type and not expr.call.args:
            self.fail('"super" used outside class', expr)
            return
        expr.info = self.type
        for arg in expr.call.args:
            arg.accept(self)

    def visit_tuple_expr(self, expr: TupleExpr) -> None:
        for item in expr.items:
            if isinstance(item, StarExpr):
                item.valid = True
            item.accept(self)

    def visit_list_expr(self, expr: ListExpr) -> None:
        for item in expr.items:
            if isinstance(item, StarExpr):
                item.valid = True
            item.accept(self)

    def visit_set_expr(self, expr: SetExpr) -> None:
        for item in expr.items:
            if isinstance(item, StarExpr):
                item.valid = True
            item.accept(self)

    def visit_dict_expr(self, expr: DictExpr) -> None:
        for key, value in expr.items:
            if key is not None:
                key.accept(self)
            value.accept(self)

    def visit_star_expr(self, expr: StarExpr) -> None:
        if not expr.valid:
            # XXX TODO Change this error message
            self.fail('Can use starred expression only as assignment target', expr)
        else:
            expr.expr.accept(self)

    def visit_yield_from_expr(self, e: YieldFromExpr) -> None:
        if not self.is_func_scope():  # not sure
            self.fail('"yield from" outside function', e, serious=True, blocker=True)
        else:
            if self.function_stack[-1].is_coroutine:
                self.fail('"yield from" in async function', e, serious=True, blocker=True)
            else:
                self.function_stack[-1].is_generator = True
        if e.expr:
            e.expr.accept(self)

    def visit_call_expr(self, expr: CallExpr) -> None:
        """Analyze a call expression.

        Some call expressions are recognized as special forms, including
        cast(...).
        """
        expr.callee.accept(self)
        if refers_to_fullname(expr.callee, 'typing.cast'):
            # Special form cast(...).
            if not self.check_fixed_args(expr, 2, 'cast'):
                return
            # Translate first argument to an unanalyzed type.
            try:
                target = self.expr_to_unanalyzed_type(expr.args[0])
            except TypeTranslationError:
                self.fail('Cast target is not a type', expr)
                return
            # Piggyback CastExpr object to the CallExpr object; it takes
            # precedence over the CallExpr semantics.
            expr.analyzed = CastExpr(expr.args[1], target)
            expr.analyzed.line = expr.line
            expr.analyzed.column = expr.column
            expr.analyzed.accept(self)
        elif refers_to_fullname(expr.callee, 'builtins.reveal_type'):
            if not self.check_fixed_args(expr, 1, 'reveal_type'):
                return
            expr.analyzed = RevealExpr(kind=REVEAL_TYPE, expr=expr.args[0])
            expr.analyzed.line = expr.line
            expr.analyzed.column = expr.column
            expr.analyzed.accept(self)
        elif refers_to_fullname(expr.callee, 'builtins.reveal_locals'):
            # Store the local variable names into the RevealExpr for use in the
            # type checking pass
            local_nodes: List[Var] = []
            if self.is_module_scope():
                # try to determine just the variable declarations in module scope
                # self.globals.values() contains SymbolTableNode's
                # Each SymbolTableNode has an attribute node that is nodes.Var
                # look for variable nodes that marked as is_inferred
                # Each symboltable node has a Var node as .node
                local_nodes = [n.node
                               for name, n in self.globals.items()
                               if getattr(n.node, 'is_inferred', False)
                               and isinstance(n.node, Var)]
            elif self.is_class_scope():
                # type = None  # type: Optional[TypeInfo]
                if self.type is not None:
                    local_nodes = [st.node
                                   for st in self.type.names.values()
                                   if isinstance(st.node, Var)]
            elif self.is_func_scope():
                # locals = None  # type: List[Optional[SymbolTable]]
                if self.locals is not None:
                    symbol_table = self.locals[-1]
                    if symbol_table is not None:
                        local_nodes = [st.node
                                       for st in symbol_table.values()
                                       if isinstance(st.node, Var)]
            expr.analyzed = RevealExpr(kind=REVEAL_LOCALS, local_nodes=local_nodes)
            expr.analyzed.line = expr.line
            expr.analyzed.column = expr.column
            expr.analyzed.accept(self)
        elif refers_to_fullname(expr.callee, 'typing.Any'):
            # Special form Any(...) no longer supported.
            self.fail('Any(...) is no longer supported. Use cast(Any, ...) instead', expr)
        elif refers_to_fullname(expr.callee, 'typing._promote'):
            # Special form _promote(...).
            if not self.check_fixed_args(expr, 1, '_promote'):
                return
            # Translate first argument to an unanalyzed type.
            try:
                target = self.expr_to_unanalyzed_type(expr.args[0])
            except TypeTranslationError:
                self.fail('Argument 1 to _promote is not a type', expr)
                return
            expr.analyzed = PromoteExpr(target)
            expr.analyzed.line = expr.line
            expr.analyzed.accept(self)
        elif refers_to_fullname(expr.callee, 'builtins.dict'):
            expr.analyzed = self.translate_dict_call(expr)
        elif refers_to_fullname(expr.callee, 'builtins.divmod'):
            if not self.check_fixed_args(expr, 2, 'divmod'):
                return
            expr.analyzed = OpExpr('divmod', expr.args[0], expr.args[1])
            expr.analyzed.line = expr.line
            expr.analyzed.accept(self)
        else:
            # Normal call expression.
            for a in expr.args:
                a.accept(self)

            if (isinstance(expr.callee, MemberExpr) and
                    isinstance(expr.callee.expr, NameExpr) and
                    expr.callee.expr.name == '__all__' and
                    expr.callee.expr.kind == GDEF and
                    expr.callee.name in ('append', 'extend')):
                if expr.callee.name == 'append' and expr.args:
                    self.add_exports(expr.args[0])
                elif (expr.callee.name == 'extend' and expr.args and
                        isinstance(expr.args[0], (ListExpr, TupleExpr))):
                    self.add_exports(expr.args[0].items)

    def translate_dict_call(self, call: CallExpr) -> Optional[DictExpr]:
        """Translate 'dict(x=y, ...)' to {'x': y, ...} and 'dict()' to {}.

        For other variants of dict(...), return None.
        """
        if not all(kind == ARG_NAMED for kind in call.arg_kinds):
            # Must still accept those args.
            for a in call.args:
                a.accept(self)
            return None
        expr = DictExpr([(StrExpr(cast(str, key)), value)  # since they are all ARG_NAMED
                         for key, value in zip(call.arg_names, call.args)])
        expr.set_line(call)
        expr.accept(self)
        return expr

    def check_fixed_args(self, expr: CallExpr, numargs: int,
                         name: str) -> bool:
        """Verify that expr has specified number of positional args.

        Return True if the arguments are valid.
        """
        s = 's'
        if numargs == 1:
            s = ''
        if len(expr.args) != numargs:
            self.fail('"%s" expects %d argument%s' % (name, numargs, s),
                      expr)
            return False
        if expr.arg_kinds != [ARG_POS] * numargs:
            self.fail('"%s" must be called with %s positional argument%s' %
                      (name, numargs, s), expr)
            return False
        return True

    def visit_member_expr(self, expr: MemberExpr) -> None:
        base = expr.expr
        base.accept(self)
        if isinstance(base, RefExpr) and isinstance(base.node, MypyFile):
            # Handle module attribute.
            sym = self.get_module_symbol(base.node, expr.name)
            if sym:
                if isinstance(sym.node, PlaceholderNode):
                    self.process_placeholder(expr.name, 'attribute', expr)
                    return
                expr.kind = sym.kind
                expr.fullname = sym.fullname
                expr.node = sym.node
        elif isinstance(base, RefExpr):
            # This branch handles the case C.bar (or cls.bar or self.bar inside
            # a classmethod/method), where C is a class and bar is a type
            # definition or a module resulting from `import bar` (or a module
            # assignment) inside class C. We look up bar in the class' TypeInfo
            # namespace.  This is done only when bar is a module or a type;
            # other things (e.g. methods) are handled by other code in
            # checkmember.
            type_info = None
            if isinstance(base.node, TypeInfo):
                # C.bar where C is a class
                type_info = base.node
            elif isinstance(base.node, Var) and self.type and self.function_stack:
                # check for self.bar or cls.bar in method/classmethod
                func_def = self.function_stack[-1]
                if not func_def.is_static and isinstance(func_def.type, CallableType):
                    formal_arg = func_def.type.argument_by_name(base.node.name)
                    if formal_arg and formal_arg.pos == 0:
                        type_info = self.type
            elif isinstance(base.node, TypeAlias) and base.node.no_args:
                assert isinstance(base.node.target, ProperType)
                if isinstance(base.node.target, Instance):
                    type_info = base.node.target.type

            if type_info:
                n = type_info.names.get(expr.name)
                if n is not None and isinstance(n.node, (MypyFile, TypeInfo, TypeAlias)):
                    if not n:
                        return
                    expr.kind = n.kind
                    expr.fullname = n.fullname
                    expr.node = n.node

    def visit_op_expr(self, expr: OpExpr) -> None:
        expr.left.accept(self)

        if expr.op in ('and', 'or'):
            inferred = infer_condition_value(expr.left, self.options)
            if ((inferred in (ALWAYS_FALSE, MYPY_FALSE) and expr.op == 'and') or
                    (inferred in (ALWAYS_TRUE, MYPY_TRUE) and expr.op == 'or')):
                expr.right_unreachable = True
                return
            elif ((inferred in (ALWAYS_TRUE, MYPY_TRUE) and expr.op == 'and') or
                    (inferred in (ALWAYS_FALSE, MYPY_FALSE) and expr.op == 'or')):
                expr.right_always = True

        expr.right.accept(self)

    def visit_comparison_expr(self, expr: ComparisonExpr) -> None:
        for operand in expr.operands:
            operand.accept(self)

    def visit_unary_expr(self, expr: UnaryExpr) -> None:
        expr.expr.accept(self)

    def visit_index_expr(self, expr: IndexExpr) -> None:
        base = expr.base
        base.accept(self)
        if (isinstance(base, RefExpr)
                and isinstance(base.node, TypeInfo)
                and not base.node.is_generic()):
            expr.index.accept(self)
        elif ((isinstance(base, RefExpr) and isinstance(base.node, TypeAlias))
              or refers_to_class_or_function(base)):
            # We need to do full processing on every iteration, since some type
            # arguments may contain placeholder types.
            self.analyze_type_application(expr)
        else:
            expr.index.accept(self)

    def analyze_type_application(self, expr: IndexExpr) -> None:
        """Analyze special form -- type application (either direct or via type aliasing)."""
        types = self.analyze_type_application_args(expr)
        if types is None:
            return
        base = expr.base
        expr.analyzed = TypeApplication(base, types)
        expr.analyzed.line = expr.line
        expr.analyzed.column = expr.column
        # Types list, dict, set are not subscriptable, prohibit this if
        # subscripted either via type alias...
        if isinstance(base, RefExpr) and isinstance(base.node, TypeAlias):
            alias = base.node
            target = get_proper_type(alias.target)
            if isinstance(target, Instance):
                name = target.type.fullname
                if (alias.no_args and  # this avoids bogus errors for already reported aliases
                        name in get_nongen_builtins(self.options.python_version) and
                        not self.is_stub_file and
                        not alias.normalized):
                    self.fail(no_subscript_builtin_alias(name, propose_alt=False), expr)
        # ...or directly.
        else:
            n = self.lookup_type_node(base)
            if (n and n.fullname in get_nongen_builtins(self.options.python_version) and
                    not self.is_stub_file):
                self.fail(no_subscript_builtin_alias(n.fullname, propose_alt=False), expr)

    def analyze_type_application_args(self, expr: IndexExpr) -> Optional[List[Type]]:
        """Analyze type arguments (index) in a type application.

        Return None if anything was incomplete.
        """
        index = expr.index
        tag = self.track_incomplete_refs()
        self.analyze_type_expr(index)
        if self.found_incomplete_ref(tag):
            return None
        types: List[Type] = []
        if isinstance(index, TupleExpr):
            items = index.items
            is_tuple = isinstance(expr.base, RefExpr) and expr.base.fullname == 'builtins.tuple'
            if is_tuple and len(items) == 2 and isinstance(items[-1], EllipsisExpr):
                items = items[:-1]
        else:
            items = [index]
        for item in items:
            try:
                typearg = self.expr_to_unanalyzed_type(item)
            except TypeTranslationError:
                self.fail('Type expected within [...]', expr)
                return None
            # We always allow unbound type variables in IndexExpr, since we
            # may be analysing a type alias definition rvalue. The error will be
            # reported elsewhere if it is not the case.
            analyzed = self.anal_type(typearg, allow_unbound_tvars=True,
                                      allow_placeholder=True)
            if analyzed is None:
                return None
            types.append(analyzed)
        return types

    def visit_slice_expr(self, expr: SliceExpr) -> None:
        if expr.begin_index:
            expr.begin_index.accept(self)
        if expr.end_index:
            expr.end_index.accept(self)
        if expr.stride:
            expr.stride.accept(self)

    def visit_cast_expr(self, expr: CastExpr) -> None:
        expr.expr.accept(self)
        analyzed = self.anal_type(expr.type)
        if analyzed is not None:
            expr.type = analyzed

    def visit_reveal_expr(self, expr: RevealExpr) -> None:
        if expr.kind == REVEAL_TYPE:
            if expr.expr is not None:
                expr.expr.accept(self)
        else:
            # Reveal locals doesn't have an inner expression, there's no
            # need to traverse inside it
            pass

    def visit_type_application(self, expr: TypeApplication) -> None:
        expr.expr.accept(self)
        for i in range(len(expr.types)):
            analyzed = self.anal_type(expr.types[i])
            if analyzed is not None:
                expr.types[i] = analyzed

    def visit_list_comprehension(self, expr: ListComprehension) -> None:
        expr.generator.accept(self)

    def visit_set_comprehension(self, expr: SetComprehension) -> None:
        expr.generator.accept(self)

    def visit_dictionary_comprehension(self, expr: DictionaryComprehension) -> None:
        with self.enter(expr):
            self.analyze_comp_for(expr)
            expr.key.accept(self)
            expr.value.accept(self)
        self.analyze_comp_for_2(expr)

    def visit_generator_expr(self, expr: GeneratorExpr) -> None:
        with self.enter(expr):
            self.analyze_comp_for(expr)
            expr.left_expr.accept(self)
        self.analyze_comp_for_2(expr)

    def analyze_comp_for(self, expr: Union[GeneratorExpr,
                                           DictionaryComprehension]) -> None:
        """Analyses the 'comp_for' part of comprehensions (part 1).

        That is the part after 'for' in (x for x in l if p). This analyzes
        variables and conditions which are analyzed in a local scope.
        """
        for i, (index, sequence, conditions) in enumerate(zip(expr.indices,
                                                              expr.sequences,
                                                              expr.condlists)):
            if i > 0:
                sequence.accept(self)
            # Bind index variables.
            self.analyze_lvalue(index)
            for cond in conditions:
                cond.accept(self)

    def analyze_comp_for_2(self, expr: Union[GeneratorExpr,
                                             DictionaryComprehension]) -> None:
        """Analyses the 'comp_for' part of comprehensions (part 2).

        That is the part after 'for' in (x for x in l if p). This analyzes
        the 'l' part which is analyzed in the surrounding scope.
        """
        expr.sequences[0].accept(self)

    def visit_lambda_expr(self, expr: LambdaExpr) -> None:
        self.analyze_arg_initializers(expr)
        self.analyze_function_body(expr)

    def visit_conditional_expr(self, expr: ConditionalExpr) -> None:
        expr.if_expr.accept(self)
        expr.cond.accept(self)
        expr.else_expr.accept(self)

    def visit_backquote_expr(self, expr: BackquoteExpr) -> None:
        expr.expr.accept(self)

    def visit__promote_expr(self, expr: PromoteExpr) -> None:
        analyzed = self.anal_type(expr.type)
        if analyzed is not None:
            expr.type = analyzed

    def visit_yield_expr(self, expr: YieldExpr) -> None:
        if not self.is_func_scope():
            self.fail('"yield" outside function', expr, serious=True, blocker=True)
        else:
            if self.function_stack[-1].is_coroutine:
                if self.options.python_version < (3, 6):
                    self.fail('"yield" in async function', expr, serious=True, blocker=True)
                else:
                    self.function_stack[-1].is_generator = True
                    self.function_stack[-1].is_async_generator = True
            else:
                self.function_stack[-1].is_generator = True
        if expr.expr:
            expr.expr.accept(self)

    def visit_await_expr(self, expr: AwaitExpr) -> None:
        if not self.is_func_scope():
            self.fail('"await" outside function', expr)
        elif not self.function_stack[-1].is_coroutine:
            self.fail('"await" outside coroutine ("async def")', expr)
        expr.expr.accept(self)

    #
    # Lookup functions
    #

    def lookup(self, name: str, ctx: Context,
               suppress_errors: bool = False) -> Optional[SymbolTableNode]:
        """Look up an unqualified (no dots) name in all active namespaces.

        Note that the result may contain a PlaceholderNode. The caller may
        want to defer in that case.

        Generate an error if the name is not defined unless suppress_errors
        is true or the current namespace is incomplete. In the latter case
        defer.
        """
        implicit_name = False
        # 1a. Name declared using 'global x' takes precedence
        if name in self.global_decls[-1]:
            if name in self.globals:
                return self.globals[name]
            if not suppress_errors:
                self.name_not_defined(name, ctx)
            return None
        # 1b. Name declared using 'nonlocal x' takes precedence
        if name in self.nonlocal_decls[-1]:
            for table in reversed(self.locals[:-1]):
                if table is not None and name in table:
                    return table[name]
            else:
                if not suppress_errors:
                    self.name_not_defined(name, ctx)
                return None
        # 2. Class attributes (if within class definition)
        if self.type and not self.is_func_scope() and name in self.type.names:
            node = self.type.names[name]
            if not node.implicit:
                if self.is_active_symbol_in_class_body(node.node):
                    return node
            else:
                # Defined through self.x assignment
                implicit_name = True
                implicit_node = node
        # 3. Local (function) scopes
        for table in reversed(self.locals):
            if table is not None and name in table:
                return table[name]
        # 4. Current file global scope
        if name in self.globals:
            return self.globals[name]
        # 5. Builtins
        b = self.globals.get('__builtins__', None)
        if b:
            assert isinstance(b.node, MypyFile)
            table = b.node.names
            if name in table:
                if name[0] == "_" and name[1] != "_":
                    if not suppress_errors:
                        self.name_not_defined(name, ctx)
                    return None
                node = table[name]
                return node
        # Give up.
        if not implicit_name and not suppress_errors:
            self.name_not_defined(name, ctx)
        else:
            if implicit_name:
                return implicit_node
        return None

    def is_active_symbol_in_class_body(self, node: Optional[SymbolNode]) -> bool:
        """Can a symbol defined in class body accessed at current statement?

        Only allow access to class attributes textually after
        the definition, so that it's possible to fall back to the
        outer scope. Example:

            class X: ...

            class C:
                X = X  # Initializer refers to outer scope

        Nested classes are an exception, since we want to support
        arbitrary forward references in type annotations.
        """
        # TODO: Forward reference to name imported in class body is not
        #       caught.
        assert self.statement  # we are at class scope
        return (node is None
                or self.is_textually_before_statement(node)
                or not self.is_defined_in_current_module(node.fullname)
                or isinstance(node, TypeInfo)
                or (isinstance(node, PlaceholderNode) and node.becomes_typeinfo))

    def is_textually_before_statement(self, node: SymbolNode) -> bool:
        """Check if a node is defined textually before the current statement

        Note that decorated functions' line number are the same as
        the top decorator.
        """
        assert self.statement
        line_diff = self.statement.line - node.line

        # The first branch handles reference an overloaded function variant inside itself,
        # this is a corner case where mypy technically deviates from runtime name resolution,
        # but it is fine because we want an overloaded function to be treated as a single unit.
        if self.is_overloaded_item(node, self.statement):
            return False
        elif isinstance(node, Decorator) and not node.is_overload:
            return line_diff > len(node.original_decorators)
        else:
            return line_diff > 0

    def is_overloaded_item(self, node: SymbolNode, statement: Statement) -> bool:
        """Check whether the function belongs to the overloaded variants"""
        if isinstance(node, OverloadedFuncDef) and isinstance(statement, FuncDef):
            in_items = statement in {item.func if isinstance(item, Decorator)
                                     else item for item in node.items}
            in_impl = (node.impl is not None and
                      ((isinstance(node.impl, Decorator) and statement is node.impl.func)
                       or statement is node.impl))
            return in_items or in_impl
        return False

    def is_defined_in_current_module(self, fullname: Optional[str]) -> bool:
        if fullname is None:
            return False
        return module_prefix(self.modules, fullname) == self.cur_mod_id

    def lookup_qualified(self, name: str, ctx: Context,
                         suppress_errors: bool = False) -> Optional[SymbolTableNode]:
        """Lookup a qualified name in all activate namespaces.

        Note that the result may contain a PlaceholderNode. The caller may
        want to defer in that case.

        Generate an error if the name is not defined unless suppress_errors
        is true or the current namespace is incomplete. In the latter case
        defer.
        """
        if '.' not in name:
            # Simple case: look up a short name.
            return self.lookup(name, ctx, suppress_errors=suppress_errors)
        parts = name.split('.')
        namespace = self.cur_mod_id
        sym = self.lookup(parts[0], ctx, suppress_errors=suppress_errors)
        if sym:
            for i in range(1, len(parts)):
                node = sym.node
                part = parts[i]
                if isinstance(node, TypeInfo):
                    nextsym = node.get(part)
                elif isinstance(node, MypyFile):
                    nextsym = self.get_module_symbol(node, part)
                    namespace = node.fullname
                elif isinstance(node, PlaceholderNode):
                    return sym
                elif isinstance(node, TypeAlias) and node.no_args:
                    assert isinstance(node.target, ProperType)
                    if isinstance(node.target, Instance):
                        nextsym = node.target.type.get(part)
                else:
                    if isinstance(node, Var):
                        typ = get_proper_type(node.type)
                        if isinstance(typ, AnyType):
                            # Allow access through Var with Any type without error.
                            return self.implicit_symbol(sym, name, parts[i:], typ)
                    # Lookup through invalid node, such as variable or function
                    nextsym = None
                if not nextsym or nextsym.module_hidden:
                    if not suppress_errors:
                        self.name_not_defined(name, ctx, namespace=namespace)
                    return None
                sym = nextsym
        return sym

    def lookup_type_node(self, expr: Expression) -> Optional[SymbolTableNode]:
        try:
            t = self.expr_to_unanalyzed_type(expr)
        except TypeTranslationError:
            return None
        if isinstance(t, UnboundType):
            n = self.lookup_qualified(t.name, expr, suppress_errors=True)
            return n
        return None

    def get_module_symbol(self, node: MypyFile, name: str) -> Optional[SymbolTableNode]:
        """Look up a symbol from a module.

        Return None if no matching symbol could be bound.
        """
        module = node.fullname
        names = node.names
        sym = names.get(name)
        if not sym:
            fullname = module + '.' + name
            if fullname in self.modules:
                sym = SymbolTableNode(GDEF, self.modules[fullname])
            elif self.is_incomplete_namespace(module):
                self.record_incomplete_ref()
            elif ('__getattr__' in names
                    and (node.is_stub
                         or self.options.python_version >= (3, 7))):
                gvar = self.create_getattr_var(names['__getattr__'], name, fullname)
                if gvar:
                    sym = SymbolTableNode(GDEF, gvar)
            elif self.is_missing_module(fullname):
                # We use the fullname of the original definition so that we can
                # detect whether two names refer to the same thing.
                var_type = AnyType(TypeOfAny.from_unimported_type)
                v = Var(name, type=var_type)
                v._fullname = fullname
                sym = SymbolTableNode(GDEF, v)
        elif sym.module_hidden:
            sym = None
        return sym

    def is_missing_module(self, module: str) -> bool:
        return module in self.missing_modules

    def implicit_symbol(self, sym: SymbolTableNode, name: str, parts: List[str],
                        source_type: AnyType) -> SymbolTableNode:
        """Create symbol for a qualified name reference through Any type."""
        if sym.node is None:
            basename = None
        else:
            basename = sym.node.fullname
        if basename is None:
            fullname = name
        else:
            fullname = basename + '.' + '.'.join(parts)
        var_type = AnyType(TypeOfAny.from_another_any, source_type)
        var = Var(parts[-1], var_type)
        var._fullname = fullname
        return SymbolTableNode(GDEF, var)

    def create_getattr_var(self, getattr_defn: SymbolTableNode,
                           name: str, fullname: str) -> Optional[Var]:
        """Create a dummy variable using module-level __getattr__ return type.

        If not possible, return None.

        Note that multiple Var nodes can be created for a single name. We
        can use the from_module_getattr and the fullname attributes to
        check if two dummy Var nodes refer to the same thing. Reusing Var
        nodes would require non-local mutable state, which we prefer to
        avoid.
        """
        if isinstance(getattr_defn.node, (FuncDef, Var)):
            node_type = get_proper_type(getattr_defn.node.type)
            if isinstance(node_type, CallableType):
                typ = node_type.ret_type
            else:
                typ = AnyType(TypeOfAny.from_error)
            v = Var(name, type=typ)
            v._fullname = fullname
            v.from_module_getattr = True
            return v
        return None

    def lookup_fully_qualified(self, fullname: str) -> SymbolTableNode:
        ret = self.lookup_fully_qualified_or_none(fullname)
        assert ret is not None
        return ret

    def lookup_fully_qualified_or_none(self, fullname: str) -> Optional[SymbolTableNode]:
        """Lookup a fully qualified name that refers to a module-level definition.

        Don't assume that the name is defined. This happens in the global namespace --
        the local module namespace is ignored. This does not dereference indirect
        refs.

        Note that this can't be used for names nested in class namespaces.
        """
        # TODO: unify/clean-up/simplify lookup methods, see #4157.
        # TODO: support nested classes (but consider performance impact,
        #       we might keep the module level only lookup for thing like 'builtins.int').
        assert '.' in fullname
        module, name = fullname.rsplit('.', maxsplit=1)
        if module not in self.modules:
            return None
        filenode = self.modules[module]
        result = filenode.names.get(name)
        if result is None and self.is_incomplete_namespace(module):
            # TODO: More explicit handling of incomplete refs?
            self.record_incomplete_ref()
        return result

    def object_type(self) -> Instance:
        return self.named_type('builtins.object')

    def str_type(self) -> Instance:
        return self.named_type('builtins.str')

    def named_type(self, fullname: str, args: Optional[List[Type]] = None) -> Instance:
        sym = self.lookup_fully_qualified(fullname)
        assert sym, "Internal error: attempted to construct unknown type"
        node = sym.node
        assert isinstance(node, TypeInfo)
        if args:
            # TODO: assert len(args) == len(node.defn.type_vars)
            return Instance(node, args)
        return Instance(node, [AnyType(TypeOfAny.special_form)] * len(node.defn.type_vars))

    def named_type_or_none(self, fullname: str,
                           args: Optional[List[Type]] = None) -> Optional[Instance]:
        sym = self.lookup_fully_qualified_or_none(fullname)
        if not sym or isinstance(sym.node, PlaceholderNode):
            return None
        node = sym.node
        if isinstance(node, TypeAlias):
            assert isinstance(node.target, Instance)  # type: ignore
            node = node.target.type
        assert isinstance(node, TypeInfo), node
        if args is not None:
            # TODO: assert len(args) == len(node.defn.type_vars)
            return Instance(node, args)
        return Instance(node, [AnyType(TypeOfAny.unannotated)] * len(node.defn.type_vars))

    def builtin_type(self, fully_qualified_name: str) -> Instance:
        """Legacy function -- use named_type() instead."""
        return self.named_type(fully_qualified_name)

    def lookup_current_scope(self, name: str) -> Optional[SymbolTableNode]:
        if self.locals[-1] is not None:
            return self.locals[-1].get(name)
        elif self.type is not None:
            return self.type.names.get(name)
        else:
            return self.globals.get(name)

    #
    # Adding symbols
    #

    def add_symbol(self,
                   name: str,
                   node: SymbolNode,
                   context: Context,
                   module_public: bool = True,
                   module_hidden: bool = False,
                   can_defer: bool = True,
                   escape_comprehensions: bool = False) -> bool:
        """Add symbol to the currently active symbol table.

        Generally additions to symbol table should go through this method or
        one of the methods below so that kinds, redefinitions, conditional
        definitions, and skipped names are handled consistently.

        Return True if we actually added the symbol, or False if we refused to do so
        (because something is not ready).

        If can_defer is True, defer current target if adding a placeholder.
        """
        if self.is_func_scope():
            kind = LDEF
        elif self.type is not None:
            kind = MDEF
        else:
            kind = GDEF
        symbol = SymbolTableNode(kind,
                                 node,
                                 module_public=module_public,
                                 module_hidden=module_hidden)
        return self.add_symbol_table_node(name, symbol, context, can_defer, escape_comprehensions)

    def add_symbol_skip_local(self, name: str, node: SymbolNode) -> None:
        """Same as above, but skipping the local namespace.

        This doesn't check for previous definition and is only used
        for serialization of method-level classes.

        Classes defined within methods can be exposed through an
        attribute type, but method-level symbol tables aren't serialized.
        This method can be used to add such classes to an enclosing,
        serialized symbol table.
        """
        # TODO: currently this is only used by named tuples. Use this method
        # also by typed dicts and normal classes, see issue #6422.
        if self.type is not None:
            names = self.type.names
            kind = MDEF
        else:
            names = self.globals
            kind = GDEF
        symbol = SymbolTableNode(kind, node)
        names[name] = symbol

    def add_symbol_table_node(self,
                              name: str,
                              symbol: SymbolTableNode,
                              context: Optional[Context] = None,
                              can_defer: bool = True,
                              escape_comprehensions: bool = False) -> bool:
        """Add symbol table node to the currently active symbol table.

        Return True if we actually added the symbol, or False if we refused
        to do so (because something is not ready or it was a no-op).

        Generate an error if there is an invalid redefinition.

        If context is None, unconditionally add node, since we can't report
        an error. Note that this is used by plugins to forcibly replace nodes!

        TODO: Prevent plugins from replacing nodes, as it could cause problems?

        Args:
            name: short name of symbol
            symbol: Node to add
            can_defer: if True, defer current target if adding a placeholder
            context: error context (see above about None value)
        """
        names = self.current_symbol_table(escape_comprehensions=escape_comprehensions)
        existing = names.get(name)
        if isinstance(symbol.node, PlaceholderNode) and can_defer:
            if context is not None:
                self.process_placeholder(name, 'name', context)
            else:
                # see note in docstring describing None contexts
                self.defer()
        if (existing is not None
                and context is not None
                and not is_valid_replacement(existing, symbol)):
            # There is an existing node, so this may be a redefinition.
            # If the new node points to the same node as the old one,
            # or if both old and new nodes are placeholders, we don't
            # need to do anything.
            old = existing.node
            new = symbol.node
            if isinstance(new, PlaceholderNode):
                # We don't know whether this is okay. Let's wait until the next iteration.
                return False
            if not is_same_symbol(old, new):
                if isinstance(new, (FuncDef, Decorator, OverloadedFuncDef, TypeInfo)):
                    self.add_redefinition(names, name, symbol)
                if not (isinstance(new, (FuncDef, Decorator))
                        and self.set_original_def(old, new)):
                    self.name_already_defined(name, context, existing)
        elif (name not in self.missing_names[-1] and '*' not in self.missing_names[-1]):
            names[name] = symbol
            self.progress = True
            return True
        return False

    def add_redefinition(self,
                         names: SymbolTable,
                         name: str,
                         symbol: SymbolTableNode) -> None:
        """Add a symbol table node that reflects a redefinition as a function or a class.

        Redefinitions need to be added to the symbol table so that they can be found
        through AST traversal, but they have dummy names of form 'name-redefinition[N]',
        where N ranges over 2, 3, ... (omitted for the first redefinition).

        Note: we always store redefinitions independently of whether they are valid or not
        (so they will be semantically analyzed), the caller should give an error for invalid
        redefinitions (such as e.g. variable redefined as a class).
        """
        i = 1
        # Don't serialize redefined nodes. They are likely to have
        # busted internal references which can cause problems with
        # serialization and they can't have any external references to
        # them.
        symbol.no_serialize = True
        while True:
            if i == 1:
                new_name = '{}-redefinition'.format(name)
            else:
                new_name = '{}-redefinition{}'.format(name, i)
            existing = names.get(new_name)
            if existing is None:
                names[new_name] = symbol
                return
            elif existing.node is symbol.node:
                # Already there
                return
            i += 1

    def add_local(self, node: Union[Var, FuncDef, OverloadedFuncDef], context: Context) -> None:
        """Add local variable or function."""
        assert self.is_func_scope()
        name = node.name
        node._fullname = name
        self.add_symbol(name, node, context)

    def add_module_symbol(self,
                          id: str,
                          as_id: str,
                          context: Context,
                          module_public: bool,
                          module_hidden: bool) -> None:
        """Add symbol that is a reference to a module object."""
        if id in self.modules:
            node = self.modules[id]
            self.add_symbol(as_id, node, context,
                            module_public=module_public,
                            module_hidden=module_hidden)
        else:
            self.add_unknown_imported_symbol(
                as_id, context, target_name=id, module_public=module_public,
                module_hidden=module_hidden
            )

    def add_imported_symbol(self,
                            name: str,
                            node: SymbolTableNode,
                            context: Context,
                            module_public: bool,
                            module_hidden: bool) -> None:
        """Add an alias to an existing symbol through import."""
        assert not module_hidden or not module_public
        symbol = SymbolTableNode(node.kind, node.node,
                                 module_public=module_public,
                                 module_hidden=module_hidden)
        self.add_symbol_table_node(name, symbol, context)

    def add_unknown_imported_symbol(self,
                                    name: str,
                                    context: Context,
                                    target_name: Optional[str],
                                    module_public: bool,
                                    module_hidden: bool) -> None:
        """Add symbol that we don't know what it points to because resolving an import failed.

        This can happen if a module is missing, or it is present, but doesn't have
        the imported attribute. The `target_name` is the name of symbol in the namespace
        it is imported from. For example, for 'from mod import x as y' the target_name is
        'mod.x'. This is currently used only to track logical dependencies.
        """
        existing = self.current_symbol_table().get(name)
        if existing and isinstance(existing.node, Var) and existing.node.is_suppressed_import:
            # This missing import was already added -- nothing to do here.
            return
        var = Var(name)
        if self.options.logical_deps and target_name is not None:
            # This makes it possible to add logical fine-grained dependencies
            # from a missing module. We can't use this by default, since in a
            # few places we assume that the full name points to a real
            # definition, but this name may point to nothing.
            var._fullname = target_name
        elif self.type:
            var._fullname = self.type.fullname + "." + name
            var.info = self.type
        else:
            var._fullname = self.qualified_name(name)
        var.is_ready = True
        any_type = AnyType(TypeOfAny.from_unimported_type, missing_import_name=var._fullname)
        var.type = any_type
        var.is_suppressed_import = True
        self.add_symbol(
            name, var, context, module_public=module_public, module_hidden=module_hidden
        )

    #
    # Other helpers
    #

    @contextmanager
    def tvar_scope_frame(self, frame: TypeVarLikeScope) -> Iterator[None]:
        old_scope = self.tvar_scope
        self.tvar_scope = frame
        yield
        self.tvar_scope = old_scope

    def defer(self, debug_context: Optional[Context] = None) -> None:
        """Defer current analysis target to be analyzed again.

        This must be called if something in the current target is
        incomplete or has a placeholder node. However, this must *not*
        be called during the final analysis iteration! Instead, an error
        should be generated. Often 'process_placeholder' is a good
        way to either defer or generate an error.

        NOTE: Some methods, such as 'anal_type', 'mark_incomplete' and
              'record_incomplete_ref', call this implicitly, or when needed.
              They are usually preferable to a direct defer() call.
        """
        assert not self.final_iteration, 'Must not defer during final iteration'
        self.deferred = True
        # Store debug info for this deferral.
        line = (debug_context.line if debug_context else
                self.statement.line if self.statement else -1)
        self.deferral_debug_context.append((self.cur_mod_id, line))

    def track_incomplete_refs(self) -> Tag:
        """Return tag that can be used for tracking references to incomplete names."""
        return self.num_incomplete_refs

    def found_incomplete_ref(self, tag: Tag) -> bool:
        """Have we encountered an incomplete reference since starting tracking?"""
        return self.num_incomplete_refs != tag

    def record_incomplete_ref(self) -> None:
        """Record the encounter of an incomplete reference and defer current analysis target."""
        self.defer()
        self.num_incomplete_refs += 1

    def mark_incomplete(self, name: str, node: Node,
                        becomes_typeinfo: bool = False,
                        module_public: bool = True,
                        module_hidden: bool = False) -> None:
        """Mark a definition as incomplete (and defer current analysis target).

        Also potentially mark the current namespace as incomplete.

        Args:
            name: The name that we weren't able to define (or '*' if the name is unknown)
            node: The node that refers to the name (definition or lvalue)
            becomes_typeinfo: Pass this to PlaceholderNode (used by special forms like
                named tuples that will create TypeInfos).
        """
        self.defer(node)
        if name == '*':
            self.incomplete = True
        elif not self.is_global_or_nonlocal(name):
            fullname = self.qualified_name(name)
            assert self.statement
            placeholder = PlaceholderNode(fullname, node, self.statement.line,
                                          becomes_typeinfo=becomes_typeinfo)
            self.add_symbol(name, placeholder,
                            module_public=module_public, module_hidden=module_hidden,
                            context=dummy_context())
        self.missing_names[-1].add(name)

    def is_incomplete_namespace(self, fullname: str) -> bool:
        """Is a module or class namespace potentially missing some definitions?

        If a name is missing from an incomplete namespace, we'll need to defer the
        current analysis target.
        """
        return fullname in self.incomplete_namespaces

    def process_placeholder(self, name: str, kind: str, ctx: Context) -> None:
        """Process a reference targeting placeholder node.

        If this is not a final iteration, defer current node,
        otherwise report an error.

        The 'kind' argument indicates if this a name or attribute expression
        (used for better error message).
        """
        if self.final_iteration:
            self.cannot_resolve_name(name, kind, ctx)
        else:
            self.defer(ctx)

    def cannot_resolve_name(self, name: str, kind: str, ctx: Context) -> None:
        self.fail('Cannot resolve {} "{}" (possible cyclic definition)'.format(kind, name), ctx)

    def qualified_name(self, name: str) -> str:
        if self.type is not None:
            return self.type._fullname + '.' + name
        elif self.is_func_scope():
            return name
        else:
            return self.cur_mod_id + '.' + name

    @contextmanager
    def enter(self,
              function: Union[FuncItem, GeneratorExpr, DictionaryComprehension]) -> Iterator[None]:
        """Enter a function, generator or comprehension scope."""
        names = self.saved_locals.setdefault(function, SymbolTable())
        self.locals.append(names)
        is_comprehension = isinstance(function, (GeneratorExpr, DictionaryComprehension))
        self.is_comprehension_stack.append(is_comprehension)
        self.global_decls.append(set())
        self.nonlocal_decls.append(set())
        # -1 since entering block will increment this to 0.
        self.block_depth.append(-1)
        self.missing_names.append(set())
        try:
            yield
        finally:
            self.locals.pop()
            self.is_comprehension_stack.pop()
            self.global_decls.pop()
            self.nonlocal_decls.pop()
            self.block_depth.pop()
            self.missing_names.pop()

    def is_func_scope(self) -> bool:
        return self.locals[-1] is not None

    def is_nested_within_func_scope(self) -> bool:
        """Are we underneath a function scope, even if we are in a nested class also?"""
        return any(l is not None for l in self.locals)

    def is_class_scope(self) -> bool:
        return self.type is not None and not self.is_func_scope()

    def is_module_scope(self) -> bool:
        return not (self.is_class_scope() or self.is_func_scope())

    def current_symbol_kind(self) -> int:
        if self.is_class_scope():
            kind = MDEF
        elif self.is_func_scope():
            kind = LDEF
        else:
            kind = GDEF
        return kind

    def current_symbol_table(self, escape_comprehensions: bool = False) -> SymbolTable:
        if self.is_func_scope():
            assert self.locals[-1] is not None
            if escape_comprehensions:
                assert len(self.locals) == len(self.is_comprehension_stack)
                # Retrieve the symbol table from the enclosing non-comprehension scope.
                for i, is_comprehension in enumerate(reversed(self.is_comprehension_stack)):
                    if not is_comprehension:
                        if i == len(self.locals) - 1:  # The last iteration.
                            # The caller of the comprehension is in the global space.
                            names = self.globals
                        else:
                            names_candidate = self.locals[-1 - i]
                            assert names_candidate is not None, \
                                "Escaping comprehension from invalid scope"
                            names = names_candidate
                        break
                else:
                    assert False, "Should have at least one non-comprehension scope"
            else:
                names = self.locals[-1]
            assert names is not None
        elif self.type is not None:
            names = self.type.names
        else:
            names = self.globals
        return names

    def is_global_or_nonlocal(self, name: str) -> bool:
        return (self.is_func_scope()
                and (name in self.global_decls[-1]
                     or name in self.nonlocal_decls[-1]))

    def add_exports(self, exp_or_exps: Union[Iterable[Expression], Expression]) -> None:
        exps = [exp_or_exps] if isinstance(exp_or_exps, Expression) else exp_or_exps
        for exp in exps:
            if isinstance(exp, StrExpr):
                self.all_exports.append(exp.value)

    def check_no_global(self,
                        name: str,
                        ctx: Context,
                        is_overloaded_func: bool = False) -> None:
        if name in self.globals:
            prev_is_overloaded = isinstance(self.globals[name], OverloadedFuncDef)
            if is_overloaded_func and prev_is_overloaded:
                self.fail("Nonconsecutive overload {} found".format(name), ctx)
            elif prev_is_overloaded:
                self.fail("Definition of '{}' missing 'overload'".format(name), ctx)
            else:
                self.name_already_defined(name, ctx, self.globals[name])

    def name_not_defined(self, name: str, ctx: Context, namespace: Optional[str] = None) -> None:
        incomplete = self.is_incomplete_namespace(namespace or self.cur_mod_id)
        if (namespace is None
                and self.type
                and not self.is_func_scope()
                and self.incomplete_type_stack[-1]
                and not self.final_iteration):
            # We are processing a class body for the first time, so it is incomplete.
            incomplete = True
        if incomplete:
            # Target namespace is incomplete, so it's possible that the name will be defined
            # later on. Defer current target.
            self.record_incomplete_ref()
            return
        message = 'Name "{}" is not defined'.format(name)
        self.fail(message, ctx, code=codes.NAME_DEFINED)

        if 'builtins.{}'.format(name) in SUGGESTED_TEST_FIXTURES:
            # The user probably has a missing definition in a test fixture. Let's verify.
            fullname = 'builtins.{}'.format(name)
            if self.lookup_fully_qualified_or_none(fullname) is None:
                # Yes. Generate a helpful note.
                self.msg.add_fixture_note(fullname, ctx)

        modules_with_unimported_hints = {
            name.split('.', 1)[0]
            for name in TYPES_FOR_UNIMPORTED_HINTS
        }
        lowercased = {
            name.lower(): name
            for name in TYPES_FOR_UNIMPORTED_HINTS
        }
        for module in modules_with_unimported_hints:
            fullname = '{}.{}'.format(module, name).lower()
            if fullname not in lowercased:
                continue
            # User probably forgot to import these types.
            hint = (
                'Did you forget to import it from "{module}"?'
                ' (Suggestion: "from {module} import {name}")'
            ).format(module=module, name=lowercased[fullname].rsplit('.', 1)[-1])
            self.note(hint, ctx, code=codes.NAME_DEFINED)

    def already_defined(self,
                        name: str,
                        ctx: Context,
                        original_ctx: Optional[Union[SymbolTableNode, SymbolNode]],
                        noun: str) -> None:
        if isinstance(original_ctx, SymbolTableNode):
            node: Optional[SymbolNode] = original_ctx.node
        elif isinstance(original_ctx, SymbolNode):
            node = original_ctx
        else:
            node = None

        if isinstance(original_ctx, SymbolTableNode) and isinstance(original_ctx.node, MypyFile):
            # Since this is an import, original_ctx.node points to the module definition.
            # Therefore its line number is always 1, which is not useful for this
            # error message.
            extra_msg = ' (by an import)'
        elif node and node.line != -1 and self.is_local_name(node.fullname):
            # TODO: Using previous symbol node may give wrong line. We should use
            #       the line number where the binding was established instead.
            extra_msg = ' on line {}'.format(node.line)
        else:
            extra_msg = ' (possibly by an import)'
        self.fail('{} "{}" already defined{}'.format(noun, unmangle(name), extra_msg), ctx,
                  code=codes.NO_REDEF)

    def name_already_defined(self,
                             name: str,
                             ctx: Context,
                             original_ctx: Optional[Union[SymbolTableNode, SymbolNode]] = None
                             ) -> None:
        self.already_defined(name, ctx, original_ctx, noun='Name')

    def attribute_already_defined(self,
                                  name: str,
                                  ctx: Context,
                                  original_ctx: Optional[Union[SymbolTableNode, SymbolNode]] = None
                                  ) -> None:
        self.already_defined(name, ctx, original_ctx, noun='Attribute')

    def is_local_name(self, name: str) -> bool:
        """Does name look like reference to a definition in the current module?"""
        return self.is_defined_in_current_module(name) or '.' not in name

    def in_checked_function(self) -> bool:
        """Should we type-check the current function?

        - Yes if --check-untyped-defs is set.
        - Yes outside functions.
        - Yes in annotated functions.
        - No otherwise.
        """
        if self.options.check_untyped_defs or not self.function_stack:
            return True

        current_index = len(self.function_stack) - 1
        while current_index >= 0:
            current_func = self.function_stack[current_index]
            if (
                isinstance(current_func, FuncItem)
                and not isinstance(current_func, LambdaExpr)
            ):
                return not current_func.is_dynamic()

            # Special case, `lambda` inherits the "checked" state from its parent.
            # Because `lambda` itself cannot be annotated.
            # `lambdas` can be deeply nested, so we try to find at least one other parent.
            current_index -= 1

        # This means that we only have a stack of `lambda` functions,
        # no regular functions.
        return True

    def fail(self,
             msg: str,
             ctx: Context,
             serious: bool = False,
             *,
             code: Optional[ErrorCode] = None,
             blocker: bool = False) -> None:
        if not serious and not self.in_checked_function():
            return
        # In case it's a bug and we don't really have context
        assert ctx is not None, msg
        self.errors.report(ctx.get_line(), ctx.get_column(), msg, blocker=blocker, code=code)

    def fail_blocker(self, msg: str, ctx: Context) -> None:
        self.fail(msg, ctx, blocker=True)

    def note(self, msg: str, ctx: Context, code: Optional[ErrorCode] = None) -> None:
        if not self.in_checked_function():
            return
        self.errors.report(ctx.get_line(), ctx.get_column(), msg, severity='note', code=code)

    def accept(self, node: Node) -> None:
        try:
            node.accept(self)
        except Exception as err:
            report_internal_error(err, self.errors.file, node.line, self.errors, self.options)

    def expr_to_analyzed_type(self,
                              expr: Expression,
                              report_invalid_types: bool = True,
                              allow_placeholder: bool = False) -> Optional[Type]:
        if isinstance(expr, CallExpr):
            expr.accept(self)
            internal_name, info = self.named_tuple_analyzer.check_namedtuple(expr, None,
                                                                             self.is_func_scope())
            if internal_name is None:
                # Some form of namedtuple is the only valid type that looks like a call
                # expression. This isn't a valid type.
                raise TypeTranslationError()
            elif not info:
                self.defer(expr)
                return None
            assert info.tuple_type, "NamedTuple without tuple type"
            fallback = Instance(info, [])
            return TupleType(info.tuple_type.items, fallback=fallback)
        typ = self.expr_to_unanalyzed_type(expr)
        return self.anal_type(typ, report_invalid_types=report_invalid_types,
                              allow_placeholder=allow_placeholder)

    def analyze_type_expr(self, expr: Expression) -> None:
        # There are certain expressions that mypy does not need to semantically analyze,
        # since they analyzed solely as type. (For example, indexes in type alias definitions
        # and base classes in class defs). External consumers of the mypy AST may need
        # them semantically analyzed, however, if they need to treat it as an expression
        # and not a type. (Which is to say, mypyc needs to do this.) Do the analysis
        # in a fresh tvar scope in order to suppress any errors about using type variables.
        with self.tvar_scope_frame(TypeVarLikeScope()):
            expr.accept(self)

    def type_analyzer(self, *,
                      tvar_scope: Optional[TypeVarLikeScope] = None,
                      allow_tuple_literal: bool = False,
                      allow_unbound_tvars: bool = False,
                      allow_placeholder: bool = False,
                      allow_required: bool = False,
                      report_invalid_types: bool = True) -> TypeAnalyser:
        if tvar_scope is None:
            tvar_scope = self.tvar_scope
        tpan = TypeAnalyser(self,
                            tvar_scope,
                            self.plugin,
                            self.options,
                            self.is_typeshed_stub_file,
                            allow_unbound_tvars=allow_unbound_tvars,
                            allow_tuple_literal=allow_tuple_literal,
                            report_invalid_types=report_invalid_types,
                            allow_placeholder=allow_placeholder,
                            allow_required=allow_required)
        tpan.in_dynamic_func = bool(self.function_stack and self.function_stack[-1].is_dynamic())
        tpan.global_scope = not self.type and not self.function_stack
        return tpan

    def expr_to_unanalyzed_type(self, node: Expression) -> ProperType:
        return expr_to_unanalyzed_type(node, self.options, self.is_stub_file)

    def anal_type(self,
                  typ: Type, *,
                  tvar_scope: Optional[TypeVarLikeScope] = None,
                  allow_tuple_literal: bool = False,
                  allow_unbound_tvars: bool = False,
                  allow_placeholder: bool = False,
                  allow_required: bool = False,
                  report_invalid_types: bool = True,
                  third_pass: bool = False) -> Optional[Type]:
        """Semantically analyze a type.

        Args:
            typ: Type to analyze (if already analyzed, this is a no-op)
            allow_placeholder: If True, may return PlaceholderType if
                encountering an incomplete definition
            third_pass: Unused; only for compatibility with old semantic
                analyzer

        Return None only if some part of the type couldn't be bound *and* it
        referred to an incomplete namespace or definition. In this case also
        defer as needed. During a final iteration this won't return None;
        instead report an error if the type can't be analyzed and return
        AnyType.

        In case of other errors, report an error message and return AnyType.

        NOTE: The caller shouldn't defer even if this returns None or a
              placeholder type.
        """
        a = self.type_analyzer(tvar_scope=tvar_scope,
                               allow_unbound_tvars=allow_unbound_tvars,
                               allow_tuple_literal=allow_tuple_literal,
                               allow_placeholder=allow_placeholder,
                               allow_required=allow_required,
                               report_invalid_types=report_invalid_types)
        tag = self.track_incomplete_refs()
        typ = typ.accept(a)
        if self.found_incomplete_ref(tag):
            # Something could not be bound yet.
            return None
        self.add_type_alias_deps(a.aliases_used)
        return typ

    def class_type(self, self_type: Type) -> Type:
        return TypeType.make_normalized(self_type)

    def schedule_patch(self, priority: int, patch: Callable[[], None]) -> None:
        self.patches.append((priority, patch))

    def report_hang(self) -> None:
        print('Deferral trace:')
        for mod, line in self.deferral_debug_context:
            print('    {}:{}'.format(mod, line))
        self.errors.report(-1, -1,
                           'INTERNAL ERROR: maximum semantic analysis iteration count reached',
                           blocker=True)

    def add_plugin_dependency(self, trigger: str, target: Optional[str] = None) -> None:
        """Add dependency from trigger to a target.

        If the target is not given explicitly, use the current target.
        """
        if target is None:
            target = self.scope.current_target()
        self.cur_mod_node.plugin_deps.setdefault(trigger, set()).add(target)

    def add_type_alias_deps(self,
                            aliases_used: Iterable[str],
                            target: Optional[str] = None) -> None:
        """Add full names of type aliases on which the current node depends.

        This is used by fine-grained incremental mode to re-check the corresponding nodes.
        If `target` is None, then the target node used will be the current scope.
        """
        if not aliases_used:
            # A basic optimization to avoid adding targets with no dependencies to
            # the `alias_deps` dict.
            return
        if target is None:
            target = self.scope.current_target()
        self.cur_mod_node.alias_deps[target].update(aliases_used)

    def is_mangled_global(self, name: str) -> bool:
        # A global is mangled if there exists at least one renamed variant.
        return unmangle(name) + "'" in self.globals

    def is_initial_mangled_global(self, name: str) -> bool:
        # If there are renamed definitions for a global, the first one has exactly one prime.
        return name == unmangle(name) + "'"

    def parse_bool(self, expr: Expression) -> Optional[bool]:
        if isinstance(expr, NameExpr):
            if expr.fullname == 'builtins.True':
                return True
            if expr.fullname == 'builtins.False':
                return False
        return None

    def set_future_import_flags(self, module_name: str) -> None:
        if module_name in FUTURE_IMPORTS:
            self.future_import_flags.add(FUTURE_IMPORTS[module_name])

    def is_future_flag_set(self, flag: str) -> bool:
        return flag in self.future_import_flags


class HasPlaceholders(TypeQuery[bool]):
    def __init__(self) -> None:
        super().__init__(any)

    def visit_placeholder_type(self, t: PlaceholderType) -> bool:
        return True


def has_placeholder(typ: Type) -> bool:
    """Check if a type contains any placeholder types (recursively)."""
    return typ.accept(HasPlaceholders())


def replace_implicit_first_type(sig: FunctionLike, new: Type) -> FunctionLike:
    if isinstance(sig, CallableType):
        if len(sig.arg_types) == 0:
            return sig
        return sig.copy_modified(arg_types=[new] + sig.arg_types[1:])
    elif isinstance(sig, Overloaded):
        return Overloaded([cast(CallableType, replace_implicit_first_type(i, new))
                           for i in sig.items])
    else:
        assert False


def refers_to_fullname(node: Expression, fullname: str) -> bool:
    """Is node a name or member expression with the given full name?"""
    if not isinstance(node, RefExpr):
        return False
    if node.fullname == fullname:
        return True
    if isinstance(node.node, TypeAlias):
        target = get_proper_type(node.node.target)
        if isinstance(target, Instance) and target.type.fullname == fullname:
            return True
    return False


def refers_to_class_or_function(node: Expression) -> bool:
    """Does semantically analyzed node refer to a class?"""
    return (isinstance(node, RefExpr) and
            isinstance(node.node, (TypeInfo, FuncDef, OverloadedFuncDef)))


def find_duplicate(list: List[T]) -> Optional[T]:
    """If the list has duplicates, return one of the duplicates.

    Otherwise, return None.
    """
    for i in range(1, len(list)):
        if list[i] in list[:i]:
            return list[i]
    return None


def remove_imported_names_from_symtable(names: SymbolTable,
                                        module: str) -> None:
    """Remove all imported names from the symbol table of a module."""
    removed: List[str] = []
    for name, node in names.items():
        if node.node is None:
            continue
        fullname = node.node.fullname
        prefix = fullname[:fullname.rfind('.')]
        if prefix != module:
            removed.append(name)
    for name in removed:
        del names[name]


def make_any_non_explicit(t: Type) -> Type:
    """Replace all Any types within in with Any that has attribute 'explicit' set to False"""
    return t.accept(MakeAnyNonExplicit())


class MakeAnyNonExplicit(TypeTranslator):
    def visit_any(self, t: AnyType) -> Type:
        if t.type_of_any == TypeOfAny.explicit:
            return t.copy_modified(TypeOfAny.special_form)
        return t

    def visit_type_alias_type(self, t: TypeAliasType) -> Type:
        return t.copy_modified(args=[a.accept(self) for a in t.args])


def apply_semantic_analyzer_patches(patches: List[Tuple[int, Callable[[], None]]]) -> None:
    """Call patch callbacks in the right order.

    This should happen after semantic analyzer pass 3.
    """
    patches_by_priority = sorted(patches, key=lambda x: x[0])
    for priority, patch_func in patches_by_priority:
        patch_func()


def names_modified_by_assignment(s: AssignmentStmt) -> List[NameExpr]:
    """Return all unqualified (short) names assigned to in an assignment statement."""
    result: List[NameExpr] = []
    for lvalue in s.lvalues:
        result += names_modified_in_lvalue(lvalue)
    return result


def names_modified_in_lvalue(lvalue: Lvalue) -> List[NameExpr]:
    """Return all NameExpr assignment targets in an Lvalue."""
    if isinstance(lvalue, NameExpr):
        return [lvalue]
    elif isinstance(lvalue, StarExpr):
        return names_modified_in_lvalue(lvalue.expr)
    elif isinstance(lvalue, (ListExpr, TupleExpr)):
        result: List[NameExpr] = []
        for item in lvalue.items:
            result += names_modified_in_lvalue(item)
        return result
    return []


def is_same_var_from_getattr(n1: Optional[SymbolNode], n2: Optional[SymbolNode]) -> bool:
    """Do n1 and n2 refer to the same Var derived from module-level __getattr__?"""
    return (isinstance(n1, Var)
            and n1.from_module_getattr
            and isinstance(n2, Var)
            and n2.from_module_getattr
            and n1.fullname == n2.fullname)


def dummy_context() -> Context:
    return TempNode(AnyType(TypeOfAny.special_form))


def is_valid_replacement(old: SymbolTableNode, new: SymbolTableNode) -> bool:
    """Can symbol table node replace an existing one?

    These are the only valid cases:

    1. Placeholder gets replaced with a non-placeholder
    2. Placeholder that isn't known to become type replaced with a
       placeholder that can become a type
    """
    if isinstance(old.node, PlaceholderNode):
        if isinstance(new.node, PlaceholderNode):
            return not old.node.becomes_typeinfo and new.node.becomes_typeinfo
        else:
            return True
    return False


def is_same_symbol(a: Optional[SymbolNode], b: Optional[SymbolNode]) -> bool:
    return (a == b
            or (isinstance(a, PlaceholderNode)
                and isinstance(b, PlaceholderNode))
            or is_same_var_from_getattr(a, b))
