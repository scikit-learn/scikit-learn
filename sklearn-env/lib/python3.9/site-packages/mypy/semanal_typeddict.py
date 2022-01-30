"""Semantic analysis of TypedDict definitions."""

from mypy.backports import OrderedDict
from typing import Optional, List, Set, Tuple
from typing_extensions import Final

from mypy.types import (
    Type, AnyType, TypeOfAny, TypedDictType, TPDICT_NAMES, RequiredType,
)
from mypy.nodes import (
    CallExpr, TypedDictExpr, Expression, NameExpr, Context, StrExpr, BytesExpr, UnicodeExpr,
    ClassDef, RefExpr, TypeInfo, AssignmentStmt, PassStmt, ExpressionStmt, EllipsisExpr, TempNode,
    DictExpr, ARG_POS, ARG_NAMED
)
from mypy.semanal_shared import SemanticAnalyzerInterface
from mypy.exprtotype import expr_to_unanalyzed_type, TypeTranslationError
from mypy.options import Options
from mypy.typeanal import check_for_explicit_any, has_any_from_unimported_type
from mypy.messages import MessageBuilder
from mypy.errorcodes import ErrorCode
from mypy import errorcodes as codes

TPDICT_CLASS_ERROR: Final = (
    "Invalid statement in TypedDict definition; " 'expected "field_name: field_type"'
)


class TypedDictAnalyzer:
    def __init__(self,
                 options: Options,
                 api: SemanticAnalyzerInterface,
                 msg: MessageBuilder) -> None:
        self.options = options
        self.api = api
        self.msg = msg

    def analyze_typeddict_classdef(self, defn: ClassDef) -> Tuple[bool, Optional[TypeInfo]]:
        """Analyze a class that may define a TypedDict.

        Assume that base classes have been analyzed already.

        Note: Unlike normal classes, we won't create a TypeInfo until
        the whole definition of the TypeDict (including the body and all
        key names and types) is complete.  This is mostly because we
        store the corresponding TypedDictType in the TypeInfo.

        Return (is this a TypedDict, new TypeInfo). Specifics:
         * If we couldn't finish due to incomplete reference anywhere in
           the definition, return (True, None).
         * If this is not a TypedDict, return (False, None).
        """
        possible = False
        for base_expr in defn.base_type_exprs:
            if isinstance(base_expr, RefExpr):
                self.api.accept(base_expr)
                if base_expr.fullname in TPDICT_NAMES or self.is_typeddict(base_expr):
                    possible = True
        if possible:
            if (len(defn.base_type_exprs) == 1 and
                    isinstance(defn.base_type_exprs[0], RefExpr) and
                    defn.base_type_exprs[0].fullname in TPDICT_NAMES):
                # Building a new TypedDict
                fields, types, required_keys = self.analyze_typeddict_classdef_fields(defn)
                if fields is None:
                    return True, None  # Defer
                info = self.build_typeddict_typeinfo(defn.name, fields, types, required_keys,
                                                     defn.line)
                defn.analyzed = TypedDictExpr(info)
                defn.analyzed.line = defn.line
                defn.analyzed.column = defn.column
                return True, info

            # Extending/merging existing TypedDicts
            typeddict_bases = []
            typeddict_bases_set = set()
            for expr in defn.base_type_exprs:
                if isinstance(expr, RefExpr) and expr.fullname in TPDICT_NAMES:
                    if 'TypedDict' not in typeddict_bases_set:
                        typeddict_bases_set.add('TypedDict')
                    else:
                        self.fail('Duplicate base class "TypedDict"', defn)
                elif isinstance(expr, RefExpr) and self.is_typeddict(expr):
                    assert expr.fullname
                    if expr.fullname not in typeddict_bases_set:
                        typeddict_bases_set.add(expr.fullname)
                        typeddict_bases.append(expr)
                    else:
                        assert isinstance(expr.node, TypeInfo)
                        self.fail('Duplicate base class "%s"' % expr.node.name, defn)
                else:
                    self.fail("All bases of a new TypedDict must be TypedDict types", defn)

            keys: List[str] = []
            types = []
            required_keys = set()
            # Iterate over bases in reverse order so that leftmost base class' keys take precedence
            for base in reversed(typeddict_bases):
                assert isinstance(base, RefExpr)
                assert isinstance(base.node, TypeInfo)
                assert isinstance(base.node.typeddict_type, TypedDictType)
                base_typed_dict = base.node.typeddict_type
                base_items = base_typed_dict.items
                valid_items = base_items.copy()
                for key in base_items:
                    if key in keys:
                        self.fail('Overwriting TypedDict field "{}" while merging'
                                  .format(key), defn)
                keys.extend(valid_items.keys())
                types.extend(valid_items.values())
                required_keys.update(base_typed_dict.required_keys)
            new_keys, new_types, new_required_keys = self.analyze_typeddict_classdef_fields(defn,
                                                                                            keys)
            if new_keys is None:
                return True, None  # Defer
            keys.extend(new_keys)
            types.extend(new_types)
            required_keys.update(new_required_keys)
            info = self.build_typeddict_typeinfo(defn.name, keys, types, required_keys, defn.line)
            defn.analyzed = TypedDictExpr(info)
            defn.analyzed.line = defn.line
            defn.analyzed.column = defn.column
            return True, info
        return False, None

    def analyze_typeddict_classdef_fields(
            self,
            defn: ClassDef,
            oldfields: Optional[List[str]] = None) -> Tuple[Optional[List[str]],
                                                            List[Type],
                                                            Set[str]]:
        """Analyze fields defined in a TypedDict class definition.

        This doesn't consider inherited fields (if any). Also consider totality,
        if given.

        Return tuple with these items:
         * List of keys (or None if found an incomplete reference --> deferral)
         * List of types for each key
         * Set of required keys
        """
        fields: List[str] = []
        types: List[Type] = []
        for stmt in defn.defs.body:
            if not isinstance(stmt, AssignmentStmt):
                # Still allow pass or ... (for empty TypedDict's).
                if (not isinstance(stmt, PassStmt) and
                    not (isinstance(stmt, ExpressionStmt) and
                         isinstance(stmt.expr, (EllipsisExpr, StrExpr)))):
                    self.fail(TPDICT_CLASS_ERROR, stmt)
            elif len(stmt.lvalues) > 1 or not isinstance(stmt.lvalues[0], NameExpr):
                # An assignment, but an invalid one.
                self.fail(TPDICT_CLASS_ERROR, stmt)
            else:
                name = stmt.lvalues[0].name
                if name in (oldfields or []):
                    self.fail('Overwriting TypedDict field "{}" while extending'
                              .format(name), stmt)
                if name in fields:
                    self.fail('Duplicate TypedDict key "{}"'.format(name), stmt)
                    continue
                # Append name and type in this case...
                fields.append(name)
                if stmt.type is None:
                    types.append(AnyType(TypeOfAny.unannotated))
                else:
                    analyzed = self.api.anal_type(stmt.type, allow_required=True)
                    if analyzed is None:
                        return None, [], set()  # Need to defer
                    types.append(analyzed)
                # ...despite possible minor failures that allow further analyzis.
                if stmt.type is None or hasattr(stmt, 'new_syntax') and not stmt.new_syntax:
                    self.fail(TPDICT_CLASS_ERROR, stmt)
                elif not isinstance(stmt.rvalue, TempNode):
                    # x: int assigns rvalue to TempNode(AnyType())
                    self.fail('Right hand side values are not supported in TypedDict', stmt)
        total: Optional[bool] = True
        if 'total' in defn.keywords:
            total = self.api.parse_bool(defn.keywords['total'])
            if total is None:
                self.fail('Value of "total" must be True or False', defn)
                total = True
        required_keys = {
            field
            for (field, t) in zip(fields, types)
            if (total or (
                isinstance(t, RequiredType) and  # type: ignore[misc]
                t.required
            )) and not (
                isinstance(t, RequiredType) and  # type: ignore[misc]
                not t.required
            )
        }
        types = [  # unwrap Required[T] to just T
            t.item if isinstance(t, RequiredType) else t  # type: ignore[misc]
            for t in types
        ]

        return fields, types, required_keys

    def check_typeddict(self,
                        node: Expression,
                        var_name: Optional[str],
                        is_func_scope: bool) -> Tuple[bool, Optional[TypeInfo]]:
        """Check if a call defines a TypedDict.

        The optional var_name argument is the name of the variable to
        which this is assigned, if any.

        Return a pair (is it a typed dict, corresponding TypeInfo).

        If the definition is invalid but looks like a TypedDict,
        report errors but return (some) TypeInfo. If some type is not ready,
        return (True, None).
        """
        if not isinstance(node, CallExpr):
            return False, None
        call = node
        callee = call.callee
        if not isinstance(callee, RefExpr):
            return False, None
        fullname = callee.fullname
        if fullname not in TPDICT_NAMES:
            return False, None
        res = self.parse_typeddict_args(call)
        if res is None:
            # This is a valid typed dict, but some type is not ready.
            # The caller should defer this until next iteration.
            return True, None
        name, items, types, total, ok = res
        if not ok:
            # Error. Construct dummy return value.
            info = self.build_typeddict_typeinfo('TypedDict', [], [], set(), call.line)
        else:
            if var_name is not None and name != var_name:
                self.fail(
                    'First argument "{}" to TypedDict() does not match variable name "{}"'.format(
                        name, var_name), node, code=codes.NAME_MATCH)
            if name != var_name or is_func_scope:
                # Give it a unique name derived from the line number.
                name += '@' + str(call.line)
            required_keys = {
                field
                for (field, t) in zip(items, types)
                if (total or (
                    isinstance(t, RequiredType) and  # type: ignore[misc]
                    t.required
                )) and not (
                    isinstance(t, RequiredType) and  # type: ignore[misc]
                    not t.required
                )
            }
            types = [  # unwrap Required[T] to just T
                t.item if isinstance(t, RequiredType) else t  # type: ignore[misc]
                for t in types
            ]
            info = self.build_typeddict_typeinfo(name, items, types, required_keys, call.line)
            info.line = node.line
            # Store generated TypeInfo under both names, see semanal_namedtuple for more details.
            if name != var_name or is_func_scope:
                self.api.add_symbol_skip_local(name, info)
        if var_name:
            self.api.add_symbol(var_name, info, node)
        call.analyzed = TypedDictExpr(info)
        call.analyzed.set_line(call.line, call.column)
        return True, info

    def parse_typeddict_args(
            self, call: CallExpr) -> Optional[Tuple[str, List[str], List[Type], bool, bool]]:
        """Parse typed dict call expression.

        Return names, types, totality, was there an error during parsing.
        If some type is not ready, return None.
        """
        # TODO: Share code with check_argument_count in checkexpr.py?
        args = call.args
        if len(args) < 2:
            return self.fail_typeddict_arg("Too few arguments for TypedDict()", call)
        if len(args) > 3:
            return self.fail_typeddict_arg("Too many arguments for TypedDict()", call)
        # TODO: Support keyword arguments
        if call.arg_kinds not in ([ARG_POS, ARG_POS], [ARG_POS, ARG_POS, ARG_NAMED]):
            return self.fail_typeddict_arg("Unexpected arguments to TypedDict()", call)
        if len(args) == 3 and call.arg_names[2] != 'total':
            return self.fail_typeddict_arg(
                'Unexpected keyword argument "{}" for "TypedDict"'.format(call.arg_names[2]), call)
        if not isinstance(args[0], (StrExpr, BytesExpr, UnicodeExpr)):
            return self.fail_typeddict_arg(
                "TypedDict() expects a string literal as the first argument", call)
        if not isinstance(args[1], DictExpr):
            return self.fail_typeddict_arg(
                "TypedDict() expects a dictionary literal as the second argument", call)
        total: Optional[bool] = True
        if len(args) == 3:
            total = self.api.parse_bool(call.args[2])
            if total is None:
                return self.fail_typeddict_arg(
                    'TypedDict() "total" argument must be True or False', call)
        dictexpr = args[1]
        res = self.parse_typeddict_fields_with_types(dictexpr.items, call)
        if res is None:
            # One of the types is not ready, defer.
            return None
        items, types, ok = res
        for t in types:
            check_for_explicit_any(t, self.options, self.api.is_typeshed_stub_file, self.msg,
                                   context=call)

        if self.options.disallow_any_unimported:
            for t in types:
                if has_any_from_unimported_type(t):
                    self.msg.unimported_type_becomes_any("Type of a TypedDict key", t, dictexpr)
        assert total is not None
        return args[0].value, items, types, total, ok

    def parse_typeddict_fields_with_types(
            self,
            dict_items: List[Tuple[Optional[Expression], Expression]],
            context: Context) -> Optional[Tuple[List[str], List[Type], bool]]:
        """Parse typed dict items passed as pairs (name expression, type expression).

        Return names, types, was there an error. If some type is not ready, return None.
        """
        seen_keys = set()
        items: List[str] = []
        types: List[Type] = []
        for (field_name_expr, field_type_expr) in dict_items:
            if isinstance(field_name_expr, (StrExpr, BytesExpr, UnicodeExpr)):
                key = field_name_expr.value
                items.append(key)
                if key in seen_keys:
                    self.fail('Duplicate TypedDict key "{}"'.format(key), field_name_expr)
                seen_keys.add(key)
            else:
                name_context = field_name_expr or field_type_expr
                self.fail_typeddict_arg("Invalid TypedDict() field name", name_context)
                return [], [], False
            try:
                type = expr_to_unanalyzed_type(field_type_expr, self.options,
                                               self.api.is_stub_file)
            except TypeTranslationError:
                if (isinstance(field_type_expr, CallExpr) and
                        isinstance(field_type_expr.callee, RefExpr) and
                        field_type_expr.callee.fullname in TPDICT_NAMES):
                    self.fail_typeddict_arg(
                        'Inline TypedDict types not supported; use assignment to define TypedDict',
                        field_type_expr)
                else:
                    self.fail_typeddict_arg('Invalid field type', field_type_expr)
                return [], [], False
            analyzed = self.api.anal_type(type, allow_required=True)
            if analyzed is None:
                return None
            types.append(analyzed)
        return items, types, True

    def fail_typeddict_arg(self, message: str,
                           context: Context) -> Tuple[str, List[str], List[Type], bool, bool]:
        self.fail(message, context)
        return '', [], [], True, False

    def build_typeddict_typeinfo(self, name: str, items: List[str],
                                 types: List[Type],
                                 required_keys: Set[str],
                                 line: int) -> TypeInfo:
        # Prefer typing then typing_extensions if available.
        fallback = (self.api.named_type_or_none('typing._TypedDict', []) or
                    self.api.named_type_or_none('typing_extensions._TypedDict', []) or
                    self.api.named_type_or_none('mypy_extensions._TypedDict', []))
        assert fallback is not None
        info = self.api.basic_new_typeinfo(name, fallback, line)
        info.typeddict_type = TypedDictType(OrderedDict(zip(items, types)), required_keys,
                                            fallback)
        return info

    # Helpers

    def is_typeddict(self, expr: Expression) -> bool:
        return (isinstance(expr, RefExpr) and isinstance(expr.node, TypeInfo) and
                expr.node.typeddict_type is not None)

    def fail(self, msg: str, ctx: Context, *, code: Optional[ErrorCode] = None) -> None:
        self.api.fail(msg, ctx, code=code)

    def note(self, msg: str, ctx: Context) -> None:
        self.api.note(msg, ctx)
