# SPDX-License-Identifier: Apache-2.0
# Copyright 2016-2017 The Meson development team

# This class contains the basic functionality needed to run any interpreter
# or an interpreter-based tool.
from __future__ import annotations

from .. import environment, mparser, mesonlib

from .baseobjects import (
    InterpreterObject,
    MesonInterpreterObject,
    MutableInterpreterObject,
    ObjectHolder,
    IterableObject,
    ContextManagerObject,

    HoldableTypes,
)

from .exceptions import (
    BreakRequest,
    ContinueRequest,
    InterpreterException,
    InvalidArguments,
    InvalidCode,
    SubdirDoneRequest,
)

from .decorators import FeatureNew
from .disabler import Disabler, is_disabled
from .helpers import default_resolve_key, flatten, resolve_second_level_holders, stringifyUserArguments
from .operator import MesonOperator
from ._unholder import _unholder

import os, copy, re, pathlib
import typing as T
import textwrap

if T.TYPE_CHECKING:
    from .baseobjects import InterpreterObjectTypeVar, SubProject, TYPE_kwargs, TYPE_var
    from ..interpreter import Interpreter

    HolderMapType = T.Dict[
        T.Union[
            T.Type[mesonlib.HoldableObject],
            T.Type[int],
            T.Type[bool],
            T.Type[str],
            T.Type[list],
            T.Type[dict],
        ],
        # For some reason, this has to be a callable and can't just be ObjectHolder[InterpreterObjectTypeVar]
        T.Callable[[InterpreterObjectTypeVar, 'Interpreter'], ObjectHolder[InterpreterObjectTypeVar]]
    ]

    FunctionType = T.Dict[
        str,
        T.Callable[[mparser.BaseNode, T.List[TYPE_var], T.Dict[str, TYPE_var]], TYPE_var]
    ]


class InvalidCodeOnVoid(InvalidCode):

    def __init__(self, op_type: str) -> None:
        super().__init__(f'Cannot perform {op_type!r} operation on void statement.')


class InterpreterBase:
    def __init__(self, source_root: str, subdir: str, subproject: 'SubProject'):
        self.source_root = source_root
        self.funcs: FunctionType = {}
        self.builtin: T.Dict[str, InterpreterObject] = {}
        # Holder maps store a mapping from an HoldableObject to a class ObjectHolder
        self.holder_map: HolderMapType = {}
        self.bound_holder_map: HolderMapType = {}
        self.subdir = subdir
        self.root_subdir = subdir
        self.subproject = subproject
        self.variables: T.Dict[str, InterpreterObject] = {}
        self.argument_depth = 0
        self.current_lineno = -1
        # Current node set during a function call. This can be used as location
        # when printing a warning message during a method call.
        self.current_node = mparser.BaseNode(-1, -1, 'sentinel')
        # This is set to `version_string` when this statement is evaluated:
        # meson.version().compare_version(version_string)
        # If it was part of a if-clause, it is used to temporally override the
        # current meson version target within that if-block.
        self.tmp_meson_version: T.Optional[str] = None

    def handle_meson_version_from_ast(self, strict: bool = True) -> None:
        # do nothing in an AST interpreter
        return

    def read_buildfile(self, fname: str, errname: str) -> str:
        try:
            with open(fname, encoding='utf-8') as f:
                return f.read()
        except UnicodeDecodeError as e:
            node = mparser.BaseNode(1, 1, errname)
            raise InvalidCode.from_node(f'Build file failed to parse as unicode: {e}', node=node)

    def load_root_meson_file(self) -> None:
        mesonfile = os.path.join(self.source_root, self.subdir, environment.build_filename)
        if not os.path.isfile(mesonfile):
            raise InvalidArguments(f'Missing Meson file in {mesonfile}')
        code = self.read_buildfile(mesonfile, mesonfile)
        if code.isspace():
            raise InvalidCode('Builder file is empty.')
        assert isinstance(code, str)
        try:
            self.ast = mparser.Parser(code, mesonfile).parse()
            self.handle_meson_version_from_ast()
        except mparser.ParseException as me:
            me.file = mesonfile
            if me.ast:
                # try to detect parser errors from new syntax added by future
                # meson versions, and just tell the user to update meson
                self.ast = me.ast
                self.handle_meson_version_from_ast()
            raise me

    def parse_project(self) -> None:
        """
        Parses project() and initializes languages, compilers etc. Do this
        early because we need this before we parse the rest of the AST.
        """
        self.evaluate_codeblock(self.ast, end=1)

    def sanity_check_ast(self) -> None:
        def _is_project(ast: mparser.CodeBlockNode) -> object:
            if not isinstance(ast, mparser.CodeBlockNode):
                raise InvalidCode('AST is of invalid type. Possibly a bug in the parser.')
            if not ast.lines:
                raise InvalidCode('No statements in code.')
            first = ast.lines[0]
            return isinstance(first, mparser.FunctionNode) and first.func_name.value == 'project'

        if not _is_project(self.ast):
            p = pathlib.Path(self.source_root).resolve()
            found = p
            for parent in p.parents:
                if (parent / 'meson.build').is_file():
                    with open(parent / 'meson.build', encoding='utf-8') as f:
                        code = f.read()

                    try:
                        ast = mparser.Parser(code, 'empty').parse()
                    except mparser.ParseException:
                        continue

                    if _is_project(ast):
                        found = parent
                        break
                else:
                    break

            error = 'first statement must be a call to project()'
            if found != p:
                raise InvalidCode(f'Not the project root: {error}\n\nDid you mean to run meson from the directory: "{found}"?')
            else:
                raise InvalidCode(f'Invalid source tree: {error}')

    def run(self) -> None:
        # Evaluate everything after the first line, which is project() because
        # we already parsed that in self.parse_project()
        try:
            self.evaluate_codeblock(self.ast, start=1)
        except SubdirDoneRequest:
            pass

    def evaluate_codeblock(self, node: mparser.CodeBlockNode, start: int = 0, end: T.Optional[int] = None) -> None:
        if node is None:
            return
        if not isinstance(node, mparser.CodeBlockNode):
            e = InvalidCode('Tried to execute a non-codeblock. Possibly a bug in the parser.')
            e.lineno = node.lineno
            e.colno = node.colno
            raise e
        statements = node.lines[start:end]
        i = 0
        while i < len(statements):
            cur = statements[i]
            try:
                self.evaluate_statement(cur)
            except Exception as e:
                if getattr(e, 'lineno', None) is None:
                    # We are doing the equivalent to setattr here and mypy does not like it
                    # NOTE: self.current_node is continually updated during processing
                    e.lineno = self.current_node.lineno                                               # type: ignore
                    e.colno = self.current_node.colno                                                 # type: ignore
                    e.file = os.path.join(self.source_root, self.subdir, environment.build_filename)  # type: ignore
                raise e
            i += 1 # In THE FUTURE jump over blocks and stuff.

    def evaluate_statement(self, cur: mparser.BaseNode) -> T.Optional[InterpreterObject]:
        self.current_node = cur
        if isinstance(cur, mparser.FunctionNode):
            return self.function_call(cur)
        elif isinstance(cur, mparser.PlusAssignmentNode):
            self.evaluate_plusassign(cur)
        elif isinstance(cur, mparser.AssignmentNode):
            self.assignment(cur)
        elif isinstance(cur, mparser.MethodNode):
            return self.method_call(cur)
        elif isinstance(cur, mparser.StringNode):
            if cur.is_fstring:
                if cur.is_multiline:
                    return self.evaluate_multiline_fstring(cur)
                else:
                    return self.evaluate_fstring(cur)
            else:
                return self._holderify(cur.value)
        elif isinstance(cur, mparser.BooleanNode):
            return self._holderify(cur.value)
        elif isinstance(cur, mparser.IfClauseNode):
            return self.evaluate_if(cur)
        elif isinstance(cur, mparser.IdNode):
            return self.get_variable(cur.value)
        elif isinstance(cur, mparser.ComparisonNode):
            return self.evaluate_comparison(cur)
        elif isinstance(cur, mparser.ArrayNode):
            return self.evaluate_arraystatement(cur)
        elif isinstance(cur, mparser.DictNode):
            return self.evaluate_dictstatement(cur)
        elif isinstance(cur, mparser.NumberNode):
            return self._holderify(cur.value)
        elif isinstance(cur, mparser.AndNode):
            return self.evaluate_andstatement(cur)
        elif isinstance(cur, mparser.OrNode):
            return self.evaluate_orstatement(cur)
        elif isinstance(cur, mparser.NotNode):
            return self.evaluate_notstatement(cur)
        elif isinstance(cur, mparser.UMinusNode):
            return self.evaluate_uminusstatement(cur)
        elif isinstance(cur, mparser.ArithmeticNode):
            return self.evaluate_arithmeticstatement(cur)
        elif isinstance(cur, mparser.ForeachClauseNode):
            self.evaluate_foreach(cur)
        elif isinstance(cur, mparser.IndexNode):
            return self.evaluate_indexing(cur)
        elif isinstance(cur, mparser.TernaryNode):
            return self.evaluate_ternary(cur)
        elif isinstance(cur, mparser.ContinueNode):
            raise ContinueRequest()
        elif isinstance(cur, mparser.BreakNode):
            raise BreakRequest()
        elif isinstance(cur, mparser.ParenthesizedNode):
            return self.evaluate_statement(cur.inner)
        elif isinstance(cur, mparser.TestCaseClauseNode):
            return self.evaluate_testcase(cur)
        else:
            raise InvalidCode("Unknown statement.")
        return None

    def evaluate_arraystatement(self, cur: mparser.ArrayNode) -> InterpreterObject:
        (arguments, kwargs) = self.reduce_arguments(cur.args)
        if len(kwargs) > 0:
            raise InvalidCode('Keyword arguments are invalid in array construction.')
        return self._holderify([_unholder(x) for x in arguments])

    @FeatureNew('dict', '0.47.0')
    def evaluate_dictstatement(self, cur: mparser.DictNode) -> InterpreterObject:
        def resolve_key(key: mparser.BaseNode) -> str:
            if not isinstance(key, mparser.StringNode):
                FeatureNew.single_use('Dictionary entry using non literal key', '0.53.0', self.subproject)
            key_holder = self.evaluate_statement(key)
            if key_holder is None:
                raise InvalidArguments('Key cannot be void.')
            str_key = _unholder(key_holder)
            if not isinstance(str_key, str):
                raise InvalidArguments('Key must be a string')
            return str_key
        arguments, kwargs = self.reduce_arguments(cur.args, key_resolver=resolve_key, duplicate_key_error='Duplicate dictionary key: {}')
        assert not arguments
        return self._holderify({k: _unholder(v) for k, v in kwargs.items()})

    def evaluate_notstatement(self, cur: mparser.NotNode) -> InterpreterObject:
        v = self.evaluate_statement(cur.value)
        if v is None:
            raise InvalidCodeOnVoid('not')
        if isinstance(v, Disabler):
            return v
        return self._holderify(v.operator_call(MesonOperator.NOT, None))

    def evaluate_if(self, node: mparser.IfClauseNode) -> T.Optional[Disabler]:
        assert isinstance(node, mparser.IfClauseNode)
        for i in node.ifs:
            # Reset self.tmp_meson_version to know if it gets set during this
            # statement evaluation.
            self.tmp_meson_version = None
            result = self.evaluate_statement(i.condition)
            if result is None:
                raise InvalidCodeOnVoid('if')
            if isinstance(result, Disabler):
                return result
            if not isinstance(result, InterpreterObject):
                raise mesonlib.MesonBugException(f'Argument to if ({result}) is not an InterpreterObject but {type(result).__name__}.')
            res = result.operator_call(MesonOperator.BOOL, None)
            if not isinstance(res, bool):
                raise InvalidCode(f'If clause {result!r} does not evaluate to true or false.')
            if res:
                prev_meson_version = mesonlib.project_meson_versions[self.subproject]
                if self.tmp_meson_version:
                    mesonlib.project_meson_versions[self.subproject] = self.tmp_meson_version
                try:
                    self.evaluate_codeblock(i.block)
                finally:
                    mesonlib.project_meson_versions[self.subproject] = prev_meson_version
                return None
        if not isinstance(node.elseblock, mparser.EmptyNode):
            self.evaluate_codeblock(node.elseblock.block)
        return None

    def evaluate_testcase(self, node: mparser.TestCaseClauseNode) -> T.Optional[Disabler]:
        result = self.evaluate_statement(node.condition)
        if isinstance(result, Disabler):
            return result
        if not isinstance(result, ContextManagerObject):
            raise InvalidCode(f'testcase clause {result!r} does not evaluate to a context manager.')
        with result:
            self.evaluate_codeblock(node.block)
        return None

    def evaluate_comparison(self, node: mparser.ComparisonNode) -> InterpreterObject:
        val1 = self.evaluate_statement(node.left)
        if val1 is None:
            raise mesonlib.MesonException('Cannot compare a void statement on the left-hand side')
        if isinstance(val1, Disabler):
            return val1
        val2 = self.evaluate_statement(node.right)
        if val2 is None:
            raise mesonlib.MesonException('Cannot compare a void statement on the right-hand side')
        if isinstance(val2, Disabler):
            return val2

        # New code based on InterpreterObjects
        operator = {
            'in': MesonOperator.IN,
            'notin': MesonOperator.NOT_IN,
            '==': MesonOperator.EQUALS,
            '!=': MesonOperator.NOT_EQUALS,
            '>': MesonOperator.GREATER,
            '<': MesonOperator.LESS,
            '>=': MesonOperator.GREATER_EQUALS,
            '<=': MesonOperator.LESS_EQUALS,
        }[node.ctype]

        # Check if the arguments should be reversed for simplicity (this essentially converts `in` to `contains`)
        if operator in (MesonOperator.IN, MesonOperator.NOT_IN):
            val1, val2 = val2, val1

        val1.current_node = node
        return self._holderify(val1.operator_call(operator, _unholder(val2)))

    def evaluate_andstatement(self, cur: mparser.AndNode) -> InterpreterObject:
        l = self.evaluate_statement(cur.left)
        if l is None:
            raise mesonlib.MesonException('Cannot compare a void statement on the left-hand side')
        if isinstance(l, Disabler):
            return l
        l_bool = l.operator_call(MesonOperator.BOOL, None)
        if not l_bool:
            return self._holderify(l_bool)
        r = self.evaluate_statement(cur.right)
        if r is None:
            raise mesonlib.MesonException('Cannot compare a void statement on the right-hand side')
        if isinstance(r, Disabler):
            return r
        return self._holderify(r.operator_call(MesonOperator.BOOL, None))

    def evaluate_orstatement(self, cur: mparser.OrNode) -> InterpreterObject:
        l = self.evaluate_statement(cur.left)
        if l is None:
            raise mesonlib.MesonException('Cannot compare a void statement on the left-hand side')
        if isinstance(l, Disabler):
            return l
        l_bool = l.operator_call(MesonOperator.BOOL, None)
        if l_bool:
            return self._holderify(l_bool)
        r = self.evaluate_statement(cur.right)
        if r is None:
            raise mesonlib.MesonException('Cannot compare a void statement on the right-hand side')
        if isinstance(r, Disabler):
            return r
        return self._holderify(r.operator_call(MesonOperator.BOOL, None))

    def evaluate_uminusstatement(self, cur: mparser.UMinusNode) -> InterpreterObject:
        v = self.evaluate_statement(cur.value)
        if v is None:
            raise InvalidCodeOnVoid('unary minus')
        if isinstance(v, Disabler):
            return v
        v.current_node = cur
        return self._holderify(v.operator_call(MesonOperator.UMINUS, None))

    def evaluate_arithmeticstatement(self, cur: mparser.ArithmeticNode) -> InterpreterObject:
        l = self.evaluate_statement(cur.left)
        if isinstance(l, Disabler):
            return l
        r = self.evaluate_statement(cur.right)
        if isinstance(r, Disabler):
            return r
        if l is None or r is None:
            raise InvalidCodeOnVoid(cur.operation)

        mapping: T.Dict[str, MesonOperator] = {
            'add': MesonOperator.PLUS,
            'sub': MesonOperator.MINUS,
            'mul': MesonOperator.TIMES,
            'div': MesonOperator.DIV,
            'mod': MesonOperator.MOD,
        }
        l.current_node = cur
        res = l.operator_call(mapping[cur.operation], _unholder(r))
        return self._holderify(res)

    def evaluate_ternary(self, node: mparser.TernaryNode) -> T.Optional[InterpreterObject]:
        assert isinstance(node, mparser.TernaryNode)
        result = self.evaluate_statement(node.condition)
        if result is None:
            raise mesonlib.MesonException('Cannot use a void statement as condition for ternary operator.')
        if isinstance(result, Disabler):
            return result
        result.current_node = node
        result_bool = result.operator_call(MesonOperator.BOOL, None)
        if result_bool:
            return self.evaluate_statement(node.trueblock)
        else:
            return self.evaluate_statement(node.falseblock)

    @FeatureNew('multiline format strings', '0.63.0')
    def evaluate_multiline_fstring(self, node: mparser.StringNode) -> InterpreterObject:
        return self.evaluate_fstring(node)

    @FeatureNew('format strings', '0.58.0')
    def evaluate_fstring(self, node: mparser.StringNode) -> InterpreterObject:
        def replace(match: T.Match[str]) -> str:
            var = str(match.group(1))
            try:
                val = _unholder(self.variables[var])
                if isinstance(val, (list, dict)):
                    FeatureNew.single_use('List or dictionary in f-string', '1.3.0', self.subproject, location=self.current_node)
                try:
                    return stringifyUserArguments(val, self.subproject)
                except InvalidArguments as e:
                    raise InvalidArguments(f'f-string: {str(e)}')
            except KeyError:
                raise InvalidCode(f'Identifier "{var}" does not name a variable.')

        res = re.sub(r'@([_a-zA-Z][_0-9a-zA-Z]*)@', replace, node.value)
        return self._holderify(res)

    def evaluate_foreach(self, node: mparser.ForeachClauseNode) -> None:
        assert isinstance(node, mparser.ForeachClauseNode)
        items = self.evaluate_statement(node.items)
        if not isinstance(items, IterableObject):
            raise InvalidArguments('Items of foreach loop do not support iterating')

        tsize = items.iter_tuple_size()
        if len(node.varnames) != (tsize or 1):
            raise InvalidArguments(f'Foreach expects exactly {tsize or 1} variables for iterating over objects of type {items.display_name()}')

        for i in items.iter_self():
            if tsize is None:
                if isinstance(i, tuple):
                    raise mesonlib.MesonBugException(f'Iteration of {items} returned a tuple even though iter_tuple_size() is None')
                self.set_variable(node.varnames[0].value, self._holderify(i))
            else:
                if not isinstance(i, tuple):
                    raise mesonlib.MesonBugException(f'Iteration of {items} did not return a tuple even though iter_tuple_size() is {tsize}')
                if len(i) != tsize:
                    raise mesonlib.MesonBugException(f'Iteration of {items} did not return a tuple even though iter_tuple_size() is {tsize}')
                for j in range(tsize):
                    self.set_variable(node.varnames[j].value, self._holderify(i[j]))
            try:
                self.evaluate_codeblock(node.block)
            except ContinueRequest:
                continue
            except BreakRequest:
                break

    def evaluate_plusassign(self, node: mparser.PlusAssignmentNode) -> None:
        assert isinstance(node, mparser.PlusAssignmentNode)
        varname = node.var_name.value
        addition = self.evaluate_statement(node.value)
        if addition is None:
            raise InvalidCodeOnVoid('plus assign')

        # Remember that all variables are immutable. We must always create a
        # full new variable and then assign it.
        old_variable = self.get_variable(varname)
        old_variable.current_node = node
        new_value = self._holderify(old_variable.operator_call(MesonOperator.PLUS, _unholder(addition)))
        self.set_variable(varname, new_value)

    def evaluate_indexing(self, node: mparser.IndexNode) -> InterpreterObject:
        assert isinstance(node, mparser.IndexNode)
        iobject = self.evaluate_statement(node.iobject)
        if iobject is None:
            raise InterpreterException('Tried to evaluate indexing on void.')
        if isinstance(iobject, Disabler):
            return iobject
        index_holder = self.evaluate_statement(node.index)
        if index_holder is None:
            raise InvalidArguments('Cannot use void statement as index.')
        index = _unholder(index_holder)

        iobject.current_node = node
        return self._holderify(iobject.operator_call(MesonOperator.INDEX, index))

    def function_call(self, node: mparser.FunctionNode) -> T.Optional[InterpreterObject]:
        func_name = node.func_name.value
        (h_posargs, h_kwargs) = self.reduce_arguments(node.args)
        (posargs, kwargs) = self._unholder_args(h_posargs, h_kwargs)
        if is_disabled(posargs, kwargs) and func_name not in {'get_variable', 'set_variable', 'unset_variable', 'is_disabler'}:
            return Disabler()
        if func_name in self.funcs:
            func = self.funcs[func_name]
            func_args = posargs
            if not getattr(func, 'no-args-flattening', False):
                func_args = flatten(posargs)
            if not getattr(func, 'no-second-level-holder-flattening', False):
                func_args, kwargs = resolve_second_level_holders(func_args, kwargs)
            self.current_node = node
            res = func(node, func_args, kwargs)
            return self._holderify(res) if res is not None else None
        else:
            self.unknown_function_called(func_name)
            return None

    def method_call(self, node: mparser.MethodNode) -> T.Optional[InterpreterObject]:
        invocable = node.source_object
        obj: T.Optional[InterpreterObject]
        if isinstance(invocable, mparser.IdNode):
            object_display_name = f'variable "{invocable.value}"'
            obj = self.get_variable(invocable.value)
        else:
            object_display_name = invocable.__class__.__name__
            obj = self.evaluate_statement(invocable)
        method_name = node.name.value
        (h_args, h_kwargs) = self.reduce_arguments(node.args)
        (args, kwargs) = self._unholder_args(h_args, h_kwargs)
        if is_disabled(args, kwargs):
            return Disabler()
        if not isinstance(obj, InterpreterObject):
            raise InvalidArguments(f'{object_display_name} is not callable.')
        # TODO: InterpreterBase **really** shouldn't be in charge of checking this
        if method_name == 'extract_objects':
            if isinstance(obj, ObjectHolder):
                self.validate_extraction(obj.held_object)
            elif not isinstance(obj, Disabler):
                raise InvalidArguments(f'Invalid operation "extract_objects" on {object_display_name} of type {type(obj).__name__}')
        obj.current_node = self.current_node = node
        res = obj.method_call(method_name, args, kwargs)
        return self._holderify(res) if res is not None else None

    def _holderify(self, res: T.Union[TYPE_var, InterpreterObject]) -> InterpreterObject:
        if isinstance(res, HoldableTypes):
            # Always check for an exact match first.
            cls = self.holder_map.get(type(res), None)
            if cls is not None:
                # Casts to Interpreter are required here since an assertion would
                # not work for the `ast` module.
                return cls(res, T.cast('Interpreter', self))
            # Try the boundary types next.
            for typ, cls in self.bound_holder_map.items():
                if isinstance(res, typ):
                    return cls(res, T.cast('Interpreter', self))
            raise mesonlib.MesonBugException(f'Object {res} of type {type(res).__name__} is neither in self.holder_map nor self.bound_holder_map.')
        elif isinstance(res, ObjectHolder):
            raise mesonlib.MesonBugException(f'Returned object {res} of type {type(res).__name__} is an object holder.')
        elif isinstance(res, MesonInterpreterObject):
            return res
        raise mesonlib.MesonBugException(f'Unknown returned object {res} of type {type(res).__name__} in the parameters.')

    def _unholder_args(self,
                       args: T.List[InterpreterObject],
                       kwargs: T.Dict[str, InterpreterObject]) -> T.Tuple[T.List[TYPE_var], TYPE_kwargs]:
        return [_unholder(x) for x in args], {k: _unholder(v) for k, v in kwargs.items()}

    def unknown_function_called(self, func_name: str) -> None:
        raise InvalidCode(f'Unknown function "{func_name}".')

    def reduce_arguments(
                self,
                args: mparser.ArgumentNode,
                key_resolver: T.Callable[[mparser.BaseNode], str] = default_resolve_key,
                duplicate_key_error: T.Optional[str] = None,
            ) -> T.Tuple[
                T.List[InterpreterObject],
                T.Dict[str, InterpreterObject]
            ]:
        assert isinstance(args, mparser.ArgumentNode)
        if args.incorrect_order():
            raise InvalidArguments('All keyword arguments must be after positional arguments.')
        self.argument_depth += 1
        reduced_pos = [self.evaluate_statement(arg) for arg in args.arguments]
        if any(x is None for x in reduced_pos):
            raise InvalidArguments('At least one value in the arguments is void.')
        reduced_kw: T.Dict[str, InterpreterObject] = {}
        for key, val in args.kwargs.items():
            reduced_key = key_resolver(key)
            assert isinstance(val, mparser.BaseNode)
            reduced_val = self.evaluate_statement(val)
            if reduced_val is None:
                raise InvalidArguments(f'Value of key {reduced_key} is void.')
            self.current_node = key
            if duplicate_key_error and reduced_key in reduced_kw:
                raise InvalidArguments(duplicate_key_error.format(reduced_key))
            reduced_kw[reduced_key] = reduced_val
        self.argument_depth -= 1
        final_kw = self.expand_default_kwargs(reduced_kw)
        return reduced_pos, final_kw

    def expand_default_kwargs(self, kwargs: T.Dict[str, T.Optional[InterpreterObject]]) -> T.Dict[str, T.Optional[InterpreterObject]]:
        if 'kwargs' not in kwargs:
            return kwargs
        to_expand = _unholder(kwargs.pop('kwargs'))
        if not isinstance(to_expand, dict):
            raise InterpreterException('Value of "kwargs" must be dictionary.')
        if 'kwargs' in to_expand:
            raise InterpreterException('Kwargs argument must not contain a "kwargs" entry. Points for thinking meta, though. :P')
        for k, v in to_expand.items():
            if k in kwargs:
                raise InterpreterException(f'Entry "{k}" defined both as a keyword argument and in a "kwarg" entry.')
            kwargs[k] = self._holderify(v)
        return kwargs

    def assignment(self, node: mparser.AssignmentNode) -> None:
        assert isinstance(node, mparser.AssignmentNode)
        if self.argument_depth != 0:
            raise InvalidArguments(textwrap.dedent('''\
                Tried to assign values inside an argument list.
                To specify a keyword argument, use : instead of =.
            '''))
        var_name = node.var_name.value
        if not isinstance(var_name, str):
            raise InvalidArguments('Tried to assign value to a non-variable.')
        value = self.evaluate_statement(node.value)
        # For mutable objects we need to make a copy on assignment
        if isinstance(value, MutableInterpreterObject):
            value = copy.deepcopy(value)
        self.set_variable(var_name, value)

    def set_variable(self, varname: str, variable: T.Union[TYPE_var, InterpreterObject], *, holderify: bool = False) -> None:
        if variable is None:
            raise InvalidCode('Can not assign void to variable.')
        if holderify:
            variable = self._holderify(variable)
        else:
            # Ensure that we are always storing ObjectHolders
            if not isinstance(variable, InterpreterObject):
                raise mesonlib.MesonBugException(f'set_variable in InterpreterBase called with a non InterpreterObject {variable} of type {type(variable).__name__}')
        if not isinstance(varname, str):
            raise InvalidCode('First argument to set_variable must be a string.')
        if varname in self.builtin:
            raise InvalidCode(f'Tried to overwrite internal variable "{varname}"')
        self.variables[varname] = variable

    def get_variable(self, varname: str) -> InterpreterObject:
        if varname in self.builtin:
            return self.builtin[varname]
        if varname in self.variables:
            return self.variables[varname]
        raise InvalidCode(f'Unknown variable "{varname}".')

    def validate_extraction(self, buildtarget: mesonlib.HoldableObject) -> None:
        raise InterpreterException('validate_extraction is not implemented in this context (please file a bug)')
