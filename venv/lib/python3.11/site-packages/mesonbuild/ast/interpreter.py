# SPDX-License-Identifier: Apache-2.0
# Copyright 2016 The Meson development team

# This class contains the basic functionality needed to run any interpreter
# or an interpreter-based tool.
from __future__ import annotations

import os
import sys
import typing as T

from .. import mparser, mesonlib
from .. import environment

from ..interpreterbase import (
    MesonInterpreterObject,
    InterpreterBase,
    InvalidArguments,
    BreakRequest,
    ContinueRequest,
    Disabler,
    default_resolve_key,
)

from ..interpreter import (
    StringHolder,
    BooleanHolder,
    IntegerHolder,
    ArrayHolder,
    DictHolder,
)

from ..mparser import (
    ArgumentNode,
    ArithmeticNode,
    ArrayNode,
    AssignmentNode,
    BaseNode,
    ElementaryNode,
    EmptyNode,
    IdNode,
    MethodNode,
    NotNode,
    PlusAssignmentNode,
    TernaryNode,
)

if T.TYPE_CHECKING:
    from .visitor import AstVisitor
    from ..interpreter import Interpreter
    from ..interpreterbase import SubProject, TYPE_nkwargs, TYPE_var
    from ..mparser import (
        AndNode,
        ComparisonNode,
        ForeachClauseNode,
        IfClauseNode,
        IndexNode,
        OrNode,
        TestCaseClauseNode,
        UMinusNode,
    )

class DontCareObject(MesonInterpreterObject):
    pass

class MockExecutable(MesonInterpreterObject):
    pass

class MockStaticLibrary(MesonInterpreterObject):
    pass

class MockSharedLibrary(MesonInterpreterObject):
    pass

class MockCustomTarget(MesonInterpreterObject):
    pass

class MockRunTarget(MesonInterpreterObject):
    pass

ADD_SOURCE = 0
REMOVE_SOURCE = 1

_T = T.TypeVar('_T')
_V = T.TypeVar('_V')


class AstInterpreter(InterpreterBase):
    def __init__(self, source_root: str, subdir: str, subproject: SubProject, visitors: T.Optional[T.List[AstVisitor]] = None):
        super().__init__(source_root, subdir, subproject)
        self.visitors = visitors if visitors is not None else []
        self.processed_buildfiles: T.Set[str] = set()
        self.assignments: T.Dict[str, BaseNode] = {}
        self.assign_vals: T.Dict[str, T.Any] = {}
        self.reverse_assignment: T.Dict[str, BaseNode] = {}
        self.funcs.update({'project': self.func_do_nothing,
                           'test': self.func_do_nothing,
                           'benchmark': self.func_do_nothing,
                           'install_headers': self.func_do_nothing,
                           'install_man': self.func_do_nothing,
                           'install_data': self.func_do_nothing,
                           'install_subdir': self.func_do_nothing,
                           'install_symlink': self.func_do_nothing,
                           'install_emptydir': self.func_do_nothing,
                           'configuration_data': self.func_do_nothing,
                           'configure_file': self.func_do_nothing,
                           'find_program': self.func_do_nothing,
                           'include_directories': self.func_do_nothing,
                           'add_global_arguments': self.func_do_nothing,
                           'add_global_link_arguments': self.func_do_nothing,
                           'add_project_arguments': self.func_do_nothing,
                           'add_project_dependencies': self.func_do_nothing,
                           'add_project_link_arguments': self.func_do_nothing,
                           'message': self.func_do_nothing,
                           'generator': self.func_do_nothing,
                           'error': self.func_do_nothing,
                           'run_command': self.func_do_nothing,
                           'assert': self.func_do_nothing,
                           'subproject': self.func_do_nothing,
                           'dependency': self.func_do_nothing,
                           'get_option': self.func_do_nothing,
                           'join_paths': self.func_do_nothing,
                           'environment': self.func_do_nothing,
                           'import': self.func_do_nothing,
                           'vcs_tag': self.func_do_nothing,
                           'add_languages': self.func_do_nothing,
                           'declare_dependency': self.func_do_nothing,
                           'files': self.func_do_nothing,
                           'executable': self.func_do_nothing,
                           'static_library': self.func_do_nothing,
                           'shared_library': self.func_do_nothing,
                           'library': self.func_do_nothing,
                           'build_target': self.func_do_nothing,
                           'custom_target': self.func_do_nothing,
                           'run_target': self.func_do_nothing,
                           'subdir': self.func_subdir,
                           'set_variable': self.func_do_nothing,
                           'get_variable': self.func_do_nothing,
                           'unset_variable': self.func_do_nothing,
                           'is_disabler': self.func_do_nothing,
                           'is_variable': self.func_do_nothing,
                           'disabler': self.func_do_nothing,
                           'jar': self.func_do_nothing,
                           'warning': self.func_do_nothing,
                           'shared_module': self.func_do_nothing,
                           'option': self.func_do_nothing,
                           'both_libraries': self.func_do_nothing,
                           'add_test_setup': self.func_do_nothing,
                           'subdir_done': self.func_do_nothing,
                           'alias_target': self.func_do_nothing,
                           'summary': self.func_do_nothing,
                           'range': self.func_do_nothing,
                           'structured_sources': self.func_do_nothing,
                           'debug': self.func_do_nothing,
                           })

    def _unholder_args(self, args: _T, kwargs: _V) -> T.Tuple[_T, _V]:
        return args, kwargs

    def _holderify(self, res: _T) -> _T:
        return res

    def func_do_nothing(self, node: BaseNode, args: T.List[TYPE_var], kwargs: T.Dict[str, TYPE_var]) -> bool:
        return True

    def load_root_meson_file(self) -> None:
        super().load_root_meson_file()
        for i in self.visitors:
            self.ast.accept(i)

    def func_subdir(self, node: BaseNode, args: T.List[TYPE_var], kwargs: T.Dict[str, TYPE_var]) -> None:
        args = self.flatten_args(args)
        if len(args) != 1 or not isinstance(args[0], str):
            sys.stderr.write(f'Unable to evaluate subdir({args}) in AstInterpreter --> Skipping\n')
            return

        prev_subdir = self.subdir
        subdir = os.path.join(prev_subdir, args[0])
        absdir = os.path.join(self.source_root, subdir)
        buildfilename = os.path.join(subdir, environment.build_filename)
        absname = os.path.join(self.source_root, buildfilename)
        symlinkless_dir = os.path.realpath(absdir)
        build_file = os.path.join(symlinkless_dir, 'meson.build')
        if build_file in self.processed_buildfiles:
            sys.stderr.write('Trying to enter {} which has already been visited --> Skipping\n'.format(args[0]))
            return
        self.processed_buildfiles.add(build_file)

        if not os.path.isfile(absname):
            sys.stderr.write(f'Unable to find build file {buildfilename} --> Skipping\n')
            return
        code = self.read_buildfile(absname, buildfilename)
        try:
            codeblock = mparser.Parser(code, absname).parse()
        except mesonlib.MesonException as me:
            me.file = absname
            raise me

        self.subdir = subdir
        for i in self.visitors:
            codeblock.accept(i)
        self.evaluate_codeblock(codeblock)
        self.subdir = prev_subdir

    def method_call(self, node: BaseNode) -> bool:
        return True

    def evaluate_fstring(self, node: mparser.StringNode) -> str:
        assert isinstance(node, mparser.StringNode)
        return node.value

    def evaluate_arraystatement(self, cur: mparser.ArrayNode) -> TYPE_var:
        return self.reduce_arguments(cur.args)[0]

    def evaluate_arithmeticstatement(self, cur: ArithmeticNode) -> int:
        self.evaluate_statement(cur.left)
        self.evaluate_statement(cur.right)
        return 0

    def evaluate_uminusstatement(self, cur: UMinusNode) -> int:
        self.evaluate_statement(cur.value)
        return 0

    def evaluate_ternary(self, node: TernaryNode) -> None:
        assert isinstance(node, TernaryNode)
        self.evaluate_statement(node.condition)
        self.evaluate_statement(node.trueblock)
        self.evaluate_statement(node.falseblock)

    def evaluate_dictstatement(self, node: mparser.DictNode) -> TYPE_nkwargs:
        def resolve_key(node: mparser.BaseNode) -> str:
            if isinstance(node, mparser.StringNode):
                return node.value
            return '__AST_UNKNOWN__'
        arguments, kwargs = self.reduce_arguments(node.args, key_resolver=resolve_key)
        assert not arguments
        self.argument_depth += 1
        for key, value in kwargs.items():
            if isinstance(key, BaseNode):
                self.evaluate_statement(key)
        self.argument_depth -= 1
        return {}

    def evaluate_plusassign(self, node: PlusAssignmentNode) -> None:
        assert isinstance(node, PlusAssignmentNode)
        # Cheat by doing a reassignment
        self.assignments[node.var_name.value] = node.value  # Save a reference to the value node
        if node.value.ast_id:
            self.reverse_assignment[node.value.ast_id] = node
        self.assign_vals[node.var_name.value] = self.evaluate_statement(node.value)

    def evaluate_indexing(self, node: IndexNode) -> int:
        return 0

    def unknown_function_called(self, func_name: str) -> None:
        pass

    def reduce_arguments(
                self,
                args: mparser.ArgumentNode,
                key_resolver: T.Callable[[mparser.BaseNode], str] = default_resolve_key,
                duplicate_key_error: T.Optional[str] = None,
            ) -> T.Tuple[T.List[TYPE_var], TYPE_nkwargs]:
        if isinstance(args, ArgumentNode):
            kwargs: T.Dict[str, TYPE_var] = {}
            for key, val in args.kwargs.items():
                kwargs[key_resolver(key)] = val
            if args.incorrect_order():
                raise InvalidArguments('All keyword arguments must be after positional arguments.')
            return self.flatten_args(args.arguments), kwargs
        else:
            return self.flatten_args(args), {}

    def evaluate_comparison(self, node: ComparisonNode) -> bool:
        self.evaluate_statement(node.left)
        self.evaluate_statement(node.right)
        return False

    def evaluate_andstatement(self, cur: AndNode) -> bool:
        self.evaluate_statement(cur.left)
        self.evaluate_statement(cur.right)
        return False

    def evaluate_orstatement(self, cur: OrNode) -> bool:
        self.evaluate_statement(cur.left)
        self.evaluate_statement(cur.right)
        return False

    def evaluate_notstatement(self, cur: NotNode) -> bool:
        self.evaluate_statement(cur.value)
        return False

    def evaluate_foreach(self, node: ForeachClauseNode) -> None:
        try:
            self.evaluate_codeblock(node.block)
        except ContinueRequest:
            pass
        except BreakRequest:
            pass

    def evaluate_if(self, node: IfClauseNode) -> None:
        for i in node.ifs:
            self.evaluate_codeblock(i.block)
        if not isinstance(node.elseblock, EmptyNode):
            self.evaluate_codeblock(node.elseblock.block)

    def get_variable(self, varname: str) -> int:
        return 0

    def assignment(self, node: AssignmentNode) -> None:
        assert isinstance(node, AssignmentNode)
        self.assignments[node.var_name.value] = node.value # Save a reference to the value node
        if node.value.ast_id:
            self.reverse_assignment[node.value.ast_id] = node
        self.assign_vals[node.var_name.value] = self.evaluate_statement(node.value) # Evaluate the value just in case

    def resolve_node(self, node: BaseNode, include_unknown_args: bool = False, id_loop_detect: T.Optional[T.List[str]] = None) -> T.Optional[T.Any]:
        def quick_resolve(n: BaseNode, loop_detect: T.Optional[T.List[str]] = None) -> T.Any:
            if loop_detect is None:
                loop_detect = []
            if isinstance(n, IdNode):
                assert isinstance(n.value, str)
                if n.value in loop_detect or n.value not in self.assignments:
                    return []
                return quick_resolve(self.assignments[n.value], loop_detect = loop_detect + [n.value])
            elif isinstance(n, ElementaryNode):
                return n.value
            else:
                return n

        if id_loop_detect is None:
            id_loop_detect = []
        result = None

        if not isinstance(node, BaseNode):
            return None

        assert node.ast_id
        if node.ast_id in id_loop_detect:
            return None # Loop detected
        id_loop_detect += [node.ast_id]

        # Try to evaluate the value of the node
        if isinstance(node, IdNode):
            result = quick_resolve(node)

        elif isinstance(node, ElementaryNode):
            result = node.value

        elif isinstance(node, NotNode):
            result = self.resolve_node(node.value, include_unknown_args, id_loop_detect)
            if isinstance(result, bool):
                result = not result

        elif isinstance(node, ArrayNode):
            result = node.args.arguments.copy()

        elif isinstance(node, ArgumentNode):
            result = node.arguments.copy()

        elif isinstance(node, ArithmeticNode):
            if node.operation != 'add':
                return None # Only handle string and array concats
            l = self.resolve_node(node.left, include_unknown_args, id_loop_detect)
            r = self.resolve_node(node.right, include_unknown_args, id_loop_detect)
            if isinstance(l, str) and isinstance(r, str):
                result = l + r # String concatenation detected
            else:
                result = self.flatten_args(l, include_unknown_args, id_loop_detect) + self.flatten_args(r, include_unknown_args, id_loop_detect)

        elif isinstance(node, MethodNode):
            src = quick_resolve(node.source_object)
            margs = self.flatten_args(node.args.arguments, include_unknown_args, id_loop_detect)
            mkwargs: T.Dict[str, TYPE_var] = {}
            method_name = node.name.value
            try:
                if isinstance(src, str):
                    result = StringHolder(src, T.cast('Interpreter', self)).method_call(method_name, margs, mkwargs)
                elif isinstance(src, bool):
                    result = BooleanHolder(src, T.cast('Interpreter', self)).method_call(method_name, margs, mkwargs)
                elif isinstance(src, int):
                    result = IntegerHolder(src, T.cast('Interpreter', self)).method_call(method_name, margs, mkwargs)
                elif isinstance(src, list):
                    result = ArrayHolder(src, T.cast('Interpreter', self)).method_call(method_name, margs, mkwargs)
                elif isinstance(src, dict):
                    result = DictHolder(src, T.cast('Interpreter', self)).method_call(method_name, margs, mkwargs)
            except mesonlib.MesonException:
                return None

        # Ensure that the result is fully resolved (no more nodes)
        if isinstance(result, BaseNode):
            result = self.resolve_node(result, include_unknown_args, id_loop_detect)
        elif isinstance(result, list):
            new_res: T.List[TYPE_var] = []
            for i in result:
                if isinstance(i, BaseNode):
                    resolved = self.resolve_node(i, include_unknown_args, id_loop_detect)
                    if resolved is not None:
                        new_res += self.flatten_args(resolved, include_unknown_args, id_loop_detect)
                else:
                    new_res += [i]
            result = new_res

        return result

    def flatten_args(self, args_raw: T.Union[TYPE_var, T.Sequence[TYPE_var]], include_unknown_args: bool = False, id_loop_detect: T.Optional[T.List[str]] = None) -> T.List[TYPE_var]:
        # Make sure we are always dealing with lists
        if isinstance(args_raw, list):
            args = args_raw
        else:
            args = [args_raw]

        flattened_args: T.List[TYPE_var] = []

        # Resolve the contents of args
        for i in args:
            if isinstance(i, BaseNode):
                resolved = self.resolve_node(i, include_unknown_args, id_loop_detect)
                if resolved is not None:
                    if not isinstance(resolved, list):
                        resolved = [resolved]
                    flattened_args += resolved
            elif isinstance(i, (str, bool, int, float)) or include_unknown_args:
                flattened_args += [i]
        return flattened_args

    def flatten_kwargs(self, kwargs: T.Dict[str, TYPE_var], include_unknown_args: bool = False) -> T.Dict[str, TYPE_var]:
        flattened_kwargs = {}
        for key, val in kwargs.items():
            if isinstance(val, BaseNode):
                resolved = self.resolve_node(val, include_unknown_args)
                if resolved is not None:
                    flattened_kwargs[key] = resolved
            elif isinstance(val, (str, bool, int, float)) or include_unknown_args:
                flattened_kwargs[key] = val
        return flattened_kwargs

    def evaluate_testcase(self, node: TestCaseClauseNode) -> Disabler | None:
        return Disabler(subproject=self.subproject)
