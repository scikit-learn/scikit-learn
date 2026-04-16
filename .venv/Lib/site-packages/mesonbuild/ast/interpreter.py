# SPDX-License-Identifier: Apache-2.0
# Copyright 2016 The Meson development team

# This class contains the basic functionality needed to run any interpreter
# or an interpreter-based tool.
from __future__ import annotations

import os
import sys
import typing as T
from collections import defaultdict
from dataclasses import dataclass
import itertools
from pathlib import Path

from .. import mparser, mesonlib, mlog
from .. import environment

from ..interpreterbase import (
    MesonInterpreterObject,
    InterpreterBase,
    InvalidArguments,
    BreakRequest,
    ContinueRequest,
    Disabler,
    default_resolve_key,
    is_disabled,
    UnknownValue,
    UndefinedVariable,
    InterpreterObject,
)

from ..interpreterbase.helpers import flatten

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
    EmptyNode,
    IdNode,
    MethodNode,
    NotNode,
    PlusAssignmentNode,
    TernaryNode,
    SymbolNode,
    Token,
    FunctionNode,
)

if T.TYPE_CHECKING:
    from .visitor import AstVisitor
    from ..interpreter import Interpreter
    from ..interpreterbase import SubProject, TYPE_var, TYPE_nvar
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

_T = T.TypeVar('_T')
_V = T.TypeVar('_V')

def _symbol(val: str) -> SymbolNode:
    return SymbolNode(Token('', '', 0, 0, 0, (0, 0), val))

# `IntrospectionFile` is to the `IntrospectionInterpreter` what `File` is to the normal `Interpreter`.
@dataclass
class IntrospectionFile:
    subdir: str
    rel: str

    def to_abs_path(self, root_dir: Path) -> Path:
        return (root_dir / self.subdir / self.rel).resolve()

    def __hash__(self) -> int:
        return hash((self.__class__.__name__, self.subdir, self.rel))

# `IntrospectionDependency` is to the `IntrospectionInterpreter` what `Dependency` is to the normal `Interpreter`.
@dataclass
class IntrospectionDependency(MesonInterpreterObject):
    name: T.Union[str, UnknownValue]
    required: T.Union[bool, UnknownValue]
    version: T.Union[T.List[str], UnknownValue]
    has_fallback: bool
    conditional: bool
    node: FunctionNode

# `IntrospectionBuildTarget` is to the `IntrospectionInterpreter` what `BuildTarget` is to the normal `Interpreter`.
@dataclass
class IntrospectionBuildTarget(MesonInterpreterObject):
    name: str
    machine: str
    id: str
    typename: str
    defined_in: str
    subdir: str
    build_by_default: T.Union[bool, UnknownValue]
    installed: T.Union[bool, UnknownValue]
    outputs: T.List[str]
    source_nodes: T.List[BaseNode]
    extra_files: BaseNode
    kwargs: T.Dict[str, TYPE_var]
    node: FunctionNode

def is_ignored_edge(src: T.Union[BaseNode, UnknownValue]) -> bool:
    return (isinstance(src, FunctionNode) and src.func_name.value not in {'files', 'get_variable'}) or isinstance(src, MethodNode)

class DataflowDAG:
    src_to_tgts: T.DefaultDict[T.Union[BaseNode, UnknownValue], T.Set[T.Union[BaseNode, UnknownValue]]]
    tgt_to_srcs: T.DefaultDict[T.Union[BaseNode, UnknownValue], T.Set[T.Union[BaseNode, UnknownValue]]]

    def __init__(self) -> None:
        self.src_to_tgts = defaultdict(set)
        self.tgt_to_srcs = defaultdict(set)

    def add_edge(self, source: T.Union[BaseNode, UnknownValue], target: T.Union[BaseNode, UnknownValue]) -> None:
        self.src_to_tgts[source].add(target)
        self.tgt_to_srcs[target].add(source)

    # Returns all nodes in the DAG that are reachable from a node in `srcs`.
    # In other words, A node `a` is part of the returned set exactly if data
    # from `srcs` flows into `a`, directly or indirectly.
    # Certain edges are ignored.
    def reachable(self, srcs: T.Set[T.Union[BaseNode, UnknownValue]], reverse: bool) -> T.Set[T.Union[BaseNode, UnknownValue]]:
        reachable = srcs.copy()
        active = srcs.copy()
        while active:
            new: T.Set[T.Union[BaseNode, UnknownValue]] = set()
            if reverse:
                for tgt in active:
                    new.update(src for src in self.tgt_to_srcs[tgt] if not is_ignored_edge(src))
            else:
                for src in active:
                    if is_ignored_edge(src):
                        continue
                    new.update(tgt for tgt in self.src_to_tgts[src])
            reachable.update(new)
            active = new
        return reachable

    # Returns all paths from src to target.
    # Certain edges are ignored.
    def find_all_paths(self, src: T.Union[BaseNode, UnknownValue], target: T.Union[BaseNode, UnknownValue]) -> T.List[T.List[T.Union[BaseNode, UnknownValue]]]:
        queue = [(src, [src])]
        paths = []
        while queue:
            cur, path = queue.pop()
            if cur == target:
                paths.append(path)
            if is_ignored_edge(cur):
                continue
            queue.extend((tgt, path + [tgt]) for tgt in self.src_to_tgts[cur])
        return paths

class AstInterpreter(InterpreterBase):
    def __init__(self, source_root: str, subdir: str, subproject: SubProject, subproject_dir: str, env: environment.Environment, visitors: T.Optional[T.List[AstVisitor]] = None):
        super().__init__(source_root, subdir, subproject, subproject_dir, env)
        self.visitors = visitors if visitors is not None else []
        self.nesting: T.List[int] = []
        self.cur_assignments: T.DefaultDict[str, T.List[T.Tuple[T.List[int], T.Union[BaseNode, UnknownValue]]]] = defaultdict(list)
        self.all_assignment_nodes: T.DefaultDict[str, T.List[AssignmentNode]] = defaultdict(list)
        # dataflow_dag is an acyclic directed graph that contains an edge
        # from one instance of `BaseNode` to another instance of `BaseNode` if
        # data flows directly from one to the other. Example: If meson.build
        # contains this:
        # var = 'foo' + '123'
        # executable(var, 'src.c')
        # var = 'bar'
        # dataflow_dag will contain an edge from the IdNode corresponding to
        # 'var' in line 2 to the ArithmeticNode corresponding to 'foo' + '123'.
        # This graph is crucial for e.g. node_to_runtime_value because we have
        # to know that 'var' in line2 is 'foo123' and not 'bar'.
        self.dataflow_dag = DataflowDAG()
        self.funcvals: T.Dict[BaseNode, T.Any] = {}
        self.tainted = False
        self.predefined_vars = {
            'meson': UnknownValue(),
            'host_machine': UnknownValue(),
            'build_machine': UnknownValue(),
            'target_machine': UnknownValue()
        }
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
                           'files': self.func_files,
                           'executable': self.func_do_nothing,
                           'static_library': self.func_do_nothing,
                           'shared_library': self.func_do_nothing,
                           'library': self.func_do_nothing,
                           'build_target': self.func_do_nothing,
                           'custom_target': self.func_do_nothing,
                           'run_target': self.func_do_nothing,
                           'subdir': self.func_subdir,
                           'set_variable': self.func_set_variable,
                           'get_variable': self.func_get_variable,
                           'unset_variable': self.func_unset_variable,
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

    def _unholder_args(self, args: T.Any, kwargs: T.Any) -> T.Tuple[T.Any, T.Any]:
        return args, kwargs

    def _holderify(self, res: T.Any) -> T.Any:
        return res

    def func_do_nothing(self, node: BaseNode, args: T.List[TYPE_var], kwargs: T.Dict[str, TYPE_var]) -> UnknownValue:
        return UnknownValue()

    def load_root_meson_file(self) -> None:
        super().load_root_meson_file()
        for i in self.visitors:
            self.ast.accept(i)

    def func_subdir(self, node: BaseNode, args: T.List[TYPE_var], kwargs: T.Dict[str, TYPE_var]) -> None:
        args = self.flatten_args(args)
        if len(args) != 1 or not isinstance(args[0], str):
            sys.stderr.write(f'Unable to evaluate subdir({args}) in AstInterpreter --> Skipping\n')
            return

        subdir, is_new = self._resolve_subdir(self.source_root, args[0])
        if not is_new:
            sys.stderr.write('Trying to enter {} which has already been visited --> Skipping\n'.format(args[0]))
            return

        if not self._evaluate_subdir(self.source_root, subdir, self.visitors):
            buildfilename = os.path.join(subdir, environment.build_filename)
            sys.stderr.write(f'Unable to find build file {buildfilename} --> Skipping\n')

    def inner_method_call(self, obj: BaseNode, method_name: str, args: T.List[TYPE_var], kwargs: T.Dict[str, TYPE_var]) -> T.Any:
        for arg in itertools.chain(args, kwargs.values()):
            if isinstance(arg, UnknownValue):
                return UnknownValue()

        if isinstance(obj, str):
            result = StringHolder(obj, T.cast('Interpreter', self)).method_call(method_name, args, kwargs)
        elif isinstance(obj, bool):
            result = BooleanHolder(obj, T.cast('Interpreter', self)).method_call(method_name, args, kwargs)
        elif isinstance(obj, int):
            result = IntegerHolder(obj, T.cast('Interpreter', self)).method_call(method_name, args, kwargs)
        elif isinstance(obj, list):
            result = ArrayHolder(obj, T.cast('Interpreter', self)).method_call(method_name, args, kwargs)
        elif isinstance(obj, dict):
            result = DictHolder(obj, T.cast('Interpreter', self)).method_call(method_name, args, kwargs)
        else:
            return UnknownValue()
        return result

    def method_call(self, node: mparser.MethodNode) -> None:
        invocable = node.source_object
        self.evaluate_statement(invocable)
        obj = self.node_to_runtime_value(invocable)
        method_name = node.name.value
        (args, kwargs) = self.reduce_arguments(node.args)
        if is_disabled(args, kwargs):
            res = Disabler()
        else:
            res = self.inner_method_call(obj, method_name, args, kwargs)
        self.funcvals[node] = res

    def evaluate_fstring(self, node: mparser.StringNode) -> None:
        pass

    def evaluate_arraystatement(self, cur: mparser.ArrayNode) -> None:
        for arg in cur.args.arguments:
            self.evaluate_statement(arg)

    def evaluate_arithmeticstatement(self, cur: ArithmeticNode) -> None:
        self.evaluate_statement(cur.left)
        self.evaluate_statement(cur.right)

    def evaluate_uminusstatement(self, cur: UMinusNode) -> None:
        self.evaluate_statement(cur.value)

    def evaluate_ternary(self, node: TernaryNode) -> None:
        assert isinstance(node, TernaryNode)
        self.evaluate_statement(node.condition)
        self.evaluate_statement(node.trueblock)
        self.evaluate_statement(node.falseblock)

    def evaluate_dictstatement(self, node: mparser.DictNode) -> None:
        for k, v in node.args.kwargs.items():
            self.evaluate_statement(k)
            self.evaluate_statement(v)

    def evaluate_indexing(self, node: IndexNode) -> None:
        self.evaluate_statement(node.iobject)
        self.evaluate_statement(node.index)

    def reduce_arguments(
                self,
                args: mparser.ArgumentNode,
                key_resolver: T.Callable[[mparser.BaseNode], str] = default_resolve_key,
                duplicate_key_error: T.Optional[str] = None,
            ) -> T.Tuple[T.List[T.Any], T.Any]:
        for arg in args.arguments:
            self.evaluate_statement(arg)
        for value in args.kwargs.values():
            self.evaluate_statement(value)
        if isinstance(args, ArgumentNode):
            kwargs = {}
            for key, val in args.kwargs.items():
                kwargs[key_resolver(key)] = val
            if args.incorrect_order():
                raise InvalidArguments('All keyword arguments must be after positional arguments.')
            return self.flatten_args(args.arguments), kwargs
        else:
            return self.flatten_args(args), {}

    def evaluate_comparison(self, node: ComparisonNode) -> None:
        self.evaluate_statement(node.left)
        self.evaluate_statement(node.right)

    def evaluate_andstatement(self, cur: AndNode) -> None:
        self.evaluate_statement(cur.left)
        self.evaluate_statement(cur.right)

    def evaluate_orstatement(self, cur: OrNode) -> None:
        self.evaluate_statement(cur.left)
        self.evaluate_statement(cur.right)

    def evaluate_notstatement(self, cur: NotNode) -> None:
        self.evaluate_statement(cur.value)

    def find_potential_writes(self, node: BaseNode) -> T.Set[str]:
        if isinstance(node, mparser.ForeachClauseNode):
            return {el.value for el in node.varnames} | self.find_potential_writes(node.block)
        elif isinstance(node, mparser.CodeBlockNode):
            ret = set()
            for line in node.lines:
                ret.update(self.find_potential_writes(line))
            return ret
        elif isinstance(node, (AssignmentNode, PlusAssignmentNode)):
            return set([node.var_name.value]) | self.find_potential_writes(node.value)
        elif isinstance(node, IdNode):
            return set()
        elif isinstance(node, ArrayNode):
            ret = set()
            for arg in node.args.arguments:
                ret.update(self.find_potential_writes(arg))
            return ret
        elif isinstance(node, mparser.DictNode):
            ret = set()
            for k, v in node.args.kwargs.items():
                ret.update(self.find_potential_writes(k))
                ret.update(self.find_potential_writes(v))
            return ret
        elif isinstance(node, FunctionNode):
            ret = set()
            for arg in node.args.arguments:
                ret.update(self.find_potential_writes(arg))
            for arg in node.args.kwargs.values():
                ret.update(self.find_potential_writes(arg))
            return ret
        elif isinstance(node, MethodNode):
            ret = self.find_potential_writes(node.source_object)
            for arg in node.args.arguments:
                ret.update(self.find_potential_writes(arg))
            for arg in node.args.kwargs.values():
                ret.update(self.find_potential_writes(arg))
            return ret
        elif isinstance(node, ArithmeticNode):
            return self.find_potential_writes(node.left) | self.find_potential_writes(node.right)
        elif isinstance(node, (mparser.NumberNode, mparser.StringNode, mparser.BreakNode, mparser.BooleanNode, mparser.ContinueNode)):
            return set()
        elif isinstance(node, mparser.IfClauseNode):
            if isinstance(node.elseblock, EmptyNode):
                ret = set()
            else:
                ret = self.find_potential_writes(node.elseblock.block)
            for i in node.ifs:
                ret.update(self.find_potential_writes(i))
            return ret
        elif isinstance(node, mparser.IndexNode):
            return self.find_potential_writes(node.iobject) | self.find_potential_writes(node.index)
        elif isinstance(node, mparser.IfNode):
            return self.find_potential_writes(node.condition) | self.find_potential_writes(node.block)
        elif isinstance(node, (mparser.ComparisonNode, mparser.OrNode, mparser.AndNode)):
            return self.find_potential_writes(node.left) | self.find_potential_writes(node.right)
        elif isinstance(node, mparser.NotNode):
            return self.find_potential_writes(node.value)
        elif isinstance(node, mparser.TernaryNode):
            return self.find_potential_writes(node.condition) | self.find_potential_writes(node.trueblock) | self.find_potential_writes(node.falseblock)
        elif isinstance(node, mparser.UMinusNode):
            return self.find_potential_writes(node.value)
        elif isinstance(node, mparser.ParenthesizedNode):
            return self.find_potential_writes(node.inner)
        raise mesonlib.MesonBugException('Unhandled node type')

    def evaluate_foreach(self, node: ForeachClauseNode) -> None:
        asses = self.find_potential_writes(node)
        for ass in asses:
            self.cur_assignments[ass].append((self.nesting.copy(), UnknownValue()))
        try:
            self.evaluate_codeblock(node.block)
        except ContinueRequest:
            pass
        except BreakRequest:
            pass
        for ass in asses:
            self.cur_assignments[ass].append((self.nesting.copy(), UnknownValue())) # In case the foreach loops 0 times.

    def evaluate_if(self, node: IfClauseNode) -> None:
        self.nesting.append(0)
        for i in node.ifs:
            self.evaluate_codeblock(i.block)
            self.nesting[-1] += 1
        if not isinstance(node.elseblock, EmptyNode):
            self.evaluate_codeblock(node.elseblock.block)
        self.nesting.pop()
        for var_name in self.cur_assignments:
            potential_values = []
            oldval = self.get_cur_value_if_defined(var_name)
            if not isinstance(oldval, UndefinedVariable):
                potential_values.append(oldval)
            for nesting, value in self.cur_assignments[var_name]:
                if len(nesting) > len(self.nesting):
                    potential_values.append(value)
            self.cur_assignments[var_name] = [(nesting, v) for (nesting, v) in self.cur_assignments[var_name] if len(nesting) <= len(self.nesting)]
            if len(potential_values) > 1 or (len(potential_values) > 0 and isinstance(oldval, UndefinedVariable)):
                uv = UnknownValue()
                for pv in potential_values:
                    self.dataflow_dag.add_edge(pv, uv)
                self.cur_assignments[var_name].append((self.nesting.copy(), uv))

    def func_files(self, node: BaseNode, args: T.List[TYPE_var], kwargs: T.Dict[str, TYPE_var]) -> T.Any:
        ret: T.List[T.Union[IntrospectionFile, UnknownValue]] = []
        for arg in args:
            if isinstance(arg, str):
                ret.append(IntrospectionFile(self.subdir, arg))
            elif isinstance(arg, UnknownValue):
                ret.append(UnknownValue())
            else:
                raise TypeError
        return ret

    def get_cur_value_if_defined(self, var_name: str) -> T.Union[BaseNode, UnknownValue, UndefinedVariable]:
        if var_name in self.predefined_vars:
            return self.predefined_vars[var_name]
        ret: T.Union[BaseNode, UnknownValue, UndefinedVariable] = UndefinedVariable()
        for nesting, value in reversed(self.cur_assignments[var_name]):
            if len(self.nesting) >= len(nesting) and self.nesting[:len(nesting)] == nesting:
                ret = value
                break
        if isinstance(ret, UndefinedVariable) and self.tainted:
            return UnknownValue()
        return ret

    def get_cur_value(self, var_name: str) -> T.Union[BaseNode, UnknownValue]:
        ret = self.get_cur_value_if_defined(var_name)
        if isinstance(ret, UndefinedVariable):
            path = mlog.get_relative_path(Path(self.current_node.filename), Path(os.getcwd()))
            mlog.warning(f"{path}:{self.current_node.lineno}:{self.current_node.colno} will always crash if executed, since a variable named `{var_name}` is not defined")
            # We could add more advanced analysis of code referencing undefined
            # variables, but it is probably not worth the effort and the
            # complexity. So we do the simplest thing, returning an
            # UnknownValue.
            return UnknownValue()
        return ret

    # The function `node_to_runtime_value` takes a node of the ast as an
    # argument and tries to return the same thing that would be passed to e.g.
    # `func_message` if you put `message(node)` in your `meson.build` file and
    # run `meson setup`. If this is not possible, `UnknownValue()` is returned.
    # There are 3 Reasons why this is sometimes impossible:
    #     1. Because the meson rewriter is imperfect and has not implemented everything yet
    #     2. Because the value is different on different machines, example:
    #     ```meson
    #     node = somedep.found()
    #     message(node)
    #     ```
    #     will print `true` on some machines and `false` on others, so
    #     `node_to_runtime_value` does not know whether to return `true` or
    #     `false` and will return `UnknownValue()`.
    #     3. Here:
    #     ```meson
    #     foreach x : [1, 2]
    #         node = x
    #         message(node)
    #     endforeach
    #     ```
    #     `node_to_runtime_value` does not know whether to return `1` or `2` and
    #     will return `UnknownValue()`.
    #
    # If you have something like
    # ```
    # node = [123, somedep.found()]
    # ```
    # `node_to_runtime_value` will return `[123, UnknownValue()]`.
    def node_to_runtime_value(self, node: T.Union[UnknownValue, BaseNode, TYPE_var]) -> T.Any:
        if isinstance(node, (mparser.StringNode, mparser.BooleanNode, mparser.NumberNode)):
            return node.value
        elif isinstance(node, mparser.StringNode):
            if node.is_fstring:
                return UnknownValue()
            else:
                return node.value
        elif isinstance(node, list):
            return [self.node_to_runtime_value(x) for x in node]
        elif isinstance(node, ArrayNode):
            return [self.node_to_runtime_value(x) for x in node.args.arguments]
        elif isinstance(node, mparser.DictNode):
            return {self.node_to_runtime_value(k): self.node_to_runtime_value(v) for k, v in node.args.kwargs.items()}
        elif isinstance(node, IdNode):
            assert len(self.dataflow_dag.tgt_to_srcs[node]) == 1
            val = next(iter(self.dataflow_dag.tgt_to_srcs[node]))
            return self.node_to_runtime_value(val)
        elif isinstance(node, (MethodNode, FunctionNode)):
            funcval = self.funcvals[node]
            if isinstance(funcval, (dict, str)):
                return funcval
            else:
                return self.node_to_runtime_value(funcval)
        elif isinstance(node, ArithmeticNode):
            left = self.node_to_runtime_value(node.left)
            right = self.node_to_runtime_value(node.right)
            if isinstance(left, list) and isinstance(right, UnknownValue):
                return left + [right]
            if isinstance(right, list) and isinstance(left, UnknownValue):
                return [left] + right
            if isinstance(left, UnknownValue) or isinstance(right, UnknownValue):
                return UnknownValue()
            if node.operation == '+':
                if isinstance(left, dict) and isinstance(right, dict):
                    ret = left.copy()
                    for k, v in right.items():
                        ret[k] = v
                    return ret
                if isinstance(left, list):
                    if not isinstance(right, list):
                        right = [right]
                    return left + right
                return left + right
            elif node.operation == '-':
                return left - right
            elif node.operation == '*':
                return left * right
            elif node.operation == '/':
                if isinstance(left, int) and isinstance(right, int):
                    return left // right
                elif isinstance(left, str) and isinstance(right, str):
                    return os.path.join(left, right).replace('\\', '/')
            elif node.operation == '%':
                if isinstance(left, int) and isinstance(right, int):
                    return left % right
        elif isinstance(node, (UnknownValue, IntrospectionBuildTarget, IntrospectionFile, IntrospectionDependency, str, bool, int)):
            return node
        elif isinstance(node, mparser.IndexNode):
            iobject = self.node_to_runtime_value(node.iobject)
            index = self.node_to_runtime_value(node.index)
            if isinstance(iobject, UnknownValue) or isinstance(index, UnknownValue):
                return UnknownValue()
            return iobject[index]
        elif isinstance(node, mparser.ComparisonNode):
            left = self.node_to_runtime_value(node.left)
            right = self.node_to_runtime_value(node.right)
            if isinstance(left, UnknownValue) or isinstance(right, UnknownValue):
                return UnknownValue()
            if node.ctype == '==':
                return left == right
            elif node.ctype == '!=':
                return left != right
            elif node.ctype == 'in':
                return left in right
            elif node.ctype == 'not in':
                return left not in right
        elif isinstance(node, mparser.TernaryNode):
            cond = self.node_to_runtime_value(node.condition)
            if isinstance(cond, UnknownValue):
                return UnknownValue()
            if cond is True:
                return self.node_to_runtime_value(node.trueblock)
            if cond is False:
                return self.node_to_runtime_value(node.falseblock)
        elif isinstance(node, mparser.OrNode):
            left = self.node_to_runtime_value(node.left)
            right = self.node_to_runtime_value(node.right)
            if isinstance(left, UnknownValue) or isinstance(right, UnknownValue):
                return UnknownValue()
            return left or right
        elif isinstance(node, mparser.AndNode):
            left = self.node_to_runtime_value(node.left)
            right = self.node_to_runtime_value(node.right)
            if isinstance(left, UnknownValue) or isinstance(right, UnknownValue):
                return UnknownValue()
            return left and right
        elif isinstance(node, mparser.UMinusNode):
            val = self.node_to_runtime_value(node.value)
            if isinstance(val, UnknownValue):
                return val
            if isinstance(val, (int, float)):
                return -val
        elif isinstance(node, mparser.NotNode):
            val = self.node_to_runtime_value(node.value)
            if isinstance(val, UnknownValue):
                return val
            if isinstance(val, bool):
                return not val
        elif isinstance(node, mparser.ParenthesizedNode):
            return self.node_to_runtime_value(node.inner)
        raise mesonlib.MesonBugException('Unhandled node type')

    def assignment(self, node: AssignmentNode) -> None:
        assert isinstance(node, AssignmentNode)
        self.evaluate_statement(node.value)
        self.cur_assignments[node.var_name.value].append((self.nesting.copy(), node.value))
        self.all_assignment_nodes[node.var_name.value].append(node)

    def evaluate_plusassign(self, node: PlusAssignmentNode) -> None:
        assert isinstance(node, PlusAssignmentNode)
        self.evaluate_statement(node.value)
        lhs = self.get_cur_value(node.var_name.value)
        newval: T.Union[UnknownValue, ArithmeticNode]
        if isinstance(lhs, UnknownValue):
            newval = UnknownValue()
        else:
            newval = mparser.ArithmeticNode(operation='+', left=lhs, operator=_symbol('+'), right=node.value)
        self.cur_assignments[node.var_name.value].append((self.nesting.copy(), newval))
        self.all_assignment_nodes[node.var_name.value].append(node)

        self.dataflow_dag.add_edge(lhs, newval)
        self.dataflow_dag.add_edge(node.value, newval)

    def func_set_variable(self, node: BaseNode, args: T.List[TYPE_var], kwargs: T.Dict[str, TYPE_var]) -> None:
        assert isinstance(node, FunctionNode)
        if bool(node.args.kwargs):
            raise InvalidArguments('set_variable accepts no keyword arguments')
        if len(node.args.arguments) != 2:
            raise InvalidArguments('set_variable requires exactly two positional arguments')
        var_name = args[0]
        value = node.args.arguments[1]
        if isinstance(var_name, UnknownValue):
            self.evaluate_statement(value)
            self.tainted = True
            return
        assert isinstance(var_name, str)
        equiv = AssignmentNode(var_name=IdNode(Token('', '', 0, 0, 0, (0, 0), var_name)), value=value, operator=_symbol('='))
        equiv.ast_id = str(id(equiv))
        self.evaluate_statement(equiv)

    def func_get_variable(self, node: BaseNode, args: T.List[TYPE_var], kwargs: T.Dict[str, TYPE_var]) -> T.Any:
        assert isinstance(node, FunctionNode)
        var_name = args[0]
        if isinstance(var_name, UnknownValue):
            val: T.Union[UnknownValue, BaseNode] = UnknownValue()
        else:
            assert isinstance(var_name, str)
            val = self.get_cur_value(var_name)
        self.dataflow_dag.add_edge(val, node)
        return val

    def func_unset_variable(self, node: BaseNode, args: T.List[TYPE_var], kwargs: T.Dict[str, TYPE_var]) -> None:
        assert isinstance(node, FunctionNode)
        if bool(node.args.kwargs):
            raise InvalidArguments('unset_variable accepts no keyword arguments')
        if len(node.args.arguments) != 1:
            raise InvalidArguments('unset_variable requires exactly one positional arguments')
        var_name = args[0]
        assert isinstance(var_name, str)
        self.cur_assignments[var_name].append((self.nesting.copy(), node))

    def nodes_to_pretty_filelist(self, root_path: Path, subdir: str, nodes: T.List[BaseNode]) -> T.List[T.Union[str, UnknownValue]]:
        def src_to_abs(src: T.Union[str, IntrospectionFile, UnknownValue]) -> T.Union[str, UnknownValue]:
            if isinstance(src, str):
                return os.path.normpath(os.path.join(root_path, subdir, src))
            elif isinstance(src, IntrospectionFile):
                return str(src.to_abs_path(root_path))
            elif isinstance(src, UnknownValue):
                return src
            else:
                raise TypeError

        rtvals: T.List[T.Any] = flatten([self.node_to_runtime_value(sn) for sn in nodes])
        return [src_to_abs(x) for x in rtvals]

    def flatten_args(self, args_raw: T.Union[TYPE_nvar, T.Sequence[TYPE_nvar]], include_unknown_args: bool = False) -> T.List[TYPE_var]:
        # Make sure we are always dealing with lists
        if isinstance(args_raw, list):
            args = args_raw
        else:
            args = [args_raw]

        # BaseNode resolves to Any. :/
        flattened_args: T.List[T.Union[TYPE_var, T.Any]] = []

        # Resolve the contents of args
        for i in args:
            if isinstance(i, BaseNode):
                resolved = self.node_to_runtime_value(i)
                if resolved is not None:
                    if not isinstance(resolved, list):
                        resolved = [resolved]
                    flattened_args += resolved
            elif isinstance(i, (str, bool, int, float, UnknownValue, IntrospectionFile)) or include_unknown_args:
                flattened_args += [i]
            else:
                raise NotImplementedError
        return flattened_args

    def evaluate_testcase(self, node: TestCaseClauseNode) -> Disabler | None:
        return Disabler(subproject=self.subproject)

    def evaluate_statement(self, cur: mparser.BaseNode) -> T.Optional[InterpreterObject]:
        if hasattr(cur, 'args'):
            for arg in cur.args.arguments:
                self.dataflow_dag.add_edge(arg, cur)
            for k, v in cur.args.kwargs.items():
                self.dataflow_dag.add_edge(v, cur)
        for attr in ['source_object', 'left', 'right', 'items', 'iobject', 'index', 'condition']:
            if hasattr(cur, attr):
                assert isinstance(getattr(cur, attr), mparser.BaseNode)
                self.dataflow_dag.add_edge(getattr(cur, attr), cur)
        if isinstance(cur, mparser.IdNode):
            self.dataflow_dag.add_edge(self.get_cur_value(cur.value), cur)
            return None
        else:
            return super().evaluate_statement(cur)

    def function_call(self, node: mparser.FunctionNode) -> T.Any:
        ret = super().function_call(node)
        if ret is not None:
            self.funcvals[node] = ret
        return ret
