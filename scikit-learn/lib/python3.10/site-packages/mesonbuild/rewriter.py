# SPDX-License-Identifier: Apache-2.0
# Copyright 2016 The Meson development team

# This tool is used to manipulate an existing Meson build definition.
#
# - add a file to a target
# - remove files from a target
# - move targets
# - reindent?
from __future__ import annotations

from .ast import IntrospectionInterpreter, BUILD_TARGET_FUNCTIONS, AstConditionLevel, AstIDGenerator, AstIndentationGenerator, AstPrinter
from .ast.interpreter import IntrospectionBuildTarget, IntrospectionDependency, _symbol
from .interpreterbase import UnknownValue, TV_func
from .interpreterbase.helpers import flatten
from mesonbuild.mesonlib import MesonException, setup_vsenv, relpath
from . import mlog, environment
from functools import wraps
from .mparser import Token, ArrayNode, ArgumentNode, ArithmeticNode, AssignmentNode, BaseNode, StringNode, BooleanNode, ElementaryNode, IdNode, FunctionNode, PlusAssignmentNode
from .mintro import IntrospectionEncoder
import json, os, re, sys, codecs
import typing as T
from pathlib import Path

if T.TYPE_CHECKING:
    import argparse
    from argparse import ArgumentParser, _FormatterClass
    from .mlog import AnsiDecorator

class RewriterException(MesonException):
    pass

# Note: when adding arguments, please also add them to the completion
# scripts in $MESONSRC/data/shell-completions/
def add_arguments(parser: ArgumentParser, formatter: _FormatterClass) -> None:
    parser.add_argument('-s', '--sourcedir', type=str, default='.', metavar='SRCDIR', help='Path to source directory.')
    parser.add_argument('-V', '--verbose', action='store_true', default=False, help='Enable verbose output')
    parser.add_argument('-S', '--skip-errors', dest='skip', action='store_true', default=False, help='Skip errors instead of aborting')
    subparsers = parser.add_subparsers(dest='type', title='Rewriter commands', description='Rewrite command to execute')

    # Target
    tgt_parser = subparsers.add_parser('target', aliases=['tgt'], help='Modify a target', formatter_class=formatter)
    tgt_parser.add_argument('-s', '--subdir', default='', dest='subdir', help='Subdirectory of the new target (only for the "add_target" action)')
    tgt_parser.add_argument('--type', dest='tgt_type', choices=rewriter_keys['target']['target_type'][2], default='executable',
                            help='Type of the target to add (only for the "add_target" action)')
    tgt_parser.add_argument('target', help='Name or ID of the target')
    tgt_parser.add_argument('operation', choices=['add', 'rm', 'add_target', 'rm_target', 'add_extra_files', 'rm_extra_files', 'info'],
                            help='Action to execute')
    tgt_parser.add_argument('sources', nargs='*', help='Sources to add/remove')

    # KWARGS
    kw_parser = subparsers.add_parser('kwargs', help='Modify keyword arguments', formatter_class=formatter)
    kw_parser.add_argument('operation', choices=rewriter_keys['kwargs']['operation'][2],
                           help='Action to execute')
    kw_parser.add_argument('function', choices=list(rewriter_func_kwargs.keys()),
                           help='Function type to modify')
    kw_parser.add_argument('id', help='ID of the function to modify (can be anything for "project")')
    kw_parser.add_argument('kwargs', nargs='*', help='Pairs of keyword and value')

    # Default options
    def_parser = subparsers.add_parser('default-options', aliases=['def'], help='Modify the project default options', formatter_class=formatter)
    def_parser.add_argument('operation', choices=rewriter_keys['default_options']['operation'][2],
                            help='Action to execute')
    def_parser.add_argument('options', nargs='*', help='Key, value pairs of configuration option')

    # JSON file/command
    cmd_parser = subparsers.add_parser('command', aliases=['cmd'], help='Execute a JSON array of commands', formatter_class=formatter)
    cmd_parser.add_argument('json', help='JSON string or file to execute')

class RequiredKeys:
    keys: T.Dict[str, T.Any]

    def __init__(self, keys: T.Dict[str, T.Any]):
        self.keys = keys

    def __call__(self, f: TV_func) -> TV_func:
        @wraps(f)
        def wrapped(*wrapped_args: T.Any, **wrapped_kwargs: T.Any) -> T.Any:
            assert len(wrapped_args) >= 2
            cmd = wrapped_args[1]
            for key, val in self.keys.items():
                typ = val[0] # The type of the value
                default = val[1] # The default value -- None is required
                choices = val[2] # Valid choices -- None is for everything
                if key not in cmd:
                    if default is not None:
                        cmd[key] = default
                    else:
                        raise RewriterException('Key "{}" is missing in object for {}'
                                                .format(key, f.__name__))
                if not isinstance(cmd[key], typ):
                    raise RewriterException('Invalid type of "{}". Required is {} but provided was {}'
                                            .format(key, typ.__name__, type(cmd[key]).__name__))
                if choices is not None:
                    assert isinstance(choices, list)
                    if cmd[key] not in choices:
                        raise RewriterException('Invalid value of "{}": Possible values are {} but provided was "{}"'
                                                .format(key, choices, cmd[key]))
            return f(*wrapped_args, **wrapped_kwargs)

        return T.cast('TV_func', wrapped)

class MTypeBase:
    node: BaseNode

    def __init__(self, node: T.Optional[BaseNode] = None):
        if node is None:
            self.node = self.new_node()
        else:
            self.node = node
        self.node_type = None
        for i in self.supported_nodes():
            if isinstance(self.node, i):
                self.node_type = i

    @classmethod
    def new_node(cls, value: T.Any = None) -> BaseNode:
        # Overwrite in derived class
        raise RewriterException('Internal error: new_node of MTypeBase was called')

    @classmethod
    def supported_nodes(cls) -> T.List[type]:
        # Overwrite in derived class
        return []

    def can_modify(self) -> bool:
        return self.node_type is not None

    def get_node(self) -> BaseNode:
        return self.node

    def add_value(self, value: T.Any) -> None:
        # Overwrite in derived class
        mlog.warning('Cannot add a value of type', mlog.bold(type(self).__name__), '--> skipping')

    def remove_value(self, value: T.Any) -> None:
        # Overwrite in derived class
        mlog.warning('Cannot remove a value of type', mlog.bold(type(self).__name__), '--> skipping')

    def remove_regex(self, value: T.Any) -> None:
        # Overwrite in derived class
        mlog.warning('Cannot remove a regex in type', mlog.bold(type(self).__name__), '--> skipping')

class MTypeStr(MTypeBase):
    def __init__(self, node: T.Optional[BaseNode] = None):
        super().__init__(node)

    @classmethod
    def new_node(cls, value: T.Optional[str] = None) -> BaseNode:
        if value is None:
            value = ''
        return StringNode(Token('string', '', 0, 0, 0, None, str(value)))

    @classmethod
    def supported_nodes(cls) -> T.List[type]:
        return [StringNode]

class MTypeBool(MTypeBase):
    def __init__(self, node: T.Optional[BaseNode] = None):
        super().__init__(node)

    @classmethod
    def new_node(cls, value: T.Optional[str] = None) -> BaseNode:
        return BooleanNode(Token('', '', 0, 0, 0, None, bool(value)))

    @classmethod
    def supported_nodes(cls) -> T.List[type]:
        return [BooleanNode]

class MTypeID(MTypeBase):
    def __init__(self, node: T.Optional[BaseNode] = None):
        super().__init__(node)

    @classmethod
    def new_node(cls, value: T.Optional[str] = None) -> BaseNode:
        if value is None:
            value = ''
        return IdNode(Token('', '', 0, 0, 0, None, str(value)))

    @classmethod
    def supported_nodes(cls) -> T.List[type]:
        return [IdNode]

class MTypeList(MTypeBase):
    node: ArrayNode

    def __init__(self, node: T.Optional[BaseNode] = None):
        super().__init__(node)

    @classmethod
    def new_node(cls, value: T.Optional[T.List[T.Any]] = None) -> ArrayNode:
        if value is None:
            value = []
        elif not isinstance(value, list):
            return cls._new_element_node(value)
        args = ArgumentNode(Token('', '', 0, 0, 0, None, ''))
        args.arguments = [cls._new_element_node(i) for i in value]
        return ArrayNode(_symbol('['), args, _symbol(']'))

    @classmethod
    def _new_element_node(cls, value: T.Any) -> BaseNode:
        # Overwrite in derived class
        raise RewriterException('Internal error: _new_element_node of MTypeList was called')

    def _ensure_array_node(self) -> None:
        if not isinstance(self.node, ArrayNode):
            tmp = self.node
            self.node = self.new_node()
            self.node.args.arguments = [tmp]

    @staticmethod
    def _check_is_equal(node: BaseNode, value: str) -> bool:
        # Overwrite in derived class
        return False

    @staticmethod
    def _check_regex_matches(node: BaseNode, regex: str) -> bool:
        # Overwrite in derived class
        return False

    def get_node(self) -> BaseNode:
        if isinstance(self.node, ArrayNode):
            if len(self.node.args.arguments) == 1:
                return self.node.args.arguments[0]
        return self.node

    @classmethod
    def supported_element_nodes(cls) -> T.List[T.Type]:
        # Overwrite in derived class
        return []

    @classmethod
    def supported_nodes(cls) -> T.List[T.Type]:
        return [ArrayNode] + cls.supported_element_nodes()

    def add_value(self, value: T.Any) -> None:
        if not isinstance(value, list):
            value = [value]
        self._ensure_array_node()
        for i in value:
            assert hasattr(self.node, 'args') # For mypy
            assert isinstance(self.node.args, ArgumentNode) # For mypy
            self.node.args.arguments += [self._new_element_node(i)]

    def _remove_helper(self, value: T.Any, equal_func: T.Callable[[T.Any, T.Any], bool]) -> None:
        def check_remove_node(node: BaseNode) -> bool:
            for j in value:
                if equal_func(i, j):
                    return True
            return False

        if not isinstance(value, list):
            value = [value]
        self._ensure_array_node()
        assert hasattr(self.node, 'args') # For mypy
        assert isinstance(self.node.args, ArgumentNode) # For mypy
        removed_list = []
        for i in self.node.args.arguments:
            if not check_remove_node(i):
                removed_list += [i]
        self.node.args.arguments = removed_list

    def remove_value(self, value: T.Any) -> None:
        self._remove_helper(value, self._check_is_equal)

    def remove_regex(self, regex: str) -> None:
        self._remove_helper(regex, self._check_regex_matches)

class MTypeStrList(MTypeList):
    def __init__(self, node: T.Optional[BaseNode] = None):
        super().__init__(node)

    @classmethod
    def _new_element_node(cls, value: str) -> StringNode:
        return StringNode(Token('string', '', 0, 0, 0, None, str(value)))

    @staticmethod
    def _check_is_equal(node: BaseNode, value: str) -> bool:
        if isinstance(node, StringNode):
            return bool(node.value == value)
        return False

    @staticmethod
    def _check_regex_matches(node: BaseNode, regex: str) -> bool:
        if isinstance(node, StringNode):
            return re.match(regex, node.value) is not None
        return False

    @classmethod
    def supported_element_nodes(cls) -> T.List[T.Type]:
        return [StringNode]

class MTypeIDList(MTypeList):
    def __init__(self, node: T.Optional[BaseNode] = None):
        super().__init__(node)

    @classmethod
    def _new_element_node(cls, value: str) -> IdNode:
        return IdNode(Token('', '', 0, 0, 0, None, str(value)))

    @staticmethod
    def _check_is_equal(node: BaseNode, value: str) -> bool:
        if isinstance(node, IdNode):
            return bool(node.value == value)
        return False

    @staticmethod
    def _check_regex_matches(node: BaseNode, regex: str) -> bool:
        if isinstance(node, StringNode):
            return re.match(regex, node.value) is not None
        return False

    @classmethod
    def supported_element_nodes(cls) -> T.List[T.Type]:
        return [IdNode]

rewriter_keys: T.Dict[str, T.Dict[str, T.Any]] = {
    'default_options': {
        'operation': (str, None, ['set', 'delete']),
        'options': (dict, {}, None)
    },
    'kwargs': {
        'function': (str, None, None),
        'id': (str, None, None),
        'operation': (str, None, ['set', 'delete', 'add', 'remove', 'remove_regex', 'info']),
        'kwargs': (dict, {}, None)
    },
    'target': {
        'target': (str, None, None),
        'operation': (str, None, ['src_add', 'src_rm', 'target_rm', 'target_add', 'extra_files_add', 'extra_files_rm', 'info']),
        'sources': (list, [], None),
        'subdir': (str, '', None),
        'target_type': (str, 'executable', ['both_libraries', 'executable', 'jar', 'library', 'shared_library', 'shared_module', 'static_library']),
    }
}

rewriter_func_kwargs = {
    'dependency': {
        'language': MTypeStr,
        'method': MTypeStr,
        'native': MTypeBool,
        'not_found_message': MTypeStr,
        'required': MTypeBool,
        'static': MTypeBool,
        'version': MTypeStrList,
        'modules': MTypeStrList
    },
    'target': {
        'build_by_default': MTypeBool,
        'build_rpath': MTypeStr,
        'dependencies': MTypeIDList,
        'gui_app': MTypeBool,
        'link_with': MTypeIDList,
        'export_dynamic': MTypeBool,
        'implib': MTypeBool,
        'install': MTypeBool,
        'install_dir': MTypeStr,
        'install_rpath': MTypeStr,
        'pie': MTypeBool
    },
    'project': {
        'default_options': MTypeStrList,
        'meson_version': MTypeStr,
        'license': MTypeStrList,
        'subproject_dir': MTypeStr,
        'version': MTypeStr
    }
}

class Rewriter:
    info_dump: T.Optional[T.Dict[str, T.Dict[str, T.Any]]]

    def __init__(self, sourcedir: str, generator: str = 'ninja', skip_errors: bool = False):
        self.sourcedir = sourcedir
        self.interpreter = IntrospectionInterpreter(sourcedir, '', generator, visitors = [AstIDGenerator(), AstIndentationGenerator(), AstConditionLevel()])
        self.skip_errors = skip_errors
        self.modified_nodes: T.List[BaseNode] = []
        self.to_remove_nodes: T.List[BaseNode] = []
        self.to_add_nodes: T.List[BaseNode] = []
        self.functions = {
            'default_options': self.process_default_options,
            'kwargs': self.process_kwargs,
            'target': self.process_target,
        }
        self.info_dump = None

    def analyze_meson(self) -> None:
        mlog.log('Analyzing meson file:', mlog.bold(os.path.join(self.sourcedir, environment.build_filename)))
        self.interpreter.analyze()
        mlog.log('  -- Project:', mlog.bold(self.interpreter.project_data['descriptive_name']))
        mlog.log('  -- Version:', mlog.cyan(self.interpreter.project_data['version']))

    def add_info(self, cmd_type: str, cmd_id: str, data: dict) -> None:
        if self.info_dump is None:
            self.info_dump = {}
        if cmd_type not in self.info_dump:
            self.info_dump[cmd_type] = {}
        self.info_dump[cmd_type][cmd_id] = data

    def print_info(self) -> None:
        if self.info_dump is None:
            return
        sys.stdout.write(json.dumps(self.info_dump, indent=2, cls=IntrospectionEncoder))

    def on_error(self) -> T.Tuple[AnsiDecorator, AnsiDecorator]:
        if self.skip_errors:
            return mlog.cyan('-->'), mlog.yellow('skipping')
        return mlog.cyan('-->'), mlog.red('aborting')

    def handle_error(self) -> None:
        if self.skip_errors:
            return None
        raise MesonException('Rewriting the meson.build failed')

    def all_assignments(self, varname: str) -> T.List[BaseNode]:
        assigned_values = []
        for ass in self.interpreter.all_assignment_nodes[varname]:
            if isinstance(ass, PlusAssignmentNode):
                continue
            assert isinstance(ass, AssignmentNode)
            assigned_values.append(ass.value)
        return assigned_values

    def find_target(self, target: str) -> T.Optional[IntrospectionBuildTarget]:
        for i in self.interpreter.targets:
            if target == i.id:
                return i

        potential_tgts = []
        for i in self.interpreter.targets:
            if target == i.name:
                potential_tgts.append(i)

        if not potential_tgts:
            potenial_tgts_1 = self.all_assignments(target)
            potenial_tgts_1 = [self.interpreter.node_to_runtime_value(el) for el in potenial_tgts_1]
            potential_tgts = [el for el in potenial_tgts_1 if isinstance(el, IntrospectionBuildTarget)]

        if not potential_tgts:
            return None
        elif len(potential_tgts) == 1:
            return potential_tgts[0]
        else:
            mlog.error('There are multiple targets matching', mlog.bold(target))
            for i in potential_tgts:
                mlog.error('  -- Target name', mlog.bold(i.name), 'with ID', mlog.bold(i.id))
            mlog.error('Please try again with the unique ID of the target', *self.on_error())
            self.handle_error()
            return None

    def find_dependency(self, dependency: str) -> T.Optional[IntrospectionDependency]:
        potential_deps = []
        for i in self.interpreter.dependencies:
            if i.name == dependency:
                potential_deps.append(i)

        checking_varnames = len(potential_deps) == 0

        if checking_varnames:
            potential_deps1 = self.all_assignments(dependency)
            potential_deps = [self.interpreter.node_to_runtime_value(el) for el in potential_deps1 if isinstance(el, FunctionNode) and el.func_name.value == 'dependency']

        if not potential_deps:
            return None
        elif len(potential_deps) == 1:
            return potential_deps[0]
        else:
            mlog.error('There are multiple dependencies matching', mlog.bold(dependency))
            for i in potential_deps:
                mlog.error('  -- Dependency name', i)
            if checking_varnames:
                mlog.error('Please try again with the name of the dependency', *self.on_error())
            self.handle_error()
            return None

    @RequiredKeys(rewriter_keys['default_options'])
    def process_default_options(self, cmd: T.Dict[str, T.Any]) -> None:
        # First, remove the old values
        kwargs_cmd: T.Dict[str, T.Any] = {
            'function': 'project',
            'id': "/",
            'operation': 'remove_regex',
            'kwargs': {
                'default_options': [f'{x}=.*' for x in cmd['options'].keys()]
            }
        }
        self.process_kwargs(kwargs_cmd)

        # Then add the new values
        if cmd['operation'] != 'set':
            return

        kwargs_cmd['operation'] = 'add'
        kwargs_cmd['kwargs']['default_options'] = []

        cdata = self.interpreter.coredata
        options = {
            **{str(k): v for k, v in cdata.optstore.items()},
            **{str(k): v for k, v in cdata.optstore.items()},
            **{str(k): v for k, v in cdata.optstore.items()},
            **{str(k): v for k, v in cdata.optstore.items()},
            **{str(k): v for k, v in cdata.optstore.items()},
        }

        for key, val in sorted(cmd['options'].items()):
            if key not in options:
                mlog.error('Unknown options', mlog.bold(key), *self.on_error())
                self.handle_error()
                continue

            try:
                val = options[key].validate_value(val)
            except MesonException as e:
                mlog.error('Unable to set', mlog.bold(key), mlog.red(str(e)), *self.on_error())
                self.handle_error()
                continue

            kwargs_cmd['kwargs']['default_options'] += [f'{key}={val}']

        self.process_kwargs(kwargs_cmd)

    @RequiredKeys(rewriter_keys['kwargs'])
    def process_kwargs(self, cmd: T.Dict[str, T.Any]) -> None:
        mlog.log('Processing function type', mlog.bold(cmd['function']), 'with id', mlog.cyan("'" + cmd['id'] + "'"))
        if cmd['function'] not in rewriter_func_kwargs:
            mlog.error('Unknown function type', cmd['function'], *self.on_error())
            return self.handle_error()
        kwargs_def = rewriter_func_kwargs[cmd['function']]

        # Find the function node to modify
        node = None
        arg_node = None
        if cmd['function'] == 'project':
            # msys bash may expand '/' to a path. It will mangle '//' to '/'
            # but in order to keep usage shell-agnostic, also allow `//` as
            # the function ID such that it will work in both msys bash and
            # other shells.
            if {'/', '//'}.isdisjoint({cmd['id']}):
                mlog.error('The ID for the function type project must be "/" or "//" not "' + cmd['id'] + '"', *self.on_error())
                return self.handle_error()
            node = self.interpreter.project_node
            arg_node = node.args
        elif cmd['function'] == 'target':
            tmp_tgt = self.find_target(cmd['id'])
            if tmp_tgt:
                node = tmp_tgt.node
                arg_node = node.args
        elif cmd['function'] == 'dependency':
            tmp_dep = self.find_dependency(cmd['id'])
            if tmp_dep:
                node = tmp_dep.node
                arg_node = node.args
        if not node:
            mlog.error('Unable to find the function node')
        assert isinstance(node, FunctionNode)
        assert isinstance(arg_node, ArgumentNode)
        # Transform the key nodes to plain strings
        kwargs = {T.cast(IdNode, k).value: v for k, v in arg_node.kwargs.items()}

        # Print kwargs info
        if cmd['operation'] == 'info':
            info_data: T.Dict[str, T.Any] = {}
            for key, val in sorted(kwargs.items()):
                info_data[key] = None
                if isinstance(val, ElementaryNode):
                    info_data[key] = val.value
                elif isinstance(val, ArrayNode):
                    data_list = []
                    for i in val.args.arguments:
                        element = None
                        if isinstance(i, ElementaryNode):
                            element = i.value
                        data_list += [element]
                    info_data[key] = data_list

            self.add_info('kwargs', '{}#{}'.format(cmd['function'], cmd['id']), info_data)
            return # Nothing else to do

        # Modify the kwargs
        num_changed = 0
        for key, val in sorted(cmd['kwargs'].items()):
            if key not in kwargs_def:
                mlog.error('Cannot modify unknown kwarg', mlog.bold(key), *self.on_error())
                self.handle_error()
                continue

            if cmd['operation'] == 'delete':
                # Remove the key from the kwargs
                if key not in kwargs:
                    mlog.log('  -- Key', mlog.bold(key), 'is already deleted')
                    continue
                mlog.log('  -- Deleting', mlog.bold(key), 'from the kwargs')
                del kwargs[key]
            elif cmd['operation'] == 'set':
                # Replace the key from the kwargs
                mlog.log('  -- Setting', mlog.bold(key), 'to', mlog.yellow(str(val)))
                kwargs[key] = kwargs_def[key].new_node(val)
            else:
                # Modify the value from the kwargs

                if key not in kwargs:
                    kwargs[key] = None
                modifier = kwargs_def[key](kwargs[key])
                if not modifier.can_modify():
                    mlog.log('  -- Skipping', mlog.bold(key), 'because it is too complex to modify')
                    continue

                # Apply the operation
                val_str = str(val)
                if cmd['operation'] == 'add':
                    mlog.log('  -- Adding', mlog.yellow(val_str), 'to', mlog.bold(key))
                    modifier.add_value(val)
                elif cmd['operation'] == 'remove':
                    mlog.log('  -- Removing', mlog.yellow(val_str), 'from', mlog.bold(key))
                    modifier.remove_value(val)
                elif cmd['operation'] == 'remove_regex':
                    mlog.log('  -- Removing all values matching', mlog.yellow(val_str), 'from', mlog.bold(key))
                    modifier.remove_regex(val)

                # Write back the result
                kwargs[key] = modifier.get_node()

            num_changed += 1

        # Convert the keys back to IdNode's
        arg_node.kwargs = {IdNode(Token('', '', 0, 0, 0, None, k)): v for k, v in kwargs.items()}
        for k, v in arg_node.kwargs.items():
            k.level = v.level
        if num_changed > 0 and node not in self.modified_nodes:
            self.modified_nodes += [node]

    def find_assignment_node(self, node: BaseNode) -> T.Optional[AssignmentNode]:
        for k, v in self.interpreter.all_assignment_nodes.items():
            for ass in v:
                if ass.value == node:
                    return ass
        return None

    def affects_no_other_targets(self, candidate: BaseNode) -> bool:
        affected = self.interpreter.dataflow_dag.reachable({candidate}, False)
        affected_targets = [x for x in affected if isinstance(x, FunctionNode) and x.func_name.value in BUILD_TARGET_FUNCTIONS]
        return len(affected_targets) == 1

    def get_relto(self, target_node: BaseNode, node: BaseNode) -> Path:
        cwd = Path(os.getcwd())
        all_paths = self.interpreter.dataflow_dag.find_all_paths(node, target_node)
        # len(all_paths) == 0 would imply that data does not flow from node to
        # target_node. This would imply that adding sources to node would not
        # add the source to the target.
        assert all_paths
        if len(all_paths) > 1:
            return None
        return (cwd / next(x for x in all_paths[0] if isinstance(x, FunctionNode)).filename).parent

    def add_src_or_extra(self, op: str, target: IntrospectionBuildTarget, newfiles: T.List[str], to_sort_nodes: T.List[T.Union[FunctionNode, ArrayNode]]) -> None:
        assert op in {'src_add', 'extra_files_add'}

        if op == 'src_add':
            old: T.Set[T.Union[BaseNode, UnknownValue]] = set(target.source_nodes)
        elif op == 'extra_files_add':
            if target.extra_files is None:
                old = set()
            else:
                old = {target.extra_files}
            tgt_function: FunctionNode = target.node

        cwd = Path(os.getcwd())
        target_dir_abs = cwd / os.path.dirname(target.node.filename)
        source_root_abs = cwd / self.interpreter.source_root

        candidates1 = self.interpreter.dataflow_dag.reachable(old, True)
        # A node is a member of the set `candidates1` exactly if data from this node
        # flow into one of the `dest` nodes. We assume that this implies that if we
        # add `foo.c` to this node, then 'foo.c' will be added to one of these
        # nodes. This assumption is not always true:
        # ar = ['a.c', 'b.c']
        # srcs = ar[1]
        # executable('name', srcs)
        # Data flows from `ar` to `srcs`, but if we add 'foo.c':
        # ar = ['a.c', 'b.c', 'foo.c']
        # srcs = ar[1]
        # executable('name', srcs)
        # this does not add 'foo.c' to `srcs`. This is a known bug/limitation of
        # the meson rewriter that could be fixed by replacing `reachable` with a
        # more advanced analysis. But this is a lot of work and I think e.g.
        # `srcs = ar[1]` is rare in real-world projects, so I will just leave
        # this for now.

        candidates2 = {x for x in candidates1 if isinstance(x, (FunctionNode, ArrayNode))}

        # If we have this meson.build file:
        # shared = ['shared.c']
        # executable('foo', shared + ['foo.c'])
        # executable('bar', shared + ['bar.c'])
        # and we are tasked with adding 'new.c' to 'foo', we should do e.g this:
        # shared = ['shared.c']
        # executable('foo', shared + ['foo.c', 'new.c'])
        # executable('bar', shared + ['bar.c'])
        # but never this:
        # shared = ['shared.c', 'new.c']
        # executable('foo', shared + ['foo.c'])
        # executable('bar', shared + ['bar.c'])
        # We do this by removing the `['shared.c']`-node from `candidates2`.
        candidates2 = {x for x in candidates2 if self.affects_no_other_targets(x)}

        def path_contains_unknowns(candidate: BaseNode) -> bool:
            all_paths = self.interpreter.dataflow_dag.find_all_paths(candidate, target.node)
            for path in all_paths:
                for el in path:
                    if isinstance(el, UnknownValue):
                        return True
            return False

        candidates2 = {x for x in candidates2 if not path_contains_unknowns(x)}

        candidates2 = {x for x in candidates2 if self.get_relto(target.node, x) is not None}

        chosen: T.Union[FunctionNode, ArrayNode] = None
        new_kwarg_flag = False
        if len(candidates2) > 0:
            # So that files(['a', 'b']) gets modified to files(['a', 'b', 'c']) instead of files(['a', 'b'], 'c')
            if len({x for x in candidates2 if isinstance(x, ArrayNode)}) > 0:
                candidates2 = {x for x in candidates2 if isinstance(x, ArrayNode)}

            # We choose one more or less arbitrary candidate
            chosen = min(candidates2, key=lambda x: (x.lineno, x.colno))
        elif op == 'src_add':
            chosen = target.node
        elif op == 'extra_files_add':
            chosen = ArrayNode(_symbol('['), ArgumentNode(Token('', tgt_function.filename, 0, 0, 0, None, '[]')), _symbol(']'))

            # this is fundamentally error prone
            self.interpreter.dataflow_dag.add_edge(chosen, target.node)

            extra_files_idnode = IdNode(Token('string', tgt_function.filename, 0, 0, 0, None, 'extra_files'))
            if tgt_function not in self.modified_nodes:
                self.modified_nodes += [tgt_function]
            new_extra_files_node: BaseNode
            if target.node.args.get_kwarg_or_default('extra_files', None) is None:
                # Target has no extra_files kwarg, create one
                new_kwarg_flag = True
                new_extra_files_node = chosen
            else:
                new_kwarg_flag = True
                old_extra_files = target.node.args.get_kwarg_or_default('extra_files', None)
                target.node.args.kwargs = {k: v for k, v in target.node.args.kwargs.items() if not (isinstance(k, IdNode) and k.value == 'extra_files')}
                new_extra_files_node = ArithmeticNode('add', old_extra_files, _symbol('+'), chosen)

            tgt_function.args.kwargs[extra_files_idnode] = new_extra_files_node

        newfiles_relto = self.get_relto(target.node, chosen)
        old_src_list: T.List[T.Any] = flatten([self.interpreter.node_to_runtime_value(sn) for sn in old])

        if op == 'src_add':
            name = 'Source'
        elif op == 'extra_files_add':
            name = 'Extra file'
        # Generate the new String nodes
        to_append = []
        added = []

        old_src_list = [(target_dir_abs / x).resolve() if isinstance(x, str) else x.to_abs_path(source_root_abs) for x in old_src_list if not isinstance(x, UnknownValue)]
        for _newf in sorted(set(newfiles)):
            newf = Path(_newf)
            if os.path.isabs(newf):
                newf = Path(newf)
            else:
                newf = source_root_abs / newf
            if newf in old_src_list:
                mlog.log('  -- ', name, mlog.green(str(newf)), 'is already defined for the target --> skipping')
                continue

            mlog.log('  -- Adding ', name.lower(), mlog.green(str(newf)), 'at',
                     mlog.yellow(f'{chosen.filename}:{chosen.lineno}'))
            added.append(newf)
            mocktarget = self.interpreter.funcvals[target.node]
            assert isinstance(mocktarget, IntrospectionBuildTarget)
            # print("adding ", str(newf), 'to', mocktarget.name) todo: should we write something to stderr?

            path = relpath(newf, newfiles_relto)
            path = codecs.encode(path, 'unicode_escape').decode() # Because the StringNode constructor does the inverse
            token = Token('string', chosen.filename, 0, 0, 0, None, path)
            to_append += [StringNode(token)]

        assert isinstance(chosen, (FunctionNode, ArrayNode))
        arg_node = chosen.args
        # Append to the AST at the right place
        arg_node.arguments += to_append

        # Mark the node as modified
        if chosen not in to_sort_nodes:
            to_sort_nodes += [chosen]
        # If the extra_files array is newly created, i.e. if new_kwarg_flag is
        # True, don't mark it as its parent function node already is, otherwise
        # this would cause double modification.
        if chosen not in self.modified_nodes and not new_kwarg_flag:
            self.modified_nodes += [chosen]

    # Utility function to get a list of the sources from a node
    def arg_list_from_node(self, n: BaseNode) -> T.List[BaseNode]:
        args = []
        if isinstance(n, FunctionNode):
            args = list(n.args.arguments)
            if n.func_name.value in BUILD_TARGET_FUNCTIONS:
                args.pop(0)
        elif isinstance(n, ArrayNode):
            args = n.args.arguments
        elif isinstance(n, ArgumentNode):
            args = n.arguments
        return args

    def rm_src_or_extra(self, op: str, target: IntrospectionBuildTarget, to_be_removed: T.List[str], to_sort_nodes: T.List[T.Union[FunctionNode, ArrayNode]]) -> None:
        assert op in {'src_rm', 'extra_files_rm'}
        cwd = Path(os.getcwd())
        source_root_abs = cwd / self.interpreter.source_root

        # Helper to find the exact string node and its parent
        def find_node(src: str) -> T.Tuple[T.Optional[BaseNode], T.Optional[StringNode]]:
            if op == 'src_rm':
                nodes = self.interpreter.dataflow_dag.reachable(set(target.source_nodes), True).union({target.node})
            elif op == 'extra_files_rm':
                nodes = self.interpreter.dataflow_dag.reachable({target.extra_files}, True)
            for i in nodes:
                if isinstance(i, UnknownValue):
                    continue
                relto = self.get_relto(target.node, i)
                if relto is not None:
                    for j in self.arg_list_from_node(i):
                        if isinstance(j, StringNode):
                            if os.path.normpath(relto / j.value) == os.path.normpath(source_root_abs / src):
                                return i, j
            return None, None

        if op == 'src_rm':
            name = 'source'
        elif op == 'extra_files_rm':
            name = 'extra file'

        for i in to_be_removed:
            # Try to find the node with the source string
            root, string_node = find_node(i)
            if root is None:
                mlog.warning('  -- Unable to find', name, mlog.green(i), 'in the target')
                continue
            if not self.affects_no_other_targets(string_node):
                mlog.warning('  -- Removing the', name, mlog.green(i), 'is too compilicated')
                continue

            if not isinstance(root, (FunctionNode, ArrayNode)):
                raise NotImplementedError # I'm lazy

            # Remove the found string node from the argument list
            arg_node = root.args
            mlog.log('  -- Removing', name, mlog.green(i), 'from',
                     mlog.yellow(f'{string_node.filename}:{string_node.lineno}'))
            arg_node.arguments.remove(string_node)

            # Mark the node as modified
            if root not in to_sort_nodes:
                to_sort_nodes += [root]
            if root not in self.modified_nodes:
                self.modified_nodes += [root]

    @RequiredKeys(rewriter_keys['target'])
    def process_target(self, cmd: T.Dict[str, T.Any]) -> None:
        mlog.log('Processing target', mlog.bold(cmd['target']), 'operation', mlog.cyan(cmd['operation']))
        target = self.find_target(cmd['target'])
        if target is None and cmd['operation'] != 'target_add':
            mlog.error('Unknown target', mlog.bold(cmd['target']), *self.on_error())
            return self.handle_error()

        # Make source paths relative to the current subdir
        def rel_source(src: str) -> str:
            subdir = os.path.abspath(os.path.join(self.sourcedir, target.subdir))
            if os.path.isabs(src):
                return os.path.relpath(src, subdir)
            elif not os.path.exists(src):
                return src # Trust the user when the source doesn't exist
            # Make sure that the path is relative to the subdir
            return os.path.relpath(os.path.abspath(src), subdir)

        if target is not None:
            cmd['sources'] = [rel_source(x) for x in cmd['sources']]

        to_sort_nodes: T.List[T.Union[FunctionNode, ArrayNode]] = []

        if cmd['operation'] in {'src_add', 'extra_files_add'}:
            self.add_src_or_extra(cmd['operation'], target, cmd['sources'], to_sort_nodes)

        elif cmd['operation'] in {'src_rm', 'extra_files_rm'}:
            self.rm_src_or_extra(cmd['operation'], target, cmd['sources'], to_sort_nodes)

        elif cmd['operation'] == 'target_add':
            if target is not None:
                mlog.error('Can not add target', mlog.bold(cmd['target']), 'because it already exists', *self.on_error())
                return self.handle_error()

            id_base = re.sub(r'[- ]', '_', cmd['target'])
            target_id = id_base + '_exe' if cmd['target_type'] == 'executable' else '_lib'
            source_id = id_base + '_sources'
            filename = os.path.join(os.getcwd(), self.interpreter.source_root, cmd['subdir'], environment.build_filename)

            # Build src list
            src_arg_node = ArgumentNode(Token('string', filename, 0, 0, 0, None, ''))
            src_arr_node = ArrayNode(_symbol('['), src_arg_node, _symbol(']'))
            src_far_node = ArgumentNode(Token('string', filename, 0, 0, 0, None, ''))
            src_fun_node = FunctionNode(IdNode(Token('id', filename, 0, 0, 0, (0, 0), 'files')), _symbol('('), src_far_node, _symbol(')'))
            src_ass_node = AssignmentNode(IdNode(Token('id', filename, 0, 0, 0, (0, 0), source_id)), _symbol('='), src_fun_node)
            src_arg_node.arguments = [StringNode(Token('string', filename, 0, 0, 0, None, x)) for x in cmd['sources']]
            src_far_node.arguments = [src_arr_node]

            # Build target
            tgt_arg_node = ArgumentNode(Token('string', filename, 0, 0, 0, None, ''))
            tgt_fun_node = FunctionNode(IdNode(Token('id', filename, 0, 0, 0, (0, 0), cmd['target_type'])), _symbol('('), tgt_arg_node, _symbol(')'))
            tgt_ass_node = AssignmentNode(IdNode(Token('id', filename, 0, 0, 0, (0, 0), target_id)), _symbol('='), tgt_fun_node)
            tgt_arg_node.arguments = [
                StringNode(Token('string', filename, 0, 0, 0, None, cmd['target'])),
                IdNode(Token('string', filename, 0, 0, 0, None, source_id))
            ]

            src_ass_node.accept(AstIndentationGenerator())
            tgt_ass_node.accept(AstIndentationGenerator())
            self.to_add_nodes += [src_ass_node, tgt_ass_node]

        elif cmd['operation'] == 'target_rm':
            to_remove: BaseNode = self.find_assignment_node(target.node)
            if to_remove is None:
                to_remove = target.node
            self.to_remove_nodes += [to_remove]
            mlog.log('  -- Removing target', mlog.green(cmd['target']), 'at',
                     mlog.yellow(f'{to_remove.filename}:{to_remove.lineno}'))

        elif cmd['operation'] == 'info':
            # T.List all sources in the target

            cwd = Path(os.getcwd())
            source_root_abs = cwd / self.interpreter.source_root

            src_list = self.interpreter.nodes_to_pretty_filelist(source_root_abs, target.subdir, target.source_nodes)
            extra_files_list = self.interpreter.nodes_to_pretty_filelist(source_root_abs, target.subdir, [target.extra_files] if target.extra_files else [])

            src_list = ['unknown' if isinstance(x, UnknownValue) else relpath(x, source_root_abs) for x in src_list]
            extra_files_list = ['unknown' if isinstance(x, UnknownValue) else relpath(x, source_root_abs) for x in extra_files_list]

            test_data = {
                'name': target.name,
                'sources': src_list,
                'extra_files': extra_files_list
            }
            self.add_info('target', target.id, test_data)

        # Sort files
        for i in to_sort_nodes:
            def convert(text: str) -> T.Union[int, str]:
                return int(text) if text.isdigit() else text.lower()

            def alphanum_key(key: str) -> T.List[T.Union[int, str]]:
                return [convert(c) for c in re.split('([0-9]+)', key)]

            def path_sorter(key: str) -> T.List[T.Tuple[bool, T.List[T.Union[int, str]]]]:
                return [(key.count('/') <= idx, alphanum_key(x)) for idx, x in enumerate(key.split('/'))]

            if isinstance(i, FunctionNode) and i.func_name.value in BUILD_TARGET_FUNCTIONS:
                src_args = i.args.arguments[1:]
                target_name = [i.args.arguments[0]]
            else:
                src_args = i.args.arguments
                target_name = []
            unknown: T.List[BaseNode] = [x for x in src_args if not isinstance(x, StringNode)]
            sources: T.List[StringNode] = [x for x in src_args if isinstance(x, StringNode)]
            sources = sorted(sources, key=lambda x: path_sorter(x.value))
            i.args.arguments = target_name + unknown + T.cast(T.List[BaseNode], sources)

    def process(self, cmd: T.Dict[str, T.Any]) -> None:
        if 'type' not in cmd:
            raise RewriterException('Command has no key "type"')
        if cmd['type'] not in self.functions:
            raise RewriterException('Unknown command "{}". Supported commands are: {}'
                                    .format(cmd['type'], list(self.functions.keys())))
        self.functions[cmd['type']](cmd)

    def apply_changes(self) -> None:
        assert all(hasattr(x, 'lineno') and hasattr(x, 'colno') and hasattr(x, 'filename') for x in self.modified_nodes)
        assert all(hasattr(x, 'lineno') and hasattr(x, 'colno') and hasattr(x, 'filename') for x in self.to_remove_nodes)
        assert all(isinstance(x, (ArrayNode, FunctionNode)) for x in self.modified_nodes)
        assert all(isinstance(x, (ArrayNode, AssignmentNode, FunctionNode)) for x in self.to_remove_nodes)
        # Sort based on line and column in reversed order
        work_nodes = [{'node': x, 'action': 'modify'} for x in self.modified_nodes]
        work_nodes += [{'node': x, 'action': 'rm'} for x in self.to_remove_nodes]
        work_nodes = sorted(work_nodes, key=lambda x: (T.cast(BaseNode, x['node']).lineno, T.cast(BaseNode, x['node']).colno), reverse=True)
        work_nodes += [{'node': x, 'action': 'add'} for x in self.to_add_nodes]

        # Generating the new replacement string
        str_list = []
        for i in work_nodes:
            new_data = ''
            if i['action'] == 'modify' or i['action'] == 'add':
                printer = AstPrinter()
                T.cast(BaseNode, i['node']).accept(printer)
                printer.post_process()
                new_data = printer.result.strip()
            data = {
                'file': T.cast(BaseNode, i['node']).filename,
                'str': new_data,
                'node': i['node'],
                'action': i['action']
            }
            str_list += [data]

        # Load build files
        files: T.Dict[str, T.Any] = {}
        for i in str_list:
            if i['file'] in files:
                continue
            fpath = os.path.realpath(T.cast(str, i['file']))
            fdata = ''
            # Create an empty file if it does not exist
            if not os.path.exists(fpath):
                with open(fpath, 'w', encoding='utf-8'):
                    pass
            with open(fpath, encoding='utf-8') as fp:
                fdata = fp.read()

            # Generate line offsets numbers
            m_lines = fdata.splitlines(True)
            offset = 0
            line_offsets = []
            for j in m_lines:
                line_offsets += [offset]
                offset += len(j)

            files[T.cast(str, i['file'])] = {
                'path': fpath,
                'raw': fdata,
                'offsets': line_offsets
            }

        # Replace in source code
        def remove_node(i: T.Dict[str, T.Any]) -> None:
            offsets = files[i['file']]['offsets']
            raw = files[i['file']]['raw']
            node = i['node']
            line = node.lineno - 1
            col = node.colno
            start = offsets[line] + col
            end = start
            if isinstance(node, (ArrayNode, FunctionNode)):
                end = offsets[node.end_lineno - 1] + node.end_colno

            # Only removal is supported for assignments
            elif isinstance(node, AssignmentNode) and i['action'] == 'rm':
                if isinstance(node.value, (ArrayNode, FunctionNode)):
                    remove_node({'file': i['file'], 'str': '', 'node': node.value, 'action': 'rm'})
                    raw = files[i['file']]['raw']
                while raw[end] != '=':
                    end += 1
                end += 1 # Handle the '='
                while raw[end] in {' ', '\n', '\t'}:
                    end += 1

            files[i['file']]['raw'] = raw[:start] + i['str'] + raw[end:]

        for i in str_list:
            if i['action'] in {'modify', 'rm'}:
                remove_node(i)
            elif i['action'] == 'add':
                files[T.cast(str, i['file'])]['raw'] += T.cast(str, i['str']) + '\n'

        # Write the files back
        for key, val in files.items():
            mlog.log('Rewriting', mlog.yellow(key))
            with open(val['path'], 'w', encoding='utf-8') as fp:
                fp.write(val['raw'])

target_operation_map = {
    'add': 'src_add',
    'rm': 'src_rm',
    'add_target': 'target_add',
    'rm_target': 'target_rm',
    'add_extra_files': 'extra_files_add',
    'rm_extra_files': 'extra_files_rm',
    'info': 'info',
}

def list_to_dict(in_list: T.List[str]) -> T.Dict[str, str]:
    result = {}
    it = iter(in_list)
    try:
        for i in it:
            # calling next(it) is not a mistake, we're taking the next element from
            # the iterator, avoiding the need to preprocess it into a sequence of
            # key value pairs.
            result[i] = next(it)
    except StopIteration:
        raise TypeError('in_list parameter of list_to_dict must have an even length.')
    return result

def generate_target(options: argparse.Namespace) -> T.List[T.Dict[str, T.Any]]:
    return [{
        'type': 'target',
        'target': options.target,
        'operation': target_operation_map[options.operation],
        'sources': options.sources,
        'subdir': options.subdir,
        'target_type': options.tgt_type,
    }]

def generate_kwargs(options: argparse.Namespace) -> T.List[T.Dict[str, T.Any]]:
    return [{
        'type': 'kwargs',
        'function': options.function,
        'id': options.id,
        'operation': options.operation,
        'kwargs': list_to_dict(options.kwargs),
    }]

def generate_def_opts(options: argparse.Namespace) -> T.List[T.Dict[str, T.Any]]:
    return [{
        'type': 'default_options',
        'operation': options.operation,
        'options': list_to_dict(options.options),
    }]

def generate_cmd(options: argparse.Namespace) -> T.List[T.Dict[str, T.Any]]:
    if os.path.exists(options.json):
        with open(options.json, encoding='utf-8') as fp:
            return T.cast(T.List[T.Dict[str, T.Any]], json.load(fp))
    else:
        return T.cast(T.List[T.Dict[str, T.Any]], json.loads(options.json))

# Map options.type to the actual type name
cli_type_map = {
    'target': generate_target,
    'tgt': generate_target,
    'kwargs': generate_kwargs,
    'default-options': generate_def_opts,
    'def': generate_def_opts,
    'command': generate_cmd,
    'cmd': generate_cmd,
}

def run(options: argparse.Namespace) -> int:
    mlog.redirect(True)
    if not options.verbose:
        mlog.set_quiet()

    try:
        setup_vsenv()
        rewriter = Rewriter(options.sourcedir, skip_errors=options.skip)
        rewriter.analyze_meson()

        if options.type is None:
            mlog.error('No command specified')
            return 1

        commands = cli_type_map[options.type](options)

        if not isinstance(commands, list):
            raise TypeError('Command is not a list')

        for i, cmd in enumerate(commands):
            if not isinstance(cmd, object):
                raise TypeError('Command is not an object')
            rewriter.process(cmd)
            rewriter.apply_changes()

            if i == len(commands) - 1: # Improves the performance, is not necessary for correctness.
                break

            rewriter.modified_nodes = []
            rewriter.to_remove_nodes = []
            rewriter.to_add_nodes = []
            # The AST changed, so we need to update every information that was derived from the AST
            rewriter.interpreter = IntrospectionInterpreter(rewriter.sourcedir, '', rewriter.interpreter.backend, visitors = [AstIDGenerator(), AstIndentationGenerator(), AstConditionLevel()])
            rewriter.analyze_meson()

        rewriter.print_info()
        return 0
    except Exception as e:
        raise e
    finally:
        mlog.set_verbose()
