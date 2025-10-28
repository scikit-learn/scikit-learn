# SPDX-License-Identifier: Apache-2.0
# Copyright 2018 The Meson development team
# Copyright © 2024-2025 Intel Corporation

# This class contains the basic functionality needed to run any interpreter
# or an interpreter-based tool

from __future__ import annotations
import os
import typing as T

from .. import compilers, environment, mesonlib, options
from ..build import Executable, Jar, SharedLibrary, SharedModule, StaticLibrary
from ..compilers import detect_compiler_for
from ..interpreterbase import InvalidArguments, SubProject, UnknownValue
from ..mesonlib import MachineChoice
from ..options import OptionKey
from ..mparser import BaseNode, ArrayNode, ElementaryNode, IdNode, FunctionNode, StringNode
from .interpreter import AstInterpreter, IntrospectionBuildTarget, IntrospectionDependency

if T.TYPE_CHECKING:
    from ..build import BuildTarget
    from ..interpreterbase import TYPE_var
    from .visitor import AstVisitor


# TODO: it would be nice to not have to duplicate this
BUILD_TARGET_FUNCTIONS = [
    'executable', 'jar', 'library', 'shared_library', 'shared_module',
    'static_library', 'both_libraries'
]

class IntrospectionHelper:
    # mimic an argparse namespace
    def __init__(self, cross_file: T.Optional[str]):
        self.cross_file = [cross_file] if cross_file is not None else []
        self.native_file: T.List[str] = []
        self.cmd_line_options: T.Dict[OptionKey, str] = {}
        self.projectoptions: T.List[str] = []

    def __eq__(self, other: object) -> bool:
        return NotImplemented

class IntrospectionInterpreter(AstInterpreter):
    # If you run `meson setup ...` the `Interpreter`-class walks over the AST.
    # If you run `meson rewrite ...` and `meson introspect meson.build ...`,
    # the `AstInterpreter`-class walks over the AST.
    # Works without a build directory.
    # Most of the code is stolen from interpreter.Interpreter .
    def __init__(self,
                 source_root: str,
                 subdir: str,
                 backend: str,
                 visitors: T.Optional[T.List[AstVisitor]] = None,
                 cross_file: T.Optional[str] = None,
                 subproject: SubProject = SubProject(''),
                 subproject_dir: str = 'subprojects',
                 env: T.Optional[environment.Environment] = None):
        options = IntrospectionHelper(cross_file)
        env_ = env or environment.Environment(source_root, None, options)
        super().__init__(source_root, subdir, subproject, subproject_dir, env_, visitors=visitors)

        self.cross_file = cross_file
        self.backend = backend
        self.project_data: T.Dict[str, T.Any] = {}
        self.targets: T.List[IntrospectionBuildTarget] = []
        self.dependencies: T.List[IntrospectionDependency] = []
        self.project_node: FunctionNode = None

        self.funcs.update({
            'add_languages': self.func_add_languages,
            'dependency': self.func_dependency,
            'executable': self.func_executable,
            'jar': self.func_jar,
            'library': self.func_library,
            'project': self.func_project,
            'shared_library': self.func_shared_lib,
            'shared_module': self.func_shared_module,
            'static_library': self.func_static_lib,
            'both_libraries': self.func_both_lib,
        })

    def func_project(self, node: BaseNode, args: T.List[TYPE_var], kwargs: T.Dict[str, TYPE_var]) -> None:
        if self.project_node:
            raise InvalidArguments('Second call to project()')
        assert isinstance(node, FunctionNode)
        self.project_node = node
        if len(args) < 1:
            raise InvalidArguments('Not enough arguments to project(). Needs at least the project name.')

        def _str_list(node: T.Any) -> T.Optional[T.List[str]]:
            if isinstance(node, ArrayNode):
                r = []
                for v in node.args.arguments:
                    if not isinstance(v, StringNode):
                        return None
                    r.append(v.value)
                return r
            if isinstance(node, StringNode):
                return [node.value]
            return None

        def create_options_dict(options: T.List[str], subproject: str = '') -> T.Mapping[OptionKey, str]:
            result: T.MutableMapping[OptionKey, str] = {}
            for o in options:
                try:
                    (key, value) = o.split('=', 1)
                except ValueError:
                    raise mesonlib.MesonException(f'Option {o!r} must have a value separated by equals sign.')
                result[OptionKey(key)] = value
            return result

        proj_name = args[0]
        proj_vers = kwargs.get('version', 'undefined')
        if isinstance(proj_vers, ElementaryNode):
            proj_vers = proj_vers.value
        if not isinstance(proj_vers, str):
            proj_vers = 'undefined'
        proj_langs = self.flatten_args(args[1:])
        # Match the value returned by ``meson.project_license()`` when
        # no ``license`` argument is specified in the ``project()`` call.
        proj_license = _str_list(kwargs.get('license', None)) or ['unknown']
        proj_license_files = _str_list(kwargs.get('license_files', None)) or []
        self.project_data = {'descriptive_name': proj_name, 'version': proj_vers, 'license': proj_license, 'license_files': proj_license_files}

        self._load_option_file()

        if not self.is_subproject() and 'subproject_dir' in kwargs:
            spdirname = kwargs['subproject_dir']
            if isinstance(spdirname, StringNode):
                assert isinstance(spdirname.value, str)
                self.subproject_dir = spdirname.value
        if not self.is_subproject():
            self.project_data['subprojects'] = []
            subprojects_dir = os.path.join(self.source_root, self.subproject_dir)
            if os.path.isdir(subprojects_dir):
                for i in os.listdir(subprojects_dir):
                    if os.path.isdir(os.path.join(subprojects_dir, i)):
                        self.do_subproject(SubProject(i))

        self.coredata.init_backend_options(self.backend)
        options = {k: v for k, v in self.environment.options.items() if self.environment.coredata.optstore.is_backend_option(k)}

        self.coredata.set_options(options)
        self._add_languages(proj_langs, True, MachineChoice.HOST)
        self._add_languages(proj_langs, True, MachineChoice.BUILD)

    def do_subproject(self, dirname: SubProject) -> None:
        subproject_dir_abs = os.path.join(self.environment.get_source_dir(), self.subproject_dir)
        subpr = os.path.join(subproject_dir_abs, dirname)
        try:
            subi = IntrospectionInterpreter(subpr, '', self.backend, cross_file=self.cross_file, subproject=dirname, subproject_dir=self.subproject_dir, env=self.environment, visitors=self.visitors)
            subi.analyze()
            subi.project_data['name'] = dirname
            self.project_data['subprojects'] += [subi.project_data]
        except (mesonlib.MesonException, RuntimeError):
            pass

    def func_add_languages(self, node: BaseNode, args: T.List[TYPE_var], kwargs: T.Dict[str, TYPE_var]) -> UnknownValue:
        kwargs = self.flatten_kwargs(kwargs)
        required = kwargs.get('required', True)
        assert isinstance(required, (bool, options.UserFeatureOption, UnknownValue)), 'for mypy'
        if isinstance(required, options.UserFeatureOption):
            required = required.is_enabled()
        if 'native' in kwargs:
            native = kwargs.get('native', False)
            self._add_languages(args, required, MachineChoice.BUILD if native else MachineChoice.HOST)
        else:
            for for_machine in [MachineChoice.BUILD, MachineChoice.HOST]:
                self._add_languages(args, required, for_machine)
        return UnknownValue()

    def _add_languages(self, raw_langs: T.List[TYPE_var], required: T.Union[bool, UnknownValue], for_machine: MachineChoice) -> None:
        langs: T.List[str] = []
        for l in self.flatten_args(raw_langs):
            if isinstance(l, str):
                langs.append(l)
            elif isinstance(l, StringNode):
                langs.append(l.value)

        for lang in sorted(langs, key=compilers.sort_clink):
            lang = lang.lower()
            if lang not in self.coredata.compilers[for_machine]:
                try:
                    comp = detect_compiler_for(self.environment, lang, for_machine, True, self.subproject)
                except mesonlib.MesonException:
                    # do we even care about introspecting this language?
                    if isinstance(required, UnknownValue) or required:
                        raise
                    else:
                        continue
                if comp:
                    self.coredata.process_compiler_options(lang, comp, self.subproject)

    def func_dependency(self, node: BaseNode, args: T.List[TYPE_var], kwargs: T.Dict[str, TYPE_var]) -> T.Optional[IntrospectionDependency]:
        assert isinstance(node, FunctionNode)
        args = self.flatten_args(args)
        kwargs = self.flatten_kwargs(kwargs)
        if not args:
            return None
        name = args[0]
        assert isinstance(name, (str, UnknownValue))
        has_fallback = 'fallback' in kwargs
        required = kwargs.get('required', True)
        version = kwargs.get('version', [])
        if not isinstance(version, list):
            version = [version]
        if any(isinstance(el, UnknownValue) for el in version):
            version = UnknownValue()
        else:
            assert all(isinstance(el, str) for el in version)
            version = T.cast(T.List[str], version)
        assert isinstance(required, (bool, UnknownValue))
        newdep = IntrospectionDependency(
            name=name,
            required=required,
            version=version,
            has_fallback=has_fallback,
            conditional=node.condition_level > 0,
            node=node)
        self.dependencies += [newdep]
        return newdep

    def build_target(self, node: BaseNode, args: T.List[TYPE_var], kwargs_raw: T.Dict[str, TYPE_var], targetclass: T.Type[BuildTarget]) -> T.Union[IntrospectionBuildTarget, UnknownValue]:
        assert isinstance(node, FunctionNode)
        args = self.flatten_args(args)
        if not args or not isinstance(args[0], str):
            return UnknownValue()
        name = args[0]
        srcqueue: T.List[BaseNode] = [node]
        extra_queue = []

        # Process the sources BEFORE flattening the kwargs, to preserve the original nodes
        if 'sources' in kwargs_raw:
            srcqueue += mesonlib.listify(kwargs_raw['sources'])

        if 'extra_files' in kwargs_raw:
            extra_queue += mesonlib.listify(kwargs_raw['extra_files'])

        kwargs = self.flatten_kwargs(kwargs_raw, True)

        oldlen = len(node.args.arguments)
        source_nodes = node.args.arguments[1:]
        for k, v in node.args.kwargs.items():
            assert isinstance(k, IdNode)
            if k.value == 'sources':
                source_nodes.append(v)
        assert oldlen == len(node.args.arguments)

        extraf_nodes = None
        for k, v in node.args.kwargs.items():
            assert isinstance(k, IdNode)
            if k.value == 'extra_files':
                assert extraf_nodes is None
                extraf_nodes = v

        # Make sure nothing can crash when creating the build class
        kwargs_reduced = {k: v for k, v in kwargs.items() if k in targetclass.known_kwargs and k in {'install', 'build_by_default', 'build_always', 'name_prefix'}}
        kwargs_reduced = {k: v.value if isinstance(v, ElementaryNode) else v for k, v in kwargs_reduced.items()}
        kwargs_reduced = {k: v for k, v in kwargs_reduced.items() if not isinstance(v, BaseNode)}
        for_machine = MachineChoice.BUILD if kwargs.get('native', False) else MachineChoice.HOST
        objects: T.List[T.Any] = []
        empty_sources: T.List[T.Any] = []
        # Passing the unresolved sources list causes errors
        kwargs_reduced['_allow_no_sources'] = True
        target = targetclass(name, self.subdir, self.subproject, for_machine, empty_sources, None, objects,
                             self.environment, self.coredata.compilers[for_machine], kwargs_reduced)
        target.process_compilers_late()

        new_target = IntrospectionBuildTarget(
            name=target.get_basename(),
            machine=target.for_machine.get_lower_case_name(),
            id=target.get_id(),
            typename=target.get_typename(),
            defined_in=os.path.normpath(os.path.join(self.source_root, self.subdir, environment.build_filename)),
            subdir=self.subdir,
            build_by_default=target.build_by_default,
            installed=target.should_install(),
            outputs=target.get_outputs(),
            source_nodes=source_nodes,
            extra_files=extraf_nodes,
            kwargs=kwargs,
            node=node)

        self.targets += [new_target]
        return new_target

    def build_library(self, node: BaseNode, args: T.List[TYPE_var], kwargs: T.Dict[str, TYPE_var]) -> T.Union[IntrospectionBuildTarget, UnknownValue]:
        default_library = self.coredata.optstore.get_value_for(OptionKey('default_library', subproject=self.subproject))
        if default_library == 'shared':
            return self.build_target(node, args, kwargs, SharedLibrary)
        elif default_library == 'static':
            return self.build_target(node, args, kwargs, StaticLibrary)
        elif default_library == 'both':
            return self.build_target(node, args, kwargs, SharedLibrary)
        return None

    def func_executable(self, node: BaseNode, args: T.List[TYPE_var], kwargs: T.Dict[str, TYPE_var]) -> T.Union[IntrospectionBuildTarget, UnknownValue]:
        return self.build_target(node, args, kwargs, Executable)

    def func_static_lib(self, node: BaseNode, args: T.List[TYPE_var], kwargs: T.Dict[str, TYPE_var]) -> T.Union[IntrospectionBuildTarget, UnknownValue]:
        return self.build_target(node, args, kwargs, StaticLibrary)

    def func_shared_lib(self, node: BaseNode, args: T.List[TYPE_var], kwargs: T.Dict[str, TYPE_var]) -> T.Union[IntrospectionBuildTarget, UnknownValue]:
        return self.build_target(node, args, kwargs, SharedLibrary)

    def func_both_lib(self, node: BaseNode, args: T.List[TYPE_var], kwargs: T.Dict[str, TYPE_var]) -> T.Union[IntrospectionBuildTarget, UnknownValue]:
        return self.build_target(node, args, kwargs, SharedLibrary)

    def func_shared_module(self, node: BaseNode, args: T.List[TYPE_var], kwargs: T.Dict[str, TYPE_var]) -> T.Union[IntrospectionBuildTarget, UnknownValue]:
        return self.build_target(node, args, kwargs, SharedModule)

    def func_library(self, node: BaseNode, args: T.List[TYPE_var], kwargs: T.Dict[str, TYPE_var]) -> T.Union[IntrospectionBuildTarget, UnknownValue]:
        return self.build_library(node, args, kwargs)

    def func_jar(self, node: BaseNode, args: T.List[TYPE_var], kwargs: T.Dict[str, TYPE_var]) -> T.Union[IntrospectionBuildTarget, UnknownValue]:
        return self.build_target(node, args, kwargs, Jar)

    def func_build_target(self, node: BaseNode, args: T.List[TYPE_var], kwargs: T.Dict[str, TYPE_var]) -> T.Union[IntrospectionBuildTarget, UnknownValue]:
        if 'target_type' not in kwargs:
            return None
        target_type = kwargs.pop('target_type')
        if isinstance(target_type, ElementaryNode):
            target_type = target_type.value
        if target_type == 'executable':
            return self.build_target(node, args, kwargs, Executable)
        elif target_type == 'shared_library':
            return self.build_target(node, args, kwargs, SharedLibrary)
        elif target_type == 'static_library':
            return self.build_target(node, args, kwargs, StaticLibrary)
        elif target_type == 'both_libraries':
            return self.build_target(node, args, kwargs, SharedLibrary)
        elif target_type == 'library':
            return self.build_library(node, args, kwargs)
        elif target_type == 'jar':
            return self.build_target(node, args, kwargs, Jar)
        return None

    def is_subproject(self) -> bool:
        return self.subproject != ''

    def analyze(self) -> None:
        self.load_root_meson_file()
        self.sanity_check_ast()
        self.parse_project()
        self.run()

    def extract_subproject_dir(self) -> T.Optional[str]:
        '''Fast path to extract subproject_dir kwarg.
           This is faster than self.parse_project() which also initialize options
           and also calls parse_project() on every subproject.
        '''
        if not self.ast.lines:
            return None
        project = self.ast.lines[0]
        # first line is always project()
        if not isinstance(project, FunctionNode):
            return None
        for kw, val in project.args.kwargs.items():
            assert isinstance(kw, IdNode), 'for mypy'
            if kw.value == 'subproject_dir':
                # mypy does not understand "and isinstance"
                if isinstance(val, StringNode):
                    return val.value
        return None

    def flatten_kwargs(self, kwargs: T.Dict[str, TYPE_var], include_unknown_args: bool = False) -> T.Dict[str, TYPE_var]:
        flattened_kwargs = {}
        for key, val in kwargs.items():
            if isinstance(val, BaseNode):
                resolved = self.node_to_runtime_value(val)
                if resolved is not None:
                    flattened_kwargs[key] = resolved
            elif isinstance(val, (str, bool, int, float)) or include_unknown_args:
                flattened_kwargs[key] = val
        return flattened_kwargs
