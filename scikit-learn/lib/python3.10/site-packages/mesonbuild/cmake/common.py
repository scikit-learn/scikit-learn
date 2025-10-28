# SPDX-License-Identifier: Apache-2.0
# Copyright 2019 The Meson development team

from __future__ import annotations

from ..mesonlib import MesonException
from ..options import OptionKey
from .. import mlog
from pathlib import Path
import typing as T

if T.TYPE_CHECKING:
    from ..environment import Environment
    from ..interpreterbase import TYPE_var

language_map = {
    'c': 'C',
    'cpp': 'CXX',
    'cuda': 'CUDA',
    'objc': 'OBJC',
    'objcpp': 'OBJCXX',
    'nasm': 'ASM_NASM',
    'cs': 'CSharp',
    'java': 'Java',
    'fortran': 'Fortran',
    'swift': 'Swift',
}

backend_generator_map = {
    'ninja': 'Ninja',
    'xcode': 'Xcode',
    'vs2010': 'Visual Studio 10 2010',
    'vs2012': 'Visual Studio 11 2012',
    'vs2013': 'Visual Studio 12 2013',
    'vs2015': 'Visual Studio 14 2015',
    'vs2017': 'Visual Studio 15 2017',
    'vs2019': 'Visual Studio 16 2019',
    'vs2022': 'Visual Studio 17 2022',
}

blacklist_cmake_defs = [
    'CMAKE_TOOLCHAIN_FILE',
    'CMAKE_PROJECT_INCLUDE',
    'MESON_PRELOAD_FILE',
    'MESON_PS_CMAKE_CURRENT_BINARY_DIR',
    'MESON_PS_CMAKE_CURRENT_SOURCE_DIR',
    'MESON_PS_DELAYED_CALLS',
    'MESON_PS_LOADED',
    'MESON_FIND_ROOT_PATH',
    'MESON_CMAKE_SYSROOT',
    'MESON_PATHS_LIST',
    'MESON_CMAKE_ROOT',
]

def cmake_is_debug(env: 'Environment') -> bool:
    if 'b_vscrt' in env.coredata.optstore:
        is_debug = env.coredata.optstore.get_value_for('buildtype') == 'debug'
        if env.coredata.optstore.get_value_for('b_vscrt') in {'mdd', 'mtd'}:
            is_debug = True
        return is_debug
    else:
        # Don't directly assign to is_debug to make mypy happy
        debug_opt = env.coredata.optstore.get_value_for('debug')
        assert isinstance(debug_opt, bool)
        return debug_opt

class CMakeException(MesonException):
    pass

class CMakeBuildFile:
    def __init__(self, file: Path, is_cmake: bool, is_temp: bool) -> None:
        self.file = file
        self.is_cmake = is_cmake
        self.is_temp = is_temp

    def __repr__(self) -> str:
        return f'<{self.__class__.__name__}: {self.file}; cmake={self.is_cmake}; temp={self.is_temp}>'

def _flags_to_list(raw: str) -> T.List[str]:
    # Convert a raw commandline string into a list of strings
    res = []
    curr = ''
    escape = False
    in_string = False
    for i in raw:
        if escape:
            # If the current char is not a quote, the '\' is probably important
            if i not in ['"', "'"]:
                curr += '\\'
            curr += i
            escape = False
        elif i == '\\':
            escape = True
        elif i in {'"', "'"}:
            in_string = not in_string
        elif i in {' ', '\n'}:
            if in_string:
                curr += i
            else:
                res += [curr]
                curr = ''
        else:
            curr += i
    res += [curr]
    res = [r for r in res if len(r) > 0]
    return res

def cmake_get_generator_args(env: 'Environment') -> T.List[str]:
    backend_name = env.coredata.optstore.get_value_for(OptionKey('backend'))
    assert isinstance(backend_name, str)
    assert backend_name in backend_generator_map
    return ['-G', backend_generator_map[backend_name]]

def cmake_defines_to_args(raw: T.List[T.Dict[str, TYPE_var]], permissive: bool = False) -> T.List[str]:
    res: T.List[str] = []

    for i in raw:
        for key, val in i.items():
            if key in blacklist_cmake_defs:
                mlog.warning('Setting', mlog.bold(key), 'is not supported. See the meson docs for cross compilation support:')
                mlog.warning('  - URL: https://mesonbuild.com/CMake-module.html#cross-compilation')
                mlog.warning('  --> Ignoring this option')
                continue
            if isinstance(val, (str, int, float)):
                res += [f'-D{key}={val}']
            elif isinstance(val, bool):
                val_str = 'ON' if val else 'OFF'
                res += [f'-D{key}={val_str}']
            else:
                raise MesonException('Type "{}" of "{}" is not supported as for a CMake define value'.format(type(val).__name__, key))

    return res

# TODO: this function will become obsolete once the `cmake_args` kwarg is dropped
def check_cmake_args(args: T.List[str]) -> T.List[str]:
    res: T.List[str] = []
    dis = ['-D' + x for x in blacklist_cmake_defs]
    assert dis  # Ensure that dis is not empty.
    for i in args:
        if any(i.startswith(x) for x in dis):
            mlog.warning('Setting', mlog.bold(i), 'is not supported. See the meson docs for cross compilation support:')
            mlog.warning('  - URL: https://mesonbuild.com/CMake-module.html#cross-compilation')
            mlog.warning('  --> Ignoring this option')
            continue
        res += [i]
    return res

class CMakeInclude:
    def __init__(self, path: Path, isSystem: bool = False):
        self.path = path
        self.isSystem = isSystem

    def __repr__(self) -> str:
        return f'<CMakeInclude: {self.path} -- isSystem = {self.isSystem}>'

class CMakeFileGroup:
    def __init__(self, data: T.Dict[str, T.Any]) -> None:
        self.defines: str = data.get('defines', '')
        self.flags = _flags_to_list(data.get('compileFlags', ''))
        self.is_generated: bool = data.get('isGenerated', False)
        self.language: str = data.get('language', 'C')
        self.sources = [Path(x) for x in data.get('sources', [])]

        # Fix the include directories
        self.includes: T.List[CMakeInclude] = []
        for i in data.get('includePath', []):
            if isinstance(i, dict) and 'path' in i:
                isSystem = i.get('isSystem', False)
                assert isinstance(isSystem, bool)
                assert isinstance(i['path'], str)
                self.includes += [CMakeInclude(Path(i['path']), isSystem)]
            elif isinstance(i, str):
                self.includes += [CMakeInclude(Path(i))]

    def log(self) -> None:
        mlog.log('flags        =', mlog.bold(', '.join(self.flags)))
        mlog.log('defines      =', mlog.bold(', '.join(self.defines)))
        mlog.log('includes     =', mlog.bold(', '.join([str(x) for x in self.includes])))
        mlog.log('is_generated =', mlog.bold('true' if self.is_generated else 'false'))
        mlog.log('language     =', mlog.bold(self.language))
        mlog.log('sources:')
        for i in self.sources:
            with mlog.nested():
                mlog.log(i.as_posix())

class CMakeTarget:
    def __init__(self, data: T.Dict[str, T.Any]) -> None:
        self.artifacts = [Path(x) for x in data.get('artifacts', [])]
        self.src_dir = Path(data.get('sourceDirectory', ''))
        self.build_dir = Path(data.get('buildDirectory', ''))
        self.name: str = data.get('name', '')
        self.full_name: str = data.get('fullName', '')
        self.install: bool = data.get('hasInstallRule', False)
        self.install_paths = [Path(x) for x in set(data.get('installPaths', []))]
        self.link_lang: str = data.get('linkerLanguage', '')
        self.link_libraries = _flags_to_list(data.get('linkLibraries', ''))
        self.link_flags = _flags_to_list(data.get('linkFlags', ''))
        self.link_lang_flags = _flags_to_list(data.get('linkLanguageFlags', ''))
        # self.link_path = Path(data.get('linkPath', ''))
        self.type: str = data.get('type', 'EXECUTABLE')
        # self.is_generator_provided: bool = data.get('isGeneratorProvided', False)
        self.files: T.List[CMakeFileGroup] = []

        for i in data.get('fileGroups', []):
            self.files += [CMakeFileGroup(i)]

    def log(self) -> None:
        mlog.log('artifacts             =', mlog.bold(', '.join([x.as_posix() for x in self.artifacts])))
        mlog.log('src_dir               =', mlog.bold(self.src_dir.as_posix()))
        mlog.log('build_dir             =', mlog.bold(self.build_dir.as_posix()))
        mlog.log('name                  =', mlog.bold(self.name))
        mlog.log('full_name             =', mlog.bold(self.full_name))
        mlog.log('install               =', mlog.bold('true' if self.install else 'false'))
        mlog.log('install_paths         =', mlog.bold(', '.join([x.as_posix() for x in self.install_paths])))
        mlog.log('link_lang             =', mlog.bold(self.link_lang))
        mlog.log('link_libraries        =', mlog.bold(', '.join(self.link_libraries)))
        mlog.log('link_flags            =', mlog.bold(', '.join(self.link_flags)))
        mlog.log('link_lang_flags       =', mlog.bold(', '.join(self.link_lang_flags)))
        # mlog.log('link_path             =', mlog.bold(self.link_path))
        mlog.log('type                  =', mlog.bold(self.type))
        # mlog.log('is_generator_provided =', mlog.bold('true' if self.is_generator_provided else 'false'))
        for idx, i in enumerate(self.files):
            mlog.log(f'Files {idx}:')
            with mlog.nested():
                i.log()

class CMakeProject:
    def __init__(self, data: T.Dict[str, T.Any]) -> None:
        self.src_dir = Path(data.get('sourceDirectory', ''))
        self.build_dir = Path(data.get('buildDirectory', ''))
        self.name: str = data.get('name', '')
        self.targets: T.List[CMakeTarget] = []

        for i in data.get('targets', []):
            self.targets += [CMakeTarget(i)]

    def log(self) -> None:
        mlog.log('src_dir   =', mlog.bold(self.src_dir.as_posix()))
        mlog.log('build_dir =', mlog.bold(self.build_dir.as_posix()))
        mlog.log('name      =', mlog.bold(self.name))
        for idx, i in enumerate(self.targets):
            mlog.log(f'Target {idx}:')
            with mlog.nested():
                i.log()

class CMakeConfiguration:
    def __init__(self, data: T.Dict[str, T.Any]) -> None:
        self.name: str = data.get('name', '')
        self.projects: T.List[CMakeProject] = []
        for i in data.get('projects', []):
            self.projects += [CMakeProject(i)]

    def log(self) -> None:
        mlog.log('name =', mlog.bold(self.name))
        for idx, i in enumerate(self.projects):
            mlog.log(f'Project {idx}:')
            with mlog.nested():
                i.log()

class SingleTargetOptions:
    def __init__(self) -> None:
        self.opts: T.Dict[str, str] = {}
        self.lang_args: T.Dict[str, T.List[str]] = {}
        self.link_args: T.List[str] = []
        self.install = 'preserve'

    def set_opt(self, opt: str, val: str) -> None:
        self.opts[opt] = val

    def append_args(self, lang: str, args: T.List[str]) -> None:
        if lang not in self.lang_args:
            self.lang_args[lang] = []
        self.lang_args[lang] += args

    def append_link_args(self, args: T.List[str]) -> None:
        self.link_args += args

    def set_install(self, install: bool) -> None:
        self.install = 'true' if install else 'false'

    def get_override_options(self, initial: T.List[str]) -> T.List[str]:
        res: T.List[str] = []
        for i in initial:
            opt = i[:i.find('=')]
            if opt not in self.opts:
                res += [i]
        res += [f'{k}={v}' for k, v in self.opts.items()]
        return res

    def get_compile_args(self, lang: str, initial: T.List[str]) -> T.List[str]:
        if lang in self.lang_args:
            return initial + self.lang_args[lang]
        return initial

    def get_link_args(self, initial: T.List[str]) -> T.List[str]:
        return initial + self.link_args

    def get_install(self, initial: bool) -> bool:
        return {'preserve': initial, 'true': True, 'false': False}[self.install]

class TargetOptions:
    def __init__(self) -> None:
        self.global_options = SingleTargetOptions()
        self.target_options: T.Dict[str, SingleTargetOptions] = {}

    def __getitem__(self, tgt: str) -> SingleTargetOptions:
        if tgt not in self.target_options:
            self.target_options[tgt] = SingleTargetOptions()
        return self.target_options[tgt]

    def get_override_options(self, tgt: str, initial: T.List[str]) -> T.List[str]:
        initial = self.global_options.get_override_options(initial)
        if tgt in self.target_options:
            initial = self.target_options[tgt].get_override_options(initial)
        return initial

    def get_compile_args(self, tgt: str, lang: str, initial: T.List[str]) -> T.List[str]:
        initial = self.global_options.get_compile_args(lang, initial)
        if tgt in self.target_options:
            initial = self.target_options[tgt].get_compile_args(lang, initial)
        return initial

    def get_link_args(self, tgt: str, initial: T.List[str]) -> T.List[str]:
        initial = self.global_options.get_link_args(initial)
        if tgt in self.target_options:
            initial = self.target_options[tgt].get_link_args(initial)
        return initial

    def get_install(self, tgt: str, initial: bool) -> bool:
        initial = self.global_options.get_install(initial)
        if tgt in self.target_options:
            initial = self.target_options[tgt].get_install(initial)
        return initial
