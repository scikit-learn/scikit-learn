# SPDX-License-Identifier: Apache-2.0
# Copyright 2013-2019 The Meson development team

from __future__ import annotations

import glob
import re
import os
import typing as T
from pathlib import Path

from .. import mesonlib
from .. import mlog
from .base import DependencyException, SystemDependency
from .detect import packages
from ..mesonlib import LibType


if T.TYPE_CHECKING:
    from ..environment import Environment
    from ..compilers import Compiler

    TV_ResultTuple = T.Tuple[T.Optional[str], T.Optional[str], bool]

class CudaDependency(SystemDependency):

    supported_languages = ['cpp', 'c', 'cuda'] # see also _default_language

    def __init__(self, environment: 'Environment', kwargs: T.Dict[str, T.Any]) -> None:
        for_machine = self.get_for_machine_from_kwargs(kwargs)
        compilers = environment.coredata.compilers[for_machine]
        machine = environment.machines[for_machine]
        language = self._detect_language(compilers)

        if language not in self.supported_languages:
            raise DependencyException(f'Language \'{language}\' is not supported by the CUDA Toolkit. Supported languages are {self.supported_languages}.')

        super().__init__('cuda', environment, kwargs, language=language)
        self.lib_modules: T.Dict[str, T.List[str]] = {}
        self.requested_modules = self.get_requested(kwargs)
        if not any(runtime in self.requested_modules for runtime in ['cudart', 'cudart_static']):
            # By default, we prefer to link the static CUDA runtime, since this is what nvcc also does by default:
            # https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#cudart-none-shared-static-cudart
            req_modules = ['cudart']
            if kwargs.get('static', True):
                req_modules = ['cudart_static']
            self.requested_modules = req_modules + self.requested_modules

        (self.cuda_path, self.version, self.is_found) = self._detect_cuda_path_and_version()
        if not self.is_found:
            return

        if not os.path.isabs(self.cuda_path):
            raise DependencyException(f'CUDA Toolkit path must be absolute, got \'{self.cuda_path}\'.')

        # Cuda target directory relative to cuda path.
        if machine.is_linux():
            # E.g. targets/x86_64-linux
            self.target_path = os.path.join('targets', f'{machine.cpu_family}-{machine.system}')
        else:
            self.target_path = '.'

        # nvcc already knows where to find the CUDA Toolkit, but if we're compiling
        # a mixed C/C++/CUDA project, we still need to make the include dir searchable
        if self.language != 'cuda' or len(compilers) > 1:
            self.incdir = os.path.join(self.cuda_path, self.target_path, 'include')
            self.compile_args += [f'-I{self.incdir}']

        arch_libdir = self._detect_arch_libdir()
        self.libdir = os.path.join(self.cuda_path, self.target_path, arch_libdir)
        mlog.debug('CUDA library directory is', mlog.bold(self.libdir))

        if 'static' not in kwargs:
            self.libtype = LibType.PREFER_STATIC

        self.is_found = self._find_requested_libraries()

    @classmethod
    def _detect_language(cls, compilers: T.Dict[str, 'Compiler']) -> str:
        for lang in cls.supported_languages:
            if lang in compilers:
                return lang
        return list(compilers.keys())[0]

    def _detect_cuda_path_and_version(self) -> TV_ResultTuple:
        self.env_var = self._default_path_env_var()
        mlog.debug('Default path env var:', mlog.bold(self.env_var))

        version_reqs = self.version_reqs
        if self.language == 'cuda':
            nvcc_version = self._strip_patch_version(self.get_compiler().version)
            mlog.debug('nvcc version:', mlog.bold(nvcc_version))
            if version_reqs:
                # make sure nvcc version satisfies specified version requirements
                (found_some, not_found, found) = mesonlib.version_compare_many(nvcc_version, version_reqs)
                if not_found:
                    msg = f'The current nvcc version {nvcc_version} does not satisfy the specified CUDA Toolkit version requirements {version_reqs}.'
                    return self._report_dependency_error(msg, (None, None, False))

            # use nvcc version to find a matching CUDA Toolkit
            version_reqs = [f'={nvcc_version}']
        else:
            nvcc_version = None

        paths = [(path, self._cuda_toolkit_version(path), default) for (path, default) in self._cuda_paths()]
        if version_reqs:
            return self._find_matching_toolkit(paths, version_reqs, nvcc_version)

        defaults = [(path, version) for (path, version, default) in paths if default]
        if defaults:
            return (defaults[0][0], defaults[0][1], True)

        platform_msg = 'set the CUDA_PATH environment variable' if self._is_windows() \
            else 'set the CUDA_PATH environment variable/create the \'/usr/local/cuda\' symbolic link'
        msg = f'Please specify the desired CUDA Toolkit version (e.g. dependency(\'cuda\', version : \'>=10.1\')) or {platform_msg} to point to the location of your desired version.'
        return self._report_dependency_error(msg, (None, None, False))

    def _find_matching_toolkit(self, paths: T.List[TV_ResultTuple], version_reqs: T.List[str], nvcc_version: T.Optional[str]) -> TV_ResultTuple:
        # keep the default paths order intact, sort the rest in the descending order
        # according to the toolkit version
        part_func: T.Callable[[TV_ResultTuple], bool] = lambda t: not t[2]
        defaults_it, rest_it = mesonlib.partition(part_func, paths)
        defaults = list(defaults_it)
        paths = defaults + sorted(rest_it, key=lambda t: mesonlib.Version(t[1]), reverse=True)
        mlog.debug(f'Search paths: {paths}')

        if nvcc_version and defaults:
            default_src = f"the {self.env_var} environment variable" if self.env_var else "the \'/usr/local/cuda\' symbolic link"
            nvcc_warning = 'The default CUDA Toolkit as designated by {} ({}) doesn\'t match the current nvcc version {} and will be ignored.'.format(default_src, os.path.realpath(defaults[0][0]), nvcc_version)
        else:
            nvcc_warning = None

        for (path, version, default) in paths:
            (found_some, not_found, found) = mesonlib.version_compare_many(version, version_reqs)
            if not not_found:
                if not default and nvcc_warning:
                    mlog.warning(nvcc_warning)
                return (path, version, True)

        if nvcc_warning:
            mlog.warning(nvcc_warning)
        return (None, None, False)

    def _default_path_env_var(self) -> T.Optional[str]:
        env_vars = ['CUDA_PATH'] if self._is_windows() else ['CUDA_PATH', 'CUDA_HOME', 'CUDA_ROOT']
        env_vars = [var for var in env_vars if var in os.environ]
        user_defaults = {os.environ[var] for var in env_vars}
        if len(user_defaults) > 1:
            mlog.warning('Environment variables {} point to conflicting toolkit locations ({}). Toolkit selection might produce unexpected results.'.format(', '.join(env_vars), ', '.join(user_defaults)))
        return env_vars[0] if env_vars else None

    def _cuda_paths(self) -> T.List[T.Tuple[str, bool]]:
        return ([(os.environ[self.env_var], True)] if self.env_var else []) \
            + (self._cuda_paths_win() if self._is_windows() else self._cuda_paths_nix())

    def _cuda_paths_win(self) -> T.List[T.Tuple[str, bool]]:
        env_vars = os.environ.keys()
        return [(os.environ[var], False) for var in env_vars if var.startswith('CUDA_PATH_')]

    def _cuda_paths_nix(self) -> T.List[T.Tuple[str, bool]]:
        # include /usr/local/cuda default only if no env_var was found
        pattern = '/usr/local/cuda-*' if self.env_var else '/usr/local/cuda*'
        return [(path, os.path.basename(path) == 'cuda') for path in glob.iglob(pattern)]

    toolkit_version_regex = re.compile(r'^CUDA Version\s+(.*)$')
    path_version_win_regex = re.compile(r'^v(.*)$')
    path_version_nix_regex = re.compile(r'^cuda-(.*)$')
    cudart_version_regex = re.compile(r'#define\s+CUDART_VERSION\s+([0-9]+)')

    def _cuda_toolkit_version(self, path: str) -> str:
        version = self._read_toolkit_version_txt(path)
        if version:
            return version
        version = self._read_cuda_runtime_api_version(path)
        if version:
            return version

        mlog.debug('Falling back to extracting version from path')
        path_version_regex = self.path_version_win_regex if self._is_windows() else self.path_version_nix_regex
        try:
            m = path_version_regex.match(os.path.basename(path))
            if m:
                return m.group(1)
            else:
                mlog.warning(f'Could not detect CUDA Toolkit version for {path}')
        except Exception as e:
            mlog.warning(f'Could not detect CUDA Toolkit version for {path}: {e!s}')

        return '0.0'

    def _read_cuda_runtime_api_version(self, path_str: str) -> T.Optional[str]:
        path = Path(path_str)
        for i in path.rglob('cuda_runtime_api.h'):
            raw = i.read_text(encoding='utf-8')
            m = self.cudart_version_regex.search(raw)
            if not m:
                continue
            try:
                vers_int = int(m.group(1))
            except ValueError:
                continue
            # use // for floor instead of / which produces a float
            major = vers_int // 1000
            minor = (vers_int - major * 1000) // 10
            return f'{major}.{minor}'
        return None

    def _read_toolkit_version_txt(self, path: str) -> T.Optional[str]:
        # Read 'version.txt' at the root of the CUDA Toolkit directory to determine the toolkit version
        version_file_path = os.path.join(path, 'version.txt')
        try:
            with open(version_file_path, encoding='utf-8') as version_file:
                version_str = version_file.readline() # e.g. 'CUDA Version 10.1.168'
                m = self.toolkit_version_regex.match(version_str)
                if m:
                    return self._strip_patch_version(m.group(1))
        except Exception as e:
            mlog.debug(f'Could not read CUDA Toolkit\'s version file {version_file_path}: {e!s}')

        return None

    @classmethod
    def _strip_patch_version(cls, version: str) -> str:
        return '.'.join(version.split('.')[:2])

    def _detect_arch_libdir(self) -> str:
        machine = self.env.machines[self.for_machine]
        arch = machine.cpu_family
        msg = '{} architecture is not supported in {} version of the CUDA Toolkit.'
        if machine.is_windows():
            libdirs = {'x86': 'Win32', 'x86_64': 'x64'}
            if arch not in libdirs:
                raise DependencyException(msg.format(arch, 'Windows'))
            return os.path.join('lib', libdirs[arch])
        elif machine.is_linux():
            return 'lib'
        elif machine.is_darwin():
            libdirs = {'x86_64': 'lib64'}
            if arch not in libdirs:
                raise DependencyException(msg.format(arch, 'macOS'))
            return libdirs[arch]
        else:
            raise DependencyException('CUDA Toolkit: unsupported platform.')

    def _find_requested_libraries(self) -> bool:
        all_found = True

        for module in self.requested_modules:
            # You should only ever link to libraries inside the cuda tree, nothing outside of it.
            # For instance, there is a
            #
            # - libnvidia-ml.so in stubs/ of the CUDA tree
            # - libnvidia-ml.so in /usr/lib/ that is provided by the nvidia drivers
            #
            # Users should never link to the latter, since its ABI may change.
            args = self.clib_compiler.find_library(module, self.env, [self.libdir, os.path.join(self.libdir, 'stubs')], self.libtype, ignore_system_dirs=True)

            if args is None:
                self._report_dependency_error(f'Couldn\'t find requested CUDA module \'{module}\'')
                all_found = False
            else:
                mlog.debug(f'Link args for CUDA module \'{module}\' are {args}')
                self.lib_modules[module] = args

        return all_found

    def _is_windows(self) -> bool:
        return self.env.machines[self.for_machine].is_windows()

    @T.overload
    def _report_dependency_error(self, msg: str) -> None: ...

    @T.overload
    def _report_dependency_error(self, msg: str, ret_val: TV_ResultTuple) -> TV_ResultTuple: ... # noqa: F811

    def _report_dependency_error(self, msg: str, ret_val: T.Optional[TV_ResultTuple] = None) -> T.Optional[TV_ResultTuple]: # noqa: F811
        if self.required:
            raise DependencyException(msg)

        mlog.debug(msg)
        return ret_val

    def log_details(self) -> str:
        module_str = ', '.join(self.requested_modules)
        return 'modules: ' + module_str

    def log_info(self) -> str:
        return self.cuda_path if self.cuda_path else ''

    def get_requested(self, kwargs: T.Dict[str, T.Any]) -> T.List[str]:
        candidates = mesonlib.extract_as_list(kwargs, 'modules')
        for c in candidates:
            if not isinstance(c, str):
                raise DependencyException('CUDA module argument is not a string.')
        return candidates

    def get_link_args(self, language: T.Optional[str] = None, raw: bool = False) -> T.List[str]:
        # when using nvcc to link, we should instead use the native driver options
        REWRITE_MODULES = {
            'cudart': ['-cudart', 'shared'],
            'cudart_static': ['-cudart', 'static'],
            'cudadevrt': ['-cudadevrt'],
        }

        args: T.List[str] = []
        for lib in self.requested_modules:
            link_args = self.lib_modules[lib]
            if language == 'cuda' and lib in REWRITE_MODULES:
                link_args = REWRITE_MODULES[lib]
                mlog.debug(f'Rewriting module \'{lib}\' to \'{link_args}\'')
            elif lib == 'cudart_static':
                machine = self.env.machines[self.for_machine]
                if machine.is_linux():
                    # extracted by running
                    #   nvcc -v foo.o
                    link_args += ['-lrt', '-lpthread', '-ldl']

            args += link_args

        return args

packages['cuda'] = CudaDependency
