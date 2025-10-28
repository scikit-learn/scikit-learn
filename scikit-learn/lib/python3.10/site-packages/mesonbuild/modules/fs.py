# SPDX-License-Identifier: Apache-2.0
# Copyright 2019 The Meson development team

from __future__ import annotations
from ntpath import sep as ntsep
from pathlib import Path
from posixpath import sep as posixsep
import hashlib
import os
import typing as T

from . import ExtensionModule, ModuleReturnValue, ModuleInfo
from .. import mlog
from ..build import BuildTarget, CustomTarget, CustomTargetIndex, InvalidArguments
from ..interpreter.type_checking import INSTALL_KW, INSTALL_MODE_KW, INSTALL_TAG_KW, NoneType
from ..interpreterbase import FeatureNew, KwargInfo, typed_kwargs, typed_pos_args, noKwargs
from ..mesonlib import File, MesonException, has_path_sep, is_windows, path_is_in_root, relpath

if T.TYPE_CHECKING:
    from . import ModuleState
    from ..build import BuildTargetTypes
    from ..interpreter import Interpreter
    from ..interpreterbase import TYPE_kwargs
    from ..mesonlib import FileOrString, FileMode

    from typing_extensions import TypedDict

    class ReadKwArgs(TypedDict):
        """Keyword Arguments for fs.read."""

        encoding: str

    class CopyKw(TypedDict):

        """Kwargs for fs.copy"""

        install: bool
        install_dir: T.Optional[str]
        install_mode: FileMode
        install_tag: T.Optional[str]


class FSModule(ExtensionModule):

    INFO = ModuleInfo('fs', '0.53.0')

    def __init__(self, interpreter: Interpreter) -> None:
        super().__init__(interpreter)
        self.methods.update({
            'as_posix': self.as_posix,
            'copyfile': self.copyfile,
            'exists': self.exists,
            'expanduser': self.expanduser,
            'hash': self.hash,
            'is_absolute': self.is_absolute,
            'is_dir': self.is_dir,
            'is_file': self.is_file,
            'is_samepath': self.is_samepath,
            'is_symlink': self.is_symlink,
            'name': self.name,
            'parent': self.parent,
            'read': self.read,
            'relative_to': self.relative_to,
            'replace_suffix': self.replace_suffix,
            'size': self.size,
            'stem': self.stem,
            'suffix': self.suffix,
        })

    def _absolute_dir(self, state: ModuleState, arg: FileOrString) -> str:
        """
        make an absolute path from a relative path, WITHOUT resolving symlinks
        """
        if isinstance(arg, File):
            return arg.absolute_path(state.source_root, state.environment.get_build_dir())
        return os.path.join(state.source_root, state.subdir, os.path.expanduser(arg))

    @staticmethod
    def _obj_to_pathstr(feature_new_prefix: str, obj: T.Union[FileOrString, BuildTargetTypes], state: ModuleState) -> str:
        if isinstance(obj, str):
            return obj

        if isinstance(obj, File):
            FeatureNew(f'{feature_new_prefix} with file', '0.59.0').use(state.subproject, location=state.current_node)
            return str(obj)

        FeatureNew(f'{feature_new_prefix} with build_tgt, custom_tgt, and custom_idx', '1.4.0').use(state.subproject, location=state.current_node)
        return state.backend.get_target_filename(obj)

    def _resolve_dir(self, state: ModuleState, arg: FileOrString) -> str:
        """
        resolves symlinks and makes absolute a directory relative to calling meson.build,
        if not already absolute
        """
        path = self._absolute_dir(state, arg)
        try:
            # accommodate unresolvable paths e.g. symlink loops
            path = os.path.realpath(path)
        except Exception:
            # return the best we could do
            pass
        return path

    @noKwargs
    @FeatureNew('fs.expanduser', '0.54.0')
    @typed_pos_args('fs.expanduser', str)
    def expanduser(self, state: ModuleState, args: T.Tuple[str], kwargs: T.Dict[str, T.Any]) -> str:
        return os.path.expanduser(args[0])

    @noKwargs
    @FeatureNew('fs.is_absolute', '0.54.0')
    @typed_pos_args('fs.is_absolute', (str, File))
    def is_absolute(self, state: ModuleState, args: T.Tuple[FileOrString], kwargs: T.Dict[str, T.Any]) -> bool:
        path = args[0]
        if isinstance(path, File):
            FeatureNew('fs.is_absolute with file', '0.59.0').use(state.subproject, location=state.current_node)
            path = str(path)
        if is_windows():
            # os.path.isabs was broken for Windows before Python 3.13, so we implement it ourselves
            path = path[:3].replace(posixsep, ntsep)
            return path.startswith(ntsep * 2) or path.startswith(':' + ntsep, 1)
        return path.startswith(posixsep)

    @noKwargs
    @FeatureNew('fs.as_posix', '0.54.0')
    @typed_pos_args('fs.as_posix', str)
    def as_posix(self, state: ModuleState, args: T.Tuple[str], kwargs: T.Dict[str, T.Any]) -> str:
        r"""
        this function assumes you are passing a Windows path, even if on a Unix-like system
        and so ALL '\' are turned to '/', even if you meant to escape a character
        """
        return args[0].replace(ntsep, posixsep)

    @noKwargs
    @typed_pos_args('fs.exists', str)
    def exists(self, state: ModuleState, args: T.Tuple[str], kwargs: T.Dict[str, T.Any]) -> bool:
        return os.path.exists(self._resolve_dir(state, args[0]))

    @noKwargs
    @typed_pos_args('fs.is_symlink', (str, File))
    def is_symlink(self, state: ModuleState, args: T.Tuple[FileOrString], kwargs: T.Dict[str, T.Any]) -> bool:
        if isinstance(args[0], File):
            FeatureNew('fs.is_symlink with file', '0.59.0').use(state.subproject, location=state.current_node)
        return os.path.islink(self._absolute_dir(state, args[0]))

    @noKwargs
    @typed_pos_args('fs.is_file', str)
    def is_file(self, state: ModuleState, args: T.Tuple[str], kwargs: T.Dict[str, T.Any]) -> bool:
        return os.path.isfile(self._resolve_dir(state, args[0]))

    @noKwargs
    @typed_pos_args('fs.is_dir', str)
    def is_dir(self, state: ModuleState, args: T.Tuple[str], kwargs: T.Dict[str, T.Any]) -> bool:
        return os.path.isdir(self._resolve_dir(state, args[0]))

    @noKwargs
    @typed_pos_args('fs.hash', (str, File), str)
    def hash(self, state: ModuleState, args: T.Tuple[FileOrString, str], kwargs: T.Dict[str, T.Any]) -> str:
        if isinstance(args[0], File):
            FeatureNew('fs.hash with file', '0.59.0').use(state.subproject, location=state.current_node)
        file = self._resolve_dir(state, args[0])
        if not os.path.isfile(file):
            raise MesonException(f'{file} is not a file and therefore cannot be hashed')
        try:
            h = hashlib.new(args[1])
        except ValueError:
            raise MesonException('hash algorithm {} is not available'.format(args[1]))
        mlog.debug('computing {} sum of {} size {} bytes'.format(args[1], file, os.stat(file).st_size))
        with open(file, mode='rb', buffering=0) as f:
            h.update(f.read())
        return h.hexdigest()

    @noKwargs
    @typed_pos_args('fs.size', (str, File))
    def size(self, state: ModuleState, args: T.Tuple[FileOrString], kwargs: T.Dict[str, T.Any]) -> int:
        if isinstance(args[0], File):
            FeatureNew('fs.size with file', '0.59.0').use(state.subproject, location=state.current_node)
        file = self._resolve_dir(state, args[0])
        if not os.path.isfile(file):
            raise MesonException(f'{file} is not a file and therefore cannot be sized')
        try:
            return os.stat(file).st_size
        except ValueError:
            raise MesonException('{} size could not be determined'.format(args[0]))

    @noKwargs
    @typed_pos_args('fs.is_samepath', (str, File), (str, File))
    def is_samepath(self, state: ModuleState, args: T.Tuple[FileOrString, FileOrString], kwargs: T.Dict[str, T.Any]) -> bool:
        if isinstance(args[0], File) or isinstance(args[1], File):
            FeatureNew('fs.is_samepath with file', '0.59.0').use(state.subproject, location=state.current_node)
        file1 = self._resolve_dir(state, args[0])
        file2 = self._resolve_dir(state, args[1])
        if not os.path.exists(file1):
            return False
        if not os.path.exists(file2):
            return False
        try:
            return os.path.samefile(file1, file2)
        except OSError:
            return False

    @noKwargs
    @typed_pos_args('fs.replace_suffix', (str, File, CustomTarget, CustomTargetIndex, BuildTarget), str)
    def replace_suffix(self, state: ModuleState, args: T.Tuple[T.Union[FileOrString, BuildTargetTypes], str], kwargs: T.Dict[str, T.Any]) -> str:
        if args[1] and not args[1].startswith('.'):
            raise ValueError(f"Invalid suffix {args[1]!r}")
        path = self._obj_to_pathstr('fs.replace_suffix', args[0], state)
        return os.path.splitext(path)[0] + args[1]

    @noKwargs
    @typed_pos_args('fs.parent', (str, File, CustomTarget, CustomTargetIndex, BuildTarget))
    def parent(self, state: ModuleState, args: T.Tuple[T.Union[FileOrString, BuildTargetTypes]], kwargs: T.Dict[str, T.Any]) -> str:
        path = self._obj_to_pathstr('fs.parent', args[0], state)
        return os.path.split(path)[0] or '.'

    @noKwargs
    @typed_pos_args('fs.name', (str, File, CustomTarget, CustomTargetIndex, BuildTarget))
    def name(self, state: ModuleState, args: T.Tuple[T.Union[FileOrString, BuildTargetTypes]], kwargs: T.Dict[str, T.Any]) -> str:
        path = self._obj_to_pathstr('fs.name', args[0], state)
        return os.path.basename(path)

    @noKwargs
    @typed_pos_args('fs.stem', (str, File, CustomTarget, CustomTargetIndex, BuildTarget))
    @FeatureNew('fs.stem', '0.54.0')
    def stem(self, state: ModuleState, args: T.Tuple[T.Union[FileOrString, BuildTargetTypes]], kwargs: T.Dict[str, T.Any]) -> str:
        path = self._obj_to_pathstr('fs.name', args[0], state)
        return os.path.splitext(os.path.basename(path))[0]

    @noKwargs
    @typed_pos_args('fs.suffix', (str, File, CustomTarget, CustomTargetIndex, BuildTarget))
    @FeatureNew('fs.suffix', '1.9.0')
    def suffix(self, state: ModuleState, args: T.Tuple[T.Union[FileOrString, BuildTargetTypes]], kwargs: T.Dict[str, T.Any]) -> str:
        path = self._obj_to_pathstr('fs.suffix', args[0], state)
        return os.path.splitext(path)[1]

    @FeatureNew('fs.read', '0.57.0')
    @typed_pos_args('fs.read', (str, File))
    @typed_kwargs('fs.read', KwargInfo('encoding', str, default='utf-8'))
    def read(self, state: ModuleState, args: T.Tuple[FileOrString], kwargs: ReadKwArgs) -> str:
        """Read a file from the source tree and return its value as a decoded
        string.

        If the encoding is not specified, the file is assumed to be utf-8
        encoded. Paths must be relative by default (to prevent accidents) and
        are forbidden to be read from the build directory (to prevent build
        loops)
        """
        path = args[0]
        encoding = kwargs['encoding']
        src_dir = state.environment.source_dir
        sub_dir = state.subdir
        build_dir = state.environment.get_build_dir()

        if isinstance(path, File):
            if path.is_built:
                raise MesonException(
                    'fs.read does not accept built files() objects')
            path = os.path.join(src_dir, path.relative_name())
        else:
            if sub_dir:
                src_dir = os.path.join(src_dir, sub_dir)
            path = os.path.join(src_dir, path)

        path = os.path.abspath(path)
        if path_is_in_root(Path(path), Path(build_dir), resolve=True):
            raise MesonException('path must not be in the build tree')
        try:
            with open(path, encoding=encoding) as f:
                data = f.read()
        except FileNotFoundError:
            raise MesonException(f'File {args[0]} does not exist.')
        except UnicodeDecodeError:
            raise MesonException(f'decoding failed for {args[0]}')
        # Reconfigure when this file changes as it can contain data used by any
        # part of the build configuration (e.g. `project(..., version:
        # fs.read_file('VERSION')` or `configure_file(...)`
        self.interpreter.add_build_def_file(path)
        return data

    @FeatureNew('fs.copyfile', '0.64.0')
    @typed_pos_args('fs.copyfile', (File, str), optargs=[str])
    @typed_kwargs(
        'fs.copyfile',
        INSTALL_KW,
        INSTALL_MODE_KW,
        INSTALL_TAG_KW,
        KwargInfo('install_dir', (str, NoneType)),
    )
    def copyfile(self, state: ModuleState, args: T.Tuple[FileOrString, T.Optional[str]],
                 kwargs: CopyKw) -> ModuleReturnValue:
        """Copy a file into the build directory at build time."""
        if kwargs['install'] and not kwargs['install_dir']:
            raise InvalidArguments('"install_dir" must be specified when "install" is true')

        src = self.interpreter.source_strings_to_files([args[0]])[0]

        # The input is allowed to have path separators, but the output may not,
        # so use the basename for the default case
        dest = args[1] if args[1] else os.path.basename(src.fname)
        if has_path_sep(dest):
            raise InvalidArguments('Destination path may not have path separators')

        ct = CustomTarget(
            dest,
            state.subdir,
            state.subproject,
            state.environment,
            state.environment.get_build_command() + ['--internal', 'copy', '@INPUT@', '@OUTPUT@'],
            [src],
            [dest],
            build_by_default=True,
            install=kwargs['install'],
            install_dir=[kwargs['install_dir']],
            install_mode=kwargs['install_mode'],
            install_tag=[kwargs['install_tag']],
            backend=state.backend,
            description='Copying file {}',
        )

        return ModuleReturnValue(ct, [ct])

    @FeatureNew('fs.relative_to', '1.3.0')
    @typed_pos_args('fs.relative_to', (str, File, CustomTarget, CustomTargetIndex, BuildTarget), (str, File, CustomTarget, CustomTargetIndex, BuildTarget))
    @noKwargs
    def relative_to(self, state: ModuleState, args: T.Tuple[T.Union[FileOrString, BuildTargetTypes], T.Union[FileOrString, BuildTargetTypes]], kwargs: TYPE_kwargs) -> str:
        def to_path(arg: T.Union[FileOrString, CustomTarget, CustomTargetIndex, BuildTarget]) -> str:
            if isinstance(arg, File):
                return arg.absolute_path(state.environment.source_dir, state.environment.build_dir)
            elif isinstance(arg, (CustomTarget, CustomTargetIndex, BuildTarget)):
                return state.backend.get_target_filename_abs(arg)
            else:
                return os.path.join(state.environment.source_dir, state.subdir, arg)

        t = to_path(args[0])
        f = to_path(args[1])

        return relpath(t, f)


def initialize(*args: T.Any, **kwargs: T.Any) -> FSModule:
    return FSModule(*args, **kwargs)
