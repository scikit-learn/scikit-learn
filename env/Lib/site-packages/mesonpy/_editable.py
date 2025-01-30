# SPDX-FileCopyrightText: 2022 The meson-python developers
#
# SPDX-License-Identifier: MIT

# This file should be standalone! It is copied during the editable hook installation.

from __future__ import annotations

import ast
import functools
import importlib.abc
import importlib.machinery
import importlib.util
import inspect
import json
import os
import pathlib
import subprocess
import sys
import typing


if typing.TYPE_CHECKING:
    from collections.abc import Sequence, Set
    from types import ModuleType
    from typing import Any, Dict, Iterator, List, Optional, Tuple, Union

    from typing_extensions import Buffer

    NodeBase = Dict[str, Union['Node', str]]
    PathStr = Union[str, os.PathLike[str]]
else:
    NodeBase = dict


if sys.version_info >= (3, 12):
    from importlib.resources.abc import Traversable, TraversableResources
elif sys.version_info >= (3, 9):
    from importlib.abc import Traversable, TraversableResources
else:
    class Traversable:
        pass
    class TraversableResources:
        pass


MARKER = 'MESONPY_EDITABLE_SKIP'
VERBOSE = 'MESONPY_EDITABLE_VERBOSE'


class MesonpyOrphan(Traversable):
    def __init__(self, name: str):
        self._name = name

    @property
    def name(self) -> str:
        return self._name

    def is_dir(self) -> bool:
        return False

    def is_file(self) -> bool:
        return False

    def iterdir(self) -> Iterator[Traversable]:
        raise FileNotFoundError()

    def open(self, *args, **kwargs):  # type: ignore
        raise FileNotFoundError()

    def joinpath(self, *descendants: PathStr) -> Traversable:
        if not descendants:
            return self
        name = os.fspath(descendants[-1]).split('/')[-1]
        return MesonpyOrphan(name)

    def __truediv__(self, child: PathStr) -> Traversable:
        return self.joinpath(child)

    def read_bytes(self) -> bytes:
        raise FileNotFoundError()

    def read_text(self, encoding: Optional[str] = None) -> str:
        raise FileNotFoundError()


class MesonpyTraversable(Traversable):
    def __init__(self, name: str, tree: Node):
        self._name = name
        self._tree = tree

    @property
    def name(self) -> str:
        return self._name

    def is_dir(self) -> bool:
        return True

    def is_file(self) -> bool:
        return False

    def iterdir(self) -> Iterator[Traversable]:
        for name, node in self._tree.items():
            yield MesonpyTraversable(name, node) if isinstance(node, dict) else pathlib.Path(node)  # type: ignore

    def open(self, *args, **kwargs):  # type: ignore
        raise IsADirectoryError()

    @staticmethod
    def _flatten(names: Tuple[PathStr, ...]) -> Iterator[str]:
        for name in names:
            yield from os.fspath(name).split('/')

    def joinpath(self, *descendants: PathStr) -> Traversable:
        if not descendants:
            return self
        names = self._flatten(descendants)
        name = next(names)
        node = self._tree.get(name)
        if isinstance(node, dict):
            return MesonpyTraversable(name, node).joinpath(*names)
        if isinstance(node, str):
            return pathlib.Path(node).joinpath(*names)
        return MesonpyOrphan(name).joinpath(*names)

    def __truediv__(self, child: PathStr) -> Traversable:
        return self.joinpath(child)

    def read_bytes(self) -> bytes:
        raise IsADirectoryError()

    def read_text(self, encoding: Optional[str] = None) -> str:
        raise IsADirectoryError()


class MesonpyReader(TraversableResources):
    def __init__(self, name: str, tree: Node):
        self._name = name
        self._tree = tree

    def files(self) -> Traversable:
        return MesonpyTraversable(self._name, self._tree)


class ExtensionFileLoader(importlib.machinery.ExtensionFileLoader):
    def __init__(self, name: str, path: str, tree: Node):
        super().__init__(name, path)
        self._tree = tree

    def get_resource_reader(self, name: str) -> TraversableResources:
        return MesonpyReader(name, self._tree)


class SourceFileLoader(importlib.machinery.SourceFileLoader):
    def __init__(self, name: str, path: str, tree: Node):
        super().__init__(name, path)
        self._tree = tree

    def set_data(self, path: Union[bytes, str], data: Buffer, *, _mode: int = ...) -> None:
        # disable saving bytecode
        pass

    def get_resource_reader(self, name: str) -> TraversableResources:
        return MesonpyReader(name, self._tree)


class SourcelessFileLoader(importlib.machinery.SourcelessFileLoader):
    def __init__(self, name: str, path: str, tree: Node):
        super().__init__(name, path)
        self._tree = tree

    def get_resource_reader(self, name: str) -> TraversableResources:
        return MesonpyReader(name, self._tree)


LOADERS = \
    [(ExtensionFileLoader, s) for s in importlib.machinery.EXTENSION_SUFFIXES] + \
    [(SourceFileLoader, s) for s in importlib.machinery.SOURCE_SUFFIXES] + \
    [(SourcelessFileLoader, s) for s in importlib.machinery.BYTECODE_SUFFIXES]


def build_module_spec(cls: type, name: str, path: str, tree: Optional[Node]) -> importlib.machinery.ModuleSpec:
    loader = cls(name, path, tree)
    spec = importlib.machinery.ModuleSpec(name, loader, origin=path)
    spec.has_location = True
    if loader.is_package(name):
        spec.submodule_search_locations = [os.path.join(__file__, name)]
    return spec


class Node(NodeBase):
    """Tree structure to store a virtual filesystem view."""

    def __missing__(self, key: str) -> Node:
        value = self[key] = Node()
        return value

    def __setitem__(self, key: Union[str, Tuple[str, ...]], value: Union[Node, str]) -> None:
        node = self
        if isinstance(key, tuple):
            for k in key[:-1]:
                node = typing.cast(Node, node[k])
            key = key[-1]
        dict.__setitem__(node, key, value)

    def __getitem__(self, key: Union[str, Tuple[str, ...]]) -> Union[Node, str]:
        node = self
        if isinstance(key, tuple):
            for k in key[:-1]:
                node = typing.cast(Node, node[k])
            key = key[-1]
        return dict.__getitem__(node, key)

    def get(self, key: Union[str, Tuple[str, ...]]) -> Optional[Union[Node, str]]:  # type: ignore[override]
        node = self
        if isinstance(key, tuple):
            for k in key[:-1]:
                v = dict.get(node, k)
                if v is None:
                    return None
                node = typing.cast(Node, v)
            key = key[-1]
        return dict.get(node, key)


def walk(src: str, exclude_files: Set[str], exclude_dirs: Set[str]) -> Iterator[str]:
    for root, dirnames, filenames in os.walk(src):
        for name in dirnames.copy():
            dirsrc = os.path.join(root, name)
            relpath = os.path.relpath(dirsrc, src)
            if relpath in exclude_dirs:
                dirnames.remove(name)
            # sort to process directories determninistically
            dirnames.sort()
        for name in sorted(filenames):
            filesrc = os.path.join(root, name)
            relpath = os.path.relpath(filesrc, src)
            if relpath in exclude_files:
                continue
            yield relpath


def collect(install_plan: Dict[str, Dict[str, Any]]) -> Node:
    tree = Node()
    for key, data in install_plan.items():
        for src, target in data.items():
            path = pathlib.Path(target['destination'])
            if path.parts[0] in {'{py_platlib}', '{py_purelib}'}:
                if key == 'install_subdirs' or key == 'targets' and os.path.isdir(src):
                    exclude_files = {os.path.normpath(x) for x in target.get('exclude_files', [])}
                    exclude_dirs = {os.path.normpath(x) for x in target.get('exclude_dirs', [])}
                    for entry in walk(src, exclude_files, exclude_dirs):
                        tree[(*path.parts[1:], *entry.split(os.sep))] = os.path.join(src, entry)
                else:
                    tree[path.parts[1:]] = src
    return tree


def find_spec(fullname: str, tree: Node) -> Optional[importlib.machinery.ModuleSpec]:
    namespace = False
    parts = fullname.split('.')

    # look for a package
    package = tree.get(tuple(parts))
    if isinstance(package, Node):
        for loader, suffix in LOADERS:
            src = package.get('__init__' + suffix)
            if isinstance(src, str):
                return build_module_spec(loader, fullname, src, package)
        else:
            namespace = True

    # look for a module
    for loader, suffix in LOADERS:
        src = tree.get((*parts[:-1], parts[-1] + suffix))
        if isinstance(src, str):
            return build_module_spec(loader, fullname, src, None)

    # namespace
    if namespace:
        spec = importlib.machinery.ModuleSpec(fullname, None, is_package=True)
        assert isinstance(spec.submodule_search_locations, list)  # make mypy happy
        spec.submodule_search_locations.append(os.path.join(__file__, fullname))
        return spec

    return None


class MesonpyMetaFinder(importlib.abc.MetaPathFinder):
    def __init__(self, package: str, names: Set[str], path: str, cmd: List[str], verbose: bool = False):
        self._name = package
        self._top_level_modules = names
        self._build_path = path
        self._build_cmd = cmd
        self._verbose = verbose
        self._loaders: List[Tuple[type, str]] = []

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self._name!r}, {self._build_path!r})'

    def find_spec(
            self,
            fullname: str,
            path: Optional[Sequence[Union[bytes, str]]] = None,
            target: Optional[ModuleType] = None
    ) -> Optional[importlib.machinery.ModuleSpec]:
        if fullname.split('.', 1)[0] not in self._top_level_modules:
            return None
        if self._build_path in os.environ.get(MARKER, '').split(os.pathsep):
            return None
        tree = self._rebuild()
        return find_spec(fullname, tree)

    def _work_to_do(self, env: dict[str, str]) -> bool:
        if sys.platform == 'win32':
            # On Windows the build command is 'meson compile' eventually with a --ninja-args= option.
            if self._build_cmd[-1].startswith('--ninja-args='):
                ninja_args = ast.literal_eval(self._build_cmd[-1].split('=', 1)[1]) + ['-n']
                dry_run_build_cmd = self._build_cmd[:-1] + [f'--ninja-args={ninja_args!r}']
            else:
                dry_run_build_cmd = self._build_cmd + ['--ninja-args=-n']
        else:
            dry_run_build_cmd = self._build_cmd + ['-n']
        # Check adapted from
        # https://github.com/mesonbuild/meson/blob/a35d4d368a21f4b70afa3195da4d6292a649cb4c/mesonbuild/mtest.py#L1635-L1636
        p = subprocess.run(dry_run_build_cmd, cwd=self._build_path, env=env, capture_output=True)
        return b'ninja: no work to do.' not in p.stdout and b'samu: nothing to do' not in p.stdout

    @functools.lru_cache(maxsize=1)
    def _rebuild(self) -> Node:
        try:
            # Skip editable wheel lookup during rebuild: during the build
            # the module we are rebuilding might be imported causing a
            # rebuild loop.
            env = os.environ.copy()
            env[MARKER] = os.pathsep.join((env.get(MARKER, ''), self._build_path))

            if self._verbose or bool(env.get(VERBOSE, '')):
                # We want to show some output only if there is some work to do.
                if self._work_to_do(env):
                    build_command = ' '.join(self._build_cmd)
                    print(f'meson-python: building {self._name}: {build_command}', flush=True)
                    subprocess.run(self._build_cmd, cwd=self._build_path, env=env, check=True)
            else:
                subprocess.run(self._build_cmd, cwd=self._build_path, env=env, stdout=subprocess.DEVNULL, check=True)
        except subprocess.CalledProcessError as exc:
            raise ImportError(f're-building the {self._name} meson-python editable wheel package failed') from exc

        install_plan_path = os.path.join(self._build_path, 'meson-info', 'intro-install_plan.json')
        with open(install_plan_path, 'r', encoding='utf8') as f:
            install_plan = json.load(f)
        return collect(install_plan)

    def _path_hook(self, path: str) -> MesonpyPathFinder:
        if os.altsep:
            path.replace(os.altsep, os.sep)
        path, _, key = path.rpartition(os.sep)
        if path == __file__:
            tree = self._rebuild()
            node = tree.get(tuple(key.split('.')))
            if isinstance(node, Node):
                return MesonpyPathFinder(node)
        raise ImportError


class MesonpyPathFinder(importlib.abc.PathEntryFinder):
    def __init__(self, tree: Node):
        self._tree = tree

    def find_spec(self, fullname: str, target: Optional[ModuleType] = None) -> Optional[importlib.machinery.ModuleSpec]:
        return find_spec(fullname, self._tree)

    def iter_modules(self, prefix: str) -> Iterator[Tuple[str, bool]]:
        yielded = set()
        for name, node in self._tree.items():
            modname = inspect.getmodulename(name)
            if modname == '__init__' or modname in yielded:
                continue
            if isinstance(node, Node):
                modname = name
                for _, suffix in LOADERS:
                    src = node.get('__init__' + suffix)
                    if isinstance(src, str):
                        yielded.add(modname)
                        yield prefix + modname, True
            elif modname and '.' not in modname:
                yielded.add(modname)
                yield prefix + modname, False


def install(package: str, names: Set[str], path: str, cmd: List[str], verbose: bool) -> None:
    finder = MesonpyMetaFinder(package, names, path, cmd, verbose)
    sys.meta_path.insert(0, finder)
    sys.path_hooks.insert(0, finder._path_hook)
