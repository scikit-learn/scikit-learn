# SPDX-License-Identifier: Apache-2.0
# Copyright 2019 The Meson development team
# Copyright © 2023-2025 Intel Corporation

from __future__ import annotations
from pathlib import Path

import abc
import os
import re
import typing as T

if T.TYPE_CHECKING:
    from ..minit import Arguments


class SampleImpl(metaclass=abc.ABCMeta):

    def __init__(self, args: Arguments):
        self.name = args.name
        self.executable_name = args.executable if args.executable else None
        self.version = args.version
        self.lowercase_token = re.sub(r'[^a-z0-9]', '_', self.name.lower())
        self.uppercase_token = self.lowercase_token.upper()
        self.capitalized_token = self.lowercase_token.capitalize()
        self.meson_version = '1.0.0'
        self.force = args.force
        self.dependencies = args.deps.split(',') if args.deps else []
        self.sources: list[Path] = []

    @abc.abstractmethod
    def create_executable(self) -> None:
        pass

    @abc.abstractmethod
    def create_library(self) -> None:
        pass

    @property
    @abc.abstractmethod
    def exe_template(self) -> str:
        pass

    @property
    @abc.abstractmethod
    def exe_meson_template(self) -> str:
        pass

    @property
    @abc.abstractmethod
    def lib_template(self) -> str:
        pass

    @property
    @abc.abstractmethod
    def lib_test_template(self) -> T.Optional[str]:
        pass

    @property
    @abc.abstractmethod
    def lib_meson_template(self) -> str:
        pass

    @property
    @abc.abstractmethod
    def source_ext(self) -> str:
        pass

    def _format_dependencies(self) -> str:
        return ''.join(f"\n  dependency('{d}')," for d in self.dependencies)

    def _format_sources(self) -> str:
        return ''.join(f"\n  '{x}'," for x in self.sources)

    def _detect_sources(self, srcfiles: T.List[Path], transform: T.Callable[[str], str]) -> T.Tuple[str, T.List[Path]]:
        # Try a source based on the executable name, fallback to one based
        # on the project name if none is found.
        expected_name = f'{transform(self.executable_name)}.{self.source_ext}'
        if any(str(x) == expected_name for x in srcfiles):
            return expected_name, srcfiles
        expected_name = f'{transform(self.lowercase_token)}.{self.source_ext}'
        if any(str(x) == expected_name for x in srcfiles):
            return expected_name, srcfiles
        return expected_name, [Path(expected_name)]

class ClassImpl(SampleImpl):

    """For Class based languages, like Java and C#"""

    def __init__(self, args: Arguments):
        super().__init__(args)
        self.source_name, self.sources = self._detect_sources(args.srcfiles, str.capitalize)

    def create_executable(self) -> None:
        if not os.path.exists(self.source_name):
            with open(self.source_name, 'w', encoding='utf-8') as f:
                f.write(self.exe_template.format(project_name=self.name,
                                                 class_name=self.capitalized_token))
        if self.force or not os.path.exists('meson.build'):
            with open('meson.build', 'w', encoding='utf-8') as f:
                f.write(self.exe_meson_template.format(project_name=self.name,
                                                       exe_name=self.executable_name,
                                                       source_name=self.source_name,
                                                       version=self.version,
                                                       meson_version=self.meson_version,
                                                       dependencies=self._format_dependencies(),
                                                       source_files=self._format_sources()))

    def create_library(self) -> None:
        test_name = f'{self.capitalized_token}_test.{self.source_ext}'
        kwargs = {
            'utoken': self.uppercase_token,
            'ltoken': self.lowercase_token,
            'class_test': f'{self.capitalized_token}_test',
            'class_name': self.capitalized_token,
            'source_file': self.source_name,
            'test_source_file': test_name,
            'test_exe_name': f'{self.lowercase_token}_test',
            'project_name': self.name,
            'lib_name': self.lowercase_token,
            'test_name': self.lowercase_token,
            'version': self.version,
            'meson_version': self.meson_version,
            'dependencies': self._format_dependencies(),
            'source_files': self._format_sources(),
        }
        if not os.path.exists(self.source_name):
            with open(self.source_name, 'w', encoding='utf-8') as f:
                f.write(self.lib_template.format(**kwargs))
        if self.lib_test_template and not os.path.exists(test_name):
            with open(test_name, 'w', encoding='utf-8') as f:
                f.write(self.lib_test_template.format(**kwargs))
        if self.force or not os.path.exists('meson.build'):
            with open('meson.build', 'w', encoding='utf-8') as f:
                f.write(self.lib_meson_template.format(**kwargs))


class FileImpl(SampleImpl):

    """File based languages without headers"""

    def __init__(self, args: Arguments):
        super().__init__(args)
        self.source_name, self.sources = self._detect_sources(args.srcfiles, str.lower)

    def create_executable(self) -> None:
        if not os.path.exists(self.source_name):
            with open(self.source_name, 'w', encoding='utf-8') as f:
                f.write(self.exe_template.format(project_name=self.name))
        if self.force or not os.path.exists('meson.build'):
            with open('meson.build', 'w', encoding='utf-8') as f:
                f.write(self.exe_meson_template.format(project_name=self.name,
                                                       exe_name=self.executable_name,
                                                       source_name=self.source_name,
                                                       version=self.version,
                                                       meson_version=self.meson_version,
                                                       dependencies=self._format_dependencies(),
                                                       source_files=self._format_sources()))

    def lib_kwargs(self) -> T.Dict[str, str]:
        """Get Language specific keyword arguments

        :return: A dictionary of key: values to fill in the templates
        """
        return {
            'utoken': self.uppercase_token,
            'ltoken': self.lowercase_token,
            'header_dir': self.lowercase_token,
            'class_name': self.capitalized_token,
            'function_name': f'{self.lowercase_token[0:3]}_func',
            'namespace': self.lowercase_token,
            'source_file': f'{self.lowercase_token}.{self.source_ext}',
            'test_source_file': f'{self.lowercase_token}_test.{self.source_ext}',
            'test_exe_name': f'{self.lowercase_token}_test',
            'project_name': self.name,
            'lib_name': self.lowercase_token,
            'test_name': self.lowercase_token,
            'version': self.version,
            'meson_version': self.meson_version,
            'dependencies': self._format_dependencies(),
            'source_files': self._format_sources(),
        }

    def create_library(self) -> None:
        test_name = f'{self.lowercase_token}_test.{self.source_ext}'
        kwargs = self.lib_kwargs()
        if not os.path.exists(self.source_name):
            with open(self.source_name, 'w', encoding='utf-8') as f:
                f.write(self.lib_template.format(**kwargs))
        if self.lib_test_template and not os.path.exists(test_name):
            with open(test_name, 'w', encoding='utf-8') as f:
                f.write(self.lib_test_template.format(**kwargs))
        if self.force or not os.path.exists('meson.build'):
            with open('meson.build', 'w', encoding='utf-8') as f:
                f.write(self.lib_meson_template.format(**kwargs))


class FileHeaderImpl(FileImpl):

    @property
    @abc.abstractmethod
    def header_ext(self) -> str:
        pass

    @property
    @abc.abstractmethod
    def lib_header_template(self) -> str:
        pass

    def lib_kwargs(self) -> T.Dict[str, str]:
        kwargs = super().lib_kwargs()
        kwargs['header_file'] = f'{self.lowercase_token}.{self.header_ext}'
        return kwargs

    def create_library(self) -> None:
        super().create_library()
        kwargs = self.lib_kwargs()
        if not os.path.exists(kwargs['header_file']):
            with open(kwargs['header_file'], 'w', encoding='utf-8') as f:
                f.write(self.lib_header_template.format_map(kwargs))
