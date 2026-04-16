# SPDX-License-Identifier: Apache-2.0
# Copyright 2013-2024 Contributors to the The Meson project
# Copyright © 2024 Intel Corporation

from __future__ import annotations
import typing as T
import configparser
import os

from . import mparser

from .cmdline import CmdLineFileParser
from .mesonlib import MesonException

if T.TYPE_CHECKING:
    from .options import ElementaryOptionValues


HOMEDIR = os.path.expanduser('~')


class MachineFileParser():
    def __init__(self, filenames: T.List[str], sourcedir: str) -> None:
        self.parser = CmdLineFileParser()
        self.constants: T.Dict[str, ElementaryOptionValues] = {'True': True, 'False': False, '~': HOMEDIR}
        self.sections: T.Dict[str, T.Dict[str, ElementaryOptionValues]] = {}

        for fname in filenames:
            try:
                with open(fname, encoding='utf-8') as f:
                    content = f.read()
            except UnicodeDecodeError as e:
                raise MesonException(f'Malformed machine file {fname!r} failed to parse as unicode: {e}')

            content = content.replace('@GLOBAL_SOURCE_ROOT@', sourcedir)
            content = content.replace('@DIRNAME@', os.path.dirname(fname))
            try:
                self.parser.read_string(content, fname)
            except configparser.Error as e:
                raise MesonException(f'Malformed machine file: {e}')

        # Parse [constants] first so they can be used in other sections
        if self.parser.has_section('constants'):
            self.constants.update(self._parse_section('constants'))

        for s in self.parser.sections():
            if s == 'constants':
                continue
            self.sections[s] = self._parse_section(s)

    def _parse_section(self, s: str) -> T.Dict[str, ElementaryOptionValues]:
        self.scope = self.constants.copy()
        section: T.Dict[str, ElementaryOptionValues] = {}
        for entry, value in self.parser.items(s):
            if ' ' in entry or '\t' in entry or "'" in entry or '"' in entry:
                raise MesonException(f'Malformed variable name {entry!r} in machine file.')
            # Windows paths...
            value = value.replace('\\', '\\\\')
            try:
                ast = mparser.Parser(value, 'machinefile', machinefile=True).parse()
                if not ast.lines:
                    raise MesonException('value cannot be empty')
                res = self._evaluate_statement(ast.lines[0])
            except MesonException as e:
                raise MesonException(f'Malformed value in machine file variable {entry!r}: {str(e)}.')
            except KeyError as e:
                raise MesonException(f'Undefined constant {e.args[0]!r} in machine file variable {entry!r}.')
            section[entry] = res
            self.scope[entry] = res
        return section

    def _evaluate_statement(self, node: mparser.BaseNode) -> ElementaryOptionValues:
        if isinstance(node, (mparser.StringNode)):
            return node.value
        elif isinstance(node, mparser.BooleanNode):
            return node.value
        elif isinstance(node, mparser.NumberNode):
            return node.value
        elif isinstance(node, mparser.ParenthesizedNode):
            return self._evaluate_statement(node.inner)
        elif isinstance(node, mparser.ArrayNode):
            a = [self._evaluate_statement(arg) for arg in node.args.arguments]
            assert all(isinstance(s, str) for s in a), 'for mypy'
            return T.cast('T.List[str]', a)
        elif isinstance(node, mparser.IdNode):
            return self.scope[node.value]
        elif isinstance(node, mparser.ArithmeticNode):
            l = self._evaluate_statement(node.left)
            r = self._evaluate_statement(node.right)
            if node.operation == '+':
                if isinstance(l, str) and isinstance(r, str):
                    return l + r
                if isinstance(l, list) and isinstance(r, list):
                    return l + r
            elif node.operation == '/':
                if isinstance(l, str) and isinstance(r, str):
                    return os.path.join(l, r)
        raise MesonException('Unsupported node type')

def parse_machine_files(filenames: T.List[str], sourcedir: str) -> T.Dict[str, T.Dict[str, ElementaryOptionValues]]:
    parser = MachineFileParser(filenames, sourcedir)
    return parser.sections


class MachineFileStore:
    def __init__(self, native_files: T.Optional[T.List[str]], cross_files: T.Optional[T.List[str]], source_dir: str):
        self.native = parse_machine_files(native_files if native_files is not None else [], source_dir)
        self.cross = parse_machine_files(cross_files if cross_files is not None else [], source_dir)
