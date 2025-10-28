# SPDX-License-Identifier: Apache-2.0
# Copyright 2020 The Meson development team
# Copyright © 2023-2024 Intel Corporation

from __future__ import annotations

import collections
import json
import os
import pathlib
import pickle
import re
import typing as T

if T.TYPE_CHECKING:
    from typing_extensions import Literal, TypedDict, NotRequired
    from ..backend.ninjabackend import TargetDependencyScannerInfo

    Require = TypedDict(
        'Require',
        {
            'logical-name': str,
            'compiled-module-path': NotRequired[str],
            'source-path': NotRequired[str],
            'unique-on-source-path': NotRequired[bool],
            'lookup-method': NotRequired[Literal['by-name', 'include-angle', 'include-quote']]
        },
    )

    Provide = TypedDict(
        'Provide',
        {
            'logical-name': str,
            'compiled-module-path': NotRequired[str],
            'source-path': NotRequired[str],
            'unique-on-source-path': NotRequired[bool],
            'is-interface': NotRequired[bool],
        },
    )

    Rule = TypedDict(
        'Rule',
        {
            'primary-output': NotRequired[str],
            'outputs': NotRequired[T.List[str]],
            'provides': NotRequired[T.List[Provide]],
            'requires': NotRequired[T.List[Require]],
        }
    )

    class Description(TypedDict):

        version: int
        revision: int
        rules: T.List[Rule]


CPP_IMPORT_RE = re.compile(r'\w*import ([a-zA-Z0-9]+);')
CPP_EXPORT_RE = re.compile(r'\w*export module ([a-zA-Z0-9]+);')

FORTRAN_INCLUDE_PAT = r"^\s*include\s*['\"](\w+\.\w+)['\"]"
FORTRAN_MODULE_PAT = r"^\s*\bmodule\b\s+(\w+)\s*(?:!+.*)*$"
FORTRAN_SUBMOD_PAT = r"^\s*\bsubmodule\b\s*\((\w+:?\w+)\)\s*(\w+)"
FORTRAN_USE_PAT = r"^\s*use,?\s*(?:non_intrinsic)?\s*(?:::)?\s*(\w+)"

FORTRAN_MODULE_RE = re.compile(FORTRAN_MODULE_PAT, re.IGNORECASE)
FORTRAN_SUBMOD_RE = re.compile(FORTRAN_SUBMOD_PAT, re.IGNORECASE)
FORTRAN_USE_RE = re.compile(FORTRAN_USE_PAT, re.IGNORECASE)

class DependencyScanner:
    def __init__(self, pickle_file: str, outfile: str):
        with open(pickle_file, 'rb') as pf:
            self.target_data: TargetDependencyScannerInfo = pickle.load(pf)
        self.outfile = outfile
        self.sources = self.target_data.sources
        self.provided_by: T.Dict[str, str] = {}
        self.exports: T.Dict[str, str] = {}
        self.imports: collections.defaultdict[str, T.List[str]] = collections.defaultdict(list)
        self.sources_with_exports: T.List[str] = []

    def scan_file(self, fname: str, lang: Literal['cpp', 'fortran']) -> None:
        if lang == 'fortran':
            self.scan_fortran_file(fname)
        else:
            self.scan_cpp_file(fname)

    def scan_fortran_file(self, fname: str) -> None:
        fpath = pathlib.Path(fname)
        modules_in_this_file = set()
        for line in fpath.read_text(encoding='utf-8', errors='ignore').split('\n'):
            import_match = FORTRAN_USE_RE.match(line)
            export_match = FORTRAN_MODULE_RE.match(line)
            submodule_export_match = FORTRAN_SUBMOD_RE.match(line)
            if import_match:
                needed = import_match.group(1).lower()
                # In Fortran you have an using declaration also for the module
                # you define in the same file. Prevent circular dependencies.
                if needed not in modules_in_this_file:
                    self.imports[fname].append(needed)
            if export_match:
                exported_module = export_match.group(1).lower()
                assert exported_module not in modules_in_this_file
                modules_in_this_file.add(exported_module)
                if exported_module in self.provided_by:
                    raise RuntimeError(f'Multiple files provide module {exported_module}.')
                self.sources_with_exports.append(fname)
                self.provided_by[exported_module] = fname
                self.exports[fname] = exported_module
            if submodule_export_match:
                # Store submodule "Foo" "Bar" as "foo:bar".
                # A submodule declaration can be both an import and an export declaration:
                #
                # submodule (a1:a2) a3
                #  - requires a1@a2.smod
                #  - produces a1@a3.smod
                parent_module_name_full = submodule_export_match.group(1).lower()
                parent_module_name = parent_module_name_full.split(':')[0]
                submodule_name = submodule_export_match.group(2).lower()
                concat_name = f'{parent_module_name}:{submodule_name}'
                self.sources_with_exports.append(fname)
                self.provided_by[concat_name] = fname
                self.exports[fname] = concat_name
                # Fortran requires that the immediate parent module must be built
                # before the current one. Thus:
                #
                # submodule (parent) parent   <- requires parent.mod (really parent.smod, but they are created at the same time)
                # submodule (a1:a2) a3        <- requires a1@a2.smod
                #
                # a3 does not depend on the a1 parent module directly, only transitively.
                self.imports[fname].append(parent_module_name_full)

    def scan_cpp_file(self, fname: str) -> None:
        fpath = pathlib.Path(fname)
        for line in fpath.read_text(encoding='utf-8', errors='ignore').split('\n'):
            import_match = CPP_IMPORT_RE.match(line)
            export_match = CPP_EXPORT_RE.match(line)
            if import_match:
                needed = import_match.group(1)
                self.imports[fname].append(needed)
            if export_match:
                exported_module = export_match.group(1)
                if exported_module in self.provided_by:
                    raise RuntimeError(f'Multiple files provide module {exported_module}.')
                self.sources_with_exports.append(fname)
                self.provided_by[exported_module] = fname
                self.exports[fname] = exported_module

    def module_name_for(self, src: str, lang: Literal['cpp', 'fortran']) -> str:
        if lang == 'fortran':
            exported = self.exports[src]
            # Module foo:bar goes to a file name foo@bar.smod
            # Module Foo goes to a file name foo.mod
            namebase = exported.replace(':', '@')
            if ':' in exported:
                extension = 'smod'
            else:
                extension = 'mod'
            return os.path.join(self.target_data.private_dir, f'{namebase}.{extension}')
        return '{}.ifc'.format(self.exports[src])

    def scan(self) -> int:
        for s, lang in self.sources:
            self.scan_file(s, lang)
        description: Description = {
            'version': 1,
            'revision': 0,
            'rules': [],
        }
        for src, lang in self.sources:
            rule: Rule = {
                'primary-output': self.target_data.source2object[src],
                'requires': [],
                'provides': [],
            }
            if src in self.sources_with_exports:
                rule['outputs'] = [self.module_name_for(src, lang)]
            if src in self.imports:
                for modname in self.imports[src]:
                    provider_src = self.provided_by.get(modname)
                    if provider_src == src:
                        continue
                    rule['requires'].append({
                        'logical-name': modname,
                    })
                    if provider_src:
                        rule['requires'][-1].update({
                            'source-path': provider_src,
                            'compiled-module-path': self.module_name_for(provider_src, lang),
                        })
            if src in self.exports:
                modname = self.exports[src]
                rule['provides'].append({
                    'logical-name': modname,
                    'source-path': src,
                    'compiled-module-path': self.module_name_for(src, lang),
                })
            description['rules'].append(rule)

        with open(self.outfile, 'w', encoding='utf-8') as f:
            json.dump(description, f)

        return 0

def run(args: T.List[str]) -> int:
    assert len(args) == 2, 'got wrong number of arguments!'
    outfile, pickle_file = args
    scanner = DependencyScanner(pickle_file, outfile)
    return scanner.scan()
