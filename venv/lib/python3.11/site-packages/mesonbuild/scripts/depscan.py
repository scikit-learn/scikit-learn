# SPDX-License-Identifier: Apache-2.0
# Copyright 2020 The Meson development team
# Copyright © 2023 Intel Corporation

from __future__ import annotations

import collections
import os
import pathlib
import pickle
import re
import typing as T

from ..backend.ninjabackend import ninja_quote

if T.TYPE_CHECKING:
    from typing_extensions import Literal
    from ..backend.ninjabackend import TargetDependencyScannerInfo

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
        self.needs: collections.defaultdict[str, T.List[str]] = collections.defaultdict(list)
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
                    self.needs[fname].append(needed)
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
                self.needs[fname].append(parent_module_name_full)

    def scan_cpp_file(self, fname: str) -> None:
        fpath = pathlib.Path(fname)
        for line in fpath.read_text(encoding='utf-8', errors='ignore').split('\n'):
            import_match = CPP_IMPORT_RE.match(line)
            export_match = CPP_EXPORT_RE.match(line)
            if import_match:
                needed = import_match.group(1)
                self.needs[fname].append(needed)
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
        with open(self.outfile, 'w', encoding='utf-8') as ofile:
            ofile.write('ninja_dyndep_version = 1\n')
            for src, lang in self.sources:
                objfilename = self.target_data.source2object[src]
                mods_and_submods_needed = []
                module_files_generated = []
                module_files_needed = []
                if src in self.sources_with_exports:
                    module_files_generated.append(self.module_name_for(src, lang))
                if src in self.needs:
                    for modname in self.needs[src]:
                        if modname not in self.provided_by:
                            # Nothing provides this module, we assume that it
                            # comes from a dependency library somewhere and is
                            # already built by the time this compilation starts.
                            pass
                        else:
                            mods_and_submods_needed.append(modname)

                for modname in mods_and_submods_needed:
                    provider_src = self.provided_by[modname]
                    provider_modfile = self.module_name_for(provider_src, lang)
                    # Prune self-dependencies
                    if provider_src != src:
                        module_files_needed.append(provider_modfile)

                quoted_objfilename = ninja_quote(objfilename, True)
                quoted_module_files_generated = [ninja_quote(x, True) for x in module_files_generated]
                quoted_module_files_needed = [ninja_quote(x, True) for x in module_files_needed]
                if quoted_module_files_generated:
                    mod_gen = '| ' + ' '.join(quoted_module_files_generated)
                else:
                    mod_gen = ''
                if quoted_module_files_needed:
                    mod_dep = '| ' + ' '.join(quoted_module_files_needed)
                else:
                    mod_dep = ''
                build_line = 'build {} {}: dyndep {}'.format(quoted_objfilename,
                                                             mod_gen,
                                                             mod_dep)
                ofile.write(build_line + '\n')
        return 0

def run(args: T.List[str]) -> int:
    assert len(args) == 2, 'got wrong number of arguments!'
    outfile, pickle_file = args
    scanner = DependencyScanner(pickle_file, outfile)
    return scanner.scan()
