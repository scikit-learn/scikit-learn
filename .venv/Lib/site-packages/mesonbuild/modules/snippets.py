# Copyright 2025 The Meson development team

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import annotations
import textwrap
import typing as T

from pathlib import Path

from . import NewExtensionModule, ModuleInfo
from ..interpreterbase import KwargInfo, typed_kwargs, typed_pos_args
from ..interpreter.type_checking import NoneType
from .. import mesonlib

if T.TYPE_CHECKING:
    from typing_extensions import TypedDict
    from . import ModuleState

    class SymbolVisibilityHeaderKW(TypedDict):
        namespace: T.Optional[str]
        api: T.Optional[str]
        compilation: T.Optional[str]
        static_compilation: T.Optional[str]
        static_only: bool


class SnippetsModule(NewExtensionModule):

    INFO = ModuleInfo('snippets', '1.10.0')

    def __init__(self) -> None:
        super().__init__()
        self.methods.update({
            'symbol_visibility_header': self.symbol_visibility_header_method,
        })

    @typed_kwargs('snippets.symbol_visibility_header',
                  KwargInfo('namespace', (str, NoneType)),
                  KwargInfo('api', (str, NoneType)),
                  KwargInfo('compilation', (str, NoneType)),
                  KwargInfo('static_compilation', (str, NoneType)),
                  KwargInfo('static_only', (bool, NoneType)))
    @typed_pos_args('snippets.symbol_visibility_header', str)
    def symbol_visibility_header_method(self, state: ModuleState, args: T.Tuple[str], kwargs: 'SymbolVisibilityHeaderKW') -> mesonlib.File:
        header_name = args[0]
        namespace = kwargs['namespace'] or state.project_name
        namespace = mesonlib.underscorify(namespace).upper()
        if namespace[0].isdigit():
            namespace = f'_{namespace}'
        api = kwargs['api'] or f'{namespace}_API'
        compilation = kwargs['compilation'] or f'{namespace}_COMPILATION'
        static_compilation = kwargs['static_compilation'] or f'{namespace}_STATIC_COMPILATION'
        static_only = kwargs['static_only']
        if static_only is None:
            default_library = state.get_option('default_library')
            static_only = default_library == 'static'
        content = textwrap.dedent('''\
            // SPDX-license-identifier: 0BSD OR CC0-1.0 OR WTFPL OR Apache-2.0 OR LGPL-2.0-or-later
            #pragma once
            ''')
        if static_only:
            content += textwrap.dedent(f'''
                #ifndef {static_compilation}
                #  define {static_compilation}
                #endif /* {static_compilation} */
                ''')
        content += textwrap.dedent(f'''
            #if (defined(_WIN32) || defined(__CYGWIN__)) && !defined({static_compilation})
            #  define {api}_EXPORT __declspec(dllexport)
            #  define {api}_IMPORT __declspec(dllimport)
            #elif defined(__OS2__) && !defined({static_compilation})
            #  define {api}_EXPORT __declspec(dllexport)
            #  define {api}_IMPORT
            #elif __GNUC__ >= 4
            #  define {api}_EXPORT __attribute__((visibility("default")))
            #  define {api}_IMPORT
            #else
            #  define {api}_EXPORT
            #  define {api}_IMPORT
            #endif

            #ifdef {compilation}
            #  define {api} {api}_EXPORT extern
            #else
            #  define {api} {api}_IMPORT extern
            #endif
            ''')
        header_path = Path(state.environment.get_build_dir(), state.subdir, header_name)
        header_path.write_text(content, encoding='utf-8')
        return mesonlib.File.from_built_file(state.subdir, header_name)

def initialize(*args: T.Any, **kwargs: T.Any) -> SnippetsModule:
    return SnippetsModule()
