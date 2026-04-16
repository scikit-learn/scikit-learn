# SPDX-License-Identifier: Apache-2.0
# Copyright 2017 The Meson development team

from __future__ import annotations

import typing as T

from .. import mlog
from .. import build
from ..compilers import Compiler
from ..interpreter.type_checking import BT_SOURCES_KW, STATIC_LIB_KWS
from ..interpreterbase.decorators import KwargInfo, permittedKwargs, typed_pos_args, typed_kwargs

from . import ExtensionModule, ModuleInfo

if T.TYPE_CHECKING:
    from typing_extensions import Literal

    from . import ModuleState
    from ..interpreter import Interpreter, kwargs as kwtypes
    from ..interpreter.type_checking import SourcesVarargsType

    class CheckKw(kwtypes.StaticLibrary):

        compiler: Compiler
        mmx: SourcesVarargsType
        sse: SourcesVarargsType
        sse2: SourcesVarargsType
        sse3: SourcesVarargsType
        ssse3: SourcesVarargsType
        sse41: SourcesVarargsType
        sse42: SourcesVarargsType
        avx: SourcesVarargsType
        avx2: SourcesVarargsType
        neon: SourcesVarargsType


# FIXME add Altivec and AVX512.
ISETS = (
    'mmx',
    'sse',
    'sse2',
    'sse3',
    'ssse3',
    'sse41',
    'sse42',
    'avx',
    'avx2',
    'neon',
)


class SimdModule(ExtensionModule):

    INFO = ModuleInfo('SIMD', '0.42.0', unstable=True)

    def __init__(self, interpreter: Interpreter):
        super().__init__(interpreter)
        self.methods.update({
            'check': self.check,
        })

    @typed_pos_args('simd.check', str)
    @typed_kwargs('simd.check',
                  KwargInfo('compiler', Compiler, required=True),
                  *[BT_SOURCES_KW.evolve(name=iset, default=None) for iset in ISETS],
                  *[a for a in STATIC_LIB_KWS if a.name != 'sources'])
    @permittedKwargs({'compiler', *ISETS, *build.known_stlib_kwargs}) # Also remove this, per above comment
    def check(self, state: ModuleState, args: T.Tuple[str], kwargs: CheckKw) -> T.List[T.Union[T.List[build.StaticLibrary], build.ConfigurationData]]:
        result: T.List[build.StaticLibrary] = []

        local_kwargs = set((*ISETS, 'compiler'))
        static_lib_kwargs = T.cast('kwtypes.StaticLibrary', {k: v for k, v in kwargs.items() if k not in local_kwargs})

        prefix = args[0]
        compiler = kwargs['compiler']
        conf = build.ConfigurationData()

        for iset in ISETS:
            sources = kwargs[iset]  # type: ignore[literal-required]
            if sources is None:
                continue

            compile_args = compiler.get_instruction_set_args(iset)
            if compile_args is None:
                mlog.log(f'Compiler supports {iset}:', mlog.red('NO'))
                continue

            if not compiler.has_multi_arguments(compile_args)[0]:
                mlog.log(f'Compiler supports {iset}:', mlog.red('NO'))
                continue
            mlog.log(f'Compiler supports {iset}:', mlog.green('YES'))
            conf.values['HAVE_' + iset.upper()] = ('1', f'Compiler supports {iset}.')

            libname = prefix + '_' + iset
            lib_kwargs = static_lib_kwargs.copy()
            lib_kwargs['sources'] = sources

            # Add compile args we derived above to those the user provided us
            # This cast is mostly a lie, but it makes mypy happy and calculating the actual key name would be difficult
            langarg_key = T.cast('Literal["c_args"]', compiler.get_language() + '_args')
            old_lang_args = lib_kwargs[langarg_key]
            all_lang_args = old_lang_args + compile_args
            lib_kwargs[langarg_key] = all_lang_args

            lib = self.interpreter.build_target(state.current_node, (libname, []), lib_kwargs, build.StaticLibrary)

            result.append(lib)

        return [result, conf]

def initialize(interp: Interpreter) -> SimdModule:
    return SimdModule(interp)
