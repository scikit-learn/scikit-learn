# SPDX-License-Identifier: Apache-2.0
# Copyright 2021 The Meson development team

from __future__ import annotations

import pathlib
import typing as T

from mesonbuild import mesonlib
from mesonbuild.build import CustomTarget, CustomTargetIndex, GeneratedList, Target
from mesonbuild.compilers import detect_compiler_for
from mesonbuild.interpreterbase.decorators import ContainerTypeInfo, FeatureDeprecated, FeatureNew, KwargInfo, typed_pos_args, typed_kwargs
from mesonbuild.mesonlib import version_compare, MachineChoice
from . import NewExtensionModule, ModuleReturnValue, ModuleInfo
from ..interpreter.type_checking import NoneType

if T.TYPE_CHECKING:
    from . import ModuleState
    from ..compilers import Compiler
    from ..interpreter import Interpreter

class JavaModule(NewExtensionModule):

    INFO = ModuleInfo('java', '0.60.0')

    def __init__(self, interpreter: Interpreter):
        super().__init__()
        self.methods.update({
            'generate_native_headers': self.generate_native_headers,
            'native_headers': self.native_headers,
        })

    def __get_java_compiler(self, state: ModuleState) -> Compiler:
        if 'java' not in state.environment.coredata.compilers[MachineChoice.BUILD]:
            detect_compiler_for(state.environment, 'java', MachineChoice.BUILD, False, state.subproject)
        return state.environment.coredata.compilers[MachineChoice.BUILD]['java']

    @FeatureNew('java.generate_native_headers', '0.62.0')
    @FeatureDeprecated('java.generate_native_headers', '1.0.0')
    @typed_pos_args(
        'java.generate_native_headers',
        varargs=(str, mesonlib.File, Target, CustomTargetIndex, GeneratedList))
    @typed_kwargs(
        'java.generate_native_headers',
        KwargInfo('classes', ContainerTypeInfo(list, str), default=[], listify=True, required=True),
        KwargInfo('package', (str, NoneType), default=None))
    def generate_native_headers(self, state: ModuleState, args: T.Tuple[T.List[mesonlib.FileOrString]],
                                kwargs: T.Dict[str, T.Optional[str]]) -> ModuleReturnValue:
        return self.__native_headers(state, args, kwargs)

    @FeatureNew('java.native_headers', '1.0.0')
    @typed_pos_args(
        'java.native_headers',
        varargs=(str, mesonlib.File, Target, CustomTargetIndex, GeneratedList))
    @typed_kwargs(
        'java.native_headers',
        KwargInfo('classes', ContainerTypeInfo(list, str), default=[], listify=True, required=True),
        KwargInfo('package', (str, NoneType), default=None))
    def native_headers(self, state: ModuleState, args: T.Tuple[T.List[mesonlib.FileOrString]],
                       kwargs: T.Dict[str, T.Optional[str]]) -> ModuleReturnValue:
        return self.__native_headers(state, args, kwargs)

    def __native_headers(self, state: ModuleState, args: T.Tuple[T.List[mesonlib.FileOrString]],
                         kwargs: T.Dict[str, T.Optional[str]]) -> ModuleReturnValue:
        classes = T.cast('T.List[str]', kwargs.get('classes'))
        package = kwargs.get('package')

        if package:
            sanitized_package = package.replace("-", "_").replace(".", "_")

        headers: T.List[str] = []
        for clazz in classes:
            sanitized_clazz = clazz.replace(".", "_")
            if package:
                headers.append(f'{sanitized_package}_{sanitized_clazz}.h')
            else:
                headers.append(f'{sanitized_clazz}.h')

        javac = self.__get_java_compiler(state)

        command = mesonlib.listify([
            javac.exelist,
            '-d',
            '@PRIVATE_DIR@',
            '-h',
            state.subdir,
            '@INPUT@',
        ])

        prefix = classes[0] if not package else package

        target = CustomTarget(f'{prefix}-native-headers',
                              state.subdir,
                              state.subproject,
                              state.environment,
                              command,
                              sources=args[0], outputs=headers, backend=state.backend)

        # It is only known that 1.8.0 won't pre-create the directory. 11 and 16
        # do not exhibit this behavior.
        if version_compare(javac.version, '1.8.0'):
            pathlib.Path(state.backend.get_target_private_dir_abs(target)).mkdir(parents=True, exist_ok=True)

        return ModuleReturnValue(target, [target])

def initialize(*args: T.Any, **kwargs: T.Any) -> JavaModule:
    return JavaModule(*args, **kwargs)
