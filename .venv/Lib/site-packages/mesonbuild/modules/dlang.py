# SPDX-License-Identifier: Apache-2.0
# Copyright 2018 The Meson development team

# This file contains the detection logic for external dependencies that
# are UI-related.
from __future__ import annotations

import json
import os
import typing as T


from . import ExtensionModule, ModuleInfo
from .. import mlog
from ..build import InvalidArguments
from ..dependencies import Dependency
from ..dependencies.dub import DubDependency
from ..interpreterbase import typed_pos_args
from ..mesonlib import Popen_safe, MesonException, listify

if T.TYPE_CHECKING:
    from typing_extensions import Literal, TypeAlias

    from . import ModuleState
    from ..interpreter.interpreter import Interpreter
    from ..interpreterbase.baseobjects import TYPE_kwargs
    from ..programs import Program

    _JSONTypes: TypeAlias = T.Union[str, int, bool, None, T.List['_JSONTypes'], T.Dict[str, '_JSONTypes']]


class DlangModule(ExtensionModule):
    class_dubbin: T.Union[Program, Literal[False], None] = None
    init_dub = False

    dubbin: T.Union[Program, Literal[False], None]

    INFO = ModuleInfo('dlang', '0.48.0')

    def __init__(self, interpreter: Interpreter):
        super().__init__(interpreter)
        self.methods.update({
            'generate_dub_file': self.generate_dub_file,
        })

    def _init_dub(self, state: ModuleState) -> None:
        if DlangModule.class_dubbin is None and DubDependency.class_dubbin is not None:
            self.dubbin = DubDependency.class_dubbin[0]
            DlangModule.class_dubbin = self.dubbin
        else:
            self.dubbin = DlangModule.class_dubbin

        if DlangModule.class_dubbin is None:
            self.dubbin = self.check_dub(state)
            DlangModule.class_dubbin = self.dubbin
        else:
            self.dubbin = DlangModule.class_dubbin

        if not self.dubbin:
            if not self.dubbin:
                raise MesonException('DUB not found.')

    @typed_pos_args('dlang.generate_dub_file', str, str)
    def generate_dub_file(self, state: ModuleState, args: T.Tuple[str, str], kwargs: TYPE_kwargs) -> None:
        if not DlangModule.init_dub:
            self._init_dub(state)

        config: T.Dict[str, _JSONTypes] = {
            'name': args[0]
        }

        config_path = os.path.join(args[1], 'dub.json')
        if os.path.exists(config_path):
            with open(config_path, encoding='utf-8') as ofile:
                try:
                    config = json.load(ofile)
                except ValueError:
                    mlog.warning('Failed to load the data in dub.json')

        warn_publishing = ['description', 'license']
        for arg in warn_publishing:
            if arg not in kwargs and \
               arg not in config:
                mlog.warning('Without', mlog.bold(arg), 'the DUB package can\'t be published')

        for key, value in kwargs.items():
            if key == 'dependencies':
                values = listify(value, flatten=False)
                data: T.Dict[str, _JSONTypes] = {}
                for dep in values:
                    if isinstance(dep, Dependency):
                        name = dep.get_name()
                        ret, res = self._call_dubbin(['describe', name])
                        if ret == 0:
                            version = dep.get_version()
                            if version is None:
                                data[name] = ''
                            else:
                                data[name] = version
                config[key] = data
            else:
                def _do_validate(v: object) -> _JSONTypes:
                    if not isinstance(v, (str, int, bool, list, dict)):
                        raise InvalidArguments('keyword arguments must be strings, numbers, booleans, arrays, or dictionaries of such')
                    if isinstance(v, list):
                        for e in v:
                            _do_validate(e)
                    if isinstance(v, dict):
                        for e in v.values():
                            _do_validate(e)
                    return T.cast('_JSONTypes', v)

                config[key] = _do_validate(value)

        with open(config_path, 'w', encoding='utf-8') as ofile:
            ofile.write(json.dumps(config, indent=4, ensure_ascii=False))

    def _call_dubbin(self, args: T.List[str], env: T.Optional[T.Mapping[str, str]] = None) -> T.Tuple[int, str]:
        assert self.dubbin is not None and self.dubbin is not False, 'for mypy'
        p, out = Popen_safe(self.dubbin.get_command() + args, env=env)[0:2]
        return p.returncode, out.strip()

    def check_dub(self, state: ModuleState) -> T.Union[Program, Literal[False]]:
        dubbin = state.find_program('dub', silent=True)
        if dubbin.found():
            try:
                p, out = Popen_safe(dubbin.get_command() + ['--version'])[0:2]
                if p.returncode != 0:
                    mlog.warning('Found dub {!r} but couldn\'t run it'
                                 ''.format(' '.join(dubbin.get_command())))
                    # Set to False instead of None to signify that we've already
                    # searched for it and not found it
                else:
                    mlog.log('Found DUB:', mlog.green('YES'), ':', mlog.bold(dubbin.get_path() or ''),
                             '({})'.format(out.strip()))
                    return dubbin
            except (FileNotFoundError, PermissionError):
                pass
        mlog.log('Found DUB:', mlog.red('NO'))
        return False

def initialize(interp: Interpreter) -> DlangModule:
    return DlangModule(interp)
