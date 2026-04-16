from __future__ import annotations

import importlib
import shutil
import json
import typing as T

from ..mesonlib import MesonException, Popen_safe
if T.TYPE_CHECKING:
    from types import ModuleType


# tomllib is present in python 3.11, before that it is a pypi module called tomli,
# we try to import tomllib, then tomli,
tomllib: T.Optional[ModuleType] = None
toml2json: T.Optional[str] = None
for t in ['tomllib', 'tomli']:
    try:
        tomllib = importlib.import_module(t)
        break
    except ImportError:
        pass
else:
    # TODO: it would be better to use an Executable here, which could be looked
    # up in the cross file or provided by a wrap. However, that will have to be
    # passed in externally, since we don't have (and I don't think we should),
    # have access to the `Environment` for that in this module.
    toml2json = shutil.which('toml2json')

class TomlImplementationMissing(MesonException):
    pass


class CargoTomlError(MesonException):
    """Exception for TOML parsing errors, keeping proper location info."""


def load_toml(filename: str) -> T.Dict[str, object]:
    if tomllib:
        try:
            with open(filename, 'rb') as f:
                raw = tomllib.load(f)
        except tomllib.TOMLDecodeError as e:
            if hasattr(e, 'msg'):
                raise CargoTomlError(e.msg, file=filename, lineno=e.lineno, colno=e.colno) from e
            else:
                raise CargoTomlError(str(e), file=filename) from e
    else:
        if toml2json is None:
            raise TomlImplementationMissing('Could not find an implementation of tomllib, nor toml2json')

        p, out, err = Popen_safe([toml2json, filename])
        if p.returncode != 0:
            error_msg = err.strip() or 'toml2json failed to decode TOML'
            raise CargoTomlError(error_msg, file=filename)

        raw = json.loads(out)

    # tomllib.load() returns T.Dict[str, T.Any] but not other implementations.
    return T.cast('T.Dict[str, object]', raw)
