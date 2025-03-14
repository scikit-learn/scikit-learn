"""setuptools.errors

Provides exceptions used by setuptools modules.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from distutils import errors as _distutils_errors

if TYPE_CHECKING:
    from typing_extensions import TypeAlias

# Re-export errors from distutils to facilitate the migration to PEP632

ByteCompileError: TypeAlias = _distutils_errors.DistutilsByteCompileError
CCompilerError: TypeAlias = _distutils_errors.CCompilerError
ClassError: TypeAlias = _distutils_errors.DistutilsClassError
CompileError: TypeAlias = _distutils_errors.CompileError
ExecError: TypeAlias = _distutils_errors.DistutilsExecError
FileError: TypeAlias = _distutils_errors.DistutilsFileError
InternalError: TypeAlias = _distutils_errors.DistutilsInternalError
LibError: TypeAlias = _distutils_errors.LibError
LinkError: TypeAlias = _distutils_errors.LinkError
ModuleError: TypeAlias = _distutils_errors.DistutilsModuleError
OptionError: TypeAlias = _distutils_errors.DistutilsOptionError
PlatformError: TypeAlias = _distutils_errors.DistutilsPlatformError
PreprocessError: TypeAlias = _distutils_errors.PreprocessError
SetupError: TypeAlias = _distutils_errors.DistutilsSetupError
TemplateError: TypeAlias = _distutils_errors.DistutilsTemplateError
UnknownFileError: TypeAlias = _distutils_errors.UnknownFileError

# The root error class in the hierarchy
BaseError: TypeAlias = _distutils_errors.DistutilsError


class InvalidConfigError(OptionError):
    """Error used for invalid configurations."""


class RemovedConfigError(OptionError):
    """Error used for configurations that were deprecated and removed."""


class RemovedCommandError(BaseError, RuntimeError):
    """Error used for commands that have been removed in setuptools.

    Since ``setuptools`` is built on ``distutils``, simply removing a command
    from ``setuptools`` will make the behavior fall back to ``distutils``; this
    error is raised if a command exists in ``distutils`` but has been actively
    removed in ``setuptools``.
    """


class PackageDiscoveryError(BaseError, RuntimeError):
    """Impossible to perform automatic discovery of packages and/or modules.

    The current project layout or given discovery options can lead to problems when
    scanning the project directory.

    Setuptools might also refuse to complete auto-discovery if an error prone condition
    is detected (e.g. when a project is organised as a flat-layout but contains
    multiple directories that can be taken as top-level packages inside a single
    distribution [*]_). In these situations the users are encouraged to be explicit
    about which packages to include or to make the discovery parameters more specific.

    .. [*] Since multi-package distributions are uncommon it is very likely that the
       developers did not intend for all the directories to be packaged, and are just
       leaving auxiliary code in the repository top-level, such as maintenance-related
       scripts.
    """
