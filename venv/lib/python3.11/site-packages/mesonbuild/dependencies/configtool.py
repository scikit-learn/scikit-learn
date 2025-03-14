# SPDX-License-Identifier: Apache-2.0
# Copyright 2013-2021 The Meson development team

from __future__ import annotations

from .base import ExternalDependency, DependencyException, DependencyTypeName
from ..mesonlib import listify, Popen_safe, Popen_safe_logged, split_args, version_compare, version_compare_many
from ..programs import find_external_program
from .. import mlog
import re
import typing as T

from mesonbuild import mesonlib

if T.TYPE_CHECKING:
    from ..environment import Environment
    from ..interpreter.type_checking import PkgConfigDefineType

class ConfigToolDependency(ExternalDependency):

    """Class representing dependencies found using a config tool.

    Takes the following extra keys in kwargs that it uses internally:
    :tools List[str]: A list of tool names to use
    :version_arg str: The argument to pass to the tool to get its version
    :skip_version str: The argument to pass to the tool to ignore its version
        (if ``version_arg`` fails, but it may start accepting it in the future)
        Because some tools are stupid and don't accept --version
    :returncode_value int: The value of the correct returncode
        Because some tools are stupid and don't return 0
    """

    tools: T.Optional[T.List[str]] = None
    tool_name: T.Optional[str] = None
    version_arg = '--version'
    skip_version: T.Optional[str] = None
    allow_default_for_cross = False
    __strip_version = re.compile(r'^[0-9][0-9.]+')

    def __init__(self, name: str, environment: 'Environment', kwargs: T.Dict[str, T.Any], language: T.Optional[str] = None, exclude_paths: T.Optional[T.List[str]] = None):
        super().__init__(DependencyTypeName('config-tool'), environment, kwargs, language=language)
        self.name = name
        # You may want to overwrite the class version in some cases
        self.tools = listify(kwargs.get('tools', self.tools))
        if not self.tool_name:
            self.tool_name = self.tools[0]
        if 'version_arg' in kwargs:
            self.version_arg = kwargs['version_arg']

        req_version_raw = kwargs.get('version', None)
        if req_version_raw is not None:
            req_version = mesonlib.stringlistify(req_version_raw)
        else:
            req_version = []
        tool, version = self.find_config(req_version, kwargs.get('returncode_value', 0), exclude_paths=exclude_paths)
        self.config = tool
        self.is_found = self.report_config(version, req_version)
        if not self.is_found:
            self.config = None
            return
        self.version = version

    def _sanitize_version(self, version: str) -> str:
        """Remove any non-numeric, non-point version suffixes."""
        m = self.__strip_version.match(version)
        if m:
            # Ensure that there isn't a trailing '.', such as an input like
            # `1.2.3.git-1234`
            return m.group(0).rstrip('.')
        return version

    def _check_and_get_version(self, tool: T.List[str], returncode: int) -> T.Tuple[bool, T.Union[str, None]]:
        """Check whether a command is valid and get its version"""
        p, out = Popen_safe(tool + [self.version_arg])[:2]
        valid = True
        if p.returncode != returncode:
            if self.skip_version:
                # maybe the executable is valid even if it doesn't support --version
                p = Popen_safe(tool + [self.skip_version])[0]
                if p.returncode != returncode:
                    valid = False
            else:
                valid = False
        version = self._sanitize_version(out.strip())
        return valid, version

    def find_config(self, versions: T.List[str], returncode: int = 0, exclude_paths: T.Optional[T.List[str]] = None) \
            -> T.Tuple[T.Optional[T.List[str]], T.Optional[str]]:
        """Helper method that searches for config tool binaries in PATH and
        returns the one that best matches the given version requirements.
        """
        exclude_paths = [] if exclude_paths is None else exclude_paths
        best_match: T.Tuple[T.Optional[T.List[str]], T.Optional[str]] = (None, None)
        for potential_bin in find_external_program(
                self.env, self.for_machine, self.tool_name,
                self.tool_name, self.tools, exclude_paths=exclude_paths,
                allow_default_for_cross=self.allow_default_for_cross):
            if not potential_bin.found():
                continue
            tool = potential_bin.get_command()
            try:
                valid, version = self._check_and_get_version(tool, returncode)
            except (FileNotFoundError, PermissionError):
                continue
            if not valid:
                continue

            # Some tools, like pcap-config don't supply a version, but also
            # don't fail with --version, in that case just assume that there is
            # only one version and return it.
            if not version:
                return (tool, None)
            if versions:
                is_found = version_compare_many(version, versions)[0]
                # This allows returning a found version without a config tool,
                # which is useful to inform the user that you found version x,
                # but y was required.
                if not is_found:
                    tool = None
            if best_match[1]:
                if version_compare(version, '> {}'.format(best_match[1])):
                    best_match = (tool, version)
            else:
                best_match = (tool, version)

        return best_match

    def report_config(self, version: T.Optional[str], req_version: T.List[str]) -> bool:
        """Helper method to print messages about the tool."""

        found_msg: T.List[T.Union[str, mlog.AnsiDecorator]] = [mlog.bold(self.tool_name), 'found:']

        if self.config is None:
            found_msg.append(mlog.red('NO'))
            if version is not None and req_version:
                found_msg.append(f'found {version!r} but need {req_version!r}')
            elif req_version:
                found_msg.append(f'need {req_version!r}')
        else:
            found_msg += [mlog.green('YES'), '({})'.format(' '.join(self.config)), version]

        mlog.log(*found_msg)

        return self.config is not None

    def get_config_value(self, args: T.List[str], stage: str) -> T.List[str]:
        p, out, err = Popen_safe_logged(self.config + args)
        if p.returncode != 0:
            if self.required:
                raise DependencyException(f'Could not generate {stage} for {self.name}.\n{err}')
            return []
        return split_args(out)

    def get_variable_args(self, variable_name: str) -> T.List[str]:
        return [f'--{variable_name}']

    @staticmethod
    def log_tried() -> str:
        return 'config-tool'

    def get_variable(self, *, cmake: T.Optional[str] = None, pkgconfig: T.Optional[str] = None,
                     configtool: T.Optional[str] = None, internal: T.Optional[str] = None,
                     system: T.Optional[str] = None, default_value: T.Optional[str] = None,
                     pkgconfig_define: PkgConfigDefineType = None) -> str:
        if configtool:
            p, out, _ = Popen_safe(self.config + self.get_variable_args(configtool))
            if p.returncode == 0:
                variable = out.strip()
                mlog.debug(f'Got config-tool variable {configtool} : {variable}')
                return variable
        if default_value is not None:
            return default_value
        raise DependencyException(f'Could not get config-tool variable and no default provided for {self!r}')
