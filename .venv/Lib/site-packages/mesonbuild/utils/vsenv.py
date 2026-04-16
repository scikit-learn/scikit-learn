# SPDX-License-Identifier: Apache-2.0
# Copyright 2012-2025 The Meson development team

from __future__ import annotations

import os
import subprocess
import json
import pathlib
import shutil
import tempfile
import pickle

from .. import mlog
from .core import MesonException
from .universal import (is_windows, windows_detect_native_arch, windows_proof_rm,
                        get_meson_command, join_args)


__all__ = [
    'setup_vsenv',
]


# If on Windows and VS is installed but not set up in the environment,
# set it to be runnable. In this way Meson can be directly invoked
# from any shell, VS Code etc.
def _setup_vsenv(force: bool) -> bool:
    if not is_windows():
        return False
    if os.environ.get('OSTYPE') == 'cygwin':
        return False
    if 'MESON_FORCE_VSENV_FOR_UNITTEST' not in os.environ:
        # VSINSTALL is set when running setvars from a Visual Studio installation
        # Tested with Visual Studio 2012 and 2017
        if 'VSINSTALLDIR' in os.environ:
            return False
        # Check explicitly for cl when on Windows
        if shutil.which('cl.exe'):
            return False
    if not force:
        if shutil.which('cc'):
            return False
        if shutil.which('gcc'):
            return False
        if shutil.which('clang'):
            return False
        if shutil.which('clang-cl'):
            return False

    root = os.environ.get("ProgramFiles(x86)") or os.environ.get("ProgramFiles")
    bat_locator_bin = pathlib.Path(root, 'Microsoft Visual Studio/Installer/vswhere.exe')
    if not bat_locator_bin.exists():
        raise MesonException(f'Could not find {bat_locator_bin}')
    bat_json = subprocess.check_output(
        [
            str(bat_locator_bin),
            '-latest',
            '-prerelease',
            '-requiresAny',
            '-requires', 'Microsoft.VisualStudio.Component.VC.Tools.x86.x64',
            '-requires', 'Microsoft.VisualStudio.Workload.WDExpress',
            '-products', '*',
            '-utf8',
            '-format',
            'json'
        ]
    )
    bat_info = json.loads(bat_json)
    if not bat_info:
        # VS installer installed but not VS itself maybe?
        raise MesonException('Could not parse vswhere.exe output')
    bat_root = pathlib.Path(bat_info[0]['installationPath'])
    if windows_detect_native_arch() == 'arm64':
        bat_path = bat_root / 'VC/Auxiliary/Build/vcvarsarm64.bat'
        if not bat_path.exists():
            bat_path = bat_root / 'VC/Auxiliary/Build/vcvarsx86_arm64.bat'
    else:
        bat_path = bat_root / 'VC/Auxiliary/Build/vcvars64.bat'
        # if VS is not found try VS Express
        if not bat_path.exists():
            bat_path = bat_root / 'VC/Auxiliary/Build/vcvarsx86_amd64.bat'
    if not bat_path.exists():
        raise MesonException(f'Could not find {bat_path}')

    mlog.log('Activating VS', bat_info[0]['catalog']['productDisplayVersion'])
    # Write a bat file that first activates VS environment, and then calls
    # a Meson script that pickles the environment into a temp file.
    with tempfile.NamedTemporaryFile(delete=False) as env_file:
        pass
    vcvars_cmd = ['call', str(bat_path)]
    pickle_cmd = get_meson_command() + ['--internal', 'pickle_env', env_file.name]
    with tempfile.NamedTemporaryFile('w', suffix='.bat', encoding='utf-8', delete=False) as bat_file:
        bat_file.write(join_args(vcvars_cmd) + '\n')
        bat_file.write(join_args(pickle_cmd))
    try:
        subprocess.check_call([bat_file.name], stdout=subprocess.DEVNULL)
        with open(env_file.name, 'rb') as f:
            vsenv = pickle.load(f)
        for k, v in vsenv.items():
            os.environ[k] = v
    finally:
        windows_proof_rm(env_file.name)
        windows_proof_rm(bat_file.name)
    return True

def setup_vsenv(force: bool = False) -> bool:
    try:
        return _setup_vsenv(force)
    except MesonException as e:
        if force:
            raise
        mlog.warning('Failed to activate VS environment:', str(e))
        return False
