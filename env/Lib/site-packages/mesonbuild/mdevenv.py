from __future__ import annotations

import os, subprocess
import argparse
import tempfile
import shutil
import itertools
import typing as T

from pathlib import Path
from . import build, minstall
from .mesonlib import (EnvironmentVariables, MesonException, join_args, is_windows, setup_vsenv,
                       get_wine_shortpath, MachineChoice, relpath)
from .options import OptionKey
from . import mlog


if T.TYPE_CHECKING:
    from .backend.backends import InstallData

POWERSHELL_EXES = {'pwsh.exe', 'powershell.exe'}

# Note: when adding arguments, please also add them to the completion
# scripts in $MESONSRC/data/shell-completions/
def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument('-C', dest='builddir', type=Path, default='.',
                        help='Path to build directory')
    parser.add_argument('--workdir', '-w', type=Path, default=None,
                        help='Directory to cd into before running (default: builddir, Since 1.0.0)')
    parser.add_argument('--dump', nargs='?', const=True,
                        help='Only print required environment (Since 0.62.0) ' +
                             'Takes an optional file path (Since 1.1.0)')
    parser.add_argument('--dump-format', default='export',
                        choices=['sh', 'export', 'vscode'],
                        help='Format used with --dump (Since 1.1.0)')
    parser.add_argument('devcmd', nargs=argparse.REMAINDER, metavar='command',
                        help='Command to run in developer environment (default: interactive shell)')

def get_windows_shell() -> T.Optional[str]:
    mesonbuild = Path(__file__).parent
    script = mesonbuild / 'scripts' / 'cmd_or_ps.ps1'
    for shell in POWERSHELL_EXES:
        try:
            command = [shell, '-noprofile', '-executionpolicy', 'bypass', '-file', str(script)]
            result = subprocess.check_output(command)
            return result.decode().strip()
        except (subprocess.CalledProcessError, OSError):
            pass
    return None

def reduce_winepath(env: T.Dict[str, str]) -> None:
    winepath = env.get('WINEPATH')
    if not winepath:
        return
    winecmd = shutil.which('wine64') or shutil.which('wine')
    if not winecmd:
        return
    env['WINEPATH'] = get_wine_shortpath([winecmd], winepath.split(';'))
    mlog.log('Meson detected wine and has set WINEPATH accordingly')

def get_env(b: build.Build, dump_fmt: T.Optional[str]) -> T.Tuple[T.Dict[str, str], T.Set[str]]:
    extra_env = EnvironmentVariables()
    extra_env.set('MESON_DEVENV', ['1'])
    extra_env.set('MESON_PROJECT_NAME', [b.project_name])

    sysroot = b.environment.properties[MachineChoice.HOST].get_sys_root()
    if sysroot:
        extra_env.set('QEMU_LD_PREFIX', [sysroot])

    env = {} if dump_fmt else os.environ.copy()
    default_fmt = '${0}' if dump_fmt in {'sh', 'export'} else None
    varnames = set()
    for i in itertools.chain(b.devenv, {extra_env}):
        env = i.get_env(env, default_fmt)
        varnames |= i.get_names()

    reduce_winepath(env)

    return env, varnames

def bash_completion_files(b: build.Build, install_data: 'InstallData') -> T.List[str]:
    from .dependencies.pkgconfig import PkgConfigDependency
    result = []
    dep = PkgConfigDependency('bash-completion', b.environment,
                              {'required': False, 'silent': True, 'version': '>=2.10'})
    if dep.found():
        prefix = b.environment.coredata.get_option(OptionKey('prefix'))
        assert isinstance(prefix, str), 'for mypy'
        datadir = b.environment.coredata.get_option(OptionKey('datadir'))
        assert isinstance(datadir, str), 'for mypy'
        datadir_abs = os.path.join(prefix, datadir)
        completionsdir = dep.get_variable(pkgconfig='completionsdir', pkgconfig_define=(('datadir', datadir_abs),))
        assert isinstance(completionsdir, str), 'for mypy'
        completionsdir_path = Path(completionsdir)
        for f in install_data.data:
            if completionsdir_path in Path(f.install_path).parents:
                result.append(f.path)
    return result

def add_gdb_auto_load(autoload_path: Path, gdb_helper: str, fname: Path) -> None:
    # Copy or symlink the GDB helper into our private directory tree
    destdir = autoload_path / fname.parent
    destdir.mkdir(parents=True, exist_ok=True)
    try:
        if is_windows():
            shutil.copy(gdb_helper, str(destdir / os.path.basename(gdb_helper)))
        else:
            os.symlink(gdb_helper, str(destdir / os.path.basename(gdb_helper)))
    except (FileExistsError, shutil.SameFileError):
        pass

def write_gdb_script(privatedir: Path, install_data: 'InstallData', workdir: Path) -> None:
    if not shutil.which('gdb'):
        return
    bdir = privatedir.parent
    autoload_basedir = privatedir / 'gdb-auto-load'
    autoload_path = Path(autoload_basedir, *bdir.parts[1:])
    have_gdb_helpers = False
    for d in install_data.data:
        if d.path.endswith('-gdb.py') or d.path.endswith('-gdb.gdb') or d.path.endswith('-gdb.scm'):
            # This GDB helper is made for a specific shared library, search if
            # we have it in our builddir.
            libname = Path(d.path).name.rsplit('-', 1)[0]
            for t in install_data.targets:
                path = Path(t.fname)
                if path.name == libname:
                    add_gdb_auto_load(autoload_path, d.path, path)
                    have_gdb_helpers = True
    if have_gdb_helpers:
        gdbinit_line = f'add-auto-load-scripts-directory {autoload_basedir}\n'
        gdbinit_path = bdir / '.gdbinit'
        first_time = False
        try:
            with gdbinit_path.open('r+', encoding='utf-8') as f:
                if gdbinit_line not in f.readlines():
                    f.write(gdbinit_line)
                    first_time = True
        except FileNotFoundError:
            gdbinit_path.write_text(gdbinit_line, encoding='utf-8')
            first_time = True
        if first_time:
            gdbinit_path = gdbinit_path.resolve()
            workdir_path = workdir.resolve()
            rel_path = Path(relpath(gdbinit_path, workdir_path))
            mlog.log('Meson detected GDB helpers and added config in', mlog.bold(str(rel_path)))
            mlog.log('To load it automatically you might need to:')
            mlog.log(' - Add', mlog.bold(f'add-auto-load-safe-path {gdbinit_path.parent}'),
                     'in', mlog.bold('~/.gdbinit'))
            if gdbinit_path.parent != workdir_path:
                mlog.log(' - Change current workdir to', mlog.bold(str(rel_path.parent)),
                         'or use', mlog.bold(f'--init-command {rel_path}'))

def dump(devenv: T.Dict[str, str], varnames: T.Set[str], dump_format: T.Optional[str], output: T.Optional[T.TextIO] = None) -> None:
    for name in varnames:
        print(f'{name}="{devenv[name]}"', file=output)
        if dump_format == 'export':
            print(f'export {name}', file=output)

def run(options: argparse.Namespace) -> int:
    privatedir = Path(options.builddir) / 'meson-private'
    buildfile = privatedir / 'build.dat'
    if not buildfile.is_file():
        raise MesonException(f'Directory {options.builddir!r} does not seem to be a Meson build directory.')
    b = build.load(options.builddir)
    workdir = options.workdir or options.builddir

    need_vsenv = T.cast('bool', b.environment.coredata.get_option(OptionKey('vsenv')))
    setup_vsenv(need_vsenv) # Call it before get_env to get vsenv vars as well
    dump_fmt = options.dump_format if options.dump else None
    devenv, varnames = get_env(b, dump_fmt)
    if options.dump:
        if options.devcmd:
            raise MesonException('--dump option does not allow running other command.')
        if options.dump is True:
            dump(devenv, varnames, dump_fmt)
        else:
            with open(options.dump, "w", encoding='utf-8') as output:
                dump(devenv, varnames, dump_fmt, output)
        return 0

    if b.environment.need_exe_wrapper():
        m = 'An executable wrapper could be required'
        exe_wrapper = b.environment.get_exe_wrapper()
        if exe_wrapper:
            cmd = ' '.join(exe_wrapper.get_command())
            m += f': {cmd}'
        mlog.log(m)

    install_data = minstall.load_install_data(str(privatedir / 'install.dat'))
    write_gdb_script(privatedir, install_data, workdir)

    args = options.devcmd
    if not args:
        prompt_prefix = f'[{b.project_name}]'
        shell_env = os.environ.get("SHELL")
        # Prefer $SHELL in a MSYS2 bash despite it being Windows
        if shell_env and os.path.exists(shell_env):
            args = [shell_env]
        elif is_windows():
            shell = get_windows_shell()
            if not shell:
                mlog.warning('Failed to determine Windows shell, fallback to cmd.exe')
            if shell in POWERSHELL_EXES:
                args = [shell, '-NoLogo', '-NoExit']
                prompt = f'function global:prompt {{  "{prompt_prefix} PS " + $PWD + "> "}}'
                args += ['-Command', prompt]
            else:
                args = [os.environ.get("COMSPEC", r"C:\WINDOWS\system32\cmd.exe")]
                args += ['/k', f'prompt {prompt_prefix} $P$G']
        else:
            args = [os.environ.get("SHELL", os.path.realpath("/bin/sh"))]
        if "bash" in args[0]:
            # Let the GC remove the tmp file
            tmprc = tempfile.NamedTemporaryFile(mode='w')
            tmprc.write('[ -e ~/.bashrc ] && . ~/.bashrc\n')
            if not os.environ.get("MESON_DISABLE_PS1_OVERRIDE"):
                tmprc.write(f'export PS1="{prompt_prefix} $PS1"\n')
            for f in bash_completion_files(b, install_data):
                tmprc.write(f'. "{f}"\n')
            tmprc.flush()
            args.append("--rcfile")
            args.append(tmprc.name)
    else:
        # Try to resolve executable using devenv's PATH
        abs_path = shutil.which(args[0], path=devenv.get('PATH', None))
        args[0] = abs_path or args[0]

    try:
        os.chdir(workdir)
        os.execvpe(args[0], args, env=devenv)
    except FileNotFoundError:
        raise MesonException(f'Command not found: {args[0]}')
    except OSError as e:
        raise MesonException(f'Command `{join_args(args)}` failed to execute: {e}')
