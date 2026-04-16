'''
This modules contains a distutils extension mechanism for Pythran
    * PythranExtension: is used as distutils's Extension
'''

import pythran.config as cfg
from pythran.tables import blas_requires
from pythran.utils import cxxid

from collections import defaultdict
try:
    from collections.abc import Iterable
except ImportError:
    from collections import Iterable
import os.path
import os
import re
import sys


from setuptools.command.build_ext import build_ext as LegacyBuildExt
from setuptools.extension import Extension



class PythranBuildExtMixIn(object):
    """Subclass of `distutils.command.build_ext.build_ext` which is required to
    build `PythranExtension` with the configured C++ compiler. It may also be
    subclassed if you want to combine with another build_ext class (NumPy,
    Cython implementations).

    """

    def build_extension(self, ext):
        StringTypes = str,

        def get_value(obj, key):
            var = getattr(obj, key)
            if isinstance(var, Iterable) and not isinstance(var, StringTypes):
                return var[0]
            else:
                return var

        def set_value(obj, key, value):
            var = getattr(obj, key)
            if isinstance(var, Iterable) and not isinstance(var, StringTypes):
                var[0] = value
            else:
                setattr(obj, key, value)

        prev = {
                # linux-like
                'preprocessor': None,
                'compiler_cxx': None,
                'compiler_so': None,
                'compiler_so_cxx': None,
                'compiler': None,
                'linker_exe': None,
                'linker_so': None,
                'linker_so_cxx': None,
                # Windows-like
                'cc': None,
        }
        # Backup compiler settings
        for key in list(prev.keys()):
            if hasattr(self.compiler, key):
                prev[key] = get_value(self.compiler, key)
            else:
                del prev[key]
        # try hard to modify the compiler
        if getattr(ext, 'cxx', None) is not None:
            for comp in prev:
                if hasattr(self.compiler, comp):
                    set_value(self.compiler, comp, ext.cxx)

        find_exe = None
        if getattr(ext, 'cc', None) is not None:
            try:
                import distutils._msvccompiler as msvc
                if not hasattr(msvc, "_find_exe"):
                    msvc = msvc.msvc
                # install hook
                find_exe = msvc._find_exe

                def _find_exe(exe, *args, **kwargs):
                    if exe == 'cl.exe':
                        exe = ext.cc
                    return find_exe(exe, *args, **kwargs)

                msvc._find_exe = _find_exe
            except (AttributeError, ImportError) as e:
                pass

        # In general, distutils uses -Wstrict-prototypes, but this option
        # is not valid for C++ code, only for C.  Remove it if it's there
        # to avoid a spurious warning on every compilation.
        for flag in cfg.cfg.get('compiler', "ignoreflags").split():
            for target in ('compiler_so', 'linker_so'):
                try:
                    while True:
                        getattr(self.compiler, target).remove(flag)
                except (AttributeError, ValueError):
                    pass

        # Remove -arch i386 if 'x86_64' is specified, otherwise incorrect
        # code is generated, at least on OSX
        if hasattr(self.compiler, 'compiler_so'):
            archs = defaultdict(list)
            for i, flag in enumerate(self.compiler.compiler_so[1:]):
                if self.compiler.compiler_so[i] == '-arch':
                    archs[flag].append(i + 1)
            if 'x86_64' in archs and 'i386' in archs:
                for i in archs['i386']:
                    self.compiler.compiler_so[i] = 'x86_64'

        try:
            return super(PythranBuildExtMixIn, self).build_extension(ext)
        finally:
            # Revert compiler settings
            for key in prev.keys():
                set_value(self.compiler, key, prev[key])

            # uninstall hook
            if find_exe is not None:
                import distutils._msvccompiler as msvc
                msvc._find_exe = find_exe


class PythranBuildExtMeta(type):

    def __getitem__(self, base):
        class PythranBuildExt(PythranBuildExtMixIn, base):
            pass

        return PythranBuildExt


class PythranBuildExt(PythranBuildExtMixIn, LegacyBuildExt, metaclass=PythranBuildExtMeta):
    pass


blas_requirements = {"/".join(map(cxxid, elem)) for elem in blas_requires}

includes_matcher = re.compile(r'^#include <pythonic/include/(.*)\.hpp>$',
                              re.MULTILINE)

def requires_blas(source):
    if not os.path.exists(source) or os.path.splitext(source)[1] != ".cpp":
        return False  # conservative

    with open(source) as fd:
        content = fd.read()

    return not blas_requirements.isdisjoint(includes_matcher.findall(content))


class PythranExtension(Extension):
    '''
    Description of a Pythran extension

    Similar to setuptools.extension.Extension except that the sources are .py files
    They must be processable by pythran, of course.

    The compilation process ends up in a native Python module.
    '''

    def __init__(self, name, sources, *args, **kwargs):
        self._kwargs = kwargs.copy()

        if all(not requires_blas(source) for source in sources):
            # Inserting at head so that user-specified config in CLI takes
            # precedence.
            kwargs['config'] = ['compiler.blas=none'] + kwargs.get('config', [])
        cfg_ext = cfg.make_extension(python=True, **kwargs)
        self.cxx = cfg_ext.pop('cxx', None)
        self.cc = cfg_ext.pop('cc', None)
        Extension.__init__(self, name, sources, *args, **cfg_ext)
        self.__dict__.pop("sources", None)

    def _update_blas_requirements(self, source):
        if requires_blas(source):
            cfg_ext = cfg.make_extension(python=True, **self._kwargs)
            for k, v in cfg_ext.items():
                setattr(self, k, v)

    @property
    def sources(self):
        import pythran.toolchain as tc
        cxx_sources = []
        for source in self._sources:
            base, ext = os.path.splitext(source)
            if ext != '.py':
                cxx_sources.append(source)
                continue
            output_file = base + '.cpp'  # target name

            if os.path.exists(source) and (not os.path.exists(output_file)
               or os.path.getmtime(output_file) < os.path.getmtime(source)):
                # get the last name in the path
                if '.' in self.name:
                    module_name = os.path.splitext(self.name)[-1][1:]
                else:
                    module_name = self.name
                tc.compile_pythranfile(source, output_file,
                                       module_name, cpponly=True)

            self._update_blas_requirements(output_file)
            cxx_sources.append(output_file)
        return cxx_sources

    @sources.setter
    def sources(self, sources):
        self._sources = sources
        for source in sources:
            self._update_blas_requirements(source)
