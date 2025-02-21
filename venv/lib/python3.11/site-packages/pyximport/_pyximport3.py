"""
Import hooks; when installed with the install() function, these hooks
allow importing .pyx files as if they were Python modules.

If you want the hook installed every time you run Python
you can add it to your Python version by adding these lines to
sitecustomize.py (which you can create from scratch in site-packages
if it doesn't exist there or somewhere else on your python path)::

    import pyximport
    pyximport.install()

For instance on the Mac with a non-system Python 2.3, you could create
sitecustomize.py with only those two lines at
/usr/local/lib/python2.3/site-packages/sitecustomize.py .

A custom distutils.core.Extension instance and setup() args
(Distribution) for for the build can be defined by a <modulename>.pyxbld
file like:

# examplemod.pyxbld
def make_ext(modname, pyxfilename):
    from distutils.extension import Extension
    return Extension(name = modname,
                     sources=[pyxfilename, 'hello.c'],
                     include_dirs=['/myinclude'] )
def make_setup_args():
    return dict(script_args=["--compiler=mingw32"])

Extra dependencies can be defined by a <modulename>.pyxdep .
See README.

Since Cython 0.11, the :mod:`pyximport` module also has experimental
compilation support for normal Python modules.  This allows you to
automatically run Cython on every .pyx and .py module that Python
imports, including parts of the standard library and installed
packages.  Cython will still fail to compile a lot of Python modules,
in which case the import mechanism will fall back to loading the
Python source modules instead.  The .py import mechanism is installed
like this::

    pyximport.install(pyimport = True)

Running this module as a top-level script will run a test and then print
the documentation.
"""

import glob
import importlib
import os
import sys
from importlib.abc import MetaPathFinder
from importlib.machinery import ExtensionFileLoader, SourceFileLoader
from importlib.util import spec_from_file_location

mod_name = "pyximport"

PY_EXT = ".py"
PYX_EXT = ".pyx"
PYXDEP_EXT = ".pyxdep"
PYXBLD_EXT = ".pyxbld"

DEBUG_IMPORT = False


def _print(message, args):
    if args:
        message = message % args
    print(message)


def _debug(message, *args):
    if DEBUG_IMPORT:
        _print(message, args)


def _info(message, *args):
    _print(message, args)


def load_source(file_path):
    import importlib.util
    from importlib.machinery import SourceFileLoader
    spec = importlib.util.spec_from_file_location("XXXX", file_path, loader=SourceFileLoader("XXXX", file_path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def get_distutils_extension(modname, pyxfilename, language_level=None):
#    try:
#        import hashlib
#    except ImportError:
#        import md5 as hashlib
#    extra = "_" + hashlib.md5(open(pyxfilename).read()).hexdigest()
#    modname = modname + extra
    extension_mod,setup_args = handle_special_build(modname, pyxfilename)
    if not extension_mod:
        if not isinstance(pyxfilename, str):
            # distutils is stupid in Py2 and requires exactly 'str'
            # => encode accidentally coerced unicode strings back to str
            pyxfilename = pyxfilename.encode(sys.getfilesystemencoding())
        from distutils.extension import Extension
        extension_mod = Extension(name = modname, sources=[pyxfilename])
        if language_level is not None:
            extension_mod.cython_directives = {'language_level': language_level}
    return extension_mod,setup_args


def handle_special_build(modname, pyxfilename):
    special_build = os.path.splitext(pyxfilename)[0] + PYXBLD_EXT
    ext = None
    setup_args={}
    if os.path.exists(special_build):
        # globls = {}
        # locs = {}
        # execfile(special_build, globls, locs)
        # ext = locs["make_ext"](modname, pyxfilename)
        mod = load_source(special_build)
        make_ext = getattr(mod,'make_ext',None)
        if make_ext:
            ext = make_ext(modname, pyxfilename)
            assert ext and ext.sources, "make_ext in %s did not return Extension" % special_build
        make_setup_args = getattr(mod, 'make_setup_args',None)
        if make_setup_args:
            setup_args = make_setup_args()
            assert isinstance(setup_args,dict), ("make_setup_args in %s did not return a dict"
                                         % special_build)
        assert ext or setup_args, ("neither make_ext nor make_setup_args %s"
                                         % special_build)
        ext.sources = [os.path.join(os.path.dirname(special_build), source)
                       for source in ext.sources]
    return ext, setup_args


def handle_dependencies(pyxfilename):
    testing = '_test_files' in globals()
    dependfile = os.path.splitext(pyxfilename)[0] + PYXDEP_EXT

    # by default let distutils decide whether to rebuild on its own
    # (it has a better idea of what the output file will be)

    # but we know more about dependencies so force a rebuild if
    # some of the dependencies are newer than the pyxfile.
    if os.path.exists(dependfile):
        with open(dependfile) as fid:
            depends = fid.readlines()
        depends = [depend.strip() for depend in depends]

        # gather dependencies in the "files" variable
        # the dependency file is itself a dependency
        files = [dependfile]
        for depend in depends:
            fullpath = os.path.join(os.path.dirname(dependfile),
                                    depend)
            files.extend(glob.glob(fullpath))

        # only for unit testing to see we did the right thing
        if testing:
            _test_files[:] = []  #$pycheck_no

        # if any file that the pyxfile depends upon is newer than
        # the pyx file, 'touch' the pyx file so that distutils will
        # be tricked into rebuilding it.
        for file in files:
            from distutils.dep_util import newer
            if newer(file, pyxfilename):
                _debug("Rebuilding %s because of %s", pyxfilename, file)
                filetime = os.path.getmtime(file)
                os.utime(pyxfilename, (filetime, filetime))
                if testing:
                    _test_files.append(file)


def build_module(name, pyxfilename, pyxbuild_dir=None, inplace=False, language_level=None):
    assert os.path.exists(pyxfilename), "Path does not exist: %s" % pyxfilename
    handle_dependencies(pyxfilename)

    extension_mod, setup_args = get_distutils_extension(name, pyxfilename, language_level)
    build_in_temp = pyxargs.build_in_temp
    sargs = pyxargs.setup_args.copy()
    sargs.update(setup_args)
    build_in_temp = sargs.pop('build_in_temp',build_in_temp)

    from . import pyxbuild
    olddir = os.getcwd()
    common = ''
    if pyxbuild_dir and sys.platform == 'win32':
        # Windows concatenates the pyxbuild_dir to the pyxfilename when
        # compiling, and then complains that the filename is too long
        common = os.path.commonprefix([pyxbuild_dir, pyxfilename])
    if len(common) > 30:
        pyxfilename = os.path.relpath(pyxfilename, common)
        pyxbuild_dir = os.path.relpath(pyxbuild_dir, common)
        os.chdir(common)
    try:
        so_path = pyxbuild.pyx_to_dll(pyxfilename, extension_mod,
                                      build_in_temp=build_in_temp,
                                      pyxbuild_dir=pyxbuild_dir,
                                      setup_args=sargs,
                                      inplace=inplace,
                                      reload_support=pyxargs.reload_support)
    finally:
        os.chdir(olddir)
    so_path = os.path.join(common, so_path)
    assert os.path.exists(so_path), "Cannot find: %s" % so_path

    junkpath = os.path.join(os.path.dirname(so_path), name+"_*")  #very dangerous with --inplace ? yes, indeed, trying to eat my files ;)
    junkstuff = glob.glob(junkpath)
    for path in junkstuff:
        if path != so_path:
            try:
                os.remove(path)
            except IOError:
                _info("Couldn't remove %s", path)

    return so_path


# import hooks

class PyxImportMetaFinder(MetaPathFinder):

    def __init__(self, extension=PYX_EXT, pyxbuild_dir=None, inplace=False, language_level=None):
        self.pyxbuild_dir = pyxbuild_dir
        self.inplace = inplace
        self.language_level = language_level
        self.extension = extension

    def find_spec(self, fullname, path, target=None):
        if not path:
            path = [os.getcwd()]  # top level import --
        if "." in fullname:
            *parents, name = fullname.split(".")
        else:
            name = fullname
        for entry in path:
            if os.path.isdir(os.path.join(entry, name)):
                # this module has child modules
                filename = os.path.join(entry, name, "__init__" + self.extension)
                submodule_locations = [os.path.join(entry, name)]
            else:
                filename = os.path.join(entry, name + self.extension)
                submodule_locations = None
            if not os.path.exists(filename):
                continue

            return spec_from_file_location(
                fullname, filename,
                loader=PyxImportLoader(filename, self.pyxbuild_dir, self.inplace, self.language_level),
                submodule_search_locations=submodule_locations)

        return None  # we don't know how to import this


class PyImportMetaFinder(MetaPathFinder):

    def __init__(self, extension=PY_EXT, pyxbuild_dir=None, inplace=False, language_level=None):
        self.pyxbuild_dir = pyxbuild_dir
        self.inplace = inplace
        self.language_level = language_level
        self.extension = extension
        self.uncompilable_modules = {}
        self.blocked_modules = ['Cython', 'pyxbuild', 'pyximport.pyxbuild',
                                'distutils', 'cython']
        self.blocked_packages = ['Cython.', 'distutils.']

    def find_spec(self, fullname, path, target=None):
        if fullname in sys.modules:
            return None
        if any([fullname.startswith(pkg) for pkg in self.blocked_packages]):
            return None
        if fullname in self.blocked_modules:
            # prevent infinite recursion
            return None

        self.blocked_modules.append(fullname)
        name = fullname
        if not path:
            path = [os.getcwd()]  # top level import --
        try:
            for entry in path:
                if os.path.isdir(os.path.join(entry, name)):
                    # this module has child modules
                    filename = os.path.join(entry, name, "__init__" + self.extension)
                    submodule_locations = [os.path.join(entry, name)]
                else:
                    filename = os.path.join(entry, name + self.extension)
                    submodule_locations = None
                if not os.path.exists(filename):
                    continue

                return spec_from_file_location(
                    fullname, filename,
                    loader=PyxImportLoader(filename, self.pyxbuild_dir, self.inplace, self.language_level),
                    submodule_search_locations=submodule_locations)
        finally:
            self.blocked_modules.pop()

        return None  # we don't know how to import this


class PyxImportLoader(ExtensionFileLoader):

    def __init__(self, filename, pyxbuild_dir, inplace, language_level):
        module_name = os.path.splitext(os.path.basename(filename))[0]
        super().__init__(module_name, filename)
        self._pyxbuild_dir = pyxbuild_dir
        self._inplace = inplace
        self._language_level = language_level

    def create_module(self, spec):
        try:
            so_path = build_module(spec.name, pyxfilename=spec.origin, pyxbuild_dir=self._pyxbuild_dir,
                                   inplace=self._inplace, language_level=self._language_level)
            self.path = so_path
            spec.origin = so_path
            return super().create_module(spec)
        except Exception as failure_exc:
            _debug("Failed to load extension module: %r" % failure_exc)
            if pyxargs.load_py_module_on_import_failure and spec.origin.endswith(PY_EXT):
                spec = importlib.util.spec_from_file_location(spec.name, spec.origin,
                                                              loader=SourceFileLoader(spec.name, spec.origin))
                mod = importlib.util.module_from_spec(spec)
                assert mod.__file__ in (spec.origin, spec.origin + 'c', spec.origin + 'o'), (mod.__file__, spec.origin)
                return mod
            else:
                tb = sys.exc_info()[2]
                import traceback
                exc = ImportError("Building module %s failed: %s" % (
                    spec.name, traceback.format_exception_only(*sys.exc_info()[:2])))
                raise exc.with_traceback(tb)

    def exec_module(self, module):
        try:
            return super().exec_module(module)
        except Exception as failure_exc:
            import traceback
            _debug("Failed to load extension module: %r" % failure_exc)
            raise ImportError("Executing module %s failed %s" % (
                    module.__file__, traceback.format_exception_only(*sys.exc_info()[:2])))


#install args
class PyxArgs(object):
    build_dir=True
    build_in_temp=True
    setup_args={}   #None


def _have_importers():
    has_py_importer = False
    has_pyx_importer = False
    for importer in sys.meta_path:
        if isinstance(importer, PyxImportMetaFinder):
            if isinstance(importer, PyImportMetaFinder):
                has_py_importer = True
            else:
                has_pyx_importer = True

    return has_py_importer, has_pyx_importer


def install(pyximport=True, pyimport=False, build_dir=None, build_in_temp=True,
            setup_args=None, reload_support=False,
            load_py_module_on_import_failure=False, inplace=False,
            language_level=None):
    """ Main entry point for pyxinstall.

    Call this to install the ``.pyx`` import hook in
    your meta-path for a single Python process.  If you want it to be
    installed whenever you use Python, add it to your ``sitecustomize``
    (as described above).

    :param pyximport: If set to False, does not try to import ``.pyx`` files.

    :param pyimport: You can pass ``pyimport=True`` to also
        install the ``.py`` import hook
        in your meta-path.  Note, however, that it is rather experimental,
        will not work at all for some ``.py`` files and packages, and will
        heavily slow down your imports due to search and compilation.
        Use at your own risk.

    :param build_dir: By default, compiled modules will end up in a ``.pyxbld``
        directory in the user's home directory.  Passing a different path
        as ``build_dir`` will override this.

    :param build_in_temp: If ``False``, will produce the C files locally. Working
        with complex dependencies and debugging becomes more easy. This
        can principally interfere with existing files of the same name.

    :param setup_args: Dict of arguments for Distribution.
        See ``distutils.core.setup()``.

    :param reload_support: Enables support for dynamic
        ``reload(my_module)``, e.g. after a change in the Cython code.
        Additional files ``<so_path>.reloadNN`` may arise on that account, when
        the previously loaded module file cannot be overwritten.

    :param load_py_module_on_import_failure: If the compilation of a ``.py``
        file succeeds, but the subsequent import fails for some reason,
        retry the import with the normal ``.py`` module instead of the
        compiled module.  Note that this may lead to unpredictable results
        for modules that change the system state during their import, as
        the second import will rerun these modifications in whatever state
        the system was left after the import of the compiled module
        failed.

    :param inplace: Install the compiled module
        (``.so`` for Linux and Mac / ``.pyd`` for Windows)
        next to the source file.

    :param language_level: The source language level to use: 2 or 3.
        The default is to use the language level of the current Python
        runtime for .py files and Py2 for ``.pyx`` files.
    """
    if setup_args is None:
        setup_args = {}
    if not build_dir:
        build_dir = os.path.join(os.path.expanduser('~'), '.pyxbld')

    global pyxargs
    pyxargs = PyxArgs()  #$pycheck_no
    pyxargs.build_dir = build_dir
    pyxargs.build_in_temp = build_in_temp
    pyxargs.setup_args = (setup_args or {}).copy()
    pyxargs.reload_support = reload_support
    pyxargs.load_py_module_on_import_failure = load_py_module_on_import_failure

    has_py_importer, has_pyx_importer = _have_importers()
    py_importer, pyx_importer = None, None

    if pyimport and not has_py_importer:
        py_importer = PyImportMetaFinder(pyxbuild_dir=build_dir, inplace=inplace,
                                         language_level=language_level)
        # make sure we import Cython before we install the import hook
        import Cython.Compiler.Main, Cython.Compiler.Pipeline, Cython.Compiler.Optimize
        sys.meta_path.insert(0, py_importer)

    if pyximport and not has_pyx_importer:
        pyx_importer = PyxImportMetaFinder(pyxbuild_dir=build_dir, inplace=inplace,
                                           language_level=language_level)
        sys.meta_path.append(pyx_importer)

    return py_importer, pyx_importer


def uninstall(py_importer, pyx_importer):
    """
    Uninstall an import hook.
    """
    try:
        sys.meta_path.remove(py_importer)
    except ValueError:
        pass

    try:
        sys.meta_path.remove(pyx_importer)
    except ValueError:
        pass


# MAIN

def show_docs():
    import __main__
    __main__.__name__ = mod_name
    for name in dir(__main__):
        item = getattr(__main__, name)
        try:
            setattr(item, "__module__", mod_name)
        except (AttributeError, TypeError):
            pass
    help(__main__)


if __name__ == '__main__':
    show_docs()
