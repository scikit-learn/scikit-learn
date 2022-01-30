"""Low-level infrastructure to find modules.

This builds on fscache.py; find_sources.py builds on top of this.
"""

import ast
import collections
import functools
import os
import re
import subprocess
import sys
from enum import Enum, unique

from typing import Dict, Iterator, List, NamedTuple, Optional, Set, Tuple, Union
from typing_extensions import Final, TypeAlias as _TypeAlias

from mypy.fscache import FileSystemCache
from mypy.options import Options
from mypy.stubinfo import is_legacy_bundled_package
from mypy import pyinfo

# Paths to be searched in find_module().
SearchPaths = NamedTuple(
    'SearchPaths',
    [('python_path', Tuple[str, ...]),  # where user code is found
     ('mypy_path', Tuple[str, ...]),  # from $MYPYPATH or config variable
     ('package_path', Tuple[str, ...]),  # from get_site_packages_dirs()
     ('typeshed_path', Tuple[str, ...]),  # paths in typeshed
     ])

# Package dirs are a two-tuple of path to search and whether to verify the module
OnePackageDir = Tuple[str, bool]
PackageDirs = List[OnePackageDir]

# Minimum and maximum Python versions for modules in stdlib as (major, minor)
StdlibVersions: _TypeAlias = Dict[str, Tuple[Tuple[int, int], Optional[Tuple[int, int]]]]

PYTHON_EXTENSIONS: Final = [".pyi", ".py"]

PYTHON2_STUB_DIR: Final = "@python2"


# TODO: Consider adding more reasons here?
# E.g. if we deduce a module would likely be found if the user were
# to set the --namespace-packages flag.
@unique
class ModuleNotFoundReason(Enum):
    # The module was not found: we found neither stubs nor a plausible code
    # implementation (with or without a py.typed file).
    NOT_FOUND = 0

    # The implementation for this module plausibly exists (e.g. we
    # found a matching folder or *.py file), but either the parent package
    # did not contain a py.typed file or we were unable to find a
    # corresponding *-stubs package.
    FOUND_WITHOUT_TYPE_HINTS = 1

    # The module was not found in the current working directory, but
    # was able to be found in the parent directory.
    WRONG_WORKING_DIRECTORY = 2

    # Stub PyPI package (typically types-pkgname) known to exist but not installed.
    APPROVED_STUBS_NOT_INSTALLED = 3

    def error_message_templates(self, daemon: bool) -> Tuple[str, List[str]]:
        doc_link = "See https://mypy.readthedocs.io/en/stable/running_mypy.html#missing-imports"
        if self is ModuleNotFoundReason.NOT_FOUND:
            msg = 'Cannot find implementation or library stub for module named "{module}"'
            notes = [doc_link]
        elif self is ModuleNotFoundReason.WRONG_WORKING_DIRECTORY:
            msg = 'Cannot find implementation or library stub for module named "{module}"'
            notes = ["You may be running mypy in a subpackage, "
                     "mypy should be run on the package root"]
        elif self is ModuleNotFoundReason.FOUND_WITHOUT_TYPE_HINTS:
            msg = (
                'Skipping analyzing "{module}": module is installed, but missing library stubs '
                'or py.typed marker'
            )
            notes = [doc_link]
        elif self is ModuleNotFoundReason.APPROVED_STUBS_NOT_INSTALLED:
            msg = (
                'Library stubs not installed for "{module}" (or incompatible with Python {pyver})'
            )
            notes = ['Hint: "python3 -m pip install {stub_dist}"']
            if not daemon:
                notes.append(
                    '(or run "mypy --install-types" to install all missing stub packages)')
            notes.append(doc_link)
        else:
            assert False
        return msg, notes


# If we found the module, returns the path to the module as a str.
# Otherwise, returns the reason why the module wasn't found.
ModuleSearchResult = Union[str, ModuleNotFoundReason]


class BuildSource:
    """A single source file."""

    def __init__(self, path: Optional[str], module: Optional[str],
                 text: Optional[str] = None, base_dir: Optional[str] = None) -> None:
        self.path = path  # File where it's found (e.g. 'xxx/yyy/foo/bar.py')
        self.module = module or '__main__'  # Module name (e.g. 'foo.bar')
        self.text = text  # Source code, if initially supplied, else None
        self.base_dir = base_dir  # Directory where the package is rooted (e.g. 'xxx/yyy')

    def __repr__(self) -> str:
        return 'BuildSource(path=%r, module=%r, has_text=%s, base_dir=%r)' % (
            self.path,
            self.module,
            self.text is not None,
            self.base_dir)


class FindModuleCache:
    """Module finder with integrated cache.

    Module locations and some intermediate results are cached internally
    and can be cleared with the clear() method.

    All file system accesses are performed through a FileSystemCache,
    which is not ever cleared by this class. If necessary it must be
    cleared by client code.
    """

    def __init__(self,
                 search_paths: SearchPaths,
                 fscache: Optional[FileSystemCache],
                 options: Optional[Options],
                 stdlib_py_versions: Optional[StdlibVersions] = None) -> None:
        self.search_paths = search_paths
        self.fscache = fscache or FileSystemCache()
        # Cache for get_toplevel_possibilities:
        # search_paths -> (toplevel_id -> list(package_dirs))
        self.initial_components: Dict[Tuple[str, ...], Dict[str, List[str]]] = {}
        # Cache find_module: id -> result
        self.results: Dict[str, ModuleSearchResult] = {}
        self.ns_ancestors: Dict[str, str] = {}
        self.options = options
        custom_typeshed_dir = None
        if options:
            custom_typeshed_dir = options.custom_typeshed_dir
        self.stdlib_py_versions = (
            stdlib_py_versions or load_stdlib_py_versions(custom_typeshed_dir)
        )
        self.python_major_ver = 3 if options is None else options.python_version[0]

    def clear(self) -> None:
        self.results.clear()
        self.initial_components.clear()
        self.ns_ancestors.clear()

    def find_lib_path_dirs(self, id: str, lib_path: Tuple[str, ...]) -> PackageDirs:
        """Find which elements of a lib_path have the directory a module needs to exist.

        This is run for the python_path, mypy_path, and typeshed_path search paths.
        """
        components = id.split('.')
        dir_chain = os.sep.join(components[:-1])  # e.g., 'foo/bar'

        dirs = []
        for pathitem in self.get_toplevel_possibilities(lib_path, components[0]):
            # e.g., '/usr/lib/python3.4/foo/bar'
            dir = os.path.normpath(os.path.join(pathitem, dir_chain))
            if self.fscache.isdir(dir):
                dirs.append((dir, True))
        return dirs

    def get_toplevel_possibilities(self, lib_path: Tuple[str, ...], id: str) -> List[str]:
        """Find which elements of lib_path could contain a particular top-level module.

        In practice, almost all modules can be routed to the correct entry in
        lib_path by looking at just the first component of the module name.

        We take advantage of this by enumerating the contents of all of the
        directories on the lib_path and building a map of which entries in
        the lib_path could contain each potential top-level module that appears.
        """

        if lib_path in self.initial_components:
            return self.initial_components[lib_path].get(id, [])

        # Enumerate all the files in the directories on lib_path and produce the map
        components: Dict[str, List[str]] = {}
        for dir in lib_path:
            try:
                contents = self.fscache.listdir(dir)
            except OSError:
                contents = []
            # False positives are fine for correctness here, since we will check
            # precisely later, so we only look at the root of every filename without
            # any concern for the exact details.
            for name in contents:
                name = os.path.splitext(name)[0]
                components.setdefault(name, []).append(dir)

        if self.python_major_ver == 2:
            components = {id: filter_redundant_py2_dirs(dirs)
                          for id, dirs in components.items()}

        self.initial_components[lib_path] = components
        return components.get(id, [])

    def find_module(self, id: str, *, fast_path: bool = False) -> ModuleSearchResult:
        """Return the path of the module source file or why it wasn't found.

        If fast_path is True, prioritize performance over generating detailed
        error descriptions.
        """
        if id not in self.results:
            top_level = id.partition('.')[0]
            use_typeshed = True
            if id in self.stdlib_py_versions:
                use_typeshed = self._typeshed_has_version(id)
            elif top_level in self.stdlib_py_versions:
                use_typeshed = self._typeshed_has_version(top_level)
            self.results[id] = self._find_module(id, use_typeshed)
            if (not fast_path
                    and self.results[id] is ModuleNotFoundReason.NOT_FOUND
                    and self._can_find_module_in_parent_dir(id)):
                self.results[id] = ModuleNotFoundReason.WRONG_WORKING_DIRECTORY
        return self.results[id]

    def _typeshed_has_version(self, module: str) -> bool:
        if not self.options:
            return True
        version = typeshed_py_version(self.options)
        min_version, max_version = self.stdlib_py_versions[module]
        return version >= min_version and (max_version is None or version <= max_version)

    def _find_module_non_stub_helper(self, components: List[str],
                                     pkg_dir: str) -> Union[OnePackageDir, ModuleNotFoundReason]:
        plausible_match = False
        dir_path = pkg_dir
        for index, component in enumerate(components):
            dir_path = os.path.join(dir_path, component)
            if self.fscache.isfile(os.path.join(dir_path, 'py.typed')):
                return os.path.join(pkg_dir, *components[:-1]), index == 0
            elif not plausible_match and (self.fscache.isdir(dir_path)
                                          or self.fscache.isfile(dir_path + ".py")):
                plausible_match = True
        if is_legacy_bundled_package(components[0], self.python_major_ver):
            if (len(components) == 1
                    or (self.find_module(components[0]) is
                        ModuleNotFoundReason.APPROVED_STUBS_NOT_INSTALLED)):
                return ModuleNotFoundReason.APPROVED_STUBS_NOT_INSTALLED
        if is_legacy_bundled_package('.'.join(components[:2]), self.python_major_ver):
            return ModuleNotFoundReason.APPROVED_STUBS_NOT_INSTALLED
        if plausible_match:
            return ModuleNotFoundReason.FOUND_WITHOUT_TYPE_HINTS
        else:
            return ModuleNotFoundReason.NOT_FOUND

    def _update_ns_ancestors(self, components: List[str], match: Tuple[str, bool]) -> None:
        path, verify = match
        for i in range(1, len(components)):
            pkg_id = '.'.join(components[:-i])
            if pkg_id not in self.ns_ancestors and self.fscache.isdir(path):
                self.ns_ancestors[pkg_id] = path
            path = os.path.dirname(path)

    def _can_find_module_in_parent_dir(self, id: str) -> bool:
        """Test if a module can be found by checking the parent directories
        of the current working directory.
        """
        working_dir = os.getcwd()
        parent_search = FindModuleCache(
            SearchPaths((), (), (), ()),
            self.fscache,
            self.options,
            stdlib_py_versions=self.stdlib_py_versions
        )
        while any(file.endswith(("__init__.py", "__init__.pyi"))
                  for file in os.listdir(working_dir)):
            working_dir = os.path.dirname(working_dir)
            parent_search.search_paths = SearchPaths((working_dir,), (), (), ())
            if not isinstance(parent_search._find_module(id, False), ModuleNotFoundReason):
                return True
        return False

    def _find_module(self, id: str, use_typeshed: bool) -> ModuleSearchResult:
        fscache = self.fscache

        # If we're looking for a module like 'foo.bar.baz', it's likely that most of the
        # many elements of lib_path don't even have a subdirectory 'foo/bar'.  Discover
        # that only once and cache it for when we look for modules like 'foo.bar.blah'
        # that will require the same subdirectory.
        components = id.split('.')
        dir_chain = os.sep.join(components[:-1])  # e.g., 'foo/bar'

        # We have two sets of folders so that we collect *all* stubs folders and
        # put them in the front of the search path
        third_party_inline_dirs: PackageDirs = []
        third_party_stubs_dirs: PackageDirs = []
        found_possible_third_party_missing_type_hints = False
        need_installed_stubs = False
        # Third-party stub/typed packages
        for pkg_dir in self.search_paths.package_path:
            stub_name = components[0] + '-stubs'
            stub_dir = os.path.join(pkg_dir, stub_name)
            if self.python_major_ver == 2:
                alt_stub_name = components[0] + '-python2-stubs'
                alt_stub_dir = os.path.join(pkg_dir, alt_stub_name)
                if fscache.isdir(alt_stub_dir):
                    stub_name = alt_stub_name
                    stub_dir = alt_stub_dir
            if fscache.isdir(stub_dir) and self._is_compatible_stub_package(stub_dir):
                stub_typed_file = os.path.join(stub_dir, 'py.typed')
                stub_components = [stub_name] + components[1:]
                path = os.path.join(pkg_dir, *stub_components[:-1])
                if fscache.isdir(path):
                    if fscache.isfile(stub_typed_file):
                        # Stub packages can have a py.typed file, which must include
                        # 'partial\n' to make the package partial
                        # Partial here means that mypy should look at the runtime
                        # package if installed.
                        if fscache.read(stub_typed_file).decode().strip() == 'partial':
                            runtime_path = os.path.join(pkg_dir, dir_chain)
                            third_party_inline_dirs.append((runtime_path, True))
                            # if the package is partial, we don't verify the module, as
                            # the partial stub package may not have a __init__.pyi
                            third_party_stubs_dirs.append((path, False))
                        else:
                            # handle the edge case where people put a py.typed file
                            # in a stub package, but it isn't partial
                            third_party_stubs_dirs.append((path, True))
                    else:
                        third_party_stubs_dirs.append((path, True))
            non_stub_match = self._find_module_non_stub_helper(components, pkg_dir)
            if isinstance(non_stub_match, ModuleNotFoundReason):
                if non_stub_match is ModuleNotFoundReason.FOUND_WITHOUT_TYPE_HINTS:
                    found_possible_third_party_missing_type_hints = True
                elif non_stub_match is ModuleNotFoundReason.APPROVED_STUBS_NOT_INSTALLED:
                    need_installed_stubs = True
            else:
                third_party_inline_dirs.append(non_stub_match)
                self._update_ns_ancestors(components, non_stub_match)
        if self.options and self.options.use_builtins_fixtures:
            # Everything should be in fixtures.
            third_party_inline_dirs.clear()
            third_party_stubs_dirs.clear()
            found_possible_third_party_missing_type_hints = False
        python_mypy_path = self.search_paths.mypy_path + self.search_paths.python_path
        candidate_base_dirs = self.find_lib_path_dirs(id, python_mypy_path)
        if use_typeshed:
            # Search for stdlib stubs in typeshed before installed
            # stubs to avoid picking up backports (dataclasses, for
            # example) when the library is included in stdlib.
            candidate_base_dirs += self.find_lib_path_dirs(id, self.search_paths.typeshed_path)
        candidate_base_dirs += third_party_stubs_dirs + third_party_inline_dirs

        # If we're looking for a module like 'foo.bar.baz', then candidate_base_dirs now
        # contains just the subdirectories 'foo/bar' that actually exist under the
        # elements of lib_path.  This is probably much shorter than lib_path itself.
        # Now just look for 'baz.pyi', 'baz/__init__.py', etc., inside those directories.
        seplast = os.sep + components[-1]  # so e.g. '/baz'
        sepinit = os.sep + '__init__'
        near_misses = []  # Collect near misses for namespace mode (see below).
        for base_dir, verify in candidate_base_dirs:
            base_path = base_dir + seplast  # so e.g. '/usr/lib/python3.4/foo/bar/baz'
            has_init = False
            dir_prefix = base_dir
            for _ in range(len(components) - 1):
                dir_prefix = os.path.dirname(dir_prefix)
            # Prefer package over module, i.e. baz/__init__.py* over baz.py*.
            for extension in PYTHON_EXTENSIONS:
                path = base_path + sepinit + extension
                suffix = '-stubs'
                if self.python_major_ver == 2:
                    if os.path.isdir(base_path + '-python2-stubs'):
                        suffix = '-python2-stubs'
                path_stubs = base_path + suffix + sepinit + extension
                if fscache.isfile_case(path, dir_prefix):
                    has_init = True
                    if verify and not verify_module(fscache, id, path, dir_prefix):
                        near_misses.append((path, dir_prefix))
                        continue
                    return path
                elif fscache.isfile_case(path_stubs, dir_prefix):
                    if verify and not verify_module(fscache, id, path_stubs, dir_prefix):
                        near_misses.append((path_stubs, dir_prefix))
                        continue
                    return path_stubs

            # In namespace mode, register a potential namespace package
            if self.options and self.options.namespace_packages:
                if fscache.exists_case(base_path, dir_prefix) and not has_init:
                    near_misses.append((base_path, dir_prefix))

            # No package, look for module.
            for extension in PYTHON_EXTENSIONS:
                path = base_path + extension
                if fscache.isfile_case(path, dir_prefix):
                    if verify and not verify_module(fscache, id, path, dir_prefix):
                        near_misses.append((path, dir_prefix))
                        continue
                    return path

        # In namespace mode, re-check those entries that had 'verify'.
        # Assume search path entries xxx, yyy and zzz, and we're
        # looking for foo.bar.baz.  Suppose near_misses has:
        #
        # - xxx/foo/bar/baz.py
        # - yyy/foo/bar/baz/__init__.py
        # - zzz/foo/bar/baz.pyi
        #
        # If any of the foo directories has __init__.py[i], it wins.
        # Else, we look for foo/bar/__init__.py[i], etc.  If there are
        # none, the first hit wins.  Note that this does not take into
        # account whether the lowest-level module is a file (baz.py),
        # a package (baz/__init__.py), or a stub file (baz.pyi) -- for
        # these the first one encountered along the search path wins.
        #
        # The helper function highest_init_level() returns an int that
        # indicates the highest level at which a __init__.py[i] file
        # is found; if no __init__ was found it returns 0, if we find
        # only foo/bar/__init__.py it returns 1, and if we have
        # foo/__init__.py it returns 2 (regardless of what's in
        # foo/bar).  It doesn't look higher than that.
        if self.options and self.options.namespace_packages and near_misses:
            levels = [highest_init_level(fscache, id, path, dir_prefix)
                      for path, dir_prefix in near_misses]
            index = levels.index(max(levels))
            return near_misses[index][0]

        # Finally, we may be asked to produce an ancestor for an
        # installed package with a py.typed marker that is a
        # subpackage of a namespace package.  We only fess up to these
        # if we would otherwise return "not found".
        ancestor = self.ns_ancestors.get(id)
        if ancestor is not None:
            return ancestor

        if need_installed_stubs:
            return ModuleNotFoundReason.APPROVED_STUBS_NOT_INSTALLED
        elif found_possible_third_party_missing_type_hints:
            return ModuleNotFoundReason.FOUND_WITHOUT_TYPE_HINTS
        else:
            return ModuleNotFoundReason.NOT_FOUND

    def _is_compatible_stub_package(self, stub_dir: str) -> bool:
        """Does a stub package support the target Python version?

        Stub packages may contain a metadata file which specifies
        whether the stubs are compatible with Python 2 and 3.
        """
        metadata_fnam = os.path.join(stub_dir, 'METADATA.toml')
        if os.path.isfile(metadata_fnam):
            # Delay import for a possible minor performance win.
            import tomli
            with open(metadata_fnam, encoding="utf-8") as f:
                metadata = tomli.loads(f.read())
            if self.python_major_ver == 2:
                return bool(metadata.get('python2', False))
            else:
                return bool(metadata.get('python3', True))
        return True

    def find_modules_recursive(self, module: str) -> List[BuildSource]:
        module_path = self.find_module(module)
        if isinstance(module_path, ModuleNotFoundReason):
            return []
        sources = [BuildSource(module_path, module, None)]

        package_path = None
        if module_path.endswith(('__init__.py', '__init__.pyi')):
            package_path = os.path.dirname(module_path)
        elif self.fscache.isdir(module_path):
            package_path = module_path
        if package_path is None:
            return sources

        # This logic closely mirrors that in find_sources. One small but important difference is
        # that we do not sort names with keyfunc. The recursive call to find_modules_recursive
        # calls find_module, which will handle the preference between packages, pyi and py.
        # Another difference is it doesn't handle nested search paths / package roots.

        seen: Set[str] = set()
        names = sorted(self.fscache.listdir(package_path))
        for name in names:
            # Skip certain names altogether
            if name in ("__pycache__", "site-packages", "node_modules") or name.startswith("."):
                continue
            subpath = os.path.join(package_path, name)

            if self.options and matches_exclude(
                subpath, self.options.exclude, self.fscache, self.options.verbosity >= 2
            ):
                continue

            if self.fscache.isdir(subpath):
                # Only recurse into packages
                if (self.options and self.options.namespace_packages) or (
                    self.fscache.isfile(os.path.join(subpath, "__init__.py"))
                    or self.fscache.isfile(os.path.join(subpath, "__init__.pyi"))
                ):
                    seen.add(name)
                    sources.extend(self.find_modules_recursive(module + '.' + name))
            else:
                stem, suffix = os.path.splitext(name)
                if stem == '__init__':
                    continue
                if stem not in seen and '.' not in stem and suffix in PYTHON_EXTENSIONS:
                    # (If we sorted names by keyfunc) we could probably just make the BuildSource
                    # ourselves, but this ensures compatibility with find_module / the cache
                    seen.add(stem)
                    sources.extend(self.find_modules_recursive(module + '.' + stem))
        return sources


def matches_exclude(subpath: str,
                    excludes: List[str],
                    fscache: FileSystemCache,
                    verbose: bool) -> bool:
    if not excludes:
        return False
    subpath_str = os.path.relpath(subpath).replace(os.sep, "/")
    if fscache.isdir(subpath):
        subpath_str += "/"
    for exclude in excludes:
        if re.search(exclude, subpath_str):
            if verbose:
                print("TRACE: Excluding {} (matches pattern {})".format(subpath_str, exclude),
                      file=sys.stderr)
            return True
    return False


def verify_module(fscache: FileSystemCache, id: str, path: str, prefix: str) -> bool:
    """Check that all packages containing id have a __init__ file."""
    if path.endswith(('__init__.py', '__init__.pyi')):
        path = os.path.dirname(path)
    for i in range(id.count('.')):
        path = os.path.dirname(path)
        if not any(fscache.isfile_case(os.path.join(path, '__init__{}'.format(extension)),
                                       prefix)
                   for extension in PYTHON_EXTENSIONS):
            return False
    return True


def highest_init_level(fscache: FileSystemCache, id: str, path: str, prefix: str) -> int:
    """Compute the highest level where an __init__ file is found."""
    if path.endswith(('__init__.py', '__init__.pyi')):
        path = os.path.dirname(path)
    level = 0
    for i in range(id.count('.')):
        path = os.path.dirname(path)
        if any(fscache.isfile_case(os.path.join(path, '__init__{}'.format(extension)),
                                   prefix)
               for extension in PYTHON_EXTENSIONS):
            level = i + 1
    return level


def mypy_path() -> List[str]:
    path_env = os.getenv('MYPYPATH')
    if not path_env:
        return []
    return path_env.split(os.pathsep)


def default_lib_path(data_dir: str,
                     pyversion: Tuple[int, int],
                     custom_typeshed_dir: Optional[str]) -> List[str]:
    """Return default standard library search paths."""
    path: List[str] = []

    if custom_typeshed_dir:
        typeshed_dir = os.path.join(custom_typeshed_dir, "stdlib")
        mypy_extensions_dir = os.path.join(custom_typeshed_dir, "stubs", "mypy-extensions")
        versions_file = os.path.join(typeshed_dir, "VERSIONS")
        if not os.path.isdir(typeshed_dir) or not os.path.isfile(versions_file):
            print("error: --custom-typeshed-dir does not point to a valid typeshed ({})".format(
                custom_typeshed_dir))
            sys.exit(2)
    else:
        auto = os.path.join(data_dir, 'stubs-auto')
        if os.path.isdir(auto):
            data_dir = auto
        typeshed_dir = os.path.join(data_dir, "typeshed", "stdlib")
        mypy_extensions_dir = os.path.join(data_dir, "typeshed", "stubs", "mypy-extensions")
    if pyversion[0] == 2:
        # Python 2 variants of certain stdlib modules are in a separate directory.
        python2_dir = os.path.join(typeshed_dir, PYTHON2_STUB_DIR)
        path.append(python2_dir)
    path.append(typeshed_dir)

    # Get mypy-extensions stubs from typeshed, since we treat it as an
    # "internal" library, similar to typing and typing-extensions.
    path.append(mypy_extensions_dir)

    # Add fallback path that can be used if we have a broken installation.
    if sys.platform != 'win32':
        path.append('/usr/local/lib/mypy')
    if not path:
        print("Could not resolve typeshed subdirectories. Your mypy install is broken.\n"
              "Python executable is located at {0}.\nMypy located at {1}".format(
                  sys.executable, data_dir), file=sys.stderr)
        sys.exit(1)
    return path


@functools.lru_cache(maxsize=None)
def get_prefixes(python_executable: Optional[str]) -> Tuple[str, str]:
    """Get the sys.base_prefix and sys.prefix for the given python.

    This runs a subprocess call to get the prefix paths of the given Python executable.
    To avoid repeatedly calling a subprocess (which can be slow!) we
    lru_cache the results.
    """
    if python_executable is None:
        return '', ''
    elif python_executable == sys.executable:
        # Use running Python's package dirs
        return pyinfo.getprefixes()
    else:
        # Use subprocess to get the package directory of given Python
        # executable
        return ast.literal_eval(
            subprocess.check_output([python_executable, pyinfo.__file__, 'getprefixes'],
            stderr=subprocess.PIPE).decode())


@functools.lru_cache(maxsize=None)
def get_site_packages_dirs(python_executable: Optional[str]) -> Tuple[List[str], List[str]]:
    """Find package directories for given python.

    This runs a subprocess call, which generates a list of the egg directories, and the site
    package directories. To avoid repeatedly calling a subprocess (which can be slow!) we
    lru_cache the results.
    """

    if python_executable is None:
        return [], []
    elif python_executable == sys.executable:
        # Use running Python's package dirs
        site_packages = pyinfo.getsitepackages()
    else:
        # Use subprocess to get the package directory of given Python
        # executable
        site_packages = ast.literal_eval(
            subprocess.check_output([python_executable, pyinfo.__file__, 'getsitepackages'],
            stderr=subprocess.PIPE).decode())
    return expand_site_packages(site_packages)


def expand_site_packages(site_packages: List[str]) -> Tuple[List[str], List[str]]:
    """Expands .pth imports in site-packages directories"""
    egg_dirs: List[str] = []
    for dir in site_packages:
        if not os.path.isdir(dir):
            continue
        pth_filenames = sorted(name for name in os.listdir(dir) if name.endswith(".pth"))
        for pth_filename in pth_filenames:
            egg_dirs.extend(_parse_pth_file(dir, pth_filename))

    return egg_dirs, site_packages


def _parse_pth_file(dir: str, pth_filename: str) -> Iterator[str]:
    """
    Mimics a subset of .pth import hook from Lib/site.py
    See https://github.com/python/cpython/blob/3.5/Lib/site.py#L146-L185
    """

    pth_file = os.path.join(dir, pth_filename)
    try:
        f = open(pth_file, "r")
    except OSError:
        return
    with f:
        for line in f.readlines():
            if line.startswith("#"):
                # Skip comment lines
                continue
            if line.startswith(("import ", "import\t")):
                # import statements in .pth files are not supported
                continue

            yield _make_abspath(line.rstrip(), dir)


def _make_abspath(path: str, root: str) -> str:
    """Take a path and make it absolute relative to root if not already absolute."""
    if os.path.isabs(path):
        return os.path.normpath(path)
    else:
        return os.path.join(root, os.path.normpath(path))


def add_py2_mypypath_entries(mypypath: List[str]) -> List[str]:
    """Add corresponding @python2 subdirectories to mypypath.

    For each path entry 'x', add 'x/@python2' before 'x' if the latter is
    a directory.
    """
    result = []
    for item in mypypath:
        python2_dir = os.path.join(item, PYTHON2_STUB_DIR)
        if os.path.isdir(python2_dir):
            # @python2 takes precedence, but we also look into the parent
            # directory.
            result.append(python2_dir)
            result.append(item)
        else:
            result.append(item)
    return result


def compute_search_paths(sources: List[BuildSource],
                         options: Options,
                         data_dir: str,
                         alt_lib_path: Optional[str] = None) -> SearchPaths:
    """Compute the search paths as specified in PEP 561.

    There are the following 4 members created:
    - User code (from `sources`)
    - MYPYPATH (set either via config or environment variable)
    - installed package directories (which will later be split into stub-only and inline)
    - typeshed
     """
    # Determine the default module search path.
    lib_path = collections.deque(
        default_lib_path(data_dir,
                         options.python_version,
                         custom_typeshed_dir=options.custom_typeshed_dir))

    if options.use_builtins_fixtures:
        # Use stub builtins (to speed up test cases and to make them easier to
        # debug).  This is a test-only feature, so assume our files are laid out
        # as in the source tree.
        # We also need to allow overriding where to look for it. Argh.
        root_dir = os.getenv('MYPY_TEST_PREFIX', None)
        if not root_dir:
            root_dir = os.path.dirname(os.path.dirname(__file__))
        lib_path.appendleft(os.path.join(root_dir, 'test-data', 'unit', 'lib-stub'))
    # alt_lib_path is used by some tests to bypass the normal lib_path mechanics.
    # If we don't have one, grab directories of source files.
    python_path: List[str] = []
    if not alt_lib_path:
        for source in sources:
            # Include directory of the program file in the module search path.
            if source.base_dir:
                dir = source.base_dir
                if dir not in python_path:
                    python_path.append(dir)

        # Do this even if running as a file, for sanity (mainly because with
        # multiple builds, there could be a mix of files/modules, so its easier
        # to just define the semantics that we always add the current director
        # to the lib_path
        # TODO: Don't do this in some cases; for motivation see see
        # https://github.com/python/mypy/issues/4195#issuecomment-341915031
        if options.bazel:
            dir = '.'
        else:
            dir = os.getcwd()
        if dir not in lib_path:
            python_path.insert(0, dir)

    # Start with a MYPYPATH environment variable at the front of the mypy_path, if defined.
    mypypath = mypy_path()

    # Add a config-defined mypy path.
    mypypath.extend(options.mypy_path)

    # If provided, insert the caller-supplied extra module path to the
    # beginning (highest priority) of the search path.
    if alt_lib_path:
        mypypath.insert(0, alt_lib_path)

    # When type checking in Python 2 module, add @python2 subdirectories of
    # path items into the search path.
    if options.python_version[0] == 2:
        mypypath = add_py2_mypypath_entries(mypypath)

    egg_dirs, site_packages = get_site_packages_dirs(options.python_executable)
    base_prefix, prefix = get_prefixes(options.python_executable)
    is_venv = base_prefix != prefix
    for site_dir in site_packages:
        assert site_dir not in lib_path
        if (site_dir in mypypath or
                any(p.startswith(site_dir + os.path.sep) for p in mypypath) or
                os.path.altsep and any(p.startswith(site_dir + os.path.altsep) for p in mypypath)):
            print("{} is in the MYPYPATH. Please remove it.".format(site_dir), file=sys.stderr)
            print("See https://mypy.readthedocs.io/en/stable/running_mypy.html"
                  "#how-mypy-handles-imports for more info", file=sys.stderr)
            sys.exit(1)
        elif site_dir in python_path and (is_venv and not site_dir.startswith(prefix)):
            print("{} is in the PYTHONPATH. Please change directory"
                  " so it is not.".format(site_dir),
                  file=sys.stderr)
            sys.exit(1)

    return SearchPaths(python_path=tuple(reversed(python_path)),
                       mypy_path=tuple(mypypath),
                       package_path=tuple(egg_dirs + site_packages),
                       typeshed_path=tuple(lib_path))


def load_stdlib_py_versions(custom_typeshed_dir: Optional[str]) -> StdlibVersions:
    """Return dict with minimum and maximum Python versions of stdlib modules.

    The contents look like
    {..., 'secrets': ((3, 6), None), 'symbol': ((2, 7), (3, 9)), ...}

    None means there is no maximum version.
    """
    typeshed_dir = custom_typeshed_dir or os.path.join(os.path.dirname(__file__), "typeshed")
    stdlib_dir = os.path.join(typeshed_dir, "stdlib")
    result = {}

    versions_path = os.path.join(stdlib_dir, "VERSIONS")
    assert os.path.isfile(versions_path), (custom_typeshed_dir, versions_path, __file__)
    with open(versions_path) as f:
        for line in f:
            line = line.split("#")[0].strip()
            if line == "":
                continue
            module, version_range = line.split(":")
            versions = version_range.split("-")
            min_version = parse_version(versions[0])
            max_version = (parse_version(versions[1])
                           if len(versions) >= 2 and versions[1].strip() else None)
            result[module] = min_version, max_version

    # Modules that are Python 2 only or have separate Python 2 stubs
    # have stubs in @python2/ and may need an override.
    python2_dir = os.path.join(stdlib_dir, PYTHON2_STUB_DIR)
    try:
        for fnam in os.listdir(python2_dir):
            fnam = fnam.replace(".pyi", "")
            max_version = result.get(fnam, ((2, 7), None))[1]
            result[fnam] = (2, 7), max_version
    except FileNotFoundError:
        # Ignore error to support installations where Python 2 stubs aren't available.
        pass

    return result


def parse_version(version: str) -> Tuple[int, int]:
    major, minor = version.strip().split(".")
    return int(major), int(minor)


def typeshed_py_version(options: Options) -> Tuple[int, int]:
    """Return Python version used for checking whether module supports typeshed."""
    # Typeshed no longer covers Python 3.x versions before 3.6, so 3.6 is
    # the earliest we can support.
    if options.python_version[0] >= 3:
        return max(options.python_version, (3, 6))
    else:
        return options.python_version


def filter_redundant_py2_dirs(dirs: List[str]) -> List[str]:
    """If dirs has <dir>/@python2 followed by <dir>, filter out the latter."""
    if len(dirs) <= 1 or not any(d.endswith(PYTHON2_STUB_DIR) for d in dirs):
        # Fast path -- nothing to do
        return dirs
    seen = []
    result = []
    for d in dirs:
        if d.endswith(PYTHON2_STUB_DIR):
            seen.append(os.path.dirname(d))
            result.append(d)
        elif d not in seen:
            result.append(d)
    return result
