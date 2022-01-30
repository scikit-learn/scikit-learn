from __future__ import absolute_import, print_function

import cython
from .. import __version__

import collections
import contextlib
import hashlib
import os
import shutil
import subprocess
import re, sys, time
import warnings
from glob import iglob
from io import open as io_open
from os.path import relpath as _relpath
from distutils.extension import Extension
from distutils.util import strtobool
import zipfile

try:
    from collections.abc import Iterable
except ImportError:
    from collections import Iterable

try:
    import gzip
    gzip_open = gzip.open
    gzip_ext = '.gz'
except ImportError:
    gzip_open = open
    gzip_ext = ''

try:
    import zlib
    zipfile_compression_mode = zipfile.ZIP_DEFLATED
except ImportError:
    zipfile_compression_mode = zipfile.ZIP_STORED

try:
    import pythran
except:
    pythran = None

from .. import Utils
from ..Utils import (cached_function, cached_method, path_exists,
    safe_makedirs, copy_file_to_dir_if_newer, is_package_dir, replace_suffix)
from ..Compiler.Main import Context, CompilationOptions, default_options

join_path = cached_function(os.path.join)
copy_once_if_newer = cached_function(copy_file_to_dir_if_newer)
safe_makedirs_once = cached_function(safe_makedirs)

if sys.version_info[0] < 3:
    # stupid Py2 distutils enforces str type in list of sources
    _fs_encoding = sys.getfilesystemencoding()
    if _fs_encoding is None:
        _fs_encoding = sys.getdefaultencoding()
    def encode_filename_in_py2(filename):
        if not isinstance(filename, bytes):
            return filename.encode(_fs_encoding)
        return filename
else:
    def encode_filename_in_py2(filename):
        return filename
    basestring = str


def _make_relative(file_paths, base=None):
    if not base:
        base = os.getcwd()
    if base[-1] != os.path.sep:
        base += os.path.sep
    return [_relpath(path, base) if path.startswith(base) else path
            for path in file_paths]


def extended_iglob(pattern):
    if '{' in pattern:
        m = re.match('(.*){([^}]+)}(.*)', pattern)
        if m:
            before, switch, after = m.groups()
            for case in switch.split(','):
                for path in extended_iglob(before + case + after):
                    yield path
            return
    if '**/' in pattern:
        seen = set()
        first, rest = pattern.split('**/', 1)
        if first:
            first = iglob(first+'/')
        else:
            first = ['']
        for root in first:
            for path in extended_iglob(join_path(root, rest)):
                if path not in seen:
                    seen.add(path)
                    yield path
            for path in extended_iglob(join_path(root, '*', '**/' + rest)):
                if path not in seen:
                    seen.add(path)
                    yield path
    else:
        for path in iglob(pattern):
            yield path


def nonempty(it, error_msg="expected non-empty iterator"):
    empty = True
    for value in it:
        empty = False
        yield value
    if empty:
        raise ValueError(error_msg)


@cached_function
def file_hash(filename):
    path = os.path.normpath(filename)
    prefix = ('%d:%s' % (len(path), path)).encode("UTF-8")
    m = hashlib.md5(prefix)
    with open(path, 'rb') as f:
        data = f.read(65000)
        while data:
            m.update(data)
            data = f.read(65000)
    return m.hexdigest()


def update_pythran_extension(ext):
    if pythran is None:
        raise RuntimeError("You first need to install Pythran to use the np_pythran directive.")
    try:
        pythran_ext = pythran.config.make_extension(python=True)
    except TypeError:  # older pythran version only
        pythran_ext = pythran.config.make_extension()

    ext.include_dirs.extend(pythran_ext['include_dirs'])
    ext.extra_compile_args.extend(pythran_ext['extra_compile_args'])
    ext.extra_link_args.extend(pythran_ext['extra_link_args'])
    ext.define_macros.extend(pythran_ext['define_macros'])
    ext.undef_macros.extend(pythran_ext['undef_macros'])
    ext.library_dirs.extend(pythran_ext['library_dirs'])
    ext.libraries.extend(pythran_ext['libraries'])
    ext.language = 'c++'

    # These options are not compatible with the way normal Cython extensions work
    for bad_option in ["-fwhole-program", "-fvisibility=hidden"]:
        try:
            ext.extra_compile_args.remove(bad_option)
        except ValueError:
            pass


def parse_list(s):
    """
    >>> parse_list("")
    []
    >>> parse_list("a")
    ['a']
    >>> parse_list("a b c")
    ['a', 'b', 'c']
    >>> parse_list("[a, b, c]")
    ['a', 'b', 'c']
    >>> parse_list('a " " b')
    ['a', ' ', 'b']
    >>> parse_list('[a, ",a", "a,", ",", ]')
    ['a', ',a', 'a,', ',']
    """
    if len(s) >= 2 and s[0] == '[' and s[-1] == ']':
        s = s[1:-1]
        delimiter = ','
    else:
        delimiter = ' '
    s, literals = strip_string_literals(s)
    def unquote(literal):
        literal = literal.strip()
        if literal[0] in "'\"":
            return literals[literal[1:-1]]
        else:
            return literal
    return [unquote(item) for item in s.split(delimiter) if item.strip()]


transitive_str = object()
transitive_list = object()
bool_or = object()

distutils_settings = {
    'name':                 str,
    'sources':              list,
    'define_macros':        list,
    'undef_macros':         list,
    'libraries':            transitive_list,
    'library_dirs':         transitive_list,
    'runtime_library_dirs': transitive_list,
    'include_dirs':         transitive_list,
    'extra_objects':        list,
    'extra_compile_args':   transitive_list,
    'extra_link_args':      transitive_list,
    'export_symbols':       list,
    'depends':              transitive_list,
    'language':             transitive_str,
    'np_pythran':           bool_or
}


@cython.locals(start=cython.Py_ssize_t, end=cython.Py_ssize_t)
def line_iter(source):
    if isinstance(source, basestring):
        start = 0
        while True:
            end = source.find('\n', start)
            if end == -1:
                yield source[start:]
                return
            yield source[start:end]
            start = end+1
    else:
        for line in source:
            yield line


class DistutilsInfo(object):

    def __init__(self, source=None, exn=None):
        self.values = {}
        if source is not None:
            for line in line_iter(source):
                line = line.lstrip()
                if not line:
                    continue
                if line[0] != '#':
                    break
                line = line[1:].lstrip()
                kind = next((k for k in ("distutils:","cython:") if line.startswith(k)), None)
                if kind is not None:
                    key, _, value = [s.strip() for s in line[len(kind):].partition('=')]
                    type = distutils_settings.get(key, None)
                    if line.startswith("cython:") and type is None: continue
                    if type in (list, transitive_list):
                        value = parse_list(value)
                        if key == 'define_macros':
                            value = [tuple(macro.split('=', 1))
                                     if '=' in macro else (macro, None)
                                     for macro in value]
                    if type is bool_or:
                        value = strtobool(value)
                    self.values[key] = value
        elif exn is not None:
            for key in distutils_settings:
                if key in ('name', 'sources','np_pythran'):
                    continue
                value = getattr(exn, key, None)
                if value:
                    self.values[key] = value

    def merge(self, other):
        if other is None:
            return self
        for key, value in other.values.items():
            type = distutils_settings[key]
            if type is transitive_str and key not in self.values:
                self.values[key] = value
            elif type is transitive_list:
                if key in self.values:
                    # Change a *copy* of the list (Trac #845)
                    all = self.values[key][:]
                    for v in value:
                        if v not in all:
                            all.append(v)
                    value = all
                self.values[key] = value
            elif type is bool_or:
                self.values[key] = self.values.get(key, False) | value
        return self

    def subs(self, aliases):
        if aliases is None:
            return self
        resolved = DistutilsInfo()
        for key, value in self.values.items():
            type = distutils_settings[key]
            if type in [list, transitive_list]:
                new_value_list = []
                for v in value:
                    if v in aliases:
                        v = aliases[v]
                    if isinstance(v, list):
                        new_value_list += v
                    else:
                        new_value_list.append(v)
                value = new_value_list
            else:
                if value in aliases:
                    value = aliases[value]
            resolved.values[key] = value
        return resolved

    def apply(self, extension):
        for key, value in self.values.items():
            type = distutils_settings[key]
            if type in [list, transitive_list]:
                value = getattr(extension, key) + list(value)
            setattr(extension, key, value)


@cython.locals(start=cython.Py_ssize_t, q=cython.Py_ssize_t,
               single_q=cython.Py_ssize_t, double_q=cython.Py_ssize_t,
               hash_mark=cython.Py_ssize_t, end=cython.Py_ssize_t,
               k=cython.Py_ssize_t, counter=cython.Py_ssize_t, quote_len=cython.Py_ssize_t)
def strip_string_literals(code, prefix='__Pyx_L'):
    """
    Normalizes every string literal to be of the form '__Pyx_Lxxx',
    returning the normalized code and a mapping of labels to
    string literals.
    """
    new_code = []
    literals = {}
    counter = 0
    start = q = 0
    in_quote = False
    hash_mark = single_q = double_q = -1
    code_len = len(code)
    quote_type = None
    quote_len = -1

    while True:
        if hash_mark < q:
            hash_mark = code.find('#', q)
        if single_q < q:
            single_q = code.find("'", q)
        if double_q < q:
            double_q = code.find('"', q)
        q = min(single_q, double_q)
        if q == -1:
            q = max(single_q, double_q)

        # We're done.
        if q == -1 and hash_mark == -1:
            new_code.append(code[start:])
            break

        # Try to close the quote.
        elif in_quote:
            if code[q-1] == u'\\':
                k = 2
                while q >= k and code[q-k] == u'\\':
                    k += 1
                if k % 2 == 0:
                    q += 1
                    continue
            if code[q] == quote_type and (
                    quote_len == 1 or (code_len > q + 2 and quote_type == code[q+1] == code[q+2])):
                counter += 1
                label = "%s%s_" % (prefix, counter)
                literals[label] = code[start+quote_len:q]
                full_quote = code[q:q+quote_len]
                new_code.append(full_quote)
                new_code.append(label)
                new_code.append(full_quote)
                q += quote_len
                in_quote = False
                start = q
            else:
                q += 1

        # Process comment.
        elif -1 != hash_mark and (hash_mark < q or q == -1):
            new_code.append(code[start:hash_mark+1])
            end = code.find('\n', hash_mark)
            counter += 1
            label = "%s%s_" % (prefix, counter)
            if end == -1:
                end_or_none = None
            else:
                end_or_none = end
            literals[label] = code[hash_mark+1:end_or_none]
            new_code.append(label)
            if end == -1:
                break
            start = q = end

        # Open the quote.
        else:
            if code_len >= q+3 and (code[q] == code[q+1] == code[q+2]):
                quote_len = 3
            else:
                quote_len = 1
            in_quote = True
            quote_type = code[q]
            new_code.append(code[start:q])
            start = q
            q += quote_len

    return "".join(new_code), literals


# We need to allow spaces to allow for conditional compilation like
# IF ...:
#     cimport ...
dependency_regex = re.compile(r"(?:^\s*from +([0-9a-zA-Z_.]+) +cimport)|"
                              r"(?:^\s*cimport +([0-9a-zA-Z_.]+(?: *, *[0-9a-zA-Z_.]+)*))|"
                              r"(?:^\s*cdef +extern +from +['\"]([^'\"]+)['\"])|"
                              r"(?:^\s*include +['\"]([^'\"]+)['\"])", re.M)
dependency_after_from_regex = re.compile(
    r"(?:^\s+\(([0-9a-zA-Z_., ]*)\)[#\n])|"
    r"(?:^\s+([0-9a-zA-Z_., ]*)[#\n])",
    re.M)


def normalize_existing(base_path, rel_paths):
    return normalize_existing0(os.path.dirname(base_path), tuple(set(rel_paths)))


@cached_function
def normalize_existing0(base_dir, rel_paths):
    """
    Given some base directory ``base_dir`` and a list of path names
    ``rel_paths``, normalize each relative path name ``rel`` by
    replacing it by ``os.path.join(base, rel)`` if that file exists.

    Return a couple ``(normalized, needed_base)`` where ``normalized``
    if the list of normalized file names and ``needed_base`` is
    ``base_dir`` if we actually needed ``base_dir``. If no paths were
    changed (for example, if all paths were already absolute), then
    ``needed_base`` is ``None``.
    """
    normalized = []
    needed_base = None
    for rel in rel_paths:
        if os.path.isabs(rel):
            normalized.append(rel)
            continue
        path = join_path(base_dir, rel)
        if path_exists(path):
            normalized.append(os.path.normpath(path))
            needed_base = base_dir
        else:
            normalized.append(rel)
    return (normalized, needed_base)


def resolve_depends(depends, include_dirs):
    include_dirs = tuple(include_dirs)
    resolved = []
    for depend in depends:
        path = resolve_depend(depend, include_dirs)
        if path is not None:
            resolved.append(path)
    return resolved


@cached_function
def resolve_depend(depend, include_dirs):
    if depend[0] == '<' and depend[-1] == '>':
        return None
    for dir in include_dirs:
        path = join_path(dir, depend)
        if path_exists(path):
            return os.path.normpath(path)
    return None


@cached_function
def package(filename):
    dir = os.path.dirname(os.path.abspath(str(filename)))
    if dir != filename and is_package_dir(dir):
        return package(dir) + (os.path.basename(dir),)
    else:
        return ()


@cached_function
def fully_qualified_name(filename):
    module = os.path.splitext(os.path.basename(filename))[0]
    return '.'.join(package(filename) + (module,))


@cached_function
def parse_dependencies(source_filename):
    # Actual parsing is way too slow, so we use regular expressions.
    # The only catch is that we must strip comments and string
    # literals ahead of time.
    with Utils.open_source_file(source_filename, error_handling='ignore') as fh:
        source = fh.read()
    distutils_info = DistutilsInfo(source)
    source, literals = strip_string_literals(source)
    source = source.replace('\\\n', ' ').replace('\t', ' ')

    # TODO: pure mode
    cimports = []
    includes = []
    externs  = []
    for m in dependency_regex.finditer(source):
        cimport_from, cimport_list, extern, include = m.groups()
        if cimport_from:
            cimports.append(cimport_from)
            m_after_from = dependency_after_from_regex.search(source, pos=m.end())
            if m_after_from:
                multiline, one_line = m_after_from.groups()
                subimports = multiline or one_line
                cimports.extend("{0}.{1}".format(cimport_from, s.strip())
                                for s in subimports.split(','))

        elif cimport_list:
            cimports.extend(x.strip() for x in cimport_list.split(","))
        elif extern:
            externs.append(literals[extern])
        else:
            includes.append(literals[include])
    return cimports, includes, externs, distutils_info


class DependencyTree(object):

    def __init__(self, context, quiet=False):
        self.context = context
        self.quiet = quiet
        self._transitive_cache = {}

    def parse_dependencies(self, source_filename):
        if path_exists(source_filename):
            source_filename = os.path.normpath(source_filename)
        return parse_dependencies(source_filename)

    @cached_method
    def included_files(self, filename):
        # This is messy because included files are textually included, resolving
        # cimports (but not includes) relative to the including file.
        all = set()
        for include in self.parse_dependencies(filename)[1]:
            include_path = join_path(os.path.dirname(filename), include)
            if not path_exists(include_path):
                include_path = self.context.find_include_file(include, None)
            if include_path:
                if '.' + os.path.sep in include_path:
                    include_path = os.path.normpath(include_path)
                all.add(include_path)
                all.update(self.included_files(include_path))
            elif not self.quiet:
                print("Unable to locate '%s' referenced from '%s'" % (filename, include))
        return all

    @cached_method
    def cimports_externs_incdirs(self, filename):
        # This is really ugly. Nested cimports are resolved with respect to the
        # includer, but includes are resolved with respect to the includee.
        cimports, includes, externs = self.parse_dependencies(filename)[:3]
        cimports = set(cimports)
        externs = set(externs)
        incdirs = set()
        for include in self.included_files(filename):
            included_cimports, included_externs, included_incdirs = self.cimports_externs_incdirs(include)
            cimports.update(included_cimports)
            externs.update(included_externs)
            incdirs.update(included_incdirs)
        externs, incdir = normalize_existing(filename, externs)
        if incdir:
            incdirs.add(incdir)
        return tuple(cimports), externs, incdirs

    def cimports(self, filename):
        return self.cimports_externs_incdirs(filename)[0]

    def package(self, filename):
        return package(filename)

    def fully_qualified_name(self, filename):
        return fully_qualified_name(filename)

    @cached_method
    def find_pxd(self, module, filename=None):
        is_relative = module[0] == '.'
        if is_relative and not filename:
            raise NotImplementedError("New relative imports.")
        if filename is not None:
            module_path = module.split('.')
            if is_relative:
                module_path.pop(0)  # just explicitly relative
            package_path = list(self.package(filename))
            while module_path and not module_path[0]:
                try:
                    package_path.pop()
                except IndexError:
                    return None   # FIXME: error?
                module_path.pop(0)
            relative = '.'.join(package_path + module_path)
            pxd = self.context.find_pxd_file(relative, None)
            if pxd:
                return pxd
        if is_relative:
            return None   # FIXME: error?
        return self.context.find_pxd_file(module, None)

    @cached_method
    def cimported_files(self, filename):
        if filename[-4:] == '.pyx' and path_exists(filename[:-4] + '.pxd'):
            pxd_list = [filename[:-4] + '.pxd']
        else:
            pxd_list = []
        # Cimports generates all possible combinations package.module
        # when imported as from package cimport module.
        for module in self.cimports(filename):
            if module[:7] == 'cython.' or module == 'cython':
                continue
            pxd_file = self.find_pxd(module, filename)
            if pxd_file is not None:
                pxd_list.append(pxd_file)
        return tuple(pxd_list)

    @cached_method
    def immediate_dependencies(self, filename):
        all = set([filename])
        all.update(self.cimported_files(filename))
        all.update(self.included_files(filename))
        return all

    def all_dependencies(self, filename):
        return self.transitive_merge(filename, self.immediate_dependencies, set.union)

    @cached_method
    def timestamp(self, filename):
        return os.path.getmtime(filename)

    def extract_timestamp(self, filename):
        return self.timestamp(filename), filename

    def newest_dependency(self, filename):
        return max([self.extract_timestamp(f) for f in self.all_dependencies(filename)])

    def transitive_fingerprint(self, filename, module, compilation_options):
        r"""
        Return a fingerprint of a cython file that is about to be cythonized.

        Fingerprints are looked up in future compilations. If the fingerprint
        is found, the cythonization can be skipped. The fingerprint must
        incorporate everything that has an influence on the generated code.
        """
        try:
            m = hashlib.md5(__version__.encode('UTF-8'))
            m.update(file_hash(filename).encode('UTF-8'))
            for x in sorted(self.all_dependencies(filename)):
                if os.path.splitext(x)[1] not in ('.c', '.cpp', '.h'):
                    m.update(file_hash(x).encode('UTF-8'))
            # Include the module attributes that change the compilation result
            # in the fingerprint. We do not iterate over module.__dict__ and
            # include almost everything here as users might extend Extension
            # with arbitrary (random) attributes that would lead to cache
            # misses.
            m.update(str((
                module.language,
                getattr(module, 'py_limited_api', False),
                getattr(module, 'np_pythran', False)
            )).encode('UTF-8'))

            m.update(compilation_options.get_fingerprint().encode('UTF-8'))
            return m.hexdigest()
        except IOError:
            return None

    def distutils_info0(self, filename):
        info = self.parse_dependencies(filename)[3]
        kwds = info.values
        cimports, externs, incdirs = self.cimports_externs_incdirs(filename)
        basedir = os.getcwd()
        # Add dependencies on "cdef extern from ..." files
        if externs:
            externs = _make_relative(externs, basedir)
            if 'depends' in kwds:
                kwds['depends'] = list(set(kwds['depends']).union(externs))
            else:
                kwds['depends'] = list(externs)
        # Add include_dirs to ensure that the C compiler will find the
        # "cdef extern from ..." files
        if incdirs:
            include_dirs = list(kwds.get('include_dirs', []))
            for inc in _make_relative(incdirs, basedir):
                if inc not in include_dirs:
                    include_dirs.append(inc)
            kwds['include_dirs'] = include_dirs
        return info

    def distutils_info(self, filename, aliases=None, base=None):
        return (self.transitive_merge(filename, self.distutils_info0, DistutilsInfo.merge)
            .subs(aliases)
            .merge(base))

    def transitive_merge(self, node, extract, merge):
        try:
            seen = self._transitive_cache[extract, merge]
        except KeyError:
            seen = self._transitive_cache[extract, merge] = {}
        return self.transitive_merge_helper(
            node, extract, merge, seen, {}, self.cimported_files)[0]

    def transitive_merge_helper(self, node, extract, merge, seen, stack, outgoing):
        if node in seen:
            return seen[node], None
        deps = extract(node)
        if node in stack:
            return deps, node
        try:
            stack[node] = len(stack)
            loop = None
            for next in outgoing(node):
                sub_deps, sub_loop = self.transitive_merge_helper(next, extract, merge, seen, stack, outgoing)
                if sub_loop is not None:
                    if loop is not None and stack[loop] < stack[sub_loop]:
                        pass
                    else:
                        loop = sub_loop
                deps = merge(deps, sub_deps)
            if loop == node:
                loop = None
            if loop is None:
                seen[node] = deps
            return deps, loop
        finally:
            del stack[node]


_dep_tree = None

def create_dependency_tree(ctx=None, quiet=False):
    global _dep_tree
    if _dep_tree is None:
        if ctx is None:
            ctx = Context(["."], CompilationOptions(default_options))
        _dep_tree = DependencyTree(ctx, quiet=quiet)
    return _dep_tree


# If this changes, change also docs/src/reference/compilation.rst
# which mentions this function
def default_create_extension(template, kwds):
    if 'depends' in kwds:
        include_dirs = kwds.get('include_dirs', []) + ["."]
        depends = resolve_depends(kwds['depends'], include_dirs)
        kwds['depends'] = sorted(set(depends + template.depends))

    t = template.__class__
    ext = t(**kwds)
    metadata = dict(distutils=kwds, module_name=kwds['name'])
    return (ext, metadata)


# This may be useful for advanced users?
def create_extension_list(patterns, exclude=None, ctx=None, aliases=None, quiet=False, language=None,
                          exclude_failures=False):
    if language is not None:
        print('Warning: passing language={0!r} to cythonize() is deprecated. '
              'Instead, put "# distutils: language={0}" in your .pyx or .pxd file(s)'.format(language))
    if exclude is None:
        exclude = []
    if patterns is None:
        return [], {}
    elif isinstance(patterns, basestring) or not isinstance(patterns, Iterable):
        patterns = [patterns]
    explicit_modules = set([m.name for m in patterns if isinstance(m, Extension)])
    seen = set()
    deps = create_dependency_tree(ctx, quiet=quiet)
    to_exclude = set()
    if not isinstance(exclude, list):
        exclude = [exclude]
    for pattern in exclude:
        to_exclude.update(map(os.path.abspath, extended_iglob(pattern)))

    module_list = []
    module_metadata = {}

    # workaround for setuptools
    if 'setuptools' in sys.modules:
        Extension_distutils = sys.modules['setuptools.extension']._Extension
        Extension_setuptools = sys.modules['setuptools'].Extension
    else:
        # dummy class, in case we do not have setuptools
        Extension_distutils = Extension
        class Extension_setuptools(Extension): pass

    # if no create_extension() function is defined, use a simple
    # default function.
    create_extension = ctx.options.create_extension or default_create_extension

    for pattern in patterns:
        if isinstance(pattern, str):
            filepattern = pattern
            template = Extension(pattern, [])  # Fake Extension without sources
            name = '*'
            base = None
            ext_language = language
        elif isinstance(pattern, (Extension_distutils, Extension_setuptools)):
            cython_sources = [s for s in pattern.sources
                              if os.path.splitext(s)[1] in ('.py', '.pyx')]
            if cython_sources:
                filepattern = cython_sources[0]
                if len(cython_sources) > 1:
                    print("Warning: Multiple cython sources found for extension '%s': %s\n"
                          "See http://cython.readthedocs.io/en/latest/src/userguide/sharing_declarations.html "
                          "for sharing declarations among Cython files." % (pattern.name, cython_sources))
            else:
                # ignore non-cython modules
                module_list.append(pattern)
                continue
            template = pattern
            name = template.name
            base = DistutilsInfo(exn=template)
            ext_language = None  # do not override whatever the Extension says
        else:
            msg = str("pattern is not of type str nor subclass of Extension (%s)"
                      " but of type %s and class %s" % (repr(Extension),
                                                        type(pattern),
                                                        pattern.__class__))
            raise TypeError(msg)

        for file in nonempty(sorted(extended_iglob(filepattern)), "'%s' doesn't match any files" % filepattern):
            if os.path.abspath(file) in to_exclude:
                continue
            module_name = deps.fully_qualified_name(file)
            if '*' in name:
                if module_name in explicit_modules:
                    continue
            elif name:
                module_name = name

            Utils.raise_error_if_module_name_forbidden(module_name)

            if module_name not in seen:
                try:
                    kwds = deps.distutils_info(file, aliases, base).values
                except Exception:
                    if exclude_failures:
                        continue
                    raise
                if base is not None:
                    for key, value in base.values.items():
                        if key not in kwds:
                            kwds[key] = value

                kwds['name'] = module_name

                sources = [file] + [m for m in template.sources if m != filepattern]
                if 'sources' in kwds:
                    # allow users to add .c files etc.
                    for source in kwds['sources']:
                        source = encode_filename_in_py2(source)
                        if source not in sources:
                            sources.append(source)
                kwds['sources'] = sources

                if ext_language and 'language' not in kwds:
                    kwds['language'] = ext_language

                np_pythran = kwds.pop('np_pythran', False)

                # Create the new extension
                m, metadata = create_extension(template, kwds)
                m.np_pythran = np_pythran or getattr(m, 'np_pythran', False)
                if m.np_pythran:
                    update_pythran_extension(m)
                module_list.append(m)

                # Store metadata (this will be written as JSON in the
                # generated C file but otherwise has no purpose)
                module_metadata[module_name] = metadata

                if file not in m.sources:
                    # Old setuptools unconditionally replaces .pyx with .c/.cpp
                    target_file = os.path.splitext(file)[0] + ('.cpp' if m.language == 'c++' else '.c')
                    try:
                        m.sources.remove(target_file)
                    except ValueError:
                        # never seen this in the wild, but probably better to warn about this unexpected case
                        print("Warning: Cython source file not found in sources list, adding %s" % file)
                    m.sources.insert(0, file)
                seen.add(name)
    return module_list, module_metadata


# This is the user-exposed entry point.
def cythonize(module_list, exclude=None, nthreads=0, aliases=None, quiet=False, force=False, language=None,
              exclude_failures=False, **options):
    """
    Compile a set of source modules into C/C++ files and return a list of distutils
    Extension objects for them.

    :param module_list: As module list, pass either a glob pattern, a list of glob
                        patterns or a list of Extension objects.  The latter
                        allows you to configure the extensions separately
                        through the normal distutils options.
                        You can also pass Extension objects that have
                        glob patterns as their sources. Then, cythonize
                        will resolve the pattern and create a
                        copy of the Extension for every matching file.

    :param exclude: When passing glob patterns as ``module_list``, you can exclude certain
                    module names explicitly by passing them into the ``exclude`` option.

    :param nthreads: The number of concurrent builds for parallel compilation
                     (requires the ``multiprocessing`` module).

    :param aliases: If you want to use compiler directives like ``# distutils: ...`` but
                    can only know at compile time (when running the ``setup.py``) which values
                    to use, you can use aliases and pass a dictionary mapping those aliases
                    to Python strings when calling :func:`cythonize`. As an example, say you
                    want to use the compiler
                    directive ``# distutils: include_dirs = ../static_libs/include/``
                    but this path isn't always fixed and you want to find it when running
                    the ``setup.py``. You can then do ``# distutils: include_dirs = MY_HEADERS``,
                    find the value of ``MY_HEADERS`` in the ``setup.py``, put it in a python
                    variable called ``foo`` as a string, and then call
                    ``cythonize(..., aliases={'MY_HEADERS': foo})``.

    :param quiet: If True, Cython won't print error, warning, or status messages during the
                  compilation.

    :param force: Forces the recompilation of the Cython modules, even if the timestamps
                  don't indicate that a recompilation is necessary.

    :param language: To globally enable C++ mode, you can pass ``language='c++'``. Otherwise, this
                     will be determined at a per-file level based on compiler directives.  This
                     affects only modules found based on file names.  Extension instances passed
                     into :func:`cythonize` will not be changed. It is recommended to rather
                     use the compiler directive ``# distutils: language = c++`` than this option.

    :param exclude_failures: For a broad 'try to compile' mode that ignores compilation
                             failures and simply excludes the failed extensions,
                             pass ``exclude_failures=True``. Note that this only
                             really makes sense for compiling ``.py`` files which can also
                             be used without compilation.

    :param annotate: If ``True``, will produce a HTML file for each of the ``.pyx`` or ``.py``
                     files compiled. The HTML file gives an indication
                     of how much Python interaction there is in
                     each of the source code lines, compared to plain C code.
                     It also allows you to see the C/C++ code
                     generated for each line of Cython code. This report is invaluable when
                     optimizing a function for speed,
                     and for determining when to :ref:`release the GIL <nogil>`:
                     in general, a ``nogil`` block may contain only "white" code.
                     See examples in :ref:`determining_where_to_add_types` or
                     :ref:`primes`.

    :param compiler_directives: Allow to set compiler directives in the ``setup.py`` like this:
                                ``compiler_directives={'embedsignature': True}``.
                                See :ref:`compiler-directives`.

    :param depfile: produce depfiles for the sources if True.
    """
    if exclude is None:
        exclude = []
    if 'include_path' not in options:
        options['include_path'] = ['.']
    if 'common_utility_include_dir' in options:
        safe_makedirs(options['common_utility_include_dir'])

    depfile = options.pop('depfile', None)

    if pythran is None:
        pythran_options = None
    else:
        pythran_options = CompilationOptions(**options)
        pythran_options.cplus = True
        pythran_options.np_pythran = True

    c_options = CompilationOptions(**options)
    cpp_options = CompilationOptions(**options); cpp_options.cplus = True
    ctx = c_options.create_context()
    options = c_options
    module_list, module_metadata = create_extension_list(
        module_list,
        exclude=exclude,
        ctx=ctx,
        quiet=quiet,
        exclude_failures=exclude_failures,
        language=language,
        aliases=aliases)
    deps = create_dependency_tree(ctx, quiet=quiet)
    build_dir = getattr(options, 'build_dir', None)

    def copy_to_build_dir(filepath, root=os.getcwd()):
        filepath_abs = os.path.abspath(filepath)
        if os.path.isabs(filepath):
            filepath = filepath_abs
        if filepath_abs.startswith(root):
            # distutil extension depends are relative to cwd
            mod_dir = join_path(build_dir,
                                os.path.dirname(_relpath(filepath, root)))
            copy_once_if_newer(filepath_abs, mod_dir)

    modules_by_cfile = collections.defaultdict(list)
    to_compile = []
    for m in module_list:
        if build_dir:
            for dep in m.depends:
                copy_to_build_dir(dep)

        cy_sources = [
            source for source in m.sources
            if os.path.splitext(source)[1] in ('.pyx', '.py')]
        if len(cy_sources) == 1:
            # normal "special" case: believe the Extension module name to allow user overrides
            full_module_name = m.name
        else:
            # infer FQMN from source files
            full_module_name = None

        new_sources = []
        for source in m.sources:
            base, ext = os.path.splitext(source)
            if ext in ('.pyx', '.py'):
                if m.np_pythran:
                    c_file = base + '.cpp'
                    options = pythran_options
                elif m.language == 'c++':
                    c_file = base + '.cpp'
                    options = cpp_options
                else:
                    c_file = base + '.c'
                    options = c_options

                # setup for out of place build directory if enabled
                if build_dir:
                    if os.path.isabs(c_file):
                      warnings.warn("build_dir has no effect for absolute source paths")
                    c_file = os.path.join(build_dir, c_file)
                    dir = os.path.dirname(c_file)
                    safe_makedirs_once(dir)

                # write out the depfile, if requested
                if depfile:
                    dependencies = deps.all_dependencies(source)
                    src_base_dir, _ = os.path.split(source)
                    if not src_base_dir.endswith(os.sep):
                        src_base_dir += os.sep
                    # paths below the base_dir are relative, otherwise absolute
                    paths = []
                    for fname in dependencies:
                        if (fname.startswith(src_base_dir) or
                            fname.startswith('.' + os.path.sep)):
                            paths.append(os.path.relpath(fname, src_base_dir))
                        else:
                            paths.append(os.path.abspath(fname))

                    depline = os.path.split(c_file)[1] + ": \\\n  "
                    depline += " \\\n  ".join(paths) + "\n"
                    with open(c_file+'.dep', 'w') as outfile:
                        outfile.write(depline)

                if os.path.exists(c_file):
                    c_timestamp = os.path.getmtime(c_file)
                else:
                    c_timestamp = -1

                # Priority goes first to modified files, second to direct
                # dependents, and finally to indirect dependents.
                if c_timestamp < deps.timestamp(source):
                    dep_timestamp, dep = deps.timestamp(source), source
                    priority = 0
                else:
                    dep_timestamp, dep = deps.newest_dependency(source)
                    priority = 2 - (dep in deps.immediate_dependencies(source))
                if force or c_timestamp < dep_timestamp:
                    if not quiet and not force:
                        if source == dep:
                            print("Compiling %s because it changed." % source)
                        else:
                            print("Compiling %s because it depends on %s." % (source, dep))
                    if not force and options.cache:
                        fingerprint = deps.transitive_fingerprint(source, m, options)
                    else:
                        fingerprint = None
                    to_compile.append((
                        priority, source, c_file, fingerprint, quiet,
                        options, not exclude_failures, module_metadata.get(m.name),
                        full_module_name))
                new_sources.append(c_file)
                modules_by_cfile[c_file].append(m)
            else:
                new_sources.append(source)
                if build_dir:
                    copy_to_build_dir(source)
        m.sources = new_sources

    if options.cache:
        if not os.path.exists(options.cache):
            os.makedirs(options.cache)
    to_compile.sort()
    # Drop "priority" component of "to_compile" entries and add a
    # simple progress indicator.
    N = len(to_compile)
    progress_fmt = "[{0:%d}/{1}] " % len(str(N))
    for i in range(N):
        progress = progress_fmt.format(i+1, N)
        to_compile[i] = to_compile[i][1:] + (progress,)

    if N <= 1:
        nthreads = 0
    if nthreads:
        # Requires multiprocessing (or Python >= 2.6)
        try:
            import multiprocessing
            pool = multiprocessing.Pool(
                nthreads, initializer=_init_multiprocessing_helper)
        except (ImportError, OSError):
            print("multiprocessing required for parallel cythonization")
            nthreads = 0
        else:
            # This is a bit more involved than it should be, because KeyboardInterrupts
            # break the multiprocessing workers when using a normal pool.map().
            # See, for example:
            # http://noswap.com/blog/python-multiprocessing-keyboardinterrupt
            try:
                result = pool.map_async(cythonize_one_helper, to_compile, chunksize=1)
                pool.close()
                while not result.ready():
                    try:
                        result.get(99999)  # seconds
                    except multiprocessing.TimeoutError:
                        pass
            except KeyboardInterrupt:
                pool.terminate()
                raise
            pool.join()
    if not nthreads:
        for args in to_compile:
            cythonize_one(*args)

    if exclude_failures:
        failed_modules = set()
        for c_file, modules in modules_by_cfile.items():
            if not os.path.exists(c_file):
                failed_modules.update(modules)
            elif os.path.getsize(c_file) < 200:
                f = io_open(c_file, 'r', encoding='iso8859-1')
                try:
                    if f.read(len('#error ')) == '#error ':
                        # dead compilation result
                        failed_modules.update(modules)
                finally:
                    f.close()
        if failed_modules:
            for module in failed_modules:
                module_list.remove(module)
            print("Failed compilations: %s" % ', '.join(sorted([
                module.name for module in failed_modules])))

    if options.cache:
        cleanup_cache(options.cache, getattr(options, 'cache_size', 1024 * 1024 * 100))
    # cythonize() is often followed by the (non-Python-buffered)
    # compiler output, flush now to avoid interleaving output.
    sys.stdout.flush()
    return module_list


if os.environ.get('XML_RESULTS'):
    compile_result_dir = os.environ['XML_RESULTS']
    def record_results(func):
        def with_record(*args):
            t = time.time()
            success = True
            try:
                try:
                    func(*args)
                except:
                    success = False
            finally:
                t = time.time() - t
                module = fully_qualified_name(args[0])
                name = "cythonize." + module
                failures = 1 - success
                if success:
                    failure_item = ""
                else:
                    failure_item = "failure"
                output = open(os.path.join(compile_result_dir, name + ".xml"), "w")
                output.write("""
                    <?xml version="1.0" ?>
                    <testsuite name="%(name)s" errors="0" failures="%(failures)s" tests="1" time="%(t)s">
                    <testcase classname="%(name)s" name="cythonize">
                    %(failure_item)s
                    </testcase>
                    </testsuite>
                """.strip() % locals())
                output.close()
        return with_record
else:
    def record_results(func):
        return func


# TODO: Share context? Issue: pyx processing leaks into pxd module
@record_results
def cythonize_one(pyx_file, c_file, fingerprint, quiet, options=None,
                  raise_on_failure=True, embedded_metadata=None, full_module_name=None,
                  progress=""):
    from ..Compiler.Main import compile_single, default_options
    from ..Compiler.Errors import CompileError, PyrexError

    if fingerprint:
        if not os.path.exists(options.cache):
            safe_makedirs(options.cache)
        # Cython-generated c files are highly compressible.
        # (E.g. a compression ratio of about 10 for Sage).
        fingerprint_file_base = join_path(
            options.cache, "%s-%s" % (os.path.basename(c_file), fingerprint))
        gz_fingerprint_file = fingerprint_file_base + gzip_ext
        zip_fingerprint_file = fingerprint_file_base + '.zip'
        if os.path.exists(gz_fingerprint_file) or os.path.exists(zip_fingerprint_file):
            if not quiet:
                print("%sFound compiled %s in cache" % (progress, pyx_file))
            if os.path.exists(gz_fingerprint_file):
                os.utime(gz_fingerprint_file, None)
                with contextlib.closing(gzip_open(gz_fingerprint_file, 'rb')) as g:
                    with contextlib.closing(open(c_file, 'wb')) as f:
                        shutil.copyfileobj(g, f)
            else:
                os.utime(zip_fingerprint_file, None)
                dirname = os.path.dirname(c_file)
                with contextlib.closing(zipfile.ZipFile(zip_fingerprint_file)) as z:
                    for artifact in z.namelist():
                        z.extract(artifact, os.path.join(dirname, artifact))
            return
    if not quiet:
        print("%sCythonizing %s" % (progress, pyx_file))
    if options is None:
        options = CompilationOptions(default_options)
    options.output_file = c_file
    options.embedded_metadata = embedded_metadata

    any_failures = 0
    try:
        result = compile_single(pyx_file, options, full_module_name=full_module_name)
        if result.num_errors > 0:
            any_failures = 1
    except (EnvironmentError, PyrexError) as e:
        sys.stderr.write('%s\n' % e)
        any_failures = 1
        # XXX
        import traceback
        traceback.print_exc()
    except Exception:
        if raise_on_failure:
            raise
        import traceback
        traceback.print_exc()
        any_failures = 1
    if any_failures:
        if raise_on_failure:
            raise CompileError(None, pyx_file)
        elif os.path.exists(c_file):
            os.remove(c_file)
    elif fingerprint:
        artifacts = list(filter(None, [
            getattr(result, attr, None)
            for attr in ('c_file', 'h_file', 'api_file', 'i_file')]))
        if len(artifacts) == 1:
            fingerprint_file = gz_fingerprint_file
            with contextlib.closing(open(c_file, 'rb')) as f:
                with contextlib.closing(gzip_open(fingerprint_file + '.tmp', 'wb')) as g:
                    shutil.copyfileobj(f, g)
        else:
            fingerprint_file = zip_fingerprint_file
            with contextlib.closing(zipfile.ZipFile(
                fingerprint_file + '.tmp', 'w', zipfile_compression_mode)) as zip:
                for artifact in artifacts:
                    zip.write(artifact, os.path.basename(artifact))
        os.rename(fingerprint_file + '.tmp', fingerprint_file)


def cythonize_one_helper(m):
    import traceback
    try:
        return cythonize_one(*m)
    except Exception:
        traceback.print_exc()
        raise


def _init_multiprocessing_helper():
    # KeyboardInterrupt kills workers, so don't let them get it
    import signal
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def cleanup_cache(cache, target_size, ratio=.85):
    try:
        p = subprocess.Popen(['du', '-s', '-k', os.path.abspath(cache)], stdout=subprocess.PIPE)
        res = p.wait()
        if res == 0:
            total_size = 1024 * int(p.stdout.read().strip().split()[0])
            if total_size < target_size:
                return
    except (OSError, ValueError):
        pass
    total_size = 0
    all = []
    for file in os.listdir(cache):
        path = join_path(cache, file)
        s = os.stat(path)
        total_size += s.st_size
        all.append((s.st_atime, s.st_size, path))
    if total_size > target_size:
        for time, size, file in reversed(sorted(all)):
            os.unlink(file)
            total_size -= size
            if total_size < target_size * ratio:
                break
