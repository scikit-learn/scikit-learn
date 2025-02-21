"""
A Cython plugin for coverage.py

Requires the coverage package at least in version 4.0 (which added the plugin API).

This plugin requires the generated C sources to be available, next to the extension module.
It parses the C file and reads the original source files from it, which are stored in C comments.
It then reports a source file to coverage.py when it hits one of its lines during line tracing.

Basically, Cython can (on request) emit explicit trace calls into the C code that it generates,
and as a general human debugging helper, it always copies the current source code line
(and its surrounding context) into the C files before it generates code for that line, e.g.

::

      /* "line_trace.pyx":147
       * def cy_add_with_nogil(a,b):
       *     cdef int z, x=a, y=b         # 1
       *     with nogil:                  # 2             # <<<<<<<<<<<<<<
       *         z = 0                    # 3
       *         z += cy_add_nogil(x, y)  # 4
       */
       __Pyx_TraceLine(147,1,__PYX_ERR(0, 147, __pyx_L4_error))
      [C code generated for file line_trace.pyx, line 147, follows here]

The crux is that multiple source files can contribute code to a single C (or C++) file
(and thus, to a single extension module) besides the main module source file (.py/.pyx),
usually shared declaration files (.pxd) but also literally included files (.pxi).

Therefore, the coverage plugin doesn't actually try to look at the file that happened
to contribute the current source line for the trace call, but simply looks up the single
.c file from which the extension was compiled (which usually lies right next to it after
the build, having the same name), and parses the code copy comments from that .c file
to recover the original source files and their code as a line-to-file mapping.

That mapping is then used to report the ``__Pyx_TraceLine()`` calls to the coverage tool.
The plugin also reports the line of source code that it found in the C file to the coverage
tool to support annotated source representations.  For this, again, it does not look at the
actual source files but only reports the source code that it found in the C code comments.

Apart from simplicity (read one file instead of finding and parsing many), part of the
reasoning here is that any line in the original sources for which there is no comment line
(and trace call) in the generated C code cannot count as executed, really, so the C code
comments are a very good source for coverage reporting.  They already filter out purely
declarative code lines that do not contribute executable code, and such (missing) lines
can then be marked as excluded from coverage analysis.
"""

from __future__ import absolute_import

import re
import os.path
import sys
from collections import defaultdict

from coverage.plugin import CoveragePlugin, FileTracer, FileReporter  # requires coverage.py 4.0+
from coverage.files import canonical_filename

from .Utils import find_root_package_dir, is_package_dir, is_cython_generated_file, open_source_file


from . import __version__


C_FILE_EXTENSIONS = ['.c', '.cpp', '.cc', '.cxx']
MODULE_FILE_EXTENSIONS = set(['.py', '.pyx', '.pxd'] + C_FILE_EXTENSIONS)


def _find_c_source(base_path):
    file_exists = os.path.exists
    for ext in C_FILE_EXTENSIONS:
        file_name = base_path + ext
        if file_exists(file_name):
            return file_name
    return None


def _find_dep_file_path(main_file, file_path, relative_path_search=False):
    abs_path = os.path.abspath(file_path)
    if not os.path.exists(abs_path) and (file_path.endswith('.pxi') or
                                         relative_path_search):
        # files are looked up relative to the main source file
        rel_file_path = os.path.join(os.path.dirname(main_file), file_path)
        if os.path.exists(rel_file_path):
            abs_path = os.path.abspath(rel_file_path)

        abs_no_ext = os.path.splitext(abs_path)[0]
        file_no_ext, extension = os.path.splitext(file_path)
        # We check if the paths match by matching the directories in reverse order.
        # pkg/module.pyx /long/absolute_path/bla/bla/site-packages/pkg/module.c should match.
        # this will match the pairs: module-module and pkg-pkg. After which there is nothing left to zip.
        abs_no_ext = os.path.normpath(abs_no_ext)
        file_no_ext = os.path.normpath(file_no_ext)
        matching_paths = zip(reversed(abs_no_ext.split(os.sep)), reversed(file_no_ext.split(os.sep)))
        for one, other in matching_paths:
            if one != other:
                break
        else:  # No mismatches detected
            matching_abs_path = os.path.splitext(main_file)[0] + extension
            if os.path.exists(matching_abs_path):
                return canonical_filename(matching_abs_path)

    # search sys.path for external locations if a valid file hasn't been found
    if not os.path.exists(abs_path):
        for sys_path in sys.path:
            test_path = os.path.realpath(os.path.join(sys_path, file_path))
            if os.path.exists(test_path):
                return canonical_filename(test_path)
    return canonical_filename(abs_path)


class Plugin(CoveragePlugin):
    # map from traced file paths to absolute file paths
    _file_path_map = None
    # map from traced file paths to corresponding C files
    _c_files_map = None
    # map from parsed C files to their content
    _parsed_c_files = None
    # map from traced files to lines that are excluded from coverage
    _excluded_lines_map = None
    # list of regex patterns for lines to exclude
    _excluded_line_patterns = ()

    def sys_info(self):
        return [('Cython version', __version__)]

    def configure(self, config):
        # Entry point for coverage "configurer".
        # Read the regular expressions from the coverage config that match lines to be excluded from coverage.
        self._excluded_line_patterns = config.get_option("report:exclude_lines")

    def file_tracer(self, filename):
        """
        Try to find a C source file for a file path found by the tracer.
        """
        if filename.startswith('<') or filename.startswith('memory:'):
            return None
        c_file = py_file = None
        filename = canonical_filename(os.path.abspath(filename))
        if self._c_files_map and filename in self._c_files_map:
            c_file = self._c_files_map[filename][0]

        if c_file is None:
            c_file, py_file = self._find_source_files(filename)
            if not c_file:
                return None  # unknown file

            # parse all source file paths and lines from C file
            # to learn about all relevant source files right away (pyx/pxi/pxd)
            # FIXME: this might already be too late if the first executed line
            #        is not from the main .pyx file but a file with a different
            #        name than the .c file (which prevents us from finding the
            #        .c file)
            _, code = self._read_source_lines(c_file, filename)
            if code is None:
                return None  # no source found

        if self._file_path_map is None:
            self._file_path_map = {}
        return CythonModuleTracer(filename, py_file, c_file, self._c_files_map, self._file_path_map)

    def file_reporter(self, filename):
        # TODO: let coverage.py handle .py files itself
        #ext = os.path.splitext(filename)[1].lower()
        #if ext == '.py':
        #    from coverage.python import PythonFileReporter
        #    return PythonFileReporter(filename)

        filename = canonical_filename(os.path.abspath(filename))
        if self._c_files_map and filename in self._c_files_map:
            c_file, rel_file_path, code = self._c_files_map[filename]
        else:
            c_file, _ = self._find_source_files(filename)
            if not c_file:
                return None  # unknown file
            rel_file_path, code = self._read_source_lines(c_file, filename)
            if code is None:
                return None  # no source found
        return CythonModuleReporter(
            c_file,
            filename,
            rel_file_path,
            code,
            self._excluded_lines_map.get(rel_file_path, frozenset())
        )

    def _find_source_files(self, filename):
        basename, ext = os.path.splitext(filename)
        ext = ext.lower()
        if ext in MODULE_FILE_EXTENSIONS:
            pass
        elif ext == '.pyd':
            # Windows extension module
            platform_suffix = re.search(r'[.]cp[0-9]+-win[_a-z0-9]*$', basename, re.I)
            if platform_suffix:
                basename = basename[:platform_suffix.start()]
        elif ext == '.so':
            # Linux/Unix/Mac extension module
            platform_suffix = re.search(r'[.](?:cpython|pypy)-[0-9]+[-_a-z0-9]*$', basename, re.I)
            if platform_suffix:
                basename = basename[:platform_suffix.start()]
        elif ext == '.pxi':
            # if we get here, it means that the first traced line of a Cython module was
            # not in the main module but in an include file, so try a little harder to
            # find the main source file
            self._find_c_source_files(os.path.dirname(filename), filename)
            if filename in self._c_files_map:
                return self._c_files_map[filename][0], None
        else:
            # none of our business
            return None, None

        c_file = filename if ext in C_FILE_EXTENSIONS else _find_c_source(basename)
        if c_file is None:
            # a module "pkg/mod.so" can have a source file "pkg/pkg.mod.c"
            package_root = find_root_package_dir.uncached(filename)
            package_path = os.path.relpath(basename, package_root).split(os.path.sep)
            if len(package_path) > 1:
                test_basepath = os.path.join(os.path.dirname(filename), '.'.join(package_path))
                c_file = _find_c_source(test_basepath)

        py_source_file = None
        if c_file:
            py_source_file = os.path.splitext(c_file)[0] + '.py'
            if not os.path.exists(py_source_file):
                py_source_file = None
            if not is_cython_generated_file(c_file, if_not_found=False):
                if py_source_file and os.path.exists(c_file):
                    # if we did not generate the C file,
                    # then we probably also shouldn't care about the .py file.
                    py_source_file = None
                c_file = None

        return c_file, py_source_file

    def _find_c_source_files(self, dir_path, source_file):
        """
        Desperately parse all C files in the directory or its package parents
        (not re-descending) to find the (included) source file in one of them.
        """
        if not os.path.isdir(dir_path):
            return
        splitext = os.path.splitext
        for filename in os.listdir(dir_path):
            ext = splitext(filename)[1].lower()
            if ext in C_FILE_EXTENSIONS:
                self._read_source_lines(os.path.join(dir_path, filename), source_file)
                if source_file in self._c_files_map:
                    return
        # not found? then try one package up
        if is_package_dir(dir_path):
            self._find_c_source_files(os.path.dirname(dir_path), source_file)

    def _read_source_lines(self, c_file, sourcefile):
        """
        Parse a Cython generated C/C++ source file and find the executable lines.
        Each executable line starts with a comment header that states source file
        and line number, as well as the surrounding range of source code lines.
        """
        if self._parsed_c_files is None:
            self._parsed_c_files = {}
        if c_file in self._parsed_c_files:
            code_lines = self._parsed_c_files[c_file]
        else:
            code_lines = self._parse_cfile_lines(c_file)
            self._parsed_c_files[c_file] = code_lines

        if self._c_files_map is None:
            self._c_files_map = {}

        for filename, code in code_lines.items():
            abs_path = _find_dep_file_path(c_file, filename,
                                           relative_path_search=True)
            self._c_files_map[abs_path] = (c_file, filename, code)

        if sourcefile not in self._c_files_map:
            return (None,) * 2  # e.g. shared library file
        return self._c_files_map[sourcefile][1:]

    def _parse_cfile_lines(self, c_file):
        """
        Parse a C file and extract all source file lines that generated executable code.
        """
        match_source_path_line = re.compile(r' */[*] +"(.*)":([0-9]+)$').match
        match_current_code_line = re.compile(r' *[*] (.*) # <<<<<<+$').match
        match_comment_end = re.compile(r' *[*]/$').match
        match_trace_line = re.compile(r' *__Pyx_TraceLine\(([0-9]+),').match
        not_executable = re.compile(
            r'\s*c(?:type)?def\s+'
            r'(?:(?:public|external)\s+)?'
            r'(?:struct|union|enum|class)'
            r'(\s+[^:]+|)\s*:'
        ).match
        if self._excluded_line_patterns:
            line_is_excluded = re.compile("|".join(["(?:%s)" % regex for regex in self._excluded_line_patterns])).search
        else:
            line_is_excluded = lambda line: False

        code_lines = defaultdict(dict)
        executable_lines = defaultdict(set)
        current_filename = None
        if self._excluded_lines_map is None:
            self._excluded_lines_map = defaultdict(set)

        with open(c_file) as lines:
            lines = iter(lines)
            for line in lines:
                match = match_source_path_line(line)
                if not match:
                    if '__Pyx_TraceLine(' in line and current_filename is not None:
                        trace_line = match_trace_line(line)
                        if trace_line:
                            executable_lines[current_filename].add(int(trace_line.group(1)))
                    continue
                filename, lineno = match.groups()
                current_filename = filename
                lineno = int(lineno)
                for comment_line in lines:
                    match = match_current_code_line(comment_line)
                    if match:
                        code_line = match.group(1).rstrip()
                        if not_executable(code_line):
                            break
                        if line_is_excluded(code_line):
                            self._excluded_lines_map[filename].add(lineno)
                            break
                        code_lines[filename][lineno] = code_line
                        break
                    elif match_comment_end(comment_line):
                        # unexpected comment format - false positive?
                        break

        # Remove lines that generated code but are not traceable.
        for filename, lines in code_lines.items():
            dead_lines = set(lines).difference(executable_lines.get(filename, ()))
            for lineno in dead_lines:
                del lines[lineno]
        return code_lines


class CythonModuleTracer(FileTracer):
    """
    Find the Python/Cython source file for a Cython module.
    """
    def __init__(self, module_file, py_file, c_file, c_files_map, file_path_map):
        super(CythonModuleTracer, self).__init__()
        self.module_file = module_file
        self.py_file = py_file
        self.c_file = c_file
        self._c_files_map = c_files_map
        self._file_path_map = file_path_map

    def has_dynamic_source_filename(self):
        return True

    def dynamic_source_filename(self, filename, frame):
        """
        Determine source file path.  Called by the function call tracer.
        """
        source_file = frame.f_code.co_filename
        try:
            return self._file_path_map[source_file]
        except KeyError:
            pass
        abs_path = _find_dep_file_path(filename, source_file)

        if self.py_file and source_file[-3:].lower() == '.py':
            # always let coverage.py handle this case itself
            self._file_path_map[source_file] = self.py_file
            return self.py_file

        assert self._c_files_map is not None
        if abs_path not in self._c_files_map:
            self._c_files_map[abs_path] = (self.c_file, source_file, None)
        self._file_path_map[source_file] = abs_path
        return abs_path


class CythonModuleReporter(FileReporter):
    """
    Provide detailed trace information for one source file to coverage.py.
    """
    def __init__(self, c_file, source_file, rel_file_path, code, excluded_lines):
        super(CythonModuleReporter, self).__init__(source_file)
        self.name = rel_file_path
        self.c_file = c_file
        self._code = code
        self._excluded_lines = excluded_lines

    def lines(self):
        """
        Return set of line numbers that are possibly executable.
        """
        return set(self._code)

    def excluded_lines(self):
        """
        Return set of line numbers that are excluded from coverage.
        """
        return self._excluded_lines

    def _iter_source_tokens(self):
        current_line = 1
        for line_no, code_line in sorted(self._code.items()):
            while line_no > current_line:
                yield []
                current_line += 1
            yield [('txt', code_line)]
            current_line += 1

    def source(self):
        """
        Return the source code of the file as a string.
        """
        if os.path.exists(self.filename):
            with open_source_file(self.filename) as f:
                return f.read()
        else:
            return '\n'.join(
                (tokens[0][1] if tokens else '')
                for tokens in self._iter_source_tokens())

    def source_token_lines(self):
        """
        Iterate over the source code tokens.
        """
        if os.path.exists(self.filename):
            with open_source_file(self.filename) as f:
                for line in f:
                    yield [('txt', line.rstrip('\n'))]
        else:
            for line in self._iter_source_tokens():
                yield [('txt', line)]


def coverage_init(reg, options):
    plugin = Plugin()
    reg.add_configurer(plugin)
    reg.add_file_tracer(plugin)
