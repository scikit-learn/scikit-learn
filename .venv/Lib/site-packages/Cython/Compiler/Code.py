#
#   Code output module
#


import cython
cython.declare(os=object, re=object, operator=object, textwrap=object,
               Template=object, Naming=object, Options=object, StringEncoding=object,
               Utils=object, SourceDescriptor=object, StringIOTree=object,
               DebugFlags=object, defaultdict=object,
               closing=object, partial=object, wraps=object,
               zlib_compress=object, bz2_compress=object, lzma_compress=object, zstd_compress=object)

import hashlib
import operator
import os
import re
import shutil
import textwrap
from dataclasses import dataclass
from string import Template
from functools import partial, wraps
from contextlib import closing, contextmanager
from collections import defaultdict

from . import Naming
from . import Options
from . import DebugFlags
from . import StringEncoding
from .. import Utils
from .Scanning import SourceDescriptor
from ..StringIOTree import StringIOTree


# Set up available compression algorithms for maximum compression.
from zlib import compress as zlib_compress
try:
    from bz2 import compress as bz2_compress
except ImportError:
    bz2_compress = None
else:
    bz2_compress = partial(bz2_compress, compresslevel=9)
#try:
#    from lzma import compress as lzma_compress
#except ImportError:
#    lzma_compress = None
try:
    from compression.zstd import (
        compress as zstd_compress,
        CompressionParameter as zstd_CompressionParameter,
        Strategy as zstd_Strategy,
    )
except ImportError:
    zstd_compress = None
else:
    zstd_compress = partial(zstd_compress, options={
        zstd_CompressionParameter.strategy: zstd_Strategy.btultra2,
        zstd_CompressionParameter.compression_level: zstd_CompressionParameter.compression_level.bounds()[1],
    })
    del zstd_CompressionParameter
    del zstd_Strategy

compression_algorithms = [
    # Note: order is important and defines values for "CYTHON_COMPRESS_STRINGS" !
    (1, 'zlib', partial(zlib_compress, level=9)),
    (2, 'bz2', bz2_compress),
    (3, 'zstd', zstd_compress),
    # LZMA is difficult to configure for efficient output from C code
    # and the default output tends to be quite large.
    #(4, 'lzma', lzma_compress),
]


renamed_py2_builtins_map = {
    # builtins that had different names in Py2 code
    'unicode'    : 'str',
    'basestring' : 'str',
    'xrange'     : 'range',
    'raw_input'  : 'input',
}

ctypedef_builtins_map = {
    # types of builtins in "ctypedef class" statements which we don't
    # import either because the names conflict with C types or because
    # the type simply is not exposed.
    'py_int'             : '&PyLong_Type',
    'py_long'            : '&PyLong_Type',
    'py_float'           : '&PyFloat_Type',
    'wrapper_descriptor' : '&PyWrapperDescr_Type',
}

basicsize_builtins_map = {
    # builtins whose type has a different tp_basicsize than sizeof(...)
    'PyTypeObject': 'PyHeapTypeObject',
}

# Builtins as of Python version ...
KNOWN_PYTHON_BUILTINS_VERSION = (3, 14, 0, 'beta', 1)
KNOWN_PYTHON_BUILTINS = frozenset([
    'ArithmeticError',
    'AssertionError',
    'AttributeError',
    'BaseException',
    'BaseExceptionGroup',
    'BlockingIOError',
    'BrokenPipeError',
    'BufferError',
    'BytesWarning',
    'ChildProcessError',
    'ConnectionAbortedError',
    'ConnectionError',
    'ConnectionRefusedError',
    'ConnectionResetError',
    'DeprecationWarning',
    'EOFError',
    'Ellipsis',
    'EncodingWarning',
    'EnvironmentError',
    'Exception',
    'ExceptionGroup',
    'False',
    'FileExistsError',
    'FileNotFoundError',
    'FloatingPointError',
    'FutureWarning',
    'GeneratorExit',
    'IOError',
    'ImportError',
    'ImportWarning',
    'IndentationError',
    'IndexError',
    'InterruptedError',
    'IsADirectoryError',
    'KeyError',
    'KeyboardInterrupt',
    'LookupError',
    'MemoryError',
    'ModuleNotFoundError',
    'NameError',
    'None',
    'NotADirectoryError',
    'NotImplemented',
    'NotImplementedError',
    'OSError',
    'OverflowError',
    'PendingDeprecationWarning',
    'PermissionError',
    'ProcessLookupError',
    'PythonFinalizationError',
    'RecursionError',
    'ReferenceError',
    'ResourceWarning',
    'RuntimeError',
    'RuntimeWarning',
    'StopAsyncIteration',
    'StopIteration',
    'SyntaxError',
    'SyntaxWarning',
    'SystemError',
    'SystemExit',
    'TabError',
    'TimeoutError',
    'True',
    'TypeError',
    'UnboundLocalError',
    'UnicodeDecodeError',
    'UnicodeEncodeError',
    'UnicodeError',
    'UnicodeTranslateError',
    'UnicodeWarning',
    'UserWarning',
    'ValueError',
    'Warning',
    'WindowsError',
    'ZeroDivisionError',
    '_IncompleteInputError',
    '__build_class__',
    '__debug__',
    '__import__',
    'abs',
    'aiter',
    'all',
    'anext',
    'any',
    'ascii',
    'bin',
    'bool',
    'breakpoint',
    'bytearray',
    'bytes',
    'callable',
    'chr',
    'classmethod',
    'compile',
    'complex',
    'copyright',
    'credits',
    'delattr',
    'dict',
    'dir',
    'divmod',
    'enumerate',
    'eval',
    'exec',
    'exit',
    'filter',
    'float',
    'format',
    'frozenset',
    'getattr',
    'globals',
    'hasattr',
    'hash',
    'help',
    'hex',
    'id',
    'input',
    'int',
    'isinstance',
    'issubclass',
    'iter',
    'len',
    'license',
    'list',
    'locals',
    'map',
    'max',
    'memoryview',
    'min',
    'next',
    'object',
    'oct',
    'open',
    'ord',
    'pow',
    'print',
    'property',
    'quit',
    'range',
    'repr',
    'reversed',
    'round',
    'set',
    'setattr',
    'slice',
    'sorted',
    'staticmethod',
    'str',
    'sum',
    'super',
    'tuple',
    'type',
    'vars',
    'zip',
])

uncachable_builtins = [
    # Global/builtin names that cannot be cached because they may or may not
    # be available at import time, for various reasons:
    ## Python 3.13+
    '_IncompleteInputError',
    'PythonFinalizationError',
    ## Python 3.11+
    'BaseExceptionGroup',
    'ExceptionGroup',
    ## - Py3.10+
    'aiter',
    'anext',
    'EncodingWarning',
    ## - Py3.7+
    'breakpoint',  # might deserve an implementation in Cython
    ## - platform specific
    'WindowsError',
    ## - others
    '_',  # e.g. used by gettext
]

special_py_methods = cython.declare(frozenset, frozenset((
    '__cinit__', '__dealloc__', '__richcmp__', '__next__',
    '__await__', '__aiter__', '__anext__',
    '__getbuffer__', '__releasebuffer__',
)))

modifier_output_mapper = {
    'inline': 'CYTHON_INLINE'
}.get

cleanup_level_for_type_prefix = cython.declare(object, {
    'ustring': None,
    'tuple': 2,
    'slice': 2,
}.get)


class IncludeCode:
    """
    An include file and/or verbatim C code to be included in the
    generated sources.
    """
    # attributes:
    #
    #  pieces    {order: unicode}: pieces of C code to be generated.
    #            For the included file, the key "order" is zero.
    #            For verbatim include code, the "order" is the "order"
    #            attribute of the original IncludeCode where this piece
    #            of C code was first added. This is needed to prevent
    #            duplication if the same include code is found through
    #            multiple cimports.
    #  location  int: where to put this include in the C sources, one
    #            of the constants INITIAL, EARLY, LATE
    #  order     int: sorting order (automatically set by increasing counter)

    # Constants for location. If the same include occurs with different
    # locations, the earliest one takes precedence.
    INITIAL = 0
    EARLY = 1
    LATE = 2

    counter = 1   # Counter for "order"

    def __init__(self, include=None, verbatim=None, late=True, initial=False):
        self.order = self.counter
        type(self).counter += 1
        self.pieces = {}

        if include:
            if include[0] == '<' and include[-1] == '>':
                self.pieces[0] = '#include {}'.format(include)
                late = False  # system include is never late
            else:
                self.pieces[0] = '#include "{}"'.format(include)

        if verbatim:
            self.pieces[self.order] = verbatim

        if initial:
            self.location = self.INITIAL
        elif late:
            self.location = self.LATE
        else:
            self.location = self.EARLY

    def dict_update(self, d, key):
        """
        Insert `self` in dict `d` with key `key`. If that key already
        exists, update the attributes of the existing value with `self`.
        """
        if key in d:
            other = d[key]
            other.location = min(self.location, other.location)
            other.pieces.update(self.pieces)
        else:
            d[key] = self

    def sortkey(self):
        return self.order

    def mainpiece(self):
        """
        Return the main piece of C code, corresponding to the include
        file. If there was no include file, return None.
        """
        return self.pieces.get(0)

    def write(self, code):
        # Write values of self.pieces dict, sorted by the keys
        for k in sorted(self.pieces):
            code.putln(self.pieces[k])


def get_utility_dir():
    # make this a function and not global variables:
    # http://trac.cython.org/cython_trac/ticket/475
    Cython_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(Cython_dir, "Utility")

read_utilities_hook = None
"""
Override the hook for reading a utilities file that contains code fragments used
by the codegen.

The hook functions takes the path of the utilities file, and returns a list
of strings, one per line.

The default behavior is to open a file relative to get_utility_dir().
"""

def read_utilities_from_utility_dir(path):
    """
    Read all lines of the file at the provided path from a path relative
    to get_utility_dir().
    """
    filename = os.path.join(get_utility_dir(), path)
    with closing(Utils.open_source_file(filename, encoding='UTF-8')) as f:
        return f.readlines()

# by default, read utilities from the utility directory.
read_utilities_hook = read_utilities_from_utility_dir


class AbstractUtilityCode:

    requires = None

    def put_code(self, globalstate: "GlobalState", used_by=None) -> None:
        pass

    def get_tree(self, **kwargs):
        return None

    def get_shared_library_scope(self, **kwargs):
        return None


class UtilityCodeBase(AbstractUtilityCode):
    """
    Support for loading utility code from a file.

    Code sections in the file can be specified as follows:

        ##### MyUtility.proto #####

        [proto declarations]

        ##### MyUtility.init #####

        [code run at module initialization]

        ##### MyUtility #####
        #@requires: MyOtherUtility
        #@substitute: naming

        [definitions]

        ##### MyUtility #####
        #@substitute: tempita

        [requires tempita substitution
         - context can't be specified here though so only
           tempita utility that requires no external context
           will benefit from this tag
         - only necessary when @required from non-tempita code]

    for prototypes and implementation respectively.  For non-python or
    -cython files backslashes should be used instead.  5 to 30 comment
    characters may be used on either side.

    If the @cname decorator is not used and this is a CythonUtilityCode,
    one should pass in the 'name' keyword argument to be used for name
    mangling of such entries.
    """

    is_cython_utility = False
    _utility_cache = {}

    match_section_title = re.compile(
        r'(.+)[.](proto(?:[.]\S+)?|impl|init|cleanup|module_state_decls|module_state_traverse|module_state_clear|export)$'
    ).match

    @staticmethod
    def get_special_comment_matcher(line_comment_char):
        return re.compile((
            # section title
            r'^%(C)s{5,30}  \s*  (?P<name> (?:\w|\.)+ )  \s*  %(C)s{5,30} |'
            # section tags and dependencies
            r'^%(C)s+  @(?P<tag> .+)'
        ) % {'C': re.escape(line_comment_char)}, re.VERBOSE).match

    @classmethod
    def _add_utility(cls, utility, name, type, lines, begin_lineno, tags=None):
        if utility is None:
            return

        code = '\n'.join(lines)
        if tags and 'substitute' in tags and 'naming' in tags['substitute']:
            try:
                new_code = Template(code).substitute(vars(Naming))
            except (KeyError, ValueError) as e:
                raise RuntimeError(
                    f"Error parsing templated utility code '{name}.{type}' at line {begin_lineno:d}: {e}")
            if new_code == code:
                raise RuntimeError(
                    f"Found useless 'substitute: naming' declaration without replacements. ({name}.{type}:{begin_lineno:d})")
            code = new_code

        # remember correct line numbers at least until after templating
        code = '\n' * begin_lineno + code

        if type == 'proto':
            utility[0] = code
        elif type == 'impl':
            utility[1] = code
        else:
            all_tags = utility[2]
            all_tags[type] = code

        if tags:
            all_tags = utility[2]
            for tag_name, tag_values in tags.items():
                all_tags.setdefault(tag_name, set()).update(tag_values)

    @classmethod
    def load_utilities_from_file(cls, path):
        utilities = cls._utility_cache.get(path)
        if utilities:
            return utilities

        _, ext = os.path.splitext(path)
        if ext in ('.pyx', '.py', '.pxd', '.pxi'):
            comment = '#'
            strip_comments = partial(re.compile(r'^\s*#(?!\s*cython\s*:).*').sub, '')
            rstrip = str.rstrip
        else:
            comment = '/'
            strip_comments = partial(re.compile(r'^\s*//.*|/\*[^*]*\*/').sub, '')
            rstrip = partial(re.compile(r'\s+(\\?)$').sub, r'\1')

        match_special = cls.get_special_comment_matcher(comment)
        match_type = cls.match_section_title

        all_lines = read_utilities_hook(path)

        utilities = defaultdict(lambda: [None, None, {}])
        lines = []
        tags = defaultdict(set)
        utility = name = type = None
        begin_lineno = 0

        for lineno, line in enumerate(all_lines):
            m = match_special(line)
            if m is None:
                lines.append(rstrip(strip_comments(line)))
            elif m.group('name'):
                cls._add_utility(utility, name, type, lines, begin_lineno, tags)

                begin_lineno = lineno + 1
                del lines[:]
                tags.clear()

                name = m.group('name')
                mtype = match_type(name)
                if mtype:
                    name, type = mtype.groups()
                else:
                    type = 'impl'
                utility = utilities[name]
            else:
                tag_value = m.group('tag')
                if ':' not in tag_value:
                    raise RuntimeError(f"Found invalid tag '{tag_value}' in utility section {name}.{type}")

                tag_name, _, tag_value = tag_value.partition(':')
                tag_name = tag_name.rstrip()
                tag_value = tag_value.strip()

                if tag_name not in ('requires', 'substitute', 'proto_block'):
                    raise RuntimeError(f"Found unknown tag name '{tag_name}' in utility section {name}.{type}")
                if not re.match(r'\S+$', tag_value):
                    raise RuntimeError(f"Found invalid tag value '{tag_value}' in utility section {name}.{type}")

                tags[tag_name].add(tag_value)
                lines.append('')  # keep line number correct

        if utility is None:
            raise ValueError("Empty utility code file")

        # Don't forget to add the last utility code
        cls._add_utility(utility, name, type, lines, begin_lineno, tags)

        utilities = dict(utilities)  # un-defaultdict-ify
        cls._utility_cache[path] = utilities
        return utilities

    @classmethod
    def load(cls, util_code_name, from_file, **kwargs):
        """
        Load utility code from a file specified by from_file (relative to
        Cython/Utility) and name util_code_name.
        """

        if '::' in util_code_name:
            from_file, util_code_name = util_code_name.rsplit('::', 1)
        assert from_file
        utilities = cls.load_utilities_from_file(from_file)
        proto, impl, tags = utilities[util_code_name]

        if tags:
            if "substitute" in tags and "tempita" in tags["substitute"]:
                if not issubclass(cls, TempitaUtilityCode):
                    return TempitaUtilityCode.load(util_code_name, from_file, **kwargs)
            orig_kwargs = kwargs.copy()
            for name, values in tags.items():
                if name in kwargs:
                    continue
                # only pass lists when we have to: most argument expect one value or None
                if name == 'requires':
                    if orig_kwargs:
                        values = [cls.load(dep, from_file, **orig_kwargs)
                                  for dep in sorted(values)]
                    else:
                        # dependencies are rarely unique, so use load_cached() when we can
                        values = [cls.load_cached(dep, from_file)
                                  for dep in sorted(values)]
                elif name == 'substitute':
                    # don't want to pass "naming" or "tempita" to the constructor
                    # since these will have been handled
                    values = values - {'naming', 'tempita'}
                    if not values:
                        continue
                elif not values:
                    values = None
                elif len(values) == 1:
                    values = list(values)[0]
                kwargs[name] = values

        if proto is not None:
            kwargs['proto'] = proto
        if impl is not None:
            kwargs['impl'] = impl

        if 'name' not in kwargs:
            kwargs['name'] = util_code_name

        if 'file' not in kwargs and from_file:
            kwargs['file'] = from_file
        return cls(**kwargs)

    @classmethod
    def load_cached(cls, utility_code_name, from_file, __cache={}):
        """
        Calls .load(), but using a per-type cache based on utility name and file name.
        """
        key = (utility_code_name, from_file, cls)
        try:
            return __cache[key]
        except KeyError:
            pass
        code = __cache[key] = cls.load(utility_code_name, from_file)
        return code

    @classmethod
    def load_as_string(cls, util_code_name, from_file, include_requires=False, **kwargs):
        """
        Load a utility code as a string. Returns (proto, implementation).

        If 'include_requires=True', concatenates all requirements before the actually
        requested utility code, separately for proto and impl part.

        In a lot of cases it may be better to use regular "load" and "CCodeWriter.put_code_here"
        since that is able to apply the code transformations to the code too.
        """
        util = cls.load(util_code_name, from_file, **kwargs)

        if not include_requires:
            return (util.format_code(util.proto),
                    util.format_code(util.impl))

        protos, impls = [], []
        def prepend(util_code):
            if util_code.requires:
                for dep in util_code.requires:
                    prepend(dep)
            if util_code.proto:
                protos.append(util_code.format_code(util_code.proto))
            if util_code.impl:
                impls.append(util_code.format_code(util_code.impl))

        prepend(util)
        return "".join(protos), "".join(impls)

    def format_code(self, code_string, replace_empty_lines=re.compile(r'\n\n+').sub):
        """
        Format a code section for output.
        """
        if code_string:
            code_string = replace_empty_lines('\n', code_string.strip()) + '\n\n'
        return code_string

    def __repr__(self):
        return "<%s(%s)>" % (type(self).__name__, self.name)

    def get_tree(self, **kwargs):
        return None

    def get_shared_library_scope(self, **kwargs):
        return None

    def __deepcopy__(self, memodict=None):
        # No need to deep-copy utility code since it's essentially immutable.
        return self

@dataclass
class SharedFunctionDecl:
    """Contains parsed declaration of shared utility function"""
    name: str
    ret: str
    params: str

class UtilityCode(UtilityCodeBase):
    """
    Stores utility code to add during code generation.

    See GlobalState.put_utility_code.

    hashes/equals by instance

    proto           C prototypes
    export          C prototypes exported from the shared utility code module
    impl            implementation code
    init            code to call on module initialization
    requires        utility code dependencies
    proto_block     the place in the resulting file where the prototype should
                    end up
    name            name of the utility code (or None)
    file            filename of the utility code file this utility was loaded
                    from (or None)
    shared_utility_functions        List of parsed declaration line of the shared utility function
    """
    code_parts = ["proto", "export", "impl", "init", "cleanup", "module_state_decls", "module_state_traverse", "module_state_clear"]

    def __init__(self, proto=None, impl=None, init=None, cleanup=None,
                 module_state_decls=None, module_state_traverse=None,
                 module_state_clear=None, requires=None,
                 proto_block='utility_code_proto', name=None, file=None, export=None):
        # proto_block: Which code block to dump prototype in. See GlobalState.
        self.proto = proto
        self.impl = impl
        self.init = init
        self.cleanup = cleanup
        self.module_state_decls = module_state_decls
        self.module_state_traverse = module_state_traverse
        self.module_state_clear = module_state_clear
        self.requires = requires
        self._cache = {}
        self.specialize_list = []
        self.proto_block = proto_block
        self.name = name
        self.file = file
        self.export = export
        self.shared_utility_functions = self.parse_export_functions(export) if export else []
        if export:
            self._validate_suitable_for_sharing()

        # cached for use in hash and eq
        self._parts_tuple = tuple(getattr(self, part, None) for part in self.code_parts)

    def parse_export_functions(self, export_proto: str) -> list:

        assert '//' not in export_proto and '/*' not in export_proto and '*/' not in export_proto, \
            f'Export block must not contain comments:\n{export_proto.strip()}\n in file {self.file}'

        parsed_protos = []
        proto_regex=r'''
            ^static\s                                         # `static` keyword
            (?P<ret_type>[^;()]+[\s*])                        # return type + modifier with optional * - e.g.: int *, float, const str *, ...
            (?P<func_name>\w+)\((?P<func_params>[^)]*)\)$     # function with params - e.g. foo(int, float, *PyObject)
        '''

        for proto in export_proto.split(';\n'):
            proto = proto.strip().replace('\n', '')
            proto = re.sub(r'\s+', ' ', proto)

            if len(proto) == 0:
                continue
            matched = re.match(proto_regex, proto, re.VERBOSE)
            assert matched is not None, \
                f"Wrong format of function definition in export block \n{proto!r}\n in {self.file}"

            ret_type, func_name, func_params = matched.groups()
            parsed_protos.append(
                SharedFunctionDecl(name=func_name.strip(), ret=ret_type.strip(), params=func_params.strip())
            )

        return parsed_protos


    def __hash__(self):
        return hash(self._parts_tuple)

    def __eq__(self, other):
        if self is other:
            return True
        self_type, other_type = type(self), type(other)
        if self_type is not other_type and not (isinstance(other, self_type) or isinstance(self, other_type)):
            return False

        return self._parts_tuple == other._parts_tuple

    def none_or_sub(self, s, context):
        """
        Format a string in this utility code with context. If None, do nothing.
        """
        if s is None:
            return None
        return s % context

    def specialize(self, pyrex_type=None, **data):
        name = self.name
        if pyrex_type is not None:
            data['type'] = pyrex_type.empty_declaration_code()
            data['type_name'] = pyrex_type.specialization_name()
            name = "%s[%s]" % (name, data['type_name'])
        # Dicts aren't hashable...
        key = tuple(sorted(data.items()))
        try:
            return self._cache[key]
        except KeyError:
            if self.requires is None:
                requires = None
            else:
                requires = [r.specialize(data) for r in self.requires]

            s = self._cache[key] = UtilityCode(
                self.none_or_sub(self.proto, data),
                self.none_or_sub(self.impl, data),
                self.none_or_sub(self.init, data),
                self.none_or_sub(self.cleanup, data),
                self.none_or_sub(self.module_state_decls, data),
                self.none_or_sub(self.module_state_traverse, data),
                self.none_or_sub(self.module_state_clear, data),
                requires,
                self.proto_block,
                name,
            )

            self.specialize_list.append(s)
            return s

    def _validate_suitable_for_sharing(self):
        code_string = getattr(self, "impl")
        if not code_string: return
        assert "NAMED_CGLOBAL(moddict_cname)" not in code_string, \
            f"moddict_cname should not be shared: {self}"

    @cython.final
    def _put_code_section(self, writer: "CCodeWriter", output: "GlobalState", code_type: str, used_by=None):
        code_string = getattr(self, code_type)
        if not code_string:
            return

        can_be_reused = code_type in ('proto', 'impl')

        code_string, result_is_module_specific = process_utility_ccode(self, output, code_string)

        used_by = f" (used by {used_by})" if used_by else ''
        name = f"{self.name}.{code_type}" if code_type != 'impl' else self.name

        writer.putln(f"/* {name}{used_by} */")

        if can_be_reused and not result_is_module_specific:
            # can be reused across modules
            writer.put_or_include(code_string, f'{self.name}_{code_type}')
        else:
            writer.put_multilines(code_string)

    def _put_init_code_section(self, output):
        if not self.init:
            return
        writer = output['init_globals']
        self._put_code_section(writer, output, 'init')
        # 'init' code can end with an 'if' statement for an error condition like:
        # if (check_ok()) ; else
        writer.putln("  " + writer.error_goto_if_PyErr(output.module_pos))
        writer.putln()

    def _put_shared_function_declarations(self, code: "CCodeWriter") -> None:
        code.putln(f'/* {self.name} */')
        for shared in self.shared_utility_functions:
            # Convert function declarations to static function pointers.
            code.putln(f'static {shared.ret}(*{shared.name})({shared.params}); /*proto*/')
        code.putln()

    def put_code(self, globalstate: "GlobalState", used_by=None) -> None:
        has_shared_utility_code = bool(
            self.shared_utility_functions and globalstate.module_node.scope.context.shared_utility_qualified_name
        )

        if self.requires and not has_shared_utility_code:
            for dependency in self.requires:
                globalstate.use_utility_code(dependency, used_by=self.name)

        if has_shared_utility_code:
            self._put_shared_function_declarations(globalstate[self.proto_block])
        globalstate.shared_utility_functions.extend(self.shared_utility_functions)

        if self.proto:
            self._put_code_section(globalstate[self.proto_block], globalstate, 'proto', used_by=used_by)
        if not has_shared_utility_code:
            self._put_code_section(globalstate[self.proto_block], globalstate, 'export')
        if self.impl and not has_shared_utility_code:
            self._put_code_section(globalstate['utility_code_def'], globalstate, 'impl', used_by=used_by)
        if self.cleanup and Options.generate_cleanup_code:
            self._put_code_section(globalstate['cleanup_globals'], globalstate, 'cleanup')
        if self.module_state_decls:
            self._put_code_section(globalstate['module_state_contents'], globalstate, 'module_state_decls')
        if self.module_state_traverse:
            self._put_code_section(globalstate['module_state_traverse_contents'], globalstate, 'module_state_traverse')
        if self.module_state_clear:
            self._put_code_section(globalstate['module_state_clear_contents'], globalstate, 'module_state_clear')

        if self.init:
            self._put_init_code_section(globalstate)


def add_macro_processor(*macro_names, regex=None, is_module_specific=False, _last_macro_processor = [None]):
    """Decorator to chain the code macro processors below.
    """
    last_processor = _last_macro_processor[0]

    def build_processor(func):
        @wraps(func)
        def process(utility_code: UtilityCode, output, code_string: str):
            # First, call the processing chain in FIFO function definition order.
            result_is_module_specific = False
            if last_processor is not None:
                code_string, result_is_module_specific = last_processor(utility_code, output, code_string)

            # Detect if we need to do something.
            if macro_names:
                for macro in macro_names:
                    if macro in code_string:
                        break
                else:
                    return code_string, result_is_module_specific

            # Process the code.
            if regex is None:
                code_string = func(utility_code, output, code_string)
            else:
                code_string = re.sub(regex, partial(func, output), code_string)

            # Make sure we found and replaced all macro occurrences.
            for macro in macro_names:
                if macro in code_string:
                    raise RuntimeError(f"Left-over utility code macro '{macro}()' found in '{utility_code.name}'")

            result_is_module_specific |= is_module_specific
            return code_string, result_is_module_specific

        _last_macro_processor[0] = process
        return process

    return build_processor


@add_macro_processor(
    'CSTRING',
    regex=r'CSTRING\(\s*"""([^"]*(?:"[^"]+)*)"""\s*\)',
)
def _wrap_c_string(_, matchobj):
    """Replace CSTRING('''xyz''') by a C compatible string, taking care of line breaks.
    """
    content = matchobj.group(1).replace('"', r'\042')
    return ''.join(
        f'"{line}\\n"\n' if not line.endswith('\\') or line.endswith('\\\\') else f'"{line[:-1]}"\n'
        for line in content.splitlines())


@add_macro_processor()
def _format_impl_code(utility_code: UtilityCode, _, impl):
    return utility_code.format_code(impl)


@add_macro_processor(
    'CALL_UNBOUND_METHOD',
    is_module_specific=True,
    regex=(
        r'CALL_UNBOUND_METHOD\('
        r'([a-zA-Z_]+),\s*'   # type cname
        r'"([^"]+)",\s*'      # method name
        r'([^),\s]+)'         # object cname
        r'((?:,[^),]+)*)'     # args*
        r'\)'
    ),
)
def _inject_unbound_method(output, matchobj):
    """Replace 'UNBOUND_METHOD(type, "name")' by a constant Python identifier cname.
    """
    type_cname, method_name, obj_cname, args = matchobj.groups()
    type_cname = '&%s' % type_cname
    args = [arg.strip() for arg in args[1:].split(',')] if args else []
    assert len(args) < 3, f"CALL_UNBOUND_METHOD() does not support {len(args):d} call arguments"
    return output.cached_unbound_method_call_code(
        f"{Naming.modulestateglobal_cname}->",
        obj_cname, type_cname, method_name, args)


@add_macro_processor(
    'PYIDENT', 'PYUNICODE',
    is_module_specific=True,
    regex=r'PY(IDENT|UNICODE)\("([^"]+)"\)',
)
def _inject_string_constant(output, matchobj):
    """Replace 'PYIDENT("xyz")' by a constant Python identifier cname.
    """
    str_type, name = matchobj.groups()
    return "%s->%s" % (
        Naming.modulestateglobal_cname,
        output.get_py_string_const(
            StringEncoding.EncodedString(name), identifier=str_type == 'IDENT').cname)


@add_macro_processor(
    'EMPTY',
    # As long as we use the same C access macros for these names, they are not module specific.
    # is_module_specific=True,
    regex=r'EMPTY\((bytes|unicode|tuple)\)',
)
def _inject_empty_collection_constant(output, matchobj):
    """Replace 'EMPTY(bytes|tuple|...)' by a constant Python identifier cname.
    """
    type_name = matchobj.group(1)
    return "%s->%s" % (
        Naming.modulestateglobal_cname,
        getattr(Naming, f'empty_{type_name}'))


@add_macro_processor(
    'CGLOBAL',  # 'NAMED_CGLOBAL',  # first is part of second and thus not needed
    is_module_specific=False,
    regex=r'(NAMED_)?CGLOBAL\(([^)]+)\)',
)
def _inject_cglobal(output, matchobj):
    is_named, name = matchobj.groups()
    if is_named:
        name = getattr(Naming, name)
    return f"{Naming.modulestateglobal_cname}->{name}"


@add_macro_processor()
def process_utility_ccode(utility_code, _, code_string):
    """Entry point for code processors, must be defined last.
    """
    return code_string


def sub_tempita(s, context, file=None, name=None, __cache={}):
    "Run tempita on string s with given context."
    if not s:
        return None

    if file:
        name = f"{file}:{name}"
    if name:
        context['__name'] = name

    try:
        template = __cache[s]
    except KeyError:
        from ..Tempita import Template
        template = __cache[s] = Template(s, name=name)

    return template.substitute(context)


class TempitaUtilityCode(UtilityCode):
    def __init__(self, name=None, proto=None, impl=None, init=None, file=None, context=None, **kwargs):
        if context is None:
            context = {}
        else:
            # prevent changes propagating back if context is shared between multiple utility codes.
            context = context.copy()
        proto = sub_tempita(proto, context, file, name)
        impl = sub_tempita(impl, context, file, name)
        init = sub_tempita(init, context, file, name)
        super().__init__(
            proto, impl, init=init, name=name, file=file, **kwargs)

    @classmethod
    def load_cached(cls, utility_code_name, from_file=None, context=None, __cache={}):
        context_key = tuple(sorted(context.items())) if context else None
        assert hash(context_key) is not None  # raise TypeError if not hashable
        key = (cls, from_file, utility_code_name, context_key)
        try:
            return __cache[key]
        except KeyError:
            pass
        code = __cache[key] = cls.load(utility_code_name, from_file, context=context)
        return code

    def none_or_sub(self, s, context):
        """
        Format a string in this utility code with context. If None, do nothing.
        """
        if s is None:
            return None
        return sub_tempita(s, context, self.file, self.name)


class LazyUtilityCode(UtilityCodeBase):
    """
    Utility code that calls a callback with the root code writer when
    available. Useful when you only have 'env' but not 'code'.
    """
    __name__ = '<lazy>'
    requires = None

    def __init__(self, callback):
        self.callback = callback

    def put_code(self, globalstate: "GlobalState", used_by=None) -> None:
        utility = self.callback(globalstate.rootwriter)
        globalstate.use_utility_code(utility, used_by=used_by)


class FunctionState:
    # return_label     string          function return point label
    # error_label      string          error catch point label
    # error_without_exception  boolean Can go to the error label without an exception (e.g. __next__ can return NULL)
    # continue_label   string          loop continue point label
    # break_label      string          loop break point label
    # return_from_error_cleanup_label string
    # label_counter    integer         counter for naming labels
    # exc_vars         (string * 3)    exception variables for reraise, or None
    # can_trace        boolean         line tracing is supported in the current context
    # scope            Scope           the scope object of the current function

    # Not used for now, perhaps later
    def __init__(self, owner, names_taken=set(), scope=None):
        self.names_taken = names_taken
        self.owner = owner
        self.scope = scope

        self.error_label = None
        self.label_counter = 0
        self.labels_used = set()
        self.return_label = self.new_label()
        self.new_error_label()
        self.continue_label = None
        self.break_label = None
        self.yield_labels = []

        self.exc_vars = None
        self.current_except = None
        self.can_trace = False
        self.gil_owned = True

        self.temps_allocated = []  # of (name, type, manage_ref, static)
        self.temps_free = {}  # (type, manage_ref) -> list of free vars with same type/managed status
        self.temps_used_type = {}  # name -> (type, manage_ref)
        self.zombie_temps = set()  # temps that must not be reused after release
        self.temp_counter = 0
        self.closure_temps = None

        # This is used to collect temporaries, useful to find out which temps
        # need to be privatized in parallel sections
        self.collect_temps_stack = []

        # This is used for the error indicator, which needs to be local to the
        # function. It used to be global, which relies on the GIL being held.
        # However, exceptions may need to be propagated through 'nogil'
        # sections, in which case we introduce a race condition.
        self.should_declare_error_indicator = False
        self.uses_error_indicator = False

        self.error_without_exception = False

        self.needs_refnanny = False

    # safety checks

    def validate_exit(self):
        # validate that all allocated temps have been freed
        if self.temps_allocated:
            leftovers = self.temps_in_use()
            if leftovers:
                msg = "TEMPGUARD: Temps left over at end of '%s': %s" % (self.scope.name, ', '.join([
                    '%s [%s]' % (name, ctype)
                    for name, ctype, is_pytemp in sorted(leftovers)]),
                )
                #print(msg)
                raise RuntimeError(msg)

    # labels

    def new_label(self, name=None):
        n: cython.size_t = self.label_counter
        self.label_counter = n + 1
        label = "%s%d" % (Naming.label_prefix, n)
        if name is not None:
            label += '_' + name
        return label

    def new_yield_label(self, expr_type='yield'):
        label = self.new_label('resume_from_%s' % expr_type)
        num_and_label = (len(self.yield_labels) + 1, label)
        self.yield_labels.append(num_and_label)
        return num_and_label

    def new_error_label(self, prefix=""):
        old_err_lbl = self.error_label
        self.error_label = self.new_label(prefix + 'error')
        return old_err_lbl

    def get_loop_labels(self):
        return (
            self.continue_label,
            self.break_label)

    def set_loop_labels(self, labels):
        (self.continue_label,
         self.break_label) = labels

    def new_loop_labels(self, prefix=""):
        old_labels = self.get_loop_labels()
        self.set_loop_labels(
            (self.new_label(prefix + "continue"),
             self.new_label(prefix + "break")))
        return old_labels

    def get_all_labels(self):
        return (
            self.continue_label,
            self.break_label,
            self.return_label,
            self.error_label)

    def set_all_labels(self, labels):
        (self.continue_label,
         self.break_label,
         self.return_label,
         self.error_label) = labels

    def all_new_labels(self):
        old_labels = self.get_all_labels()
        new_labels = []
        for old_label, name in zip(old_labels, ['continue', 'break', 'return', 'error']):
            if old_label:
                new_labels.append(self.new_label(name))
            else:
                new_labels.append(old_label)
        self.set_all_labels(new_labels)
        return old_labels

    def use_label(self, lbl):
        self.labels_used.add(lbl)

    def label_used(self, lbl):
        return lbl in self.labels_used

    # temp handling

    def allocate_temp(self, type, manage_ref, static=False, reusable=True):
        """
        Allocates a temporary (which may create a new one or get a previously
        allocated and released one of the same type). Type is simply registered
        and handed back, but will usually be a PyrexType.

        If type.needs_refcounting, manage_ref comes into play. If manage_ref is set to
        True, the temp will be decref-ed on return statements and in exception
        handling clauses. Otherwise the caller has to deal with any reference
        counting of the variable.

        If not type.needs_refcounting, then manage_ref will be ignored, but it
        still has to be passed. It is recommended to pass False by convention
        if it is known that type will never be a reference counted type.

        static=True marks the temporary declaration with "static".
        This is only used when allocating backing store for a module-level
        C array literals.

        if reusable=False, the temp will not be reused after release.

        A C string referring to the variable is returned.
        """
        if type.is_cv_qualified and not type.is_reference:
            type = type.cv_base_type
        elif type.is_reference and not type.is_fake_reference:
            type = type.ref_base_type
        elif type.is_cfunction:
            from . import PyrexTypes
            type = PyrexTypes.c_ptr_type(type)  # A function itself isn't an l-value
        elif type.is_cpp_class and not type.is_fake_reference and self.scope.directives['cpp_locals']:
            self.scope.use_utility_code(UtilityCode.load_cached("OptionalLocals", "CppSupport.cpp"))
        if not type.needs_refcounting:
            # Make manage_ref canonical, so that manage_ref will always mean
            # a decref is needed.
            manage_ref = False

        freelist = self.temps_free.get((type, manage_ref))
        if reusable and freelist is not None and freelist[0]:
            result = freelist[0].pop()
            freelist[1].remove(result)
        else:
            while True:
                self.temp_counter += 1
                result = "%s%d" % (Naming.codewriter_temp_prefix, self.temp_counter)
                if result not in self.names_taken: break
            self.temps_allocated.append((result, type, manage_ref, static))
            if not reusable:
                self.zombie_temps.add(result)
        self.temps_used_type[result] = (type, manage_ref)
        if DebugFlags.debug_temp_code_comments:
            self.owner.putln("/* %s allocated (%s)%s */" % (result, type, "" if reusable else " - zombie"))

        if self.collect_temps_stack:
            self.collect_temps_stack[-1].add((result, type))

        return result

    def release_temp(self, name):
        """
        Releases a temporary so that it can be reused by other code needing
        a temp of the same type.
        """
        type, manage_ref = self.temps_used_type[name]
        freelist = self.temps_free.get((type, manage_ref))
        if freelist is None:
            freelist = ([], set())  # keep order in list and make lookups in set fast
            self.temps_free[(type, manage_ref)] = freelist
        if name in freelist[1]:
            raise RuntimeError("Temp %s freed twice!" % name)
        if name not in self.zombie_temps:
            freelist[0].append(name)
        freelist[1].add(name)
        if DebugFlags.debug_temp_code_comments:
            self.owner.putln("/* %s released %s*/" % (
                name, " - zombie" if name in self.zombie_temps else ""))

    def temps_in_use(self):
        """Return a list of (cname,type,manage_ref) tuples of temp names and their type
        that are currently in use.
        """
        used = []
        for name, type, manage_ref, static in self.temps_allocated:
            freelist = self.temps_free.get((type, manage_ref))
            if freelist is None or name not in freelist[1]:
                used.append((name, type, manage_ref and type.needs_refcounting))
        return used

    def temps_holding_reference(self):
        """Return a list of (cname,type) tuples of temp names and their type
        that are currently in use. This includes only temps
        with a reference counted type which owns its reference.
        """
        return [(name, type)
                for name, type, manage_ref in self.temps_in_use()
                if manage_ref and type.needs_refcounting]

    def all_managed_temps(self):
        """Return a list of (cname, type) tuples of refcount-managed Python objects.
        """
        return [(cname, type)
                for cname, type, manage_ref, static in self.temps_allocated
                if manage_ref]

    def all_free_managed_temps(self):
        """Return a list of (cname, type) tuples of refcount-managed Python
        objects that are not currently in use.  This is used by
        try-except and try-finally blocks to clean up temps in the
        error case.
        """
        return sorted([  # Enforce deterministic order.
            (cname, type)
            for (type, manage_ref), freelist in self.temps_free.items() if manage_ref
            for cname in freelist[0]
        ])

    def start_collecting_temps(self):
        """
        Useful to find out which temps were used in a code block
        """
        self.collect_temps_stack.append(set())

    def stop_collecting_temps(self):
        return self.collect_temps_stack.pop()

    def init_closure_temps(self, scope):
        self.closure_temps = ClosureTempAllocator(scope)


class NumConst:
    """Global info about a Python number constant held by GlobalState.

    cname       string
    value       string
    py_type     string     int, long, float
    value_code  string     evaluation code if different from value
    """

    def __init__(self, cname, value, py_type, value_code=None):
        self.cname = cname
        self.value = value
        self.py_type = py_type
        self.value_code = value_code or value


class PyObjectConst:
    """Global info about a generic constant held by GlobalState.
    """
    # cname       string
    # type        PyrexType

    def __init__(self, cname, type):
        self.cname = cname
        self.type = type


cython.declare(possible_unicode_identifier=object, possible_bytes_identifier=object,
               replace_identifier=object, find_alphanums=object)
possible_unicode_identifier = re.compile(r"(?![0-9])\w+$", re.U).match
possible_bytes_identifier = re.compile(br"(?![0-9])\w+$").match
replace_identifier = re.compile(r'[^a-zA-Z0-9_]+').sub
find_alphanums = re.compile('([a-zA-Z0-9]+)').findall

class StringConst:
    """Global info about a C string constant held by GlobalState.
    """
    # cname            string
    # text             EncodedString or BytesLiteral
    # escaped_value    str        The string value as C code byte sequence.
    # py_strings       {(identifier, encoding) : PyStringConst}
    # c_used           boolean  Is the plain C string used (or only the Python object?)

    def __init__(self, cname, text, byte_string):
        self.cname = cname
        self.text = text
        self.escaped_value = StringEncoding.escape_byte_string(byte_string)
        self.py_strings = None
        self.c_used = False

    def get_py_string_const(self, encoding, identifier=None):
        text = self.text
        intern: cython.bint
        is_unicode: cython.bint

        if identifier or encoding is None:
            # unicode string
            encoding = encoding_key = None
            is_unicode = True
        else:
            # bytes
            is_unicode = False
            encoding = encoding.lower()
            if encoding in ('utf8', 'utf-8', 'ascii', 'usascii', 'us-ascii'):
                encoding = None
                encoding_key = None
            else:
                encoding_key = ''.join(find_alphanums(encoding))

        if identifier:
            intern = True
        elif identifier is None:
            if isinstance(text, bytes):
                intern = bool(possible_bytes_identifier(text))
            else:
                intern = bool(possible_unicode_identifier(text))
        else:
            intern = False

        key = (intern, is_unicode, encoding_key)
        if self.py_strings is None:
            self.py_strings = {}
        else:
            try:
                return self.py_strings[key]
            except KeyError:
                pass

        pystring_cname = (
            f"{Naming.interned_prefixes['str'] if intern else Naming.py_const_prefix}"
            f"{'u' if is_unicode else 'b'}"
            f"{'_' + encoding_key if encoding_key else ''}"
            f"_{self.cname[len(Naming.const_prefix):]}"
        )

        py_string = PyStringConst(pystring_cname, encoding, intern, is_unicode)
        self.py_strings[key] = py_string
        return py_string


class PyStringConst:
    """Global info about a Python string constant held by GlobalState.
    """
    # cname       string
    # encoding    string
    # intern      boolean
    # is_unicode  boolean

    def __init__(self, cname, encoding, intern=False, is_unicode=False):
        self.cname = cname
        self.encoding = encoding
        self.is_unicode = is_unicode
        self.intern = intern

    def __lt__(self, other):
        return self.cname < other.cname


class GlobalState:
    # filename_table   {string : int}  for finding filename table indexes
    # filename_list    [string]        filenames in filename table order
    # input_file_contents dict         contents (=list of lines) of any file that was used as input
    #                                  to create this output C code.  This is
    #                                  used to annotate the comments.
    #
    # utility_codes   set                IDs of used utility code (to avoid reinsertion)
    #
    # declared_cnames  {string:Entry}  used in a transition phase to merge pxd-declared
    #                                  constants etc. into the pyx-declared ones (i.e,
    #                                  check if constants are already added).
    #                                  In time, hopefully the literals etc. will be
    #                                  supplied directly instead.
    #
    # const_cnames_used  dict          global counter for unique constant identifiers
    # shared_utility_functions         List of parsed declaration lines of the shared utility functions

    # parts            {string:CCodeWriter}


    # interned_strings
    # consts
    # interned_nums

    # directives       set             Temporary variable used to track
    #                                  the current set of directives in the code generation
    #                                  process.

    directives = {}

    code_layout = [
        'h_code',
        'filename_table',
        'utility_code_proto_before_types',
        'numeric_typedefs',           # Let these detailed individual parts stay!,
        'complex_type_declarations',  # as the proper solution is to make a full DAG...
        'type_declarations',          # More coarse-grained blocks would simply hide
        'utility_code_proto',         # the ugliness, not fix it
        'module_declarations',
        'typeinfo',
        'before_global_var',
        'global_var',
        'string_decls',
        'decls',
        'late_includes',
        'module_state',
        'module_state_contents',  # can be used to inject declarations into the modulestate struct
        'module_state_end',
        'constant_name_defines',
        'module_state_clear',
        'module_state_clear_contents',
        'module_state_clear_end',
        'module_state_traverse',
        'module_state_traverse_contents',
        'module_state_traverse_end',
        'module_code',  # user code goes here
        'module_exttypes',
        'initfunc_declarations',
        'init_module',
        'pystring_table',
        'cached_builtins',
        'cached_constants',
        'init_constants',
        'init_codeobjects',
        'init_globals',  # (utility code called at init-time)
        'cleanup_globals',
        'cleanup_module',
        'main_method',
        'utility_code_pragmas',  # silence some irrelevant warnings in utility code
        'utility_code_def',
        'utility_code_pragmas_end',  # clean-up the utility_code_pragmas
        'end'
    ]

    # h files can only have a much smaller list of sections
    h_code_layout = [
        'h_code',
        'utility_code_proto_before_types',
        'type_declarations',
        'utility_code_proto',
        'end'
    ]

    def __init__(self, writer, module_node, code_config, common_utility_include_dir=None):
        self.filename_table = {}
        self.filename_list = []
        self.input_file_contents = {}
        self.utility_codes = set()
        self.declared_cnames = {}
        self.in_utility_code_generation = False
        self.code_config = code_config
        self.common_utility_include_dir = common_utility_include_dir
        self.parts = {}
        self.module_node = module_node  # because some utility code generation needs it
                                        # (generating backwards-compatible Get/ReleaseBuffer

        self.const_cnames_used = {}
        self.string_const_index = {}
        self.dedup_const_index = {}
        self.pyunicode_ptr_const_index = {}
        self.codeobject_constants = []
        self.num_const_index = {}
        self.arg_default_constants = []
        self.const_array_counters = {}  # counts of differently prefixed arrays of constants
        self.cached_cmethods = {}
        self.initialised_constants = set()
        self.shared_utility_functions = []

        writer.set_global_state(self)
        self.rootwriter = writer

    def initialize_main_c_code(self):
        rootwriter = self.rootwriter
        for i, part in enumerate(self.code_layout):
            w = self.parts[part] = rootwriter.insertion_point()
            if i > 0:
                w.putln("/* #### Code section: %s ### */" % part)

        w = self.parts['cached_builtins']
        w.start_initcfunc(
            "int __Pyx_InitCachedBuiltins("
            f"{Naming.modulestatetype_cname} *{Naming.modulestatevalue_cname})")
        w.putln(f"CYTHON_UNUSED_VAR({Naming.modulestatevalue_cname});")

        w = self.parts['cached_constants']
        w.start_initcfunc(
            "int __Pyx_InitCachedConstants("
            f"{Naming.modulestatetype_cname} *{Naming.modulestatevalue_cname})",
            refnanny=True)
        w.putln(f"CYTHON_UNUSED_VAR({Naming.modulestatevalue_cname});")
        w.put_setup_refcount_context(StringEncoding.EncodedString("__Pyx_InitCachedConstants"))

        w = self.parts['init_globals']
        w.start_initcfunc("int __Pyx_InitGlobals(void)")

        w = self.parts['init_constants']
        w.start_initcfunc(
            "int __Pyx_InitConstants("
            f"{Naming.modulestatetype_cname} *{Naming.modulestatevalue_cname})")
        w.putln(f"CYTHON_UNUSED_VAR({Naming.modulestatevalue_cname});")

        if not Options.generate_cleanup_code:
            del self.parts['cleanup_globals']
        else:
            w = self.parts['cleanup_globals']
            w.start_initcfunc(
                "void __Pyx_CleanupGlobals("
                f"{Naming.modulestatetype_cname} *{Naming.modulestatevalue_cname})")
            w.putln(f"CYTHON_UNUSED_VAR({Naming.modulestatevalue_cname});")

        code = self.parts['utility_code_proto']
        code.putln("")
        code.putln("/* --- Runtime support code (head) --- */")

        code = self.parts['utility_code_def']
        if self.code_config.emit_linenums:
            code.write('\n#line 1 "cython_utility"\n')
        code.putln("")
        code.putln("/* --- Runtime support code --- */")

    def initialize_main_h_code(self):
        rootwriter = self.rootwriter
        for part in self.h_code_layout:
            self.parts[part] = rootwriter.insertion_point()

    def finalize_main_c_code(self):
        self.close_global_decls()

        #
        # utility_code_def
        #
        code = self.parts['utility_code_def']
        util = TempitaUtilityCode.load_cached("TypeConversions", "TypeConversion.c")
        code.put(util.format_code(util.impl))
        code.putln("")

        #
        # utility code pragmas
        #
        code = self.parts['utility_code_pragmas']
        util = UtilityCode.load_cached("UtilityCodePragmas", "ModuleSetupCode.c")
        code.putln(util.format_code(util.impl))
        code.putln("")
        code = self.parts['utility_code_pragmas_end']
        util = UtilityCode.load_cached("UtilityCodePragmasEnd", "ModuleSetupCode.c")
        code.putln(util.format_code(util.impl))
        code.putln("")

    def __getitem__(self, key):
        return self.parts[key]

    #
    # Global constants, interned objects, etc.
    #
    def close_global_decls(self):
        # This is called when it is known that no more global declarations will
        # declared.
        self.generate_const_declarations()

        w = self.parts['cached_builtins']
        w.putln("return 0;")
        if w.label_used(w.error_label):
            w.put_label(w.error_label)
            w.putln("return -1;")
        w.putln("}")
        w.exit_cfunc_scope()

        w = self.parts['cached_constants']
        for const_type in ["tuple", "slice"]:
            if const_type in self.const_array_counters:
                self.immortalize_constants(
                    w.name_in_module_state(Naming.pyrex_prefix + const_type),
                    self.const_array_counters[const_type],
                    w)
        w.put_finish_refcount_context()
        w.putln("return 0;")
        if w.label_used(w.error_label):
            w.put_label(w.error_label)
            w.put_finish_refcount_context()
            w.putln("return -1;")
        w.putln("}")
        w.exit_cfunc_scope()

        for part in ['init_globals', 'init_constants']:
            w = self.parts[part]
            w.putln("return 0;")
            if w.label_used(w.error_label):
                w.put_label(w.error_label)
                w.putln("return -1;")
            w.putln("}")
            w.exit_cfunc_scope()

        if Options.generate_cleanup_code:
            w = self.parts['cleanup_globals']
            w.putln("}")
            w.exit_cfunc_scope()

        if Options.generate_cleanup_code:
            w = self.parts['cleanup_module']
            w.putln("}")
            w.exit_cfunc_scope()

    def put_pyobject_decl(self, entry):
        self['global_var'].putln("static PyObject *%s;" % entry.cname)

    # constant handling at code generation time

    def get_cached_constants_writer(self, target=None):
        if target is not None:
            if target in self.initialised_constants:
                # Return None on second/later calls to prevent duplicate creation code.
                return None
            self.initialised_constants.add(target)
        return self.parts['cached_constants']

    def get_int_const(self, str_value, longness=False):
        py_type = longness and 'long' or 'int'
        try:
            c = self.num_const_index[(str_value, py_type)]
        except KeyError:
            c = self.new_num_const(str_value, py_type)
        return c

    def get_float_const(self, str_value, value_code):
        try:
            c = self.num_const_index[(str_value, 'float')]
        except KeyError:
            c = self.new_num_const(str_value, 'float', value_code)
        return c

    def get_py_const(self, prefix, dedup_key=None):
        if dedup_key is not None:
            const = self.dedup_const_index.get(dedup_key)
            if const is not None:
                return const
        const = self.new_array_const_cname(prefix)
        if dedup_key is not None:
            self.dedup_const_index[dedup_key] = const
        return const

    def get_argument_default_const(self, type):
        cname = self.new_const_cname('')
        c = PyObjectConst(cname, type)
        self.arg_default_constants.append(c)
        # Argument default constants aren't currently cleaned up.
        # If that changes, it needs to account for the fact that they
        # aren't just Python objects
        return c

    def get_string_const(self, text, c_used=True):
        # return a C string constant, creating a new one if necessary
        if text.is_unicode:
            byte_string = text.utf8encode()
        else:
            byte_string = text.byteencode()
        try:
            c = self.string_const_index[byte_string]
        except KeyError:
            c = self.new_string_const(text, byte_string)
        if c_used:
            c.c_used = True
        return c

    def get_pyunicode_ptr_const(self, text):
        # return a Py_UNICODE[] constant, creating a new one if necessary
        assert text.is_unicode
        try:
            c = self.pyunicode_ptr_const_index[text]
        except KeyError:
            c = self.pyunicode_ptr_const_index[text] = self.new_const_cname()
        return c

    def get_py_string_const(self, text, identifier=None):
        # return a Python string constant, creating a new one if necessary
        c_string: StringConst = self.get_string_const(text, c_used=False)
        py_string = c_string.get_py_string_const(text.encoding, identifier)
        return py_string

    def get_py_codeobj_const(self, node):
        idx = len(self.codeobject_constants)
        name = f"{Naming.codeobjtab_cname}[{idx}]"
        self.codeobject_constants.append(node)
        return name

    def get_interned_identifier(self, text):
        return self.get_py_string_const(text, identifier=True)

    def new_string_const(self, text, byte_string):
        cname = self.new_string_const_cname(byte_string)
        c = StringConst(cname, text, byte_string)
        self.string_const_index[byte_string] = c
        return c

    def new_num_const(self, value, py_type, value_code=None):
        cname = self.new_num_const_cname(value, py_type)
        c = NumConst(cname, value, py_type, value_code)
        self.num_const_index[(value, py_type)] = c
        return c

    def new_string_const_cname(self, bytes_value):
        # Create a new globally-unique nice name for a C string constant.
        value = bytes_value.decode('ASCII', 'ignore')
        return self.new_const_cname(value=value)

    def unique_const_cname(self, format_str):  # type: (str) -> str
        used = self.const_cnames_used
        cname = value = format_str.format(sep='', counter='')
        while cname in used:
            counter = used[value] = used[value] + 1
            cname = format_str.format(sep='_', counter=counter)
        used[cname] = 1
        return cname

    def new_num_const_cname(self, value, py_type):  # type: (str, str) -> str
        if py_type == 'long':
            value += 'L'
            py_type = 'int'
        prefix = Naming.interned_prefixes[py_type]

        value = value.replace('.', '_').replace('+', '_').replace('-', 'neg_')
        if len(value) > 42:
            # update tests/run/large_integer_T5290.py in case the amount is changed
            cname = self.unique_const_cname(
                prefix + "large{counter}_" + value[:18] + "_xxx_" + value[-18:])
        else:
            cname = "%s%s" % (prefix, value)
        return cname

    def new_const_cname(self, prefix='', value=''):
        value = replace_identifier('_', value)[:32].strip('_')
        name_suffix = self.unique_const_cname(value + "{sep}{counter}")
        if prefix:
            prefix = Naming.interned_prefixes[prefix]
        else:
            prefix = Naming.const_prefix
        return "%s%s" % (prefix, name_suffix)

    def new_array_const_cname(self, prefix: str):
        count = self.const_array_counters.get(prefix, 0)
        self.const_array_counters[prefix] = count+1
        return f"{Naming.pyrex_prefix}{prefix}[{count}]"

    def get_cached_unbound_method(self, type_cname, method_name):
        key = (type_cname, method_name)
        try:
            cname = self.cached_cmethods[key]
        except KeyError:
            cname = self.cached_cmethods[key] = self.new_const_cname(
                'umethod', '%s_%s' % (type_cname, method_name))
        return cname

    def cached_unbound_method_call_code(self, modulestate_cname, obj_cname, type_cname, method_name, arg_cnames):
        # admittedly, not the best place to put this method, but it is reused by UtilityCode and ExprNodes ...
        utility_code_name = "CallUnboundCMethod%d" % len(arg_cnames)
        self.use_utility_code(UtilityCode.load_cached(utility_code_name, "ObjectHandling.c"))
        cache_cname = self.get_cached_unbound_method(type_cname, method_name)
        args = [obj_cname] + arg_cnames
        return "__Pyx_%s(&%s%s, %s)" % (
            utility_code_name,
            modulestate_cname,
            cache_cname,
            ', '.join(args),
        )

    def add_cached_builtin_decl(self, entry):
        if entry.is_builtin and entry.is_const:
            if self.should_declare(entry.cname, entry):
                self.put_pyobject_decl(entry)
                name = entry.name
                if name in renamed_py2_builtins_map:
                    name = renamed_py2_builtins_map[name]
                self.put_cached_builtin_init(
                    entry.pos, StringEncoding.EncodedString(name),
                    entry.cname)

    def put_cached_builtin_init(self, pos, name, cname):
        w = self.parts['cached_builtins']
        cname_in_modulestate = w.name_in_main_c_code_module_state(
            self.get_interned_identifier(name).cname)
        self.use_utility_code(
            UtilityCode.load_cached("GetBuiltinName", "ObjectHandling.c"))
        w.putln('%s = __Pyx_GetBuiltinName(%s); if (!%s) %s' % (
            cname,
            cname_in_modulestate,
            cname,
            w.error_goto(pos)))

    def generate_const_declarations(self):
        self.generate_cached_methods_decls()
        self.generate_object_constant_decls()
        self.generate_codeobject_constants()
        # generate code for string and numeric constants as late as possible
        # to allow new constants be to created by the earlier stages.
        # (although the constants themselves are written early)
        self.generate_string_constants()
        self.generate_num_constants()

    def _generate_module_array_traverse_and_clear(self, struct_attr_cname, count, may_have_refcycles=True):
        counter_type = 'int' if count < 2**15 else 'Py_ssize_t'
        visit_call = "Py_VISIT" if may_have_refcycles else "__Pyx_VISIT_CONST"

        writer = self.parts['module_state_traverse']
        writer.putln(f"for ({counter_type} i=0; i<{count}; ++i) {{ {visit_call}(traverse_module_state->{struct_attr_cname}[i]); }}")

        writer = self.parts['module_state_clear']
        writer.putln(f"for ({counter_type} i=0; i<{count}; ++i) {{ Py_CLEAR(clear_module_state->{struct_attr_cname}[i]); }}")

    def generate_object_constant_decls(self):
        consts = [(len(c.cname), c.cname, c)
                  for c in self.arg_default_constants]
        consts.sort()
        for _, cname, c in consts:
            self.parts['module_state'].putln("%s;" % c.type.declaration_code(cname))
            if not c.type.needs_refcounting:
                # Note that py_constants is used for all argument defaults
                # which aren't necessarily PyObjects, so aren't appropriate
                # to clear.
                continue

            self.parts['module_state_clear'].put_xdecref_clear(
                f"clear_module_state->{cname}",
                c.type,
                clear_before_decref=True,
                nanny=False,
            )

            if c.type.is_memoryviewslice:
                # TODO: Implement specific to type like CodeWriter.put_xdecref_clear()
                cname += "->memview"

            self.parts['module_state_traverse'].putln(
                f"Py_VISIT(traverse_module_state->{cname});")

        for prefix, count in sorted(self.const_array_counters.items()):
            struct_attr_cname = f"{Naming.pyrex_prefix}{prefix}"
            self.parts['module_state'].putln(f"PyObject *{struct_attr_cname}[{count}];")

            # The constant tuples/slices that we create can never participate in reference cycles.
            self._generate_module_array_traverse_and_clear(struct_attr_cname, count, may_have_refcycles=False)

            cleanup_level = cleanup_level_for_type_prefix(prefix)
            if cleanup_level is not None and cleanup_level <= Options.generate_cleanup_code:
                part_writer = self.parts['cleanup_globals']
                part_writer.put(f"for (size_t i=0; i<{count}; ++i) ")
                part_writer.putln(
                    "{ Py_CLEAR(%s); }" %
                        part_writer.name_in_main_c_code_module_state(f"{struct_attr_cname}[i]")
                )

    def generate_cached_methods_decls(self):
        if not self.cached_cmethods:
            return

        decl = self.parts['module_state']
        init = self.parts['cached_builtins']

        init.putln("")
        init.putln("/* Cached unbound methods */")

        cnames = []
        for (type_cname, method_name), cname in sorted(self.cached_cmethods.items()):
            cnames.append(cname)
            method_name_cname = self.get_interned_identifier(StringEncoding.EncodedString(method_name)).cname
            decl.putln('__Pyx_CachedCFunction %s;' % (
                cname))
            # split type reference storage as it might not be static
            init.putln('%s.type = (PyObject*)%s;' % (
                init.name_in_main_c_code_module_state(cname), type_cname))
            # method name string isn't static in limited api
            init.putln(
                f'{init.name_in_main_c_code_module_state(cname)}.method_name = '
                f'&{init.name_in_main_c_code_module_state(method_name_cname)};')

        if Options.generate_cleanup_code:
            cleanup = self.parts['cleanup_globals']
            for cname in cnames:
                cleanup.putln(f"Py_CLEAR({init.name_in_main_c_code_module_state(cname)}.method);")

    def generate_string_constants(self):
        c_consts = []
        py_bytes_consts = []
        py_unicode_consts = []

        # Split into buckets.
        for _, _, c in sorted([(len(c.cname), c.cname, c) for c in self.string_const_index.values()]):
            if c.c_used:
                c_consts.append((len(c.cname), c.cname, c.escaped_value))
            if c.py_strings:
                for py_string in c.py_strings.values():
                    text = c.text
                    if py_string.is_unicode and not isinstance(text, str):
                        text = StringEncoding.EncodedString(text.decode(py_string.encoding or 'UTF-8'))

                    (py_unicode_consts if py_string.is_unicode else py_bytes_consts).append((
                        py_string.intern and py_string.is_unicode,
                        py_string.cname,
                        text,
                    ))

        c_consts.sort()
        py_bytes_consts.sort()
        py_unicode_consts.sort()

        # Generate C string constants.
        decls_writer = self.parts['string_decls']
        for _, cname, escaped_value in c_consts:
            cliteral = StringEncoding.split_string_literal(escaped_value)
            decls_writer.putln(
                f'static const char {cname}[] = "{cliteral}";',
                safe=True,  # Braces in user strings are not for indentation.
            )

        # Generate legacy Py_UNICODE[] constants.
        for c, cname in sorted(self.pyunicode_ptr_const_index.items()):
            utf16_array, utf32_array = StringEncoding.encode_pyunicode_string(c)
            if utf16_array:
                # Narrow and wide representations differ
                decls_writer.putln("#ifdef Py_UNICODE_WIDE")
            decls_writer.putln("static Py_UNICODE %s[] = { %s };" % (cname, utf32_array))
            if utf16_array:
                decls_writer.putln("#else")
                decls_writer.putln("static Py_UNICODE %s[] = { %s };" % (cname, utf16_array))
                decls_writer.putln("#endif")

        # Generate stringtab and Python string constants.
        py_string_count = len(py_bytes_consts) + len(py_unicode_consts)
        self.parts['module_state'].putln(f"PyObject *{Naming.stringtab_cname}[{py_string_count}];")
        self._generate_module_array_traverse_and_clear(Naming.stringtab_cname, py_string_count, may_have_refcycles=False)

        self.generate_pystring_constants(py_unicode_consts, py_bytes_consts)

    def generate_pystring_constants(self, text_strings: list, byte_strings: list):
        # Concatenate all strings into one byte sequence and build a length index array.
        defines = self.parts['constant_name_defines']

        bytes_values = []
        first_interned: cython.Py_ssize_t = -1
        stringtab_pos: cython.Py_ssize_t = 0

        # For (Unicode) text strings, the index stores the character lengths after UTF8 decoding.
        for i, (is_interned, cname, text) in enumerate(text_strings):
            bytes_values.append(text.encode('utf-8'))
            if first_interned == -1 and is_interned:
                first_interned = i
            defines.putln(f"#define {cname} {Naming.stringtab_cname}[{stringtab_pos}]")
            stringtab_pos += 1

        stringtab_bytes_start: cython.Py_ssize_t = len(text_strings)

        # For bytes objects, the index stores the byte lengths, ignoring the initial Unicode string.
        for _, cname, text in byte_strings:
            bytes_values.append(text.byteencode() if text.encoding else text.utf8encode())
            defines.putln(f"#define {cname} {Naming.stringtab_cname}[{stringtab_pos}]")
            stringtab_pos += 1

        index = list(map(len, bytes_values))
        concat_bytes = b''.join(bytes_values)

        w = self.parts['init_constants']
        w.putln("{")  # Start code block.

        # Store the index of string lengths.
        w.putln(
            "const struct { "
            f"const unsigned int length: {max(index).bit_length()}; "
            "} "
            f"index[] = {{{','.join(['{%d}' % length for length in index])}}};",
        )

        # Store and decompress the string data.
        self.use_utility_code(UtilityCode.load_cached("DecompressString", "StringTools.c"))

        has_if = False
        for algo_number, algo_name, compress in reversed(compression_algorithms):
            if compress is None:
                continue
            compressed_bytes = compress(concat_bytes)
            if len(compressed_bytes) >= len(concat_bytes) - 10:
                continue

            if algo_name == 'zlib':
                # Use zlib as fallback if the selected compression module is not available.
                assert algo_number == 1, f"Compression algorithm no. 1 must be 'zlib' to be used as fallback."
                guard = "(CYTHON_COMPRESS_STRINGS) != 0"
            elif algo_name == 'zstd':
                # 'compression.zstd' was added in Python 3.14.
                guard = f"(CYTHON_COMPRESS_STRINGS) == {algo_number} && __PYX_LIMITED_VERSION_HEX >= 0x030e0000"
            else:
                guard = f"(CYTHON_COMPRESS_STRINGS) == {algo_number}"

            w.putln(f"#{'if' if not has_if else 'elif'} {guard} /* compression: {algo_name} ({len(compressed_bytes)} bytes) */")
            has_if = True
            escaped_bytes = StringEncoding.split_string_literal(
                StringEncoding.escape_byte_string(compressed_bytes))
            w.putln(f'const char* const cstring = "{escaped_bytes}";', safe=True)
            w.putln(f'PyObject *data = __Pyx_DecompressString(cstring, {len(compressed_bytes)}, {algo_number});')
            w.putln(w.error_goto_if_null('data', self.module_pos))

            w.putln('const char* const bytes = __Pyx_PyBytes_AsString(data);')
            w.putln("#if !CYTHON_ASSUME_SAFE_MACROS")
            w.putln(f'if (likely(bytes)); else {{ Py_DECREF(data); {w.error_goto(self.module_pos)} }}')
            w.putln('#endif')

        if has_if:
            w.putln(f"#else /* compression: none ({len(concat_bytes)} bytes) */")
        escaped_bytes = StringEncoding.split_string_literal(
            StringEncoding.escape_byte_string(concat_bytes))
        w.putln(f'const char* const bytes = "{escaped_bytes}";', safe=True)
        w.putln('PyObject *data = NULL;')  # Always allow xdecref below.
        w.putln("CYTHON_UNUSED_VAR(__Pyx_DecompressString);")
        if has_if:
            w.putln("#endif")

        # Populate stringtab.
        w.putln(f"PyObject **stringtab = {w.name_in_main_c_code_module_state(Naming.stringtab_cname)};")
        w.putln("Py_ssize_t pos = 0;")

        # Unpack Unicode strings.
        if stringtab_bytes_start > 0:
            # Note: We could decode the concatenated Unicode string in one go, but this has a drawback:
            # If most strings are ASCII/Latin-1 or at most BMP, then a single non-BMP string in the mix
            # will make all strings use 4 bytes of RAM per character during initialisation, until we finish
            # splitting the user substrings. In addition to using more memory, this might not even be faster
            # because it must copy Unicode slices between different character sizes.
            # We avoid this by repeatedly calling PyUnicode_DecodeUTF8() for each substring.
            w.putln(f"for ({'int' if stringtab_bytes_start < 2**15 else 'Py_ssize_t'} i = 0; i < {stringtab_bytes_start}; i++) {{")
            w.putln("Py_ssize_t bytes_length = index[i].length;")

            w.putln("PyObject *string = PyUnicode_DecodeUTF8(bytes + pos, bytes_length, NULL);")
            if first_interned >= 0:
                w.putln(f"if (likely(string) && i >= {first_interned}) PyUnicode_InternInPlace(&string);")
            w.putln("if (unlikely(!string)) {")
            w.putln("Py_XDECREF(data);")
            w.putln(w.error_goto(self.module_pos))
            w.putln('}')

            w.putln("stringtab[i] = string;")
            w.putln("pos += bytes_length;")
            w.putln("}")  # for()

        # Unpack byte strings.
        if stringtab_bytes_start < len(index):
            w.putln(f"for ({'int' if len(index) < 2**15 else 'Py_ssize_t'} i = {stringtab_bytes_start}; i < {len(index)}; i++) {{")
            w.putln("Py_ssize_t bytes_length = index[i].length;")

            w.putln("PyObject *string = PyBytes_FromStringAndSize(bytes + pos, bytes_length);")
            w.putln("stringtab[i] = string;")
            w.putln("pos += bytes_length;")

            w.putln("if (unlikely(!string)) {")
            w.putln("Py_XDECREF(data);")
            w.putln(w.error_goto(self.module_pos))
            w.putln('}')

            w.putln("}")  # for()

        w.putln("Py_XDECREF(data);")

        # Set up hash values.
        w.putln(f"for (Py_ssize_t i = 0; i < {len(index)}; i++) {{")
        w.putln("if (unlikely(PyObject_Hash(stringtab[i]) == -1)) {")
        w.putln(w.error_goto(self.module_pos))
        w.putln('}')
        w.putln('}')

        # Unicode strings are not trivially immortal but require certain rules.
        # See https://github.com/python/cpython/blob/920de7ccdcfa7284b6d23a124771b17c66dd3e4f/Objects/unicodeobject.c#L713-L739
        # But we can make bytes strings immortal.
        if stringtab_bytes_start < len(index):
            self.immortalize_constants(f"stringtab + {stringtab_bytes_start}", len(index) - stringtab_bytes_start, w)

        w.putln("}")  # close block

    def generate_codeobject_constants(self):
        w = self.parts['init_codeobjects']
        init_function = (
            f"int __Pyx_CreateCodeObjects({Naming.modulestatetype_cname} *{Naming.modulestatevalue_cname})"
        )

        if not self.codeobject_constants:
            w.start_initcfunc(init_function)
            w.putln(f"CYTHON_UNUSED_VAR({Naming.modulestatevalue_cname});")
            w.putln("return 0;")
            w.exit_cfunc_scope()
            w.putln("}")
            return

        # Create a downsized config struct and build code objects from it.
        max_flags = 0x3ff  # to be adapted when we start using new flags
        max_func_args = 1
        max_kwonly_args = 1
        max_posonly_args = 1
        max_vars = 1
        max_line = 1
        for node in self.codeobject_constants:
            def_node = node.def_node
            if not def_node.is_generator_expression:
                max_func_args = max(max_func_args, len(def_node.args) - def_node.num_kwonly_args)
                max_kwonly_args = max(max_kwonly_args, def_node.num_kwonly_args)
                max_posonly_args = max(max_posonly_args, def_node.num_posonly_args)
            max_vars = max(max_vars, len(node.varnames))
            max_line = max(max_line, def_node.pos[1])

        w.put(textwrap.dedent(f"""\
        typedef struct {{
            unsigned int argcount : {max_func_args.bit_length()};
            unsigned int num_posonly_args : {max_posonly_args.bit_length()};
            unsigned int num_kwonly_args : {max_kwonly_args.bit_length()};
            unsigned int nlocals : {max_vars.bit_length()};
            unsigned int flags : {max_flags.bit_length()};
            unsigned int first_line : {max_line.bit_length()};
        }} __Pyx_PyCode_New_function_description;
        """))

        self.use_utility_code(UtilityCode.load_cached("NewCodeObj", "ModuleSetupCode.c"))

        w.start_initcfunc(init_function)

        w.putln("PyObject* tuple_dedup_map = PyDict_New();")
        w.putln("if (unlikely(!tuple_dedup_map)) return -1;")

        for node in self.codeobject_constants:
            node.generate_codeobj(w, "bad")

        w.putln("Py_DECREF(tuple_dedup_map);")
        w.putln("return 0;")

        w.putln("bad:")
        w.putln("Py_DECREF(tuple_dedup_map);")
        w.putln("return -1;")
        w.exit_cfunc_scope()
        w.putln("}")

        code_object_count = len(self.codeobject_constants)
        self.parts['module_state'].putln(f"PyObject *{Naming.codeobjtab_cname}[{code_object_count}];")
        # The code objects that we generate only contain plain constants and can never participate in reference cycles.
        self._generate_module_array_traverse_and_clear(Naming.codeobjtab_cname, code_object_count, may_have_refcycles=False)

    def generate_num_constants(self):
        consts = [(c.py_type, len(c.value.lstrip('-')), c.value.lstrip('-'), c.value, c.value_code, c)
                  for c in self.num_const_index.values()]
        consts.sort()
        if not consts:
            return

        constant_count = len(consts)
        self.parts['module_state'].putln(f"PyObject *{Naming.numbertab_cname}[{constant_count}];")
        # Numeric constants can never participate in reference cycles.
        self._generate_module_array_traverse_and_clear(Naming.numbertab_cname, constant_count, may_have_refcycles=False)

        float_constants = []
        int_constants_by_bytesize = [[]]  # [[1 byte], [2 bytes], [4 bytes], [8 bytes]]
        large_constants = []
        int_constant_count = 0
        int_suffix = ''

        for py_type, _, _, value, value_code, c in consts:
            cname = c.cname
            if py_type == 'float':
                float_constants.append((cname, value_code))
            else:
                number_value = Utils.str_to_number(value)
                bit_length = number_value.bit_length()
                if bit_length <= 63:
                    while (bit_length + 8) // 8 > 1 << (len(int_constants_by_bytesize) - 1):
                        int_constants_by_bytesize.append([])
                        # Our <= 31-bit integer values pass happily as 'int32' without further modifiers,
                        # but MSVC misinterprets a negative '-2147483648' (== INT_MIN) and similar values as
                        # "that's 'uint32' just with a minus sign", where '-(2147483648)' == '2147483648'.
                        # See https://learn.microsoft.com/en-us/cpp/error-messages/compiler-warnings/compiler-warning-level-2-c4146?view=msvc-170
                        int_suffix = 'LL'[:len(int_constants_by_bytesize) - 2]
                    int_constant_count += 1
                    int_constants_by_bytesize[-1].append((cname, f"{number_value}{int_suffix}"))
                else:
                    large_constants.append((cname, number_value))

        w = self.parts['init_constants']
        defines = self.parts['constant_name_defines']

        def store_array(w, name: str, ctype: str, constants: list):
            c: tuple
            values = ','.join([c[1] for c in constants])
            w.putln(f"{ctype} const {name}[] = {{{values}}};")

        def generate_forloop_start(w, end: cython.Py_ssize_t):
            counter_type = 'int' if end < 2**15 else 'Py_ssize_t'
            w.putln(f"for ({counter_type} i = 0; i < {end}; i++) {{")

        def assign_constant(w, error_pos, rhs_code: str):
            w.putln(f"numbertab[i] = {rhs_code};")
            w.putln(w.error_goto_if_null("numbertab[i]", error_pos))

        def define_constants(defines, constants: list, start_offset: cython.Py_ssize_t = 0):
            i: cython.Py_ssize_t
            c: tuple
            numbertab_cname: str = Naming.numbertab_cname
            for i, c in enumerate(constants):
                cname: str = c[0]
                defines.putln(f"#define {cname} {numbertab_cname}[{start_offset + i}]")

        constant_offset: cython.Py_ssize_t = 0

        if float_constants:
            w.putln("{")
            w.putln(f"PyObject **numbertab = {w.name_in_main_c_code_module_state(Naming.numbertab_cname)};")

            store_array(w, "c_constants", 'double', float_constants)
            define_constants(defines, float_constants, constant_offset)

            generate_forloop_start(w, len(float_constants))
            assign_constant(w, self.module_pos, "PyFloat_FromDouble(c_constants[i])")
            w.putln("}")  # for()

            w.putln("}")
            constant_offset += len(float_constants)

        if int_constant_count > 0:
            w.putln("{")
            w.putln(f"PyObject **numbertab = {w.name_in_main_c_code_module_state(Naming.numbertab_cname)} + {constant_offset};")

            int_types = ['', 'int8_t', 'int16_t', 'int32_t', 'int64_t']
            array_access = "%s"
            int_constants_seen: cython.Py_ssize_t = 0
            byte_size: cython.int
            for byte_size, constants in enumerate(int_constants_by_bytesize, 1):
                if not constants:
                    continue

                array_name = f"cint_constants_{1 << (byte_size - 1)}"
                store_array(w, array_name, int_types[byte_size], constants)
                define_constants(defines, constants, constant_offset + int_constants_seen)

                read_item = f"{array_name}[i - {int_constants_seen}]"
                int_constants_seen += len(constants)
                array_access %= (
                    read_item if byte_size == len(int_constants_by_bytesize)  # last is simple access
                    else f"(i < {int_constants_seen} ? {read_item} : %s)"  # otherwise, access arrays step by step
                )

            generate_forloop_start(w, int_constant_count)
            capi_func = "PyLong_FromLong" if len(int_constants_by_bytesize) <= 3 else "PyLong_FromLongLong"
            assign_constant(w, self.module_pos, f"{capi_func}({array_access})")
            w.putln("}")  # for()

            w.putln("}")
            constant_offset += int_constant_count

        if large_constants:
            # We store large integer constants in a single '\0'-separated C string of base32 digits.
            def to_base32(number):
                is_neg: bool = number < 0
                if is_neg:
                    number = -number

                digits = bytearray()
                while number:
                    digit: cython.uint = number & 31
                    digit_char: cython.char = b'0123456789abcdefghijklmnopqrstuv'[digit]
                    digits.append(digit_char)
                    number >>= 5
                if not digits:
                    return b'0'

                if is_neg:
                    digits.append(ord(b'-'))
                digits.reverse()
                return digits

            w.putln("{")
            w.putln(f"PyObject **numbertab = {w.name_in_main_c_code_module_state(Naming.numbertab_cname)} + {constant_offset};")
            c_string = b'\\000'.join([to_base32(c[1]) for c in large_constants]).decode('ascii')
            w.putln(f'const char* c_constant = "{StringEncoding.split_string_literal(c_string)}";')
            define_constants(defines, large_constants, constant_offset)

            generate_forloop_start(w, len(large_constants))
            w.putln("char *end_pos;")
            assign_constant(w, self.module_pos, "PyLong_FromString(c_constant, &end_pos, 32)")
            w.putln("c_constant = end_pos + 1;")
            w.putln("}")  # for()

            w.putln("}")

        self.immortalize_constants(
            w.name_in_main_c_code_module_state(Naming.numbertab_cname),
            constant_count,
            w)

    @staticmethod
    def immortalize_constants(array_cname, constant_count, writer):
        writer.putln("#if CYTHON_IMMORTAL_CONSTANTS")
        writer.putln("{")
        writer.putln(f"PyObject **table = {array_cname};")
        writer.putln(f"for (Py_ssize_t i=0; i<{constant_count}; ++i) {{")
        writer.putln("#if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING")
        # We don't want to set the refcount on shared constants (e.g. cached integers)
        # because setting the refcount isn't thread-safe. The chances are that most of the constants
        # that this applies to are already immortal though so that isn't a great loss.
        writer.putln("#if PY_VERSION_HEX < 0x030E0000")
        writer.putln("if (_Py_IsOwnedByCurrentThread(table[i]) && Py_REFCNT(table[i]) == 1)")
        writer.putln("#else")
        writer.putln("if (PyUnstable_Object_IsUniquelyReferenced(table[i]))")
        writer.putln("#endif")
        writer.putln("{")
        writer.putln("Py_SET_REFCNT(table[i], _Py_IMMORTAL_REFCNT_LOCAL);")
        writer.putln("}")
        writer.putln("#else")
        writer.putln("Py_SET_REFCNT(table[i], _Py_IMMORTAL_INITIAL_REFCNT);")
        writer.putln("#endif")
        writer.putln("}")  # for()
        writer.putln("}")
        writer.putln("#endif")

    # The functions below are there in a transition phase only
    # and will be deprecated. They are called from Nodes.BlockNode.
    # The copy&paste duplication is intentional in order to be able
    # to see quickly how BlockNode worked, until this is replaced.

    def should_declare(self, cname, entry):
        if cname in self.declared_cnames:
            other = self.declared_cnames[cname]
            assert str(entry.type) == str(other.type)
            assert entry.init == other.init
            return False
        else:
            self.declared_cnames[cname] = entry
            return True

    #
    # File name state
    #

    def lookup_filename(self, source_desc):
        entry = source_desc.get_filenametable_entry()
        try:
            index = self.filename_table[entry]
        except KeyError:
            index = len(self.filename_list)
            self.filename_list.append(source_desc)
            self.filename_table[entry] = index
        return index

    def commented_file_contents(self, source_desc):
        try:
            return self.input_file_contents[source_desc]
        except KeyError:
            pass
        source_file = source_desc.get_lines(encoding='ASCII', error_handling='ignore')
        F = [' * ' + (
                line.replace(
                    '*/', '*[inserted by cython to avoid comment closer]/'
                ).replace(
                    '/*', '/[inserted by cython to avoid comment start]*'
                ) if '/' in line else line)
            for line in source_file
        ]
        if not F: F.append('')
        self.input_file_contents[source_desc] = F
        return F

    #
    # Utility code state
    #

    def use_utility_code(self, utility_code, used_by=None):
        """
        Adds code to the C file. utility_code should
        a) implement __eq__/__hash__ for the purpose of knowing whether the same
           code has already been included
        b) implement put_code, which takes a globalstate instance

        See UtilityCode.
        """
        if utility_code and utility_code not in self.utility_codes:
            self.utility_codes.add(utility_code)
            utility_code.put_code(self, used_by=used_by)

    def use_entry_utility_code(self, entry):
        if entry is None:
            return
        if entry.utility_code:
            self.use_utility_code(entry.utility_code)
        if entry.utility_code_definition:
            self.use_utility_code(entry.utility_code_definition)
        from . import PyrexTypes
        for tp in PyrexTypes.get_all_subtypes(entry.type):
            if hasattr(tp, "entry") and tp.entry is not entry:
                self.use_entry_utility_code(tp.entry)


def funccontext_property(func):
    name = func.__name__
    attribute_of = operator.attrgetter(name)
    def get(self):
        return attribute_of(self.funcstate)
    def set(self, value):
        setattr(self.funcstate, name, value)
    return property(get, set)


class CCodeConfig:
    # emit_linenums       boolean         write #line pragmas?
    # emit_code_comments  boolean         copy the original code into C comments?
    # c_line_in_traceback boolean         append the c file and line number to the traceback for exceptions?

    def __init__(self, emit_linenums=True, emit_code_comments=True, c_line_in_traceback=True):
        self.emit_code_comments = emit_code_comments
        self.emit_linenums = emit_linenums
        self.c_line_in_traceback = c_line_in_traceback


class CCodeWriter:
    """
    Utility class to output C code.

    When creating an insertion point one must care about the state that is
    kept:
    - formatting state (level, bol) is cloned and used in insertion points
      as well
    - labels, temps, exc_vars: One must construct a scope in which these can
      exist by calling enter_cfunc_scope/exit_cfunc_scope (these are for
      sanity checking and forward compatibility). Created insertion points
      looses this scope and cannot access it.
    - marker: Not copied to insertion point
    - filename_table, filename_list, input_file_contents: All codewriters
      coming from the same root share the same instances simultaneously.
    """

    # f                   file            output file
    # buffer              StringIOTree

    # level               int             indentation level
    # bol                 bool            beginning of line?
    # marker              string          comment to emit before next line
    # funcstate           FunctionState   contains state local to a C function used for code
    #                                     generation (labels and temps state etc.)
    # globalstate         GlobalState     contains state global for a C file (input file info,
    #                                     utility code, declared constants etc.)
    # pyclass_stack       list            used during recursive code generation to pass information
    #                                     about the current class one is in
    # code_config         CCodeConfig     configuration options for the C code writer

    @cython.locals(create_from='CCodeWriter')
    def __init__(self, create_from=None, buffer=None, copy_formatting=False):
        if buffer is None: buffer = StringIOTree()
        self.buffer = buffer
        self.last_pos = None
        self.last_marked_pos = None
        self.pyclass_stack = []

        self.funcstate = None
        self.globalstate = None
        self.code_config = None
        self.level = 0
        self.call_level = 0
        self.bol = 1

        if create_from is not None:
            # Use same global state
            self.set_global_state(create_from.globalstate)
            self.funcstate = create_from.funcstate
            # Clone formatting state
            if copy_formatting:
                self.level = create_from.level
                self.bol = create_from.bol
                self.call_level = create_from.call_level
            self.last_pos = create_from.last_pos
            self.last_marked_pos = create_from.last_marked_pos

    def create_new(self, create_from, buffer, copy_formatting):
        # polymorphic constructor -- very slightly more versatile
        # than using __class__
        result = CCodeWriter(create_from, buffer, copy_formatting)
        return result

    def set_global_state(self, global_state):
        assert self.globalstate is None  # prevent overwriting once it's set
        self.globalstate = global_state
        self.code_config = global_state.code_config

    def copyto(self, f):
        self.buffer.copyto(f)

    def getvalue(self):
        return self.buffer.getvalue()

    def write(self, s):
        if '\n' in s:
            self._write_lines(s)
        else:
            self._write_to_buffer(s)

    @cython.final
    def _write_lines(self, s):
        # Cygdb needs to know which Cython source line corresponds to which C line.
        # Therefore, we write this information into "self.buffer.markers" and then write it from there
        # into cython_debug/cython_debug_info_* (see ModuleNode._serialize_lineno_map).
        filename_line = self.last_marked_pos[:2] if self.last_marked_pos else (None, 0)
        self.buffer.markers.extend([filename_line] * s.count('\n'))

        self._write_to_buffer(s)

    def _write_to_buffer(self, s):
        self.buffer.write(s)

    def insertion_point(self):
        other = self.create_new(create_from=self, buffer=self.buffer.insertion_point(), copy_formatting=True)
        return other

    def new_writer(self):
        """
        Creates a new CCodeWriter connected to the same global state, which
        can later be inserted using insert.
        """
        return CCodeWriter(create_from=self)

    def insert(self, writer):
        """
        Inserts the contents of another code writer (created with
        the same global state) in the current location.

        It is ok to write to the inserted writer also after insertion.
        """
        assert writer.globalstate is self.globalstate
        self.buffer.insert(writer.buffer)

    # Properties delegated to function scope
    @funccontext_property
    def label_counter(self): pass
    @funccontext_property
    def return_label(self): pass
    @funccontext_property
    def error_label(self): pass
    @funccontext_property
    def labels_used(self): pass
    @funccontext_property
    def continue_label(self): pass
    @funccontext_property
    def break_label(self): pass
    @funccontext_property
    def return_from_error_cleanup_label(self): pass
    @funccontext_property
    def yield_labels(self): pass

    def label_interceptor(self, new_labels, orig_labels, skip_to_label=None, pos=None, trace=True):
        """
        Helper for generating multiple label interceptor code blocks.

        @param new_labels: the new labels that should be intercepted
        @param orig_labels: the original labels that we should dispatch to after the interception
        @param skip_to_label: a label to skip to before starting the code blocks
        @param pos: the node position to mark for each interceptor block
        @param trace: add a trace line for the pos marker or not
        """
        for label, orig_label in zip(new_labels, orig_labels):
            if not self.label_used(label):
                continue
            if skip_to_label:
                # jump over the whole interception block
                self.put_goto(skip_to_label)
                skip_to_label = None

            if pos is not None:
                self.mark_pos(pos, trace=trace)
            self.put_label(label)
            yield (label, orig_label)
            self.put_goto(orig_label)

    # Functions delegated to function scope
    def new_label(self, name=None):    return self.funcstate.new_label(name)
    def new_error_label(self, *args):  return self.funcstate.new_error_label(*args)
    def new_yield_label(self, *args):  return self.funcstate.new_yield_label(*args)
    def get_loop_labels(self):         return self.funcstate.get_loop_labels()
    def set_loop_labels(self, labels): return self.funcstate.set_loop_labels(labels)
    def new_loop_labels(self, *args):  return self.funcstate.new_loop_labels(*args)
    def get_all_labels(self):          return self.funcstate.get_all_labels()
    def set_all_labels(self, labels):  return self.funcstate.set_all_labels(labels)
    def all_new_labels(self):          return self.funcstate.all_new_labels()
    def use_label(self, lbl):          return self.funcstate.use_label(lbl)
    def label_used(self, lbl):         return self.funcstate.label_used(lbl)


    def enter_cfunc_scope(self, scope):
        self.funcstate = FunctionState(self, scope=scope)

    def exit_cfunc_scope(self):
        if self.funcstate is None:
            return
        self.funcstate.validate_exit()
        self.funcstate = None

    def start_initcfunc(self, signature, scope=None, refnanny=False):
        """
        Init code helper function to start a cfunc scope and generate
        the prototype and function header ("static SIG {") of the function.
        """
        proto = self.globalstate.parts['initfunc_declarations']
        proto.putln(f"static CYTHON_SMALL_CODE {signature}; /*proto*/")
        self.enter_cfunc_scope(scope)
        self.putln("")
        self.putln(f"static {signature} {{")
        if refnanny:
            self.put_declare_refcount_context()

    def start_slotfunc(self, class_scope, return_type, c_slot_name, args_signature, needs_funcstate=True, needs_prototype=False):
        # Slot functions currently live in the class scope as they don't have direct access to the module state.
        slotfunc_cname = class_scope.mangle_internal(c_slot_name)
        declaration = f"static {return_type.declaration_code(slotfunc_cname)}({args_signature})"

        if needs_prototype:
            self.globalstate['decls'].putln(declaration.replace("CYTHON_UNUSED ", "") + "; /*proto*/")
        if needs_funcstate:
            self.enter_cfunc_scope(class_scope)
        self.putln("")
        self.putln(declaration + " {")

    # constant handling

    def get_py_int(self, str_value, longness):
        return self.name_in_module_state(
            self.globalstate.get_int_const(str_value, longness).cname
        )

    def get_py_float(self, str_value, value_code):
        return self.name_in_module_state(
            self.globalstate.get_float_const(str_value, value_code).cname
        )

    def get_py_const(self, prefix, dedup_key=None):
        return self.name_in_module_state(
            self.globalstate.get_py_const(prefix, dedup_key)
        )

    def get_string_const(self, text):
        return self.globalstate.get_string_const(text).cname

    def get_pyunicode_ptr_const(self, text):
        return self.globalstate.get_pyunicode_ptr_const(text)

    def get_py_string_const(self, text, identifier=None):
        cname = self.globalstate.get_py_string_const(
            text, identifier).cname
        return self.name_in_module_state(cname)

    def get_py_codeobj_const(self, node):
        return self.name_in_module_state(self.globalstate.get_py_codeobj_const(node))

    def get_argument_default_const(self, type):
        return self.name_in_module_state(self.globalstate.get_argument_default_const(type).cname)

    def intern(self, text):
        return self.get_py_string_const(text)

    def intern_identifier(self, text):
        return self.get_py_string_const(text, identifier=True)

    def get_cached_constants_writer(self, target=None):
        return self.globalstate.get_cached_constants_writer(target)

    def name_in_module_state(self, cname):
        if self.funcstate.scope is None:
            # This is a mess. For example, within the codeobj generation
            # funcstate.scope is None while evaluating the strings, but not while
            # evaluating the code objects themselves. Right now it doesn't matter
            # because it all ends up going to the same place, but to actually turn
            # it into something useful this mess will need to be fixed.
            return self.name_in_main_c_code_module_state(cname)
        return self.funcstate.scope.name_in_module_state(cname)

    @staticmethod
    def name_in_main_c_code_module_state(cname):
        # The functions where this applies to have the modulestate passed
        # as an argument to them and so it's better use that argument than
        # to try to get it from a global variable.
        return f"{Naming.modulestatevalue_cname}->{cname}"

    @staticmethod
    def name_in_slot_module_state(cname):
        # TODO - eventually this will go through PyType_GetModuleByDef
        # in cases where it's supported.
        return f"{Naming.modulestateglobal_cname}->{cname}"

    def namespace_cname_in_module_state(self, scope):
        if scope.is_py_class_scope:
            return scope.namespace_cname
        else:
            return self.name_in_module_state(scope.namespace_cname)

    def typeptr_cname_in_module_state(self, type):
        if type.is_extension_type:
            return self.name_in_module_state(type.typeptr_cname)
        else:
            return type.typeptr_cname

    # code generation

    def putln(self, code="", safe=False):
        if self.last_pos and self.bol:
            self.emit_marker()
        if self.code_config.emit_linenums and self.last_marked_pos:
            source_desc, line, _ = self.last_marked_pos
            self._write_lines(f'\n#line {line} "{source_desc.get_escaped_description()}"\n')
        if code:
            if safe:
                self.put_safe(code)
            else:
                self.put(code)
        self._write_lines("\n")
        self.bol = 1

    def mark_pos(self, pos, trace=True):
        if pos is None:
            return
        if self.last_marked_pos and self.last_marked_pos[:2] == pos[:2]:
            return
        self.last_pos = (pos, trace)

    @cython.final
    def emit_marker(self):
        pos, trace = self.last_pos
        self.last_marked_pos = pos
        self.last_pos = None
        self._write_lines("\n")
        if self.code_config.emit_code_comments:
            self.indent()
            self._write_lines(self._build_marker(pos))
        if trace:
            self.write_trace_line(pos)

    @cython.final
    def write_trace_line(self, pos):
        if self.funcstate and self.funcstate.can_trace and self.globalstate.directives['linetrace']:
            self.indent()
            self._write_lines(
                f'__Pyx_TraceLine({pos[1]:d},{self.pos_to_offset(pos):d},{not self.funcstate.gil_owned:d},{self.error_goto(pos)})\n')

    @cython.final
    def _build_marker(self, pos):
        source_desc, line, col = pos
        assert isinstance(source_desc, SourceDescriptor)
        contents = self.globalstate.commented_file_contents(source_desc)
        lines = contents[max(0, line-3):line]  # line numbers start at 1
        lines[-1] += '             # <<<<<<<<<<<<<<'
        lines += contents[line:line+2]
        code = "\n".join(lines)
        return f'/* "{source_desc.get_escaped_description()}":{line:d}\n{code}\n*/\n'

    def put_safe(self, code):
        # put code, but ignore {}
        self.write(code)
        self.bol = 0

    @cython.final
    def put_or_include(self, code, name):
        include_dir = self.globalstate.common_utility_include_dir
        if include_dir and len(code) > 1024:
            hash = hashlib.sha256(code.encode('utf8')).hexdigest()
            include_file = f"{name}_{hash}.h"
            path = os.path.join(include_dir, include_file)
            if not os.path.exists(path):
                tmp_path = f'{path}.tmp{os.getpid()}'
                done = False
                try:
                    with Utils.open_new_file(tmp_path) as f:
                        f.write(code)
                    shutil.move(tmp_path, path)
                    done = True
                except (FileExistsError, PermissionError):
                    # If a different process created the file faster than us,
                    # renaming can fail on Windows.  It's ok if the file is there now.
                    if not os.path.exists(path):
                        raise
                finally:
                    if not done and os.path.exists(tmp_path):
                        os.unlink(tmp_path)
            # We use forward slashes in the include path to assure identical code generation
            # under Windows and Posix.  C/C++ compilers should still understand it.
            c_path = path.replace('\\', '/')
            code = f'#include "{c_path}"\n'
        self.put_multilines(code)

    @cython.final
    def put_multilines(self, code):
        # We assume that the code is consistently indented and just needs overall indenting.
        # We also don't need to indent the first line since "self.put()" will do it for us.
        if self.level and '\n' in code:
            code = ("  " * self.level).join(code.splitlines(keepends=True))
        self.put(code)

    def put(self, code):
        fix_indent = False
        dl: cython.Py_ssize_t = code.count("{") if "{" in code else 0
        if "}" in code:
            dl -= code.count("}")
            if dl < 0:
                self.level += dl
            elif dl == 0 and code[0] == "}":
                # special cases like "} else {" need a temporary dedent
                fix_indent = True
                self.level -= 1
        if self.bol:
            self.indent()
        self.write(code)
        self.bol = 0
        if dl > 0:
            self.level += dl
        elif fix_indent:
            self.level += 1

    def put_code_here(self, utility: UtilityCode):
        # Puts the impl section of the utility code directly to the current position.
        # Ensure we don't have a proto section (but do allow init and cleanup sections
        # because they might be useful in future).
        assert not utility.proto, utility.name
        utility._put_code_section(self, self.globalstate, "impl")
        utility._put_init_code_section(self.globalstate)
        if utility.cleanup and Options.generate_cleanup_code:
            utility._put_code_section(
                self.globalstate['cleanup_globals'], self.globalstate, "cleanup")

    @cython.final
    def increase_indent(self):
        self.level += 1

    @cython.final
    def decrease_indent(self):
        self.level -= 1

    @cython.final
    def begin_block(self):
        self.putln("{")
        self.increase_indent()

    @cython.final
    def end_block(self):
        self.decrease_indent()
        self.putln("}")

    @cython.final
    def indent(self):
        self._write_to_buffer("  " * self.level)

    def get_py_version_hex(self, pyversion):
        return "0x%02X%02X%02X%02X" % (tuple(pyversion) + (0,0,0,0))[:4]

    def put_label(self, lbl):
        if lbl in self.funcstate.labels_used:
            self.putln("%s:;" % lbl)

    def put_goto(self, lbl):
        self.funcstate.use_label(lbl)
        self.putln("goto %s;" % lbl)

    def put_var_declaration(self, entry, storage_class="",
                            dll_linkage=None, definition=True):
        #print "Code.put_var_declaration:", entry.name, "definition =", definition ###
        if entry.visibility == 'private' and not (definition or entry.defined_in_pxd):
            #print "...private and not definition, skipping", entry.cname ###
            return
        if entry.visibility == "private" and not entry.used:
            #print "...private and not used, skipping", entry.cname ###
            return
        if not entry.cf_used:
            self.put('CYTHON_UNUSED ')
        if storage_class:
            self.put("%s " % storage_class)
        if entry.is_cpp_optional:
            self.put(entry.type.cpp_optional_declaration_code(
                entry.cname, dll_linkage=dll_linkage))
        else:
            self.put(entry.type.declaration_code(
                entry.cname, dll_linkage=dll_linkage))
        if entry.init is not None:
            self.put_safe(" = %s" % entry.type.literal_code(entry.init))
        elif entry.type.is_pyobject:
            self.put(" = NULL")
        self.putln(";")
        self.globalstate.use_entry_utility_code(entry)

    def put_temp_declarations(self, func_context: FunctionState):
        for name, type, manage_ref, static in func_context.temps_allocated:
            if type.is_cpp_class and not type.is_fake_reference and func_context.scope.directives['cpp_locals']:
                decl = type.cpp_optional_declaration_code(name)
            else:
                decl = type.declaration_code(name)
            if type.is_pyobject:
                self.putln("%s = NULL;" % decl)
            elif type.is_memoryviewslice:
                self.putln("%s = %s;" % (decl, type.literal_code(type.default_value)))
            else:
                self.putln("%s%s;" % (static and "static " or "", decl))

        if func_context.should_declare_error_indicator:
            if self.funcstate.uses_error_indicator:
                unused = ''
            else:
                unused = 'CYTHON_UNUSED '
            # Initialize these variables to silence compiler warnings
            self.putln("%sint %s = 0;" % (unused, Naming.lineno_cname))
            self.putln("%sconst char *%s = NULL;" % (unused, Naming.filename_cname))
            self.putln("%sint %s = 0;" % (unused, Naming.clineno_cname))

    def put_generated_by(self):
        self.putln(Utils.GENERATED_BY_MARKER)
        self.putln("")

    def put_h_guard(self, guard):
        self.putln("#ifndef %s" % guard)
        self.putln("#define %s" % guard)

    def unlikely(self, cond):
        if Options.gcc_branch_hints:
            return 'unlikely(%s)' % cond
        else:
            return cond

    def build_function_modifiers(self, modifiers, mapper=modifier_output_mapper):
        if not modifiers:
            return ''
        return '%s ' % ' '.join([mapper(m,m) for m in modifiers])

    # Python objects and reference counting

    def entry_as_pyobject(self, entry):
        type = entry.type
        if (not entry.is_self_arg and not entry.type.is_complete()
                or entry.type.is_extension_type):
            return "(PyObject *)" + entry.cname
        else:
            return entry.cname

    def as_pyobject(self, cname, type):
        from .PyrexTypes import py_object_type, typecast
        return typecast(py_object_type, type, cname)

    def put_gotref(self, cname, type):
        type.generate_gotref(self, cname)

    def put_giveref(self, cname, type):
        type.generate_giveref(self, cname)

    def put_xgiveref(self, cname, type):
        type.generate_xgiveref(self, cname)

    def put_xgotref(self, cname, type):
        type.generate_xgotref(self, cname)

    def put_incref(self, cname, type, nanny=True):
        # Note: original put_Memslice_Incref/Decref also added in some utility code
        # this is unnecessary since the relevant utility code is loaded anyway if a memoryview is used
        # and so has been removed. However, it's potentially a feature that might be useful here
        type.generate_incref(self, cname, nanny=nanny)

    def put_xincref(self, cname, type, nanny=True):
        type.generate_xincref(self, cname, nanny=nanny)

    def put_decref(self, cname, type, nanny=True, have_gil=True):
        type.generate_decref(self, cname, nanny=nanny, have_gil=have_gil)

    def put_xdecref(self, cname, type, nanny=True, have_gil=True):
        type.generate_xdecref(self, cname, nanny=nanny, have_gil=have_gil)

    def put_decref_clear(self, cname, type, clear_before_decref=False, nanny=True, have_gil=True):
        type.generate_decref_clear(self, cname, clear_before_decref=clear_before_decref,
                              nanny=nanny, have_gil=have_gil)

    def put_xdecref_clear(self, cname, type, clear_before_decref=False, nanny=True, have_gil=True):
        type.generate_xdecref_clear(self, cname, clear_before_decref=clear_before_decref,
                              nanny=nanny, have_gil=have_gil)

    def put_decref_set(self, cname, type, rhs_cname):
        type.generate_decref_set(self, cname, rhs_cname)

    def put_xdecref_set(self, cname, type, rhs_cname):
        type.generate_xdecref_set(self, cname, rhs_cname)

    def put_incref_memoryviewslice(self, slice_cname, type, have_gil):
        # TODO ideally this would just be merged into "put_incref"
        type.generate_incref_memoryviewslice(self, slice_cname, have_gil=have_gil)

    def put_var_incref_memoryviewslice(self, entry, have_gil):
        self.put_incref_memoryviewslice(entry.cname, entry.type, have_gil=have_gil)

    def put_var_gotref(self, entry):
        self.put_gotref(entry.cname, entry.type)

    def put_var_giveref(self, entry):
        self.put_giveref(entry.cname, entry.type)

    def put_var_xgotref(self, entry):
        self.put_xgotref(entry.cname, entry.type)

    def put_var_xgiveref(self, entry):
        self.put_xgiveref(entry.cname, entry.type)

    def put_var_incref(self, entry, **kwds):
        self.put_incref(entry.cname, entry.type, **kwds)

    def put_var_xincref(self, entry, **kwds):
        self.put_xincref(entry.cname, entry.type, **kwds)

    def put_var_decref(self, entry, **kwds):
        self.put_decref(entry.cname, entry.type, **kwds)

    def put_var_xdecref(self, entry, **kwds):
        self.put_xdecref(entry.cname, entry.type, **kwds)

    def put_var_decref_clear(self, entry, **kwds):
        self.put_decref_clear(entry.cname, entry.type, clear_before_decref=entry.in_closure, **kwds)

    def put_var_decref_set(self, entry, rhs_cname, **kwds):
        self.put_decref_set(entry.cname, entry.type, rhs_cname, **kwds)

    def put_var_xdecref_set(self, entry, rhs_cname, **kwds):
        self.put_xdecref_set(entry.cname, entry.type, rhs_cname, **kwds)

    def put_var_xdecref_clear(self, entry, **kwds):
        self.put_xdecref_clear(entry.cname, entry.type, clear_before_decref=entry.in_closure, **kwds)

    def put_var_decrefs(self, entries, used_only = 0):
        for entry in entries:
            if not used_only or entry.used:
                if entry.xdecref_cleanup:
                    self.put_var_xdecref(entry)
                else:
                    self.put_var_decref(entry)

    def put_var_xdecrefs(self, entries):
        for entry in entries:
            self.put_var_xdecref(entry)

    def put_var_xdecrefs_clear(self, entries):
        for entry in entries:
            self.put_var_xdecref_clear(entry)

    def put_make_object_deferred(self, cname):
        # Deferred reference counting is probably only worthwhile on global classes
        # that we expect to be long-term accessible.  So for now exclude it if not
        # at class or module scope.
        if (self.funcstate.scope.is_module_scope or
                self.funcstate.scope.is_c_class_scope or
                self.funcstate.scope.is_py_class_scope):
            self.putln("#if CYTHON_COMPILING_IN_CPYTHON && PY_VERSION_HEX >= 0x030E0000")
            self.putln(f"PyUnstable_Object_EnableDeferredRefcount({cname});")
            self.putln("#endif")

    def put_init_to_py_none(self, cname, type, nanny=True):
        from .PyrexTypes import py_object_type, typecast
        py_none = typecast(type, py_object_type, "Py_None")
        if nanny:
            self.putln("%s = %s; __Pyx_INCREF(Py_None);" % (cname, py_none))
        else:
            self.putln("%s = %s; Py_INCREF(Py_None);" % (cname, py_none))

    def put_init_var_to_py_none(self, entry, template = "%s", nanny=True):
        code = template % entry.cname
        #if entry.type.is_extension_type:
        #    code = "((PyObject*)%s)" % code
        self.put_init_to_py_none(code, entry.type, nanny)
        if entry.in_closure:
            self.put_giveref('Py_None')

    def put_pymethoddef(self, entry, term, allow_skip=True, wrapper_code_writer=None):
        is_number_slot = False
        if entry.is_special or entry.name == '__getattribute__':
            from . import TypeSlots
            if entry.name not in special_py_methods:
                if TypeSlots.is_binop_number_slot(entry.name):
                    # It's useful if numeric binops are created with meth coexist
                    # so they can be called directly by looking up the name, skipping the
                    # dispatch wrapper that enables the reverse slots.  This is most useful
                    # when c_api_binop_methods is False, but there's no reason not to do it
                    # all the time
                    is_number_slot = True
                elif entry.name == '__getattr__' and not self.globalstate.directives['fast_getattr']:
                    pass
                # Python's typeobject.c will automatically fill in our slot
                # in add_operators() (called by PyType_Ready) with a value
                # that's better than ours.
                elif allow_skip:
                    return

        method_flags = entry.signature.method_flags()
        if not method_flags:
            return
        if entry.is_special:
            method_flags += [TypeSlots.method_coexist]
        func_ptr = wrapper_code_writer.put_pymethoddef_wrapper(entry) if wrapper_code_writer else entry.func_cname
        # Add required casts, but try not to shadow real warnings.
        cast = entry.signature.method_function_type()
        if cast != 'PyCFunction':
            func_ptr = '(void(*)(void))(%s)%s' % (cast, func_ptr)
        entry_name = entry.name.as_c_string_literal()
        if is_number_slot:
            # Unlike most special functions, binop numeric operator slots are actually generated here
            # (to ensure that they can be looked up). However, they're sometimes guarded by the preprocessor
            # so a bit of extra logic is needed
            slot = TypeSlots.get_slot_table(self.globalstate.directives).get_slot_by_method_name(entry.name)
            preproc_guard = slot.preprocessor_guard_code()
            if preproc_guard:
                self.putln(preproc_guard)
        self.putln(
            '{%s, (PyCFunction)%s, %s, %s}%s' % (
                entry_name,
                func_ptr,
                "|".join(method_flags),
                entry.doc_cname if entry.doc else '0',
                term))
        if is_number_slot and preproc_guard:
            self.putln("#endif")

    def put_pymethoddef_wrapper(self, entry):
        func_cname = entry.func_cname
        if entry.is_special:
            method_flags = entry.signature.method_flags() or []
            from .TypeSlots import method_noargs
            if method_noargs in method_flags:
                # Special NOARGS methods really take no arguments besides 'self', but PyCFunction expects one.
                func_cname = Naming.method_wrapper_prefix + func_cname
                self.putln("static PyObject *%s(PyObject *self, CYTHON_UNUSED PyObject *arg) {" % func_cname)
                func_call = "%s(self)" % entry.func_cname
                if entry.name == "__next__":
                    self.putln("PyObject *res = %s;" % func_call)
                    # tp_iternext can return NULL without an exception
                    self.putln("if (!res && !PyErr_Occurred()) { PyErr_SetNone(PyExc_StopIteration); }")
                    self.putln("return res;")
                else:
                    self.putln("return %s;" % func_call)
                self.putln("}")
        return func_cname

    # GIL methods

    def use_fast_gil_utility_code(self):
        if self.globalstate.directives['fast_gil']:
            self.globalstate.use_utility_code(UtilityCode.load_cached("FastGil", "ModuleSetupCode.c"))
        else:
            self.globalstate.use_utility_code(UtilityCode.load_cached("NoFastGil", "ModuleSetupCode.c"))

    def put_ensure_gil(self, declare_gilstate=True, variable=None):
        """
        Acquire the GIL. The generated code is safe even when no PyThreadState
        has been allocated for this thread (for threads not initialized by
        using the Python API). Additionally, the code generated by this method
        may be called recursively.
        """
        if self.globalstate.directives['subinterpreters_compatible'] != 'no':
            from .Errors import warning
            warning(
                self.last_marked_pos,
                "Acquiring the GIL is currently very unlikely to work correctly with subinterpreters.",
                2
            )
        self.globalstate.use_utility_code(
            UtilityCode.load_cached("ForceInitThreads", "ModuleSetupCode.c"))
        self.use_fast_gil_utility_code()
        if not variable:
            variable = '__pyx_gilstate_save'
            if declare_gilstate:
                self.put("PyGILState_STATE ")
        self.putln("%s = __Pyx_PyGILState_Ensure();" % variable)

    def put_release_ensured_gil(self, variable=None):
        """
        Releases the GIL, corresponds to `put_ensure_gil`.
        """
        self.use_fast_gil_utility_code()
        if not variable:
            variable = '__pyx_gilstate_save'
        self.putln("__Pyx_PyGILState_Release(%s);" % variable)

    def put_acquire_freethreading_lock(self):
        self.putln("#if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING")
        self.putln(f"PyMutex_Lock(&{Naming.parallel_freethreading_mutex});")
        self.putln("#endif")

    def put_release_freethreading_lock(self):
        self.putln("#if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING")
        self.putln(f"PyMutex_Unlock(&{Naming.parallel_freethreading_mutex});")
        self.putln("#endif")

    def put_acquire_gil(self, variable=None, unknown_gil_state=True):
        """
        Acquire the GIL. The thread's thread state must have been initialized
        by a previous `put_release_gil`
        """
        self.use_fast_gil_utility_code()
        self.putln("__Pyx_FastGIL_Forget();")
        if variable:
            self.putln('_save = %s;' % variable)
        if unknown_gil_state:
            func_name = "__Pyx_RestoreUnknownThread"
            self.globalstate.use_utility_code(
                UtilityCode.load_cached("ReleaseUnknownGil", "ModuleSetupCode.c"))
        else:
            func_name = "PyEval_RestoreThread"
        self.putln(f"{func_name}(_save);")

    def put_release_gil(self, variable=None, unknown_gil_state=True):
        "Release the GIL, corresponds to `put_acquire_gil`."
        self.use_fast_gil_utility_code()
        if unknown_gil_state:
            self.globalstate.use_utility_code(
                UtilityCode.load_cached("ReleaseUnknownGil", "ModuleSetupCode.c"))
            func_name = "__Pyx_SaveUnknownThread"
            result_type = "__Pyx_UnknownThreadState"
        else:
            func_name = "PyEval_SaveThread"
            result_type = "PyThreadState *"
        self.putln(f"{result_type} _save;")
        self.putln(f"_save = {func_name}();")
        if variable:
            self.putln('%s = _save;' % variable)
        self.putln("__Pyx_FastGIL_Remember();")

    def declare_gilstate(self):
        self.putln("PyGILState_STATE __pyx_gilstate_save;")

    # error handling

    def put_error_if_neg(self, pos, value):
        # TODO this path is almost _never_ taken, yet this macro makes is slower!
        # return self.putln("if (unlikely(%s < 0)) %s" % (value, self.error_goto(pos)))
        return self.putln("if (%s < (0)) %s" % (value, self.error_goto(pos)))

    def put_error_if_unbound(self, pos, entry, in_nogil_context=False, unbound_check_code=None):
        nogil_tag = "Nogil" if in_nogil_context else ""
        if entry.from_closure:
            func = "RaiseClosureNameError"
        elif entry.type.is_cpp_class and entry.is_cglobal:
            func = "RaiseCppGlobalNameError"
        elif entry.type.is_cpp_class and entry.is_variable and not entry.is_member and entry.scope.is_c_class_scope:
            # there doesn't seem to be a good way to detecting an instance-attribute of a C class
            # (is_member is only set for class attributes)
            func = "RaiseCppAttributeError"
        else:
            func = "RaiseUnboundLocalError"

        self.globalstate.use_utility_code(
                UtilityCode.load_cached(f"{func}{nogil_tag}", "ObjectHandling.c"))

        if not unbound_check_code:
            unbound_check_code = entry.type.check_for_null_code(entry.cname)
        self.putln('if (unlikely(!%s)) { %s(%s); %s }' % (
                                unbound_check_code,
                                f"__Pyx_{func}{nogil_tag}",
                                entry.name.as_c_string_literal(),
                                self.error_goto(pos)))

    def set_error_info(self, pos, used=False):
        self.funcstate.should_declare_error_indicator = True
        if used:
            self.funcstate.uses_error_indicator = True
        return "__PYX_MARK_ERR_POS(%s, %s)" % (
            self.lookup_filename(pos[0]),
            pos[1])

    def error_goto(self, pos, used=True):
        lbl = self.funcstate.error_label
        self.funcstate.use_label(lbl)
        if pos is None:
            return 'goto %s;' % lbl
        self.funcstate.should_declare_error_indicator = True
        if used:
            self.funcstate.uses_error_indicator = True
        return "__PYX_ERR(%s, %s, %s)" % (
            self.lookup_filename(pos[0]),
            pos[1],
            lbl)

    def error_goto_if(self, cond, pos):
        return "if (%s) %s" % (self.unlikely(cond), self.error_goto(pos))

    def error_goto_if_null(self, cname, pos):
        return self.error_goto_if("!%s" % cname, pos)

    def error_goto_if_neg(self, cname, pos):
        # Add extra parentheses to silence clang warnings about constant conditions.
        return self.error_goto_if("(%s < 0)" % cname, pos)

    def error_goto_if_PyErr(self, pos):
        return self.error_goto_if("PyErr_Occurred()", pos)

    def lookup_filename(self, filename):
        return self.globalstate.lookup_filename(filename)

    def put_declare_refcount_context(self):
        self.putln('__Pyx_RefNannyDeclarations')

    def put_setup_refcount_context(self, name, acquire_gil=False):
        name = name.as_c_string_literal()  # handle unicode names
        if acquire_gil:
            self.globalstate.use_utility_code(
                UtilityCode.load_cached("ForceInitThreads", "ModuleSetupCode.c"))
        self.putln('__Pyx_RefNannySetupContext(%s, %d);' % (name, acquire_gil and 1 or 0))

    def put_finish_refcount_context(self, nogil=False):
        self.putln("__Pyx_RefNannyFinishContextNogil()" if nogil else "__Pyx_RefNannyFinishContext();")

    def put_add_traceback(self, qualified_name, include_cline=True):
        """
        Build a Python traceback for propagating exceptions.

        qualified_name should be the qualified name of the function.
        """
        qualified_name = qualified_name.as_c_string_literal()  # handle unicode names
        format_tuple = (
            qualified_name,
            Naming.clineno_cname if include_cline else 0,
            Naming.lineno_cname,
            Naming.filename_cname,
        )

        self.funcstate.uses_error_indicator = True
        self.putln('__Pyx_AddTraceback(%s, %s, %s, %s);' % format_tuple)

    def put_unraisable(self, qualified_name, nogil=False):
        """
        Generate code to print a Python warning for an unraisable exception.

        qualified_name should be the qualified name of the function.
        """
        self.funcstate.uses_error_indicator = True
        self.putln('__Pyx_WriteUnraisable("%s", %s, %s, %s, %d, %d);' % (
            qualified_name,
            Naming.clineno_cname,
            Naming.lineno_cname,
            Naming.filename_cname,
            self.globalstate.directives['unraisable_tracebacks'],
            nogil,
        ))
        self.globalstate.use_utility_code(
            UtilityCode.load_cached("WriteUnraisableException", "Exceptions.c"))

    def is_tracing(self):
        return self.globalstate.directives['profile'] or self.globalstate.directives['linetrace']

    def pos_to_offset(self, pos):
        """
        Calculate a fake 'instruction offset' from a node position as 31 bit int (32 bit signed).
        """
        scope = self.funcstate.scope
        while scope and pos not in scope.node_positions_to_offset:
            scope = scope.parent_scope
        return scope.node_positions_to_offset[pos] if scope else 0

    def put_trace_declarations(self, is_generator=False):
        self.putln('__Pyx_TraceDeclarationsGen' if is_generator else '__Pyx_TraceDeclarationsFunc')

    def put_trace_frame_init(self, codeobj=None):
        if codeobj:
            self.putln('__Pyx_TraceFrameInit(%s)' % codeobj)

    def put_trace_start(self, name, pos, nogil=False, is_generator=False, is_cpdef_func=False):
        trace_func = "__Pyx_TraceStartGen" if is_generator else "__Pyx_TraceStartFunc"
        self.putln(
            f'{trace_func}('
            f'{name.as_c_string_literal()}, '
            f'{Naming.filetable_cname}[{self.lookup_filename(pos[0])}], '
            f'{pos[1]}, '
            f'{self.pos_to_offset(pos):d}, '
            f'{nogil:d}, '
            f'{Naming.skip_dispatch_cname if is_cpdef_func else "0"}, '
            f'{self.error_goto(pos)}'
            ');'
        )

    def put_trace_exit(self, nogil=False):
        self.putln(f"__Pyx_PyMonitoring_ExitScope({bool(nogil):d});")

    def put_trace_yield(self, retvalue_cname, pos):
        error_goto = self.error_goto(pos)
        self.putln(f"__Pyx_TraceYield({retvalue_cname}, {self.pos_to_offset(pos)}, {error_goto});")

    def put_trace_resume(self, pos):
        scope = self.funcstate.scope
        # pos[1] is probably not the first line, so try to find the first line of the generator function.
        first_line = scope.scope_class.pos[1] if scope.scope_class else pos[1]
        name = scope.name.as_c_string_literal()
        filename_index = self.lookup_filename(pos[0])
        error_goto = self.error_goto(pos)
        self.putln(
            '__Pyx_TraceResumeGen('
            f'{name}, '
            f'{Naming.filetable_cname}[{filename_index}], '
            f'{first_line}, '
            f'{self.pos_to_offset(pos)}, '
            f'{error_goto}'
            ');'
        )

    def put_trace_exception(self, pos, reraise=False, fresh=False):
        self.putln(f"__Pyx_TraceException({self.pos_to_offset(pos)}, {bool(reraise):d}, {bool(fresh):d});")

    def put_trace_exception_propagating(self):
        self.putln(f"__Pyx_TraceException({Naming.lineno_cname}, 0, 0);")

    def put_trace_exception_handled(self, pos):
        self.putln(f"__Pyx_TraceExceptionHandled({self.pos_to_offset(pos)});")

    def put_trace_unwind(self, pos, nogil=False):
        self.putln(f"__Pyx_TraceExceptionUnwind({self.pos_to_offset(pos)}, {bool(nogil):d});")

    def put_trace_stopiteration(self, pos, value):
        error_goto = self.error_goto(pos)
        self.putln(f"__Pyx_TraceStopIteration({value}, {self.pos_to_offset(pos)}, {error_goto});")

    def put_trace_return(self, retvalue_cname, pos, return_type=None, nogil=False):
        extra_arg = ""
        trace_func = "__Pyx_TraceReturnValue"

        if return_type is None:
            pass
        elif return_type.is_pyobject:
            retvalue_cname = return_type.as_pyobject(retvalue_cname)
        elif return_type.is_void:
            retvalue_cname = 'Py_None'
        elif return_type.is_string:
            # We don't know if the C string is 0-terminated, but we cannot convert if it's not.
            retvalue_cname = 'Py_None'
        elif return_type.to_py_function:
            trace_func = "__Pyx_TraceReturnCValue"
            extra_arg = f", {return_type.to_py_function}"
        else:
            # We don't have a Python visible return value but we still need to report that we returned.
            # 'None' may not be a misleading (it's false, for one), but it's hopefully better than nothing.
            retvalue_cname = 'Py_None'

        error_handling = self.error_goto(pos)
        self.putln(f"{trace_func}({retvalue_cname}{extra_arg}, {self.pos_to_offset(pos)}, {bool(nogil):d}, {error_handling});")

    def put_cpp_placement_new(self, target,
                              _utility_code=UtilityCode.load("DefaultPlacementNew", "CppSupport.cpp")):
        self.globalstate.use_utility_code(_utility_code)
        self.putln(f"__Pyx_default_placement_construct(&({target}));")

    def putln_openmp(self, string):
        self.putln("#ifdef _OPENMP")
        self.putln(string)
        self.putln("#endif /* _OPENMP */")

    def undef_builtin_expect(self, cond):
        """
        Redefine the macros likely() and unlikely to no-ops, depending on
        condition 'cond'
        """
        self.putln("#if %s" % cond)
        self.putln("    #undef likely")
        self.putln("    #undef unlikely")
        self.putln("    #define likely(x)   (x)")
        self.putln("    #define unlikely(x) (x)")
        self.putln("#endif")

    def redef_builtin_expect(self, cond):
        self.putln("#if %s" % cond)
        self.putln("    #undef likely")
        self.putln("    #undef unlikely")
        self.putln("    #define likely(x)   __builtin_expect(!!(x), 1)")
        self.putln("    #define unlikely(x) __builtin_expect(!!(x), 0)")
        self.putln("#endif")


class PyrexCodeWriter:
    # f                file      output file
    # level            int       indentation level

    def __init__(self, outfile_name):
        self.f = Utils.open_new_file(outfile_name)
        self.level = 0

    def putln(self, code):
        self.f.write("%s%s\n" % (" " * self.level, code))

    def indent(self):
        self.level += 1

    def dedent(self):
        self.level -= 1


class PyxCodeWriter:
    """
    Can be used for writing out some Cython code.
    """

    def __init__(self, buffer=None, indent_level=0, context=None, encoding='ascii'):
        self.buffer = buffer or StringIOTree()
        self.level = indent_level
        self.original_level = indent_level
        self.context = context
        self.encoding = encoding
        self._insertion_points = {}

    def indent(self, levels=1):
        self.level += levels
        return True

    def dedent(self, levels=1):
        self.level -= levels

    @contextmanager
    def indenter(self, line):
        """
        with pyx_code.indenter("for i in range(10):"):
            pyx_code.putln("print i")
        """
        self.putln(line)
        self.indent()
        yield
        self.dedent()

    def empty(self):
        return self.buffer.empty()

    def getvalue(self):
        result = self.buffer.getvalue()
        if isinstance(result, bytes):
            result = result.decode(self.encoding)
        return result

    def putln(self, line, context=None):
        if context is None:
            if self.context is not None:
                context = self.context
        if context is not None:
            line = sub_tempita(line, context)
        # Avoid indenting empty lines.
        self.buffer.write(f"{self.level * '    '}{line}\n" if line else "\n")

    def put_chunk(self, chunk, context=None):
        if context is None:
            if self.context is not None:
                context = self.context
        if context is not None:
            chunk = sub_tempita(chunk, context)

        chunk = _indent_chunk(chunk, self.level * 4)
        self.buffer.write(chunk)

    def insertion_point(self):
        return type(self)(self.buffer.insertion_point(), self.level, self.context)

    def reset(self):
        # resets the buffer so that nothing gets written. Most useful
        # for abandoning all work in a specific insertion point
        self.buffer.reset()
        self.level = self.original_level

    def named_insertion_point(self, name):
        self._insertion_points[name] = self.insertion_point()

    def __getitem__(self, name):
        return self._insertion_points[name]


@cython.final
@cython.ccall
def _indent_chunk(chunk: str, indentation_length: cython.int) -> str:
    """Normalise leading space to the intended indentation and strip empty lines.
    """
    assert '\t' not in chunk
    lines = chunk.splitlines(keepends=True)
    if not lines:
        return chunk
    last_line = lines[-1].rstrip(' ')
    if last_line:
        lines[-1] = last_line
    else:
        del lines[-1]
        if not lines:
            return '\n'

    # Count minimal (non-empty) indentation and strip empty lines.
    min_indentation: cython.int = len(chunk) + 1
    line_indentation: cython.int
    line: str
    i: cython.int
    for i, line in enumerate(lines):
        line_indentation = _count_indentation(line)
        if line_indentation + 1 == len(line):
            lines[i] = '\n'
        elif line_indentation < min_indentation:
            min_indentation = line_indentation

    if min_indentation > len(chunk):
        # All empty lines.
        min_indentation = 0

    if min_indentation < indentation_length:
        add_indent = ' ' * (indentation_length - min_indentation)
        lines = [
            add_indent + line if line != '\n' else '\n'
            for line in lines
        ]
    elif min_indentation > indentation_length:
        start: cython.int = min_indentation - indentation_length
        lines = [
            line[start:] if line != '\n' else '\n'
            for line in lines
        ]

    return ''.join(lines)


@cython.exceptval(-1)
@cython.cfunc
def _count_indentation(s: str) -> cython.int:
    i: cython.int = 0
    ch: cython.Py_UCS4
    for i, ch in enumerate(s):
        if ch != ' ':
            break
    return i


class ClosureTempAllocator:
    def __init__(self, klass):
        self.klass = klass
        self.temps_allocated = {}
        self.temps_free = {}
        self.temps_count = 0

    def reset(self):
        for type, cnames in self.temps_allocated.items():
            self.temps_free[type] = list(cnames)

    def allocate_temp(self, type):
        if type not in self.temps_allocated:
            self.temps_allocated[type] = []
            self.temps_free[type] = []
        elif self.temps_free[type]:
            return self.temps_free[type].pop(0)
        cname = '%s%d' % (Naming.codewriter_temp_prefix, self.temps_count)
        self.klass.declare_var(pos=None, name=cname, cname=cname, type=type, is_cdef=True)
        self.temps_allocated[type].append(cname)
        self.temps_count += 1
        return cname
