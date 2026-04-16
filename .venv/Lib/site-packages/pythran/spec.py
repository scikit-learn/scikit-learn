'''
This module provides a dummy parser for pythran annotations.
    * spec_parser reads the specs from a python module and returns them.
'''
from pythran.types.conversion import pytype_to_pretty_type

from collections import defaultdict
from itertools import product, chain
import re
import ply.lex as lex
import ply.yacc as yacc

from pythran.typing import List, Set, Dict, NDArray, Tuple, Pointer, Fun, Type
from pythran.syntax import PythranSyntaxError
from pythran.config import cfg


def ambiguous_types(ty0, ty1):
    from numpy import complex64, complex128
    from numpy import float32, float64
    from numpy import int8, int16, int32, int64, intp, intc
    from numpy import uint8, uint16, uint32, uint64, uintp, uintc
    try:
        from numpy import complex256, float128
    except ImportError:
        complex256 = complex128
        float128 = float64

    if isinstance(ty0, tuple):
        if len(ty0) != len(ty1):
            return False
        return all(ambiguous_types(t0, t1) for t0, t1 in zip(ty0, ty1))

    ambiguous_float_types = float, float64
    if ty0 in ambiguous_float_types and ty1 in ambiguous_float_types:
        return True

    ambiguous_cplx_types = complex, complex128
    if ty0 in ambiguous_cplx_types and ty1 in ambiguous_cplx_types:
        return True

    ambiguous_int_types = int64, int
    if ty0 in ambiguous_int_types and ty1 in ambiguous_int_types:
        return True

    if type(ty0) is not type(ty1):
        return False

    if not hasattr(ty0, '__args__'):
        return False

    if type(ty0) is NDArray:
        # no ambiguity for dtype
        return ambiguous_types(ty0.__args__[1:], ty1.__args__[1:])
    else:
        return ambiguous_types(ty0.__args__, ty1.__args__)


def istransposed(t):
    if not isinstance(t, NDArray):
        return False
    if len(t.__args__) - 1 != 2:
        return False
    return t.__args__[1] == t.__args__[2] == slice(-1, None, None)


def istransposable(t):
    if not isinstance(t, NDArray):
        return False
    if len(t.__args__) - 1 != 2:
        return False
    return all(s.step == 1 for s in t.__args__[1:])


class Spec(object):
    '''
    Result of spec parsing.

    ``functions'' is a mapping from function name to a tuple of signatures
    ``capsule'' is a mapping from function name to signature
    '''

    def __init__(self, functions, capsules=None, ufuncs=None):
        self.functions = dict(functions)
        self.capsules = capsules or dict()
        self.ufuncs = ufuncs or dict()

        # normalize function signatures
        for fname, signatures in functions.items():
            if not isinstance(signatures, tuple):
                self.functions[fname] = (signatures,)

        if not self:
            from pythran.log import logger
            logger.warning("No pythran specification, "
                           "nothing will be exported")

    def keys(self):
        return chain(self.functions.keys(), self.capsules.keys(),
                     self.ufuncs.keys())

    def __bool__(self):
        return any((self.functions, self.capsules, self.ufuncs))

    __nonzero__ = __bool__

    def to_docstrings(self, docstrings):
        for func_name, signatures in self.functions.items():
            sigdocs = signatures_to_string(func_name, signatures)
            docstring_prototypes = 'Supported prototypes:\n{}'.format(sigdocs)
            docstring_py = docstrings.get(func_name, '')
            if not docstring_py:
                docstring = docstring_prototypes
            else:
                parts = docstring_py.split('\n\n', 1)
                docstring = parts[0] + '\n\n    ' + docstring_prototypes
                if len(parts) == 2:
                    docstring += '\n\n' + parts[1]
            docstrings[func_name] = docstring


class SpecParser(object):

    """
    A parser that scans a file lurking for lines such as the one below.

    It then generates a pythran-compatible signature to inject into compile.
#pythran export a((float,(int, uint8),str list) list list)
#pythran export a(str)
#pythran export a( (str,str), int, int16 list list)
#pythran export a( {str} )
"""

    # lex part
    dtypes = {
        'bool': 'BOOL',
        'byte': 'BYTE',
        'complex': 'COMPLEX',
        'int': 'INT',
        'float': 'FLOAT',
        'uint8': 'UINT8',
        'uint16': 'UINT16',
        'uint32': 'UINT32',
        'uint64': 'UINT64',
        'uintc': 'UINTC',
        'uintp': 'UINTP',
        'int8': 'INT8',
        'int16': 'INT16',
        'int32': 'INT32',
        'int64': 'INT64',
        'intc': 'INTC',
        'intp': 'INTP',
        'float32': 'FLOAT32',
        'float64': 'FLOAT64',
        'float128': 'FLOAT128',
        'complex64': 'COMPLEX64',
        'complex128': 'COMPLEX128',
        'complex256': 'COMPLEX256',
        }

    reserved = {
        '#pythran': 'PYTHRAN',
        'export': 'EXPORT',
        'order': 'ORDER',
        'capsule': 'CAPSULE',
        'ufunc': 'UFUNC',
        'or': 'OR',
        'list': 'LIST',
        'type': 'TYPE',
        'tuple': 'TUPLE',
        'set': 'SET',
        'dict': 'DICT',
        'slice': 'SLICE',
        'str': 'STR',
        'None': 'NONE',
        }
    reserved.update(dtypes)

    tokens = ('IDENTIFIER', 'NUM', 'COLUMN', 'LPAREN', 'RPAREN', 'CRAP', 'OPT',
              'LARRAY', 'RARRAY', 'STAR', 'COMMA') + tuple(reserved.values())

    crap = [tok for tok in tokens if tok != 'PYTHRAN']
    some_crap = [tok for tok in crap if tok not in ('LPAREN', 'COMMA')]

    # token <> regexp binding
    t_CRAP = r'[^,:\(\)\[\]*?0-9]'
    t_COMMA = r','
    t_COLUMN = r':'
    t_LPAREN = r'\('
    t_RPAREN = r'\)'
    t_RARRAY = r'\]'
    t_LARRAY = r'\['
    t_STAR = r'\*'
    t_OPT = r'\?'
    t_NUM = r'[1-9][0-9]*'

    precedence = (
        ('left', 'OR'),
        ('left', 'LIST', 'DICT', 'SET'),
    )

    def t_IDENTIFER(self, t):
        t.type = SpecParser.reserved.get(t.value, 'IDENTIFIER')
        return t

    t_IDENTIFER.__doc__ = r'\#?[a-zA-Z_][a-zA-Z_0-9]*'

    # skipped characters
    t_ignore = ' \t\n\r'

    # error handling
    def t_error(self, t):
        t.lexer.skip(1)

    # yacc part

    def p_exports(self, p):
        if len(p) > 1:
            isnative = len(p) == 6
            if len(p) == 6:
                target = self.exports
            elif p[3] == "capsule":
                target = self.native_exports
            else:
                target = self.ufunc_exports
            for key, val in p[len(p)-3]:
                target[key] += val

    p_exports.__doc__ = '''exports :
                   | PYTHRAN EXPORT export_list opt_craps exports
                   | PYTHRAN EXPORT CAPSULE export_list opt_craps exports
                   | PYTHRAN EXPORT UFUNC export_list opt_craps exports'''

    def p_export_list(self, p):
        p[0] = (p[1],) if len(p) == 2 else (p[1] + (p[3],))

    p_export_list.__doc__ = '''export_list : export
                  | export_list COMMA export'''

    def p_export(self, p):
        # unlikely case: the IDENTIFIER is an otherwise reserved name
        if len(p) > 2:
            sigs = p[3] or ((),)
        else:
            sigs = ()
        p[0] = p[1], sigs
        self.export_info[p[1]] += p.lexpos(1),

    p_export.__doc__ = '''export : IDENTIFIER LPAREN opt_param_types RPAREN
                  | IDENTIFIER
                  | EXPORT LPAREN opt_param_types RPAREN
                  | ORDER LPAREN opt_param_types RPAREN'''

    def p_opt_craps(self, p):
        pass

    p_opt_craps.__doc__ = '''opt_craps :
                     | some_crap opt_all_craps'''

    def p_opt_all_craps(self, p):
        pass

    p_opt_all_craps.__doc__ = '''opt_all_craps :
                     | crap opt_all_craps'''

    def p_crap(self, p):
        pass

    p_crap.__doc__ = 'crap : ' + '\n| '.join(crap)

    def p_some_crap(self, p):
        pass

    p_some_crap.__doc__ = 'some_crap : ' + '\n| '.join(some_crap)

    def p_dtype(self, p):
        import numpy
        p[0] = eval(p[1], numpy.__dict__),

    p_dtype.__doc__ = 'dtype : ' + '\n| '.join(dtypes.values())

    def p_opt_param_types(self, p):
        p[0] = p[1] if len(p) == 2 else tuple()

    p_opt_param_types.__doc__ = '''opt_param_types :
                     | param_types'''

    def p_opt_types(self, p):
        p[0] = p[1] if len(p) == 2 else tuple()

    p_opt_types.__doc__ = '''opt_types :
                     | types'''

    def p_param_types(self, p):
        if len(p) == 2 or (len(p) == 3 and p[2] == ','):
            p[0] = tuple((t,) for t in p[1])
        elif len(p) == 3 and p[2] == '?':
            p[0] = tuple((t,) for t in p[1]) + ((),)
        elif len(p) == 4:
            p[0] = tuple((t,) + ts for t in p[1] for ts in p[3])
        else:
            p[0] = tuple((t,) + ts for t in p[1] for ts in p[4]) + ((),)

    p_param_types.__doc__ = '''param_types : type
                       | type OPT
                       | type COMMA
                       | type COMMA param_types
                       | type OPT COMMA default_types'''

    def p_default_types(self, p):
        if len(p) == 3:
            p[0] = tuple((t,) for t in p[1]) + ((),)
        else:
            p[0] = tuple((t,) + ts for t in p[1] for ts in p[4]) + ((),)
    p_default_types.__doc__ = '''default_types : type OPT
                       | type OPT COMMA default_types'''

    def p_types(self, p):
        if len(p) == 2:
            p[0] = tuple((t,) for t in p[1])
        else:
            p[0] = tuple((t,) + ts for t in p[1] for ts in p[3])

    p_types.__doc__ = '''types : type
                 | type COMMA types'''

    def p_array_type(self, p):
        if len(p) == 2:
            p[0] = p[1][0],
        elif len(p) == 5 and p[4] == ']':
            def args(t):
                return t.__args__ if isinstance(t, NDArray) else (t,)
            p[0] = tuple(NDArray[args(t) + p[3]] for t in p[1])

    p_array_type.__doc__ = '''array_type : dtype
                | array_type LARRAY array_indices RARRAY'''

    def p_type(self, p):
        if len(p) == 2:
            p[0] = p[1],
        elif len(p) == 3 and p[2] == 'list':
            p[0] = tuple(List[t] for t in p[1])
        elif len(p) == 3 and p[2] == 'set':
            p[0] = tuple(Set[t] for t in p[1])
        elif len(p) == 3 and p[2] == 'type':
            p[0] = tuple(Type[t] for t in p[1])
        elif len(p) == 3:
            if p[2] is None:
                expanded = []
                for nd in p[1]:
                    expanded.append(nd)
                    if istransposable(nd):
                        expanded.append(NDArray[nd.__args__[0], -1::, -1::])
                p[0] = tuple(expanded)
            elif p[2] == "F":
                for nd in p[1]:
                    if not istransposable(nd):
                        msg = ("Invalid Pythran spec. F order is only valid "
                               "for 2D plain arrays")
                        self.p_error(p, msg)
                p[0] = tuple(NDArray[nd.__args__[0], -1::, -1::]
                             for nd in p[1])
            else:
                p[0] = p[1]
        elif len(p) == 5 and p[4] == ')':
            p[0] = tuple(Fun[args, r]
                         for r in p[1]
                         for args in (product(*p[3])
                                      if len(p[3]) > 1 else p[3]))
        elif len(p) == 5:
            p[0] = tuple(Dict[k, v] for k in p[1] for v in p[3])
        elif len(p) == 4 and p[2] == 'or':
            p[0] = p[1] + p[3]
        elif len(p) == 4 and p[3] == ')':
            p[0] = tuple(Tuple[t] for t in p[2])
        elif len(p) == 4 and p[3] == ']':
            p[0] = p[2]
        else:
            msg = "Invalid Pythran spec. Unknown text '{0}'".format(p.value)
            self.p_error(p, msg)

    p_type.__doc__ = '''type : term
                | array_type opt_order
                | pointer_type
                | type LIST
                | type SET
                | generic_type TYPE
                | type LPAREN opt_types RPAREN
                | type COLUMN type DICT
                | LPAREN types RPAREN
                | LARRAY type RARRAY
                | type OR type
                '''

    def p_generic_type(self, p):
        if len(p) == 2:
            p[0] = {'list': list, 'set': set, 'dict': dict, 'tuple': tuple}.get(p[1], p[1]),
        elif len(p) == 3 and p[2] == 'type':
            p[0] = tuple(Type[t] for t in p[1])
        elif len(p) == 4 and p[2] == 'or':
            p[0] = p[1] + p[3]
        else:
            msg = "Invalid Pythran spec. Unknown text '{0}'".format(p.value)
            self.p_error(p, msg)

    p_generic_type.__doc__ = '''generic_type : term
                | LIST
                | SET
                | generic_type TYPE
                | DICT
                | TUPLE
                | generic_type OR generic_type
                '''

    def p_opt_order(self, p):
        if len(p) > 1:
            if p[3] not in 'CF':
                msg = "Invalid Pythran spec. Unknown order '{}'".format(p[3])
                self.p_error(p, msg)
            p[0] = p[3]
        else:
            p[0] = None

    p_opt_order.__doc__ = '''opt_order :
                     | ORDER LPAREN IDENTIFIER RPAREN'''

    def p_pointer_type(self, p):
        p[0] = Pointer[p[1][0]]

    p_pointer_type.__doc__ = '''pointer_type : dtype STAR'''

    def p_array_indices(self, p):
        if len(p) == 2:
            p[0] = p[1],
        else:
            p[0] = (p[1],) + p[3]

    p_array_indices.__doc__ = '''array_indices : array_index
                         | array_index COMMA array_indices'''

    def p_array_index(self, p):
        if len(p) == 3:
            p[0] = slice(0, -1, -1)
        elif len(p) == 1 or p[1] == ':':
            p[0] = slice(0, -1, 1)
        else:
            p[0] = slice(0, int(p[1]), 1)

    p_array_index.__doc__ = '''array_index :
                       | NUM
                       | COLUMN
                       | COLUMN COLUMN'''

    def p_term(self, p):
        if p[1] == 'str':
            p[0] = str
        elif p[1] == 'slice':
            p[0] = slice
        elif p[1] == 'None':
            p[0] = type(None)
        else:
            p[0] = p[1][0]

    p_term.__doc__ = '''term : STR
                | NONE
                | SLICE
                | dtype'''


    def PythranSpecError(self, msg, lexpos=None):
        err = PythranSyntaxError(msg)
        if lexpos is not None:
            line_start = self.input_text.rfind('\n', 0, lexpos) + 1
            err.offset = lexpos - line_start
            err.lineno = 1 + self.input_text.count('\n', 0, lexpos)
        if self.input_file:
            err.filename = self.input_file
        return err

    def p_error(self, p):
        if p.type == 'IDENTIFIER':
            alt = {'double': 'float64',
                   'void': 'None',
                   'char': 'int8',
                   'short': 'int16',
                   'long': 'int64',
                   }.get(p.value)
            if alt:
                hint = " Did you mean `{}`?".format(alt)
            else:
                hint = ''
            raise self.PythranSpecError(
                "Unsupported identifier `{}` at that point.{}".format(p.value,
                                                                     hint),
                p.lexpos)
        else:
            raise self.PythranSpecError(
                "Unexpected token `{}` at that point.".format(p.value),
                p.lexpos)

    def __init__(self):
        self.lexer = lex.lex(module=self, debug=False)
        # Do not write the table for better compatibility across ply version
        self.parser = yacc.yacc(module=self,
                                debug=False,
                                write_tables=False)

    def __call__(self, text, input_file=None):
        self.exports = defaultdict(tuple)
        self.native_exports = defaultdict(tuple)
        self.ufunc_exports = defaultdict(tuple)
        self.export_info = defaultdict(tuple)
        self.input_text = text
        self.input_file = input_file

        lines = []
        in_pythran_export = False
        for line in text.split("\n"):
            if re.match(r'\s*#\s*pythran', line):
                in_pythran_export = True
                lines.append(re.sub(r'\s*#\s*pythran', '#pythran', line))
            elif in_pythran_export:
                stripped = line.strip()
                if stripped.startswith('#'):
                    lines.append(line.replace('#', ''))
                else:
                    in_pythran_export = not stripped
                    lines.append('')
            else:
                in_pythran_export &= not line.strip()
                lines.append('')

        pythran_data = '\n'.join(lines)
        self.parser.parse(pythran_data, lexer=self.lexer, debug=False)

        for key, overloads in self.native_exports.items():
            if len(overloads) > 1:
                msg = "Overloads not supported for capsule '{}'".format(key)
                loc = self.export_info[key][-1]
                raise self.PythranSpecError(msg, loc)
            self.native_exports[key] = overloads[0]

        for key, overloads in self.exports.items():
            if len(overloads) > cfg.getint("typing", "max_export_overloads"):
                raise self.PythranSpecError(
                    "Too many overloads for function '{}', probably due to "
                    "automatic generation of C-style and Fortran-style memory "
                    "layout. Please force a layout using `order(C)` or "
                    "`order(F)` in the array signature".format(key))

            for i, ty_i in enumerate(overloads):
                sty_i = spec_to_string(key, ty_i)
                for ty_j in overloads[i+1:]:
                    sty_j = spec_to_string(key, ty_j)
                    if sty_i == sty_j:
                        msg = "Duplicate export entry {}.".format(sty_i)
                        loc = self.export_info[key][-1]
                        raise self.PythranSpecError(msg, loc)
                    if ambiguous_types(ty_i, ty_j):
                        msg = "Ambiguous overloads\n\t{}\n\t{}.".format(sty_i,
                                                                        sty_j)
                        loc = self.export_info[key][i]
                        raise self.PythranSpecError(msg, loc)

        return Spec(self.exports, self.native_exports, self.ufunc_exports)


class ExtraSpecParser(SpecParser):
    '''
    Extension of SpecParser that works on extra .pythran files
    '''

    def __call__(self, text, input_file=None):
        # make the code looks like a regular pythran file
        text = re.sub(r'^\s*export', '#pythran export', text,
                      flags=re.MULTILINE)
        return super(ExtraSpecParser, self).__call__(text, input_file)


def spec_to_string(function_name, spec):
    arguments_types = [pytype_to_pretty_type(t) for t in spec]
    return '{}({})'.format(function_name, ', '.join(arguments_types))


def signatures_to_string(func_name, signatures):
    # filter out transposed version, they are confusing for some users
    # and can generate very long docstring that break MSVC
    sigdocs = [spec_to_string(func_name, sig) for sig in signatures
               if not any(istransposed(t) for t in sig)]
    if not sigdocs:
        sigdocs = [spec_to_string(func_name, sig) for sig in signatures]
    return ''.join('\n    - ' + sigdoc for sigdoc in sigdocs)


def spec_parser(text):
    return SpecParser()(text)


def parse_pytypes(s):
    fake_def = '#pythran export fake({})'.format(s)
    specs = SpecParser()(fake_def)
    return specs.functions['fake'][0]


def load_specfile(filepath):
    with open(filepath) as fd:
        return ExtraSpecParser()(fd.read(), input_file=filepath)
