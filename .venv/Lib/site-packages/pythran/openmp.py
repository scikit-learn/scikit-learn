'''
This modules contains OpenMP-related stuff.
    * OMPDirective is used to represent OpenMP annotations in the AST
    * GatherOMPData turns OpenMP-like string annotations into metadata
'''

from pythran.passmanager import Transformation
import pythran.metadata as metadata
from pythran.spec import parse_pytypes
from pythran.types.conversion import pytype_to_ctype
from pythran.utils import isstr

from gast import AST
import gast as ast
import re


keywords = {
    'atomic',
    'barrier',
    'capture',
    'cancel',
    'collapse',
    'copyin',
    'copyprivate',
    'critical',
    'declare',
    'default',
    'final',
    'firstprivate',
    'flush',
    'for',
    'if',
    'initializer',
    'lastprivate',
    'master',
    'mergeable',
    'none',
    'nowait',
    'num_threads',
    'omp',
    'ordered',
    'parallel',
    'private',
    'read',
    'reduction',
    'schedule',
    'section',
    'sections',
    'shared',
    'simd',
    'single',
    'task',
    'taskwait',
    'taskyield',
    'threadprivate',
    'untied',
    'update',
    'write'
}

declare_keywords = {
    'omp_in',
    'omp_init',
    'omp_orig',
    'omp_out',
    'omp_priv',
}

reserved_contex = {
    'critical',
    'declare',
    'default',
    'schedule',
    'reduction',
}


def is_declare_typename(s, offset, bounds):
    start = s.rfind(':', 0, offset - 1)
    stop = s.rfind(':', offset + 1)
    if start > 0 and stop > 0:
        bounds.extend((start + 1, stop))
        return True
    else:
        return False


class OMPDirective(AST):
    '''Turn a string into a context-dependent metadata.
    >>> o = OMPDirective("omp for private(a,b) shared(c)")
    >>> o.s
    'omp for private({},{}) shared({})'
    >>> [ type(dep) for dep in o.deps ]
    [<class 'gast.gast.Name'>, <class 'gast.gast.Name'>, \
<class 'gast.gast.Name'>]
    >>> [ dep.id for dep in o.deps ]
    ['a', 'b', 'c']
    '''

    def __init__(self, *args):  # no positional argument to be deep copyable
        super(OMPDirective, self).__init__()
        if not args:
            return

        self.deps = []
        self.private_deps = []
        self.shared_deps = []


        def tokenize(s):
            '''A simple contextual "parser" for an OpenMP string'''
            # not completely satisfying if there are strings in if expressions
            out = ''
            par_count = 0
            curr_index = 0
            in_reserved_context = False
            in_declare = False
            in_shared = in_private = False
            while curr_index < len(s):
                bounds = []
                if in_declare and is_declare_typename(s, curr_index, bounds):
                    start, stop = bounds
                    pytypes = parse_pytypes(s[start:stop])
                    out += ', '.join(map(pytype_to_ctype, pytypes))
                    curr_index = stop
                    continue
                m = re.match(r'^([a-zA-Z_]\w*)', s[curr_index:])
                if m:
                    word = m.group(0)
                    curr_index += len(word)
                    if(in_reserved_context or
                         (in_declare and word in declare_keywords) or
                         (par_count == 0 and word in keywords)):
                        out += word
                        in_reserved_context = word in reserved_contex
                        in_declare |= word == 'declare'
                        in_private |= word == 'private'
                        in_shared |= word == 'shared'
                    else:
                        out += '{}'
                        self.deps.append(ast.Name(word, ast.Load(),
                                                  None, None))
                        isattr = re.match(r'^\s*(\.\s*[a-zA-Z_]\w*)', s[curr_index:])
                        if isattr:
                            attr = isattr.group(0)
                            curr_index += len(attr)
                            self.deps[-1] = ast.Attribute(self.deps[-1],
                                                          attr[1:], ast.Load())
                        if in_private:
                            self.private_deps.append(self.deps[-1])
                        if in_shared:
                            self.shared_deps.append(self.deps[-1])
                elif s[curr_index] == '(':
                    par_count += 1
                    curr_index += 1
                    out += '('
                elif s[curr_index] == ')':
                    par_count -= 1
                    curr_index += 1
                    out += ')'
                    if par_count == 0:
                        in_reserved_context = False
                        in_shared = in_private = False
                else:
                    if s[curr_index] in ',:':
                        in_reserved_context = False
                    out += s[curr_index]
                    curr_index += 1
            return out

        self.s = tokenize(args[0])
        self._fields = ('deps', 'shared_deps', 'private_deps')


##
class GatherOMPData(Transformation):
    '''Walks node and collect string comments looking for OpenMP directives.'''

    # there is a special handling for If and Expr, so not listed here
    statements = ("FunctionDef", "Return", "Delete", "Assign", "AugAssign",
                  "Print", "For", "While", "Raise", "TryExcept", "TryFinally",
                  "Assert", "Import", "ImportFrom", "Pass", "Break",)

    # these fields hold statement lists
    statement_lists = ("body", "orelse", "finalbody",)

    def __init__(self):
        Transformation.__init__(self)
        # Remap self.visit_XXXX() to self.attach_data() generic method
        for s in GatherOMPData.statements:
            setattr(self, "visit_" + s, self.attach_data)
        self.current = list()

    def isompdirective(self, node):
        return isstr(node) and node.value.startswith("omp ")

    def visit_Expr(self, node):
        if self.isompdirective(node.value):
            self.current.append(node.value.value)
            return None
        else:
            self.attach_data(node)
            return node

    def visit_If(self, node):
        if self.isompdirective(node.test):
            self.visit(ast.Expr(node.test))
            return self.visit(ast.If(ast.Constant(1, None),
                                     node.body, node.orelse))
        else:
            return self.attach_data(node)

    def attach_data(self, node):
        '''Generic method called for visit_XXXX() with XXXX in
        GatherOMPData.statements list

        '''
        if self.current:
            for curr in self.current:
                md = OMPDirective(curr)
                metadata.add(node, md)
            self.current = list()
        # add a Pass to hold some directives
        for field_name, field in ast.iter_fields(node):
            if field_name in GatherOMPData.statement_lists:
                if(field and
                   isinstance(field[-1], ast.Expr) and
                   self.isompdirective(field[-1].value)):
                    field.append(ast.Pass())
        self.generic_visit(node)

        # add an If to hold scoping OpenMP directives
        directives = metadata.get(node, OMPDirective)
        field_names = {n for n, _ in ast.iter_fields(node)}
        has_no_scope = field_names.isdisjoint(GatherOMPData.statement_lists)
        if directives and has_no_scope:
            # some directives create a scope, but the holding stmt may not
            # artificially create one here if needed
            sdirective = ''.join(d.s for d in directives)
            scoping = ('parallel', 'task', 'section')
            if any(s in sdirective for s in scoping):
                metadata.clear(node, OMPDirective)
                node = ast.If(ast.Constant(1, None), [node], [])
                for directive in directives:
                    metadata.add(node, directive)

        return node
