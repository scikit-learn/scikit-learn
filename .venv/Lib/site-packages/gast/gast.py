import sys as _sys
import ast as _ast
from ast import boolop, cmpop, excepthandler, expr, expr_context, operator
from ast import slice, stmt, unaryop, mod, AST
from ast import iter_child_nodes, walk

try:
    from ast import TypeIgnore
except ImportError:
    class TypeIgnore(AST):
        pass

try:
    from ast import pattern
except ImportError:
    class pattern(AST):
        pass


try:
    from ast import type_param
except ImportError:
    class type_param(AST):
        pass


def _make_node(Name, Fields, Attributes, Bases):

    # This constructor is used a lot during conversion from ast to gast,
    # then as the primary way to build ast nodes. So we tried to optimized it
    # for speed and not for readability.
    def create_node(self, *args, **kwargs):
        if len(args) > len(Fields):
            raise TypeError(
                "{} constructor takes at most {} positional arguments".
                format(Name, len(Fields)))

        # it's faster to iterate rather than zipping or enumerate
        for i in range(len(args)):
            setattr(self, Fields[i], args[i])
        if kwargs:  # cold branch
            self.__dict__.update(kwargs)

    setattr(_sys.modules[__name__],
            Name,
            type(Name,
                 Bases,
                 {'__init__': create_node,
                  '_fields': Fields,
                  '_attributes': Attributes}))


_nodes = (
    # mod
    ('Module', (('body', 'type_ignores'), (), (mod,))),
    ('Interactive', (('body',), (), (mod,))),
    ('Expression', (('body',), (), (mod,))),
    ('FunctionType', (('argtypes', 'returns'), (), (mod,))),
    ('Suite', (('body',), (), (mod,))),

    # stmt
    ('FunctionDef', (('name', 'args', 'body', 'decorator_list', 'returns',
                      'type_comment', 'type_params'),
                     ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
                     (stmt,))),
    ('AsyncFunctionDef', (('name', 'args', 'body', 'decorator_list', 'returns',
                           'type_comment', 'type_params',),
                          ('lineno', 'col_offset',
                           'end_lineno', 'end_col_offset',),
                          (stmt,))),
    ('ClassDef', (('name', 'bases', 'keywords', 'body', 'decorator_list',
                   'type_params',),
                  ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
                  (stmt,))),
    ('Return', (('value',),
                ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
                (stmt,))),
    ('Delete', (('targets',),
                ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
                (stmt,))),
    ('Assign', (('targets', 'value', 'type_comment'),
                ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
                (stmt,))),
    ('TypeAlias', (('name', 'type_params', 'value'),
                  ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
                  (stmt,))),
    ('AugAssign', (('target', 'op', 'value',),
                   ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
                   (stmt,))),
    ('AnnAssign', (('target', 'annotation', 'value', 'simple',),
                   ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
                   (stmt,))),
    ('Print', (('dest', 'values', 'nl',),
               ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
               (stmt,))),
    ('For', (('target', 'iter', 'body', 'orelse', 'type_comment'),
             ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
             (stmt,))),
    ('AsyncFor', (('target', 'iter', 'body', 'orelse', 'type_comment'),
                  ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
                  (stmt,))),
    ('While', (('test', 'body', 'orelse',),
               ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
               (stmt,))),
    ('If', (('test', 'body', 'orelse',),
            ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
            (stmt,))),
    ('With', (('items', 'body', 'type_comment'),
              ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
              (stmt,))),
    ('AsyncWith', (('items', 'body', 'type_comment'),
                   ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
                   (stmt,))),
    ('Match', (('subject', 'cases'),
                   ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
                   (stmt,))),
    ('Raise', (('exc', 'cause',),
               ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
               (stmt,))),
    ('Try', (('body', 'handlers', 'orelse', 'finalbody',),
             ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
             (stmt,))),
    ('TryStar', (('body', 'handlers', 'orelse', 'finalbody',),
             ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
             (stmt,))),
    ('Assert', (('test', 'msg',),
                ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
                (stmt,))),
    ('Import', (('names',),
                ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
                (stmt,))),
    ('ImportFrom', (('module', 'names', 'level',),
                    ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
                    (stmt,))),
    ('Exec', (('body', 'globals', 'locals',),
              ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
              (stmt,))),
    ('Global', (('names',),
                ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
                (stmt,))),
    ('Nonlocal', (('names',),
                  ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
                  (stmt,))),
    ('Expr', (('value',),
              ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
              (stmt,))),
    ('Pass', ((), ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
              (stmt,))),
    ('Break', ((), ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
               (stmt,))),
    ('Continue', ((),
                  ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
                  (stmt,))),

    # expr

    ('BoolOp', (('op', 'values',),
                ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
                (expr,))),
    ('NamedExpr', (('target', 'value',),
                ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
                (expr,))),
    ('BinOp', (('left', 'op', 'right',),
               ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
               (expr,))),
    ('UnaryOp', (('op', 'operand',),
                 ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
                 (expr,))),
    ('Lambda', (('args', 'body',),
                ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
                (expr,))),
    ('IfExp', (('test', 'body', 'orelse',),
               ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
               (expr,))),
    ('Dict', (('keys', 'values',),
              ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
              (expr,))),
    ('Set', (('elts',),
             ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
             (expr,))),
    ('ListComp', (('elt', 'generators',),
                  ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
                  (expr,))),
    ('SetComp', (('elt', 'generators',),
                 ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
                 (expr,))),
    ('DictComp', (('key', 'value', 'generators',),
                  ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
                  (expr,))),
    ('GeneratorExp', (('elt', 'generators',),
                      ('lineno', 'col_offset',
                       'end_lineno', 'end_col_offset',),
                      (expr,))),
    ('Await', (('value',),
               ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
               (expr,))),
    ('Yield', (('value',),
               ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
               (expr,))),
    ('YieldFrom', (('value',),
                   ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
                   (expr,))),
    ('Compare', (('left', 'ops', 'comparators',),
                 ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
                 (expr,))),
    ('Call', (('func', 'args', 'keywords',),
              ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
              (expr,))),
    ('Repr', (('value',),
              ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
              (expr,))),
    ('FormattedValue', (('value', 'conversion', 'format_spec',),
                        ('lineno', 'col_offset',
                         'end_lineno', 'end_col_offset',),
                        (expr,))),
    ('JoinedStr', (('values',),
                   ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
                   (expr,))),
    ('Constant', (('value', 'kind'),
                  ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
                  (expr,))),
    ('Attribute', (('value', 'attr', 'ctx',),
                   ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
                   (expr,))),
    ('Subscript', (('value', 'slice', 'ctx',),
                   ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
                   (expr,))),
    ('Starred', (('value', 'ctx',),
                 ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
                 (expr,))),
    ('Name', (('id', 'ctx', 'annotation', 'type_comment'),
              ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
              (expr,))),
    ('List', (('elts', 'ctx',),
              ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
              (expr,))),
    ('Tuple', (('elts', 'ctx',),
               ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
               (expr,))),

    # expr_context
    ('Load', ((), (), (expr_context,))),
    ('Store', ((), (), (expr_context,))),
    ('Del', ((), (), (expr_context,))),
    ('AugLoad', ((), (), (expr_context,))),
    ('AugStore', ((), (), (expr_context,))),
    ('Param', ((), (), (expr_context,))),

    # slice
    ('Slice', (('lower', 'upper', 'step'),
               ('lineno', 'col_offset', 'end_lineno', 'end_col_offset',),
               (slice,))),

    # boolop
    ('And', ((), (), (boolop,))),
    ('Or', ((), (), (boolop,))),

    # operator
    ('Add', ((), (), (operator,))),
    ('Sub', ((), (), (operator,))),
    ('Mult', ((), (), (operator,))),
    ('MatMult', ((), (), (operator,))),
    ('Div', ((), (), (operator,))),
    ('Mod', ((), (), (operator,))),
    ('Pow', ((), (), (operator,))),
    ('LShift', ((), (), (operator,))),
    ('RShift', ((), (), (operator,))),
    ('BitOr', ((), (), (operator,))),
    ('BitXor', ((), (), (operator,))),
    ('BitAnd', ((), (), (operator,))),
    ('FloorDiv', ((), (), (operator,))),

    # unaryop
    ('Invert', ((), (), (unaryop, AST,))),
    ('Not', ((), (), (unaryop, AST,))),
    ('UAdd', ((), (), (unaryop, AST,))),
    ('USub', ((), (), (unaryop, AST,))),

    # cmpop
    ('Eq', ((), (), (cmpop,))),
    ('NotEq', ((), (), (cmpop,))),
    ('Lt', ((), (), (cmpop,))),
    ('LtE', ((), (), (cmpop,))),
    ('Gt', ((), (), (cmpop,))),
    ('GtE', ((), (), (cmpop,))),
    ('Is', ((), (), (cmpop,))),
    ('IsNot', ((), (), (cmpop,))),
    ('In', ((), (), (cmpop,))),
    ('NotIn', ((), (), (cmpop,))),

    # comprehension
    ('comprehension', (('target', 'iter', 'ifs', 'is_async'), (), (AST,))),

    # excepthandler
    ('ExceptHandler', (('type', 'name', 'body'),
                       ('lineno', 'col_offset',
                        'end_lineno', 'end_col_offset'),
                       (excepthandler,))),

    # arguments
    ('arguments', (('args', 'posonlyargs', 'vararg', 'kwonlyargs',
                    'kw_defaults', 'kwarg', 'defaults'), (), (AST,))),

    # keyword
    ('keyword', (('arg', 'value'),
                 ('lineno', 'col_offset', 'end_lineno', 'end_col_offset'),
                 (AST,))),

    # alias
    ('alias', (('name', 'asname'),
               ('lineno', 'col_offset', 'end_lineno', 'end_col_offset'),
               (AST,))),

    # withitem
    ('withitem', (('context_expr', 'optional_vars'), (), (AST,))),

    # match_case
    ('match_case', (('pattern', 'guard', 'body'), (), (AST,))),

    # pattern
    ('MatchValue', (('value',),
                    ('lineno', 'col_offset', 'end_lineno', 'end_col_offset'),
                    (pattern,))),
    ('MatchSingleton', (('value',),
                        ('lineno', 'col_offset',
                         'end_lineno', 'end_col_offset'),
                        (pattern,))),
    ('MatchSequence', (('patterns',),
                       ('lineno', 'col_offset',
                        'end_lineno', 'end_col_offset'),
                       (pattern,))),
    ('MatchMapping', (('keys', 'patterns', 'rest'),
                      ('lineno', 'col_offset',
                       'end_lineno', 'end_col_offset'),
                      (pattern,))),
    ('MatchClass', (('cls', 'patterns', 'kwd_attrs', 'kwd_patterns'),
                    ('lineno', 'col_offset',
                     'end_lineno', 'end_col_offset'),
                    (pattern,))),
    ('MatchStar', (('name',),
                   ('lineno', 'col_offset',
                    'end_lineno', 'end_col_offset'),
                   (pattern,))),
    ('MatchAs', (('pattern', 'name'),
                   ('lineno', 'col_offset',
                    'end_lineno', 'end_col_offset'),
                   (pattern,))),
    ('MatchOr', (('patterns',),
                 ('lineno', 'col_offset',
                  'end_lineno', 'end_col_offset'),
                 (pattern,))),

    # type_ignore
    ('type_ignore', ((), ('lineno', 'tag'), (TypeIgnore,))),

    # type_param
    ('TypeVar', (('name', 'bound',),
                 ('lineno', 'col_offset',
                  'end_lineno', 'end_col_offset'),
                 (type_param,))),
    ('ParamSpec', (('name',),
                 ('lineno', 'col_offset',
                  'end_lineno', 'end_col_offset'),
                 (type_param,))),
    ('TypeVarTuple', (('name',),
                 ('lineno', 'col_offset',
                  'end_lineno', 'end_col_offset'),
                 (type_param,))),
    )




for name, descr in _nodes:
    _make_node(name, *descr)

if _sys.version_info.major == 2:
    from .ast2 import ast_to_gast, gast_to_ast
if _sys.version_info.major == 3:
    from .ast3 import ast_to_gast, gast_to_ast


def parse(*args, **kwargs):
    return ast_to_gast(_ast.parse(*args, **kwargs))


def unparse(gast_obj):
    from .unparser import unparse
    return unparse(gast_obj)


def literal_eval(node_or_string):
    if isinstance(node_or_string, AST):
        node_or_string = gast_to_ast(node_or_string)
    return _ast.literal_eval(node_or_string)


def get_docstring(node, clean=True):
    if not isinstance(node, (FunctionDef, ClassDef, Module)):
        raise TypeError("%r can't have docstrings" % node.__class__.__name__)
    if node.body and isinstance(node.body[0], Expr) and \
       isinstance(node.body[0].value, Constant):
        if clean:
            import inspect
            holder = node.body[0].value
            return inspect.cleandoc(getattr(holder, holder._fields[0]))
        return node.body[0].value.s


# the following are directly imported from python3.8's Lib/ast.py  #

def copy_location(new_node, old_node):
    """
    Copy source location (`lineno`, `col_offset`, `end_lineno`, and
    `end_col_offset` attributes) from *old_node* to *new_node* if possible,
    and return *new_node*.
    """
    for attr in 'lineno', 'col_offset', 'end_lineno', 'end_col_offset':
        if attr in old_node._attributes and attr in new_node._attributes \
           and hasattr(old_node, attr):
            setattr(new_node, attr, getattr(old_node, attr))
    return new_node


def fix_missing_locations(node):
    """
    When you compile a node tree with compile(), the compiler expects lineno
    and col_offset attributes for every node that supports them.  This is
    rather tedious to fill in for generated nodes, so this helper adds these
    attributes recursively where not already set, by setting them to the values
    of the parent node.  It works recursively starting at *node*.
    """
    def _fix(node, lineno, col_offset, end_lineno, end_col_offset):
        if 'lineno' in node._attributes:
            if not hasattr(node, 'lineno'):
                node.lineno = lineno
            else:
                lineno = node.lineno
        if 'end_lineno' in node._attributes:
            if not hasattr(node, 'end_lineno'):
                node.end_lineno = end_lineno
            else:
                end_lineno = node.end_lineno
        if 'col_offset' in node._attributes:
            if not hasattr(node, 'col_offset'):
                node.col_offset = col_offset
            else:
                col_offset = node.col_offset
        if 'end_col_offset' in node._attributes:
            if not hasattr(node, 'end_col_offset'):
                node.end_col_offset = end_col_offset
            else:
                end_col_offset = node.end_col_offset
        for child in iter_child_nodes(node):
            _fix(child, lineno, col_offset, end_lineno, end_col_offset)
    _fix(node, 1, 0, 1, 0)
    return node


def increment_lineno(node, n=1):
    """
    Increment the line number and end line number of each node in the tree
    starting at *node* by *n*. This is useful to "move code" to a different
    location in a file.
    """
    for child in walk(node):
        if 'lineno' in child._attributes:
            child.lineno = (getattr(child, 'lineno', 0) or 0) + n
        if 'end_lineno' in child._attributes:
            child.end_lineno = (getattr(child, 'end_lineno', 0) or 0) + n
    return node

if _sys.version_info.major == 3 and _sys.version_info.minor >= 13:
    dump = _ast.dump
else:
    # Code import from Lib/ast.py
    #
    # minor changes: getattr(x, y, ...) is None => getattr(x, y, 42) is None
    #
    def dump(
        node, annotate_fields=True, include_attributes=False,
        # *,  # removed for compatibility with python2 :-/
        indent=None, show_empty=False,
    ):
        """
        Return a formatted dump of the tree in node.  This is mainly useful for
        debugging purposes.  If annotate_fields is true (by default),
        the returned string will show the names and the values for fields.
        If annotate_fields is false, the result string will be more compact by
        omitting unambiguous field names.  Attributes such as line
        numbers and column offsets are not dumped by default.  If this is wanted,
        include_attributes can be set to true.  If indent is a non-negative
        integer or string, then the tree will be pretty-printed with that indent
        level. None (the default) selects the single line representation.
        If show_empty is False, then empty lists and fields that are None
        will be omitted from the output for better readability.
        """
        def _format(node, level=0):
            if indent is not None:
                level += 1
                prefix = '\n' + indent * level
                sep = ',\n' + indent * level
            else:
                prefix = ''
                sep = ', '
            if isinstance(node, AST):
                cls = type(node)
                args = []
                args_buffer = []
                allsimple = True
                keywords = annotate_fields
                for name in node._fields:
                    try:
                        value = getattr(node, name)
                    except AttributeError:
                        keywords = True
                        continue
                    if value is None and getattr(cls, name, 42) is None:
                        keywords = True
                        continue
                    if (
                        not show_empty
                        and (value is None or value == [])
                        # Special cases:
                        # `Constant(value=None)` and `MatchSingleton(value=None)`
                        and not isinstance(node, (Constant, MatchSingleton))
                    ):
                        args_buffer.append(repr(value))
                        continue
                    elif not keywords:
                        args.extend(args_buffer)
                        args_buffer = []
                    value, simple = _format(value, level)
                    allsimple = allsimple and simple
                    if keywords:
                        args.append('%s=%s' % (name, value))
                    else:
                        args.append(value)
                if include_attributes and node._attributes:
                    for name in node._attributes:
                        try:
                            value = getattr(node, name)
                        except AttributeError:
                            continue
                        if value is None and getattr(cls, name, 42) is None:
                            continue
                        value, simple = _format(value, level)
                        allsimple = allsimple and simple
                        args.append('%s=%s' % (name, value))
                if allsimple and len(args) <= 3:
                    return '%s(%s)' % (node.__class__.__name__, ', '.join(args)), not args
                return '%s(%s%s)' % (node.__class__.__name__, prefix, sep.join(args)), False
            elif isinstance(node, list):
                if not node:
                    return '[]', True
                return '[%s%s]' % (prefix, sep.join(_format(x, level)[0] for x in node)), False
            return repr(node), True

        if not isinstance(node, AST):
            raise TypeError('expected AST, got %r' % node.__class__.__name__)
        if indent is not None and not isinstance(indent, str):
            indent = ' ' * indent
        return _format(node)[0]
