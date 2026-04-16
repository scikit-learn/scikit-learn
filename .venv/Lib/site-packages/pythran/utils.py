""" Common function use for AST manipulation. """

import gast as ast
from pythran.tables import MODULES
from pythran.conversion import mangle, demangle
from functools import reduce
from contextlib import contextmanager


def isstr(node):
    return isinstance(node, ast.Constant) and isinstance(node.value, str)


def isintegral(node):
    return isinstance(node, ast.Constant) and isinstance(node.value, (int,
                                                                      bool))


def isnum(node):
    return isinstance(node, ast.Constant) and isinstance(node.value, (int,
                                                                      float,
                                                                      bool))


def isextslice(node):
    if not isinstance(node, ast.Tuple):
        return False
    return any(isinstance(elt, ast.Slice) for elt in node.elts)


def ispowi(node):
    if not isinstance(node.op, ast.Pow):
        return False
    attr = 'right' if isinstance(node, ast.BinOp) else 'value'
    if not isintegral(getattr(node, attr)):
        return False
    return getattr(node, attr).value >= 0


def attr_to_path(node):
    """ Compute path and final object for an attribute node """

    def get_intrinsic_path(modules, attr):
        """ Get function path and intrinsic from an ast.Attribute.  """
        if isinstance(attr, ast.Name):
            return modules[demangle(attr.id)], (demangle(attr.id),)
        elif isinstance(attr, ast.Attribute):
            module, path = get_intrinsic_path(modules, attr.value)
            return module[attr.attr], path + (attr.attr,)
    obj, path = get_intrinsic_path(MODULES, node)
    # hasattr check because `obj` can be a dict (for modules)
    if hasattr(obj, 'isliteral') and not obj.isliteral():
        path = path[:-1] + ('functor', path[-1])
    return obj, ('pythonic', ) + path


def path_to_attr(path):
    """
    Transform path to ast.Attribute.

    >>> import gast as ast
    >>> path = ('builtins', 'my', 'constant')
    >>> value = path_to_attr(path)
    >>> ref = ast.Attribute(
    ...     value=ast.Attribute(value=ast.Name(id="builtins",
    ...                                        ctx=ast.Load(),
    ...                                        annotation=None,
    ...                                        type_comment=None),
    ...                         attr="my", ctx=ast.Load()),
    ...     attr="constant", ctx=ast.Load())
    >>> ast.dump(ref) == ast.dump(value)
    True
    """
    return reduce(lambda hpath, last: ast.Attribute(hpath, last, ast.Load()),
                  path[1:], ast.Name(mangle(path[0]), ast.Load(), None, None))


def path_to_node(path):
    """
    Retrieve a symbol in MODULES based on its path
    >>> path = ('math', 'pi')
    >>> path_to_node(path) #doctest: +ELLIPSIS
    <pythran.intrinsic.ConstantIntr object at 0x...>
    """
    if len(path) == 1:
        return MODULES[path[0]]
    else:
        return path_to_node(path[:-1])[path[-1]]


def isattr(node):
    return (isinstance(node, ast.Call) and
            getattr(node.func, 'attr', None) == 'getattr')


def get_variable(assignable):
    """
    Return modified variable name.

    >>> import gast as ast
    >>> ref = ast.Subscript(
    ...     value=ast.Subscript(
    ...         value=ast.Name('a', ast.Load(), None, None),
    ...         slice=ast.Name('i', ast.Load(), None, None),
    ...         ctx=ast.Load()),
    ...     slice=ast.Name('j', ast.Load(), None, None),
    ...     ctx=ast.Load())
    >>> get_variable(ref).id
    'a'
    """
    msg = "Only name and subscript can be assigned."
    assert isinstance(assignable, (ast.Name, ast.Subscript)), msg
    while isinstance(assignable, ast.Subscript) or isattr(assignable):
        if isattr(assignable):
            assignable = assignable.args[0]
        else:
            assignable = assignable.value
    return assignable


def pythran_builtin(name):
    return MODULES['builtins']['pythran'][name]


def pythran_builtin_path(name):
    assert name in MODULES['builtins']['pythran']
    return ('builtins', 'pythran', name)


def pythran_builtin_attr(name):
    return path_to_attr(pythran_builtin_path(name))


def cxxid(name):
    from pythran.tables import cxx_keywords
    return name + '_' * (name in cxx_keywords)


def quote_cxxstring(s):
    subs = (('\\', '\\\\'),
            ('\n', '\\n'),
            ('\r', '\\r'),
            ('"', '\\"'),
           )
    quoted = s
    for f, t in subs:
        quoted = quoted.replace(f, t)
    return quoted


@contextmanager
def pushpop(l, v):
    l.append(v)
    yield
    l.pop()
