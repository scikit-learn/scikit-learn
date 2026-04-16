"""
A simple XPath-like language for tree traversal.

This works by creating a filter chain of generator functions.  Each
function selects a part of the expression, e.g. a child node, a
specific descendant or a node that holds an attribute.
"""


import re
import operator

path_tokenizer = re.compile(
    r"("
    r"b?'[^']*'|b?\"[^\"]*\"|"
    r"//?|"
    r"\(\)|"
    r"==?|"
    r"[/.*\[\]()@])|"
    r"([^/\[\]()@=\s]+)|"
    r"\s+"
    ).findall

def iterchildren(node, attr_name):
    # returns an iterable of all child nodes of that name
    child = getattr(node, attr_name)
    if child is not None:
        if type(child) is list:
            return child
        else:
            return [child]
    else:
        return ()

def _get_first_or_none(it):
    try:
        return next(it)
    except StopIteration:
        return None

def type_name(node):
    return node.__class__.__name__.split('.')[-1]

def parse_func(next, token):
    name = token[1]
    token = next()
    if token[0] != '(':
        raise ValueError("Expected '(' after function name '%s'" % name)
    predicate = handle_predicate(next, token)
    return name, predicate

def handle_func_not(next, token):
    """
    not(...)
    """
    name, predicate = parse_func(next, token)

    def select(result):
        for node in result:
            if _get_first_or_none(predicate([node])) is None:
                yield node
    return select

def handle_name(next, token):
    """
    /NodeName/
    or
    func(...)
    """
    name = token[1]
    if name in functions:
        return functions[name](next, token)
    def select(result):
        for node in result:
            for attr_name in node.child_attrs:
                for child in iterchildren(node, attr_name):
                    if type_name(child) == name:
                        yield child
    return select

def handle_star(next, token):
    """
    /*/
    """
    def select(result):
        for node in result:
            for name in node.child_attrs:
                yield from iterchildren(node, name)
    return select

def handle_dot(next, token):
    """
    /./
    """
    def select(result):
        return result
    return select

def handle_descendants(next, token):
    """
    //...
    """
    token = next()
    if token[0] == "*":
        def iter_recursive(node):
            for name in node.child_attrs:
                for child in iterchildren(node, name):
                    yield child
                    yield from iter_recursive(child)
    elif not token[0]:
        node_name = token[1]
        def iter_recursive(node):
            for name in node.child_attrs:
                for child in iterchildren(node, name):
                    if type_name(child) == node_name:
                        yield child
                    yield from iter_recursive(child)
    else:
        raise ValueError("Expected node name after '//'")

    def select(result):
        for node in result:
            yield from iter_recursive(node)

    return select


def handle_attribute(next, token):
    token = next()
    if token[0]:
        raise ValueError("Expected attribute name")
    name = token[1]
    value = None
    token = next.peek()
    if token[0] == '=':
        next()
        value = parse_path_value(next)

    readattr = operator.attrgetter(name)
    if value is None:
        def select(result):
            for node in result:
                try:
                    attr_value = readattr(node)
                except AttributeError:
                    continue
                if attr_value is not None:
                    yield attr_value
    else:
        def select(result):
            for node in result:
                try:
                    attr_value = readattr(node)
                except AttributeError:
                    continue
                if attr_value == value:
                    yield attr_value
                elif (isinstance(attr_value, bytes) and isinstance(value, str) and
                        attr_value == value.encode()):
                    # allow a bytes-to-string comparison too
                    yield attr_value

    return select


def parse_path_value(next):
    token = next()
    value = token[0]
    if value:
        if value[:1] == "'" or value[:1] == '"':
            assert value[-1] == value[0]
            return value[1:-1]
        if value[:2] == "b'" or value[:2] == 'b"':
            assert value[-1] == value[1]
            return value[2:-1].encode('UTF-8')
        try:
            return int(value)
        except ValueError:
            pass
    elif token[1].isdigit():
        return int(token[1])
    else:
        name = token[1].lower()
        if name == 'true':
            return True
        elif name == 'false':
            return False
    raise ValueError(f"Invalid attribute predicate: '{value}'")


def handle_predicate(next, token):
    token = next()

    and_conditions = [[]]
    or_conditions = [and_conditions]

    while token[0] not in (']', ')'):
        and_conditions[-1].append( operations[token[0]](next, token) )
        try:
            token = next()
        except StopIteration:
            break
        else:
            if token[0] == "/":
                token = next()

        if not token[0]:
            if token[1] == 'and':
                and_conditions.append([])
                token = next()
            elif token[1] == 'or':
                and_conditions = [[]]
                or_conditions.append(and_conditions)
                token = next()

    if not and_conditions[-1]:
        raise ValueError("Incomplete predicate")

    def select(result):
        for node in result:
            node_base = (node,)
            for and_conditions in or_conditions:
                for condition in and_conditions:
                    subresult = iter(node_base)
                    for select in condition:
                        subresult = select(subresult)
                    predicate_result = _get_first_or_none(subresult)
                    if predicate_result is None:
                        # Fail current 'and' condition and skip to next 'or' condition.
                        break
                else:
                    # All 'and' conditions matched, report and skip following 'or' conditions.
                    yield node
                    break

    return select


operations = {
    "@":  handle_attribute,
    "":   handle_name,
    "*":  handle_star,
    ".":  handle_dot,
    "//": handle_descendants,
    "[":  handle_predicate,
    }

functions = {
    'not' : handle_func_not
    }


class _LookAheadTokenizer:
    def __init__(self, path):
        self._tokens = [
            (special, text)
            for (special, text) in path_tokenizer(path)
            if special or text
        ]
        self._tokens.reverse()  # allow efficient .pop()

    def peek(self, default=(None, None)):
        return self._tokens[-1] if self._tokens else default

    def __call__(self):
        try:
            return self._tokens.pop()
        except IndexError:
            raise StopIteration from None


def _build_path_iterator(path):
    # parse pattern
    _next = _LookAheadTokenizer(path)
    token = _next()
    selector = []
    while 1:
        try:
            selector.append(operations[token[0]](_next, token))
        except StopIteration:
            raise ValueError("invalid path")
        try:
            token = _next()
            if token[0] == "/":
                token = _next()
        except StopIteration:
            break
    return selector

# main module API

def iterfind(node, path):
    selector_chain = _build_path_iterator(path)
    result = iter((node,))
    for select in selector_chain:
        result = select(result)
    return result

def find_first(node, path):
    return _get_first_or_none(iterfind(node, path))

def find_all(node, path):
    return list(iterfind(node, path))
