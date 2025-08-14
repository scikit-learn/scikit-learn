import re

# When sphinx (including the napoleon extension) parses the parameters
# section of a docstring, it converts the information into field lists.
# Some items in the list are for the parameter type. When the type fields
# are processed, the text is split and some tokens are turned into
# pending_xref nodes. These nodes are responsible for creating links.
#
# numpydoc does not create field lists, so the type information is
# not placed into fields that can be processed to make links. Instead,
# when parsing the type information we identify tokens that are link
# worthy and wrap them around a :obj: role.

# Note: we never split on commas that are not followed by a space
# You risk creating bad rst markup if you do so.

QUALIFIED_NAME_RE = re.compile(
    # e.g int, numpy.array, ~numpy.array, .class_in_current_module
    r"^" r"[~\.]?" r"[a-zA-Z_]\w*" r"(?:\.[a-zA-Z_]\w*)*" r"$"
)

CONTAINER_SPLIT_RE = re.compile(
    # splits dict(str, int) into
    #    ['dict', '[', 'str', ', ', 'int', ']', '']
    r"(\s*[\[\]\(\)\{\}]\s*|,\s+)"
)

CONTAINER_SPLIT_REJECT_RE = re.compile(
    # Leads to bad markup e.g.
    # {int}qualified_name
    r"[\]\)\}]\w"
)

DOUBLE_QUOTE_SPLIT_RE = re.compile(
    # splits 'callable ``f(x0, *args)`` or ``f(x0, y0, *args)``' into
    #    ['callable ', '``f(x0, *args)``', ' or ', '``f(x0, y0, *args)``', '']
    r"(``.+?``)"
)

ROLE_SPLIT_RE = re.compile(
    # splits to preserve ReST roles
    r"(:\w+:`.+?(?<!\\)`)"
)

SINGLE_QUOTE_SPLIT_RE = re.compile(
    # splits to preserve quoted expressions roles
    r"(`.+?`)"
)

TEXT_SPLIT_RE = re.compile(
    # splits on ' or ', ' | ', ', ' and ' '
    r"(\s+or\s+|\s+\|\s+|,\s+|\s+)"
)

CONTAINER_CHARS = set("[](){}")

# Save people some time and add some common standard aliases
DEFAULT_LINKS = {
    # Python
    "None": ":data:`python:None`",
    "bool": ":ref:`bool <python:bltin-boolean-values>`",
    "boolean": ":ref:`bool <python:bltin-boolean-values>`",
    "True": ":data:`python:True`",
    "False": ":data:`python:False`",
    "list": ":class:`python:list`",
    "tuple": ":class:`python:tuple`",
    "str": ":class:`python:str`",
    "string": ":class:`python:str`",
    "dict": ":class:`python:dict`",
    "float": ":class:`python:float`",
    "int": ":class:`python:int`",
    "callable": ":func:`python:callable`",
    "iterable": ":term:`python:iterable`",
    "sequence": ":term:`python:sequence`",
    "contextmanager": ":func:`python:contextlib.contextmanager`",
    "namedtuple": ":func:`python:collections.namedtuple`",
    "generator": ":term:`python:generator`",
    # NumPy
    "array": "numpy.ndarray",
    "ndarray": "numpy.ndarray",
    "np.ndarray": "numpy.ndarray",
    "array-like": ":term:`numpy:array_like`",
    "array_like": ":term:`numpy:array_like`",
    "scalar": ":ref:`scalar <numpy:arrays.scalars>`",
    "RandomState": "numpy.random.RandomState",
    "np.random.RandomState": "numpy.random.RandomState",
    "np.inf": ":data:`numpy.inf`",
    "np.nan": ":data:`numpy.nan`",
    "numpy": ":mod:`numpy`",
}


def make_xref(param_type, xref_aliases, xref_ignore):
    """Parse and apply appropriate sphinx role(s) to `param_type`.

    The :obj: role is the default.

    Parameters
    ----------
    param_type : str
        text
    xref_aliases : dict
        Mapping used to resolve common abbreviations and aliases
        to fully qualified names that can be cross-referenced.
    xref_ignore : set or "all"
        A set containing words not to cross-reference. Instead of a set, the
        string 'all' can be given to ignore all unrecognized terms.
        Unrecognized terms include those that are not in `xref_aliases` and
        are not already wrapped in a reST role.

    Returns
    -------
    out : str
        Text with fully-qualified names and terms that may be wrapped in a
        ``:obj:`` role.
    """
    ignore_set = xref_ignore
    wrap_unknown = True
    if isinstance(xref_ignore, str):
        if xref_ignore.lower() == "all":
            wrap_unknown = False
            ignore_set = set()
        else:
            raise TypeError(f"xref_ignore must be a set or 'all', got {xref_ignore}")

    if param_type in xref_aliases:
        link, title = xref_aliases[param_type], param_type
        param_type = link
    else:
        link = title = param_type

    if QUALIFIED_NAME_RE.match(link) and link not in ignore_set:
        if link != title:
            return f":obj:`{title} <{link}>`"
        if wrap_unknown:
            return f":obj:`{link}`"
        return link

    def _split_and_apply_re(s, pattern):
        """
        Split string using the regex pattern,
        apply main function to the parts that do not match the pattern,
        combine the results
        """
        results = []
        tokens = pattern.split(s)
        n = len(tokens)
        if n > 1:
            for i, tok in enumerate(tokens):
                if pattern.match(tok):
                    results.append(tok)
                else:
                    res = make_xref(tok, xref_aliases, xref_ignore)
                    # Opening brackets immediately after a role is
                    # bad markup. Detect that and add backslash.
                    # :role:`type`( to :role:`type`\(
                    if res and res[-1] == "`" and i < n - 1:
                        next_char = tokens[i + 1][0]
                        if next_char in "([{":
                            res += "\\"
                    results.append(res)

            return "".join(results)
        return s

    # The cases are dealt with in an order the prevents
    # conflict.
    # Then the strategy is:
    #   - Identify a pattern we are not interested in
    #   - split off the pattern
    #   - re-apply the function to the other parts
    #   - join the results with the pattern

    # Unsplittable literal
    if "``" in param_type:
        return _split_and_apply_re(param_type, DOUBLE_QUOTE_SPLIT_RE)

    # Any roles
    if ":`" in param_type:
        return _split_and_apply_re(param_type, ROLE_SPLIT_RE)

    # Any quoted expressions
    if "`" in param_type:
        return _split_and_apply_re(param_type, SINGLE_QUOTE_SPLIT_RE)

    # Any sort of bracket '[](){}'
    if any(c in CONTAINER_CHARS for c in param_type):
        if CONTAINER_SPLIT_REJECT_RE.search(param_type):
            return param_type
        return _split_and_apply_re(param_type, CONTAINER_SPLIT_RE)

    # Common splitter tokens
    return _split_and_apply_re(param_type, TEXT_SPLIT_RE)
