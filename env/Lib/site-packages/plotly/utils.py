import textwrap
from pprint import PrettyPrinter

from _plotly_utils.utils import *
from _plotly_utils.data_utils import *

# Pretty printing
def _list_repr_elided(v, threshold=200, edgeitems=3, indent=0, width=80):
    """
    Return a string representation for of a list where list is elided if
    it has more than n elements

    Parameters
    ----------
    v : list
        Input list
    threshold :
        Maximum number of elements to display

    Returns
    -------
    str
    """
    if isinstance(v, list):
        open_char, close_char = "[", "]"
    elif isinstance(v, tuple):
        open_char, close_char = "(", ")"
    else:
        raise ValueError("Invalid value of type: %s" % type(v))

    if len(v) <= threshold:
        disp_v = v
    else:
        disp_v = list(v[:edgeitems]) + ["..."] + list(v[-edgeitems:])

    v_str = open_char + ", ".join([str(e) for e in disp_v]) + close_char

    v_wrapped = "\n".join(
        textwrap.wrap(
            v_str,
            width=width,
            initial_indent=" " * (indent + 1),
            subsequent_indent=" " * (indent + 1),
        )
    ).strip()
    return v_wrapped


class ElidedWrapper(object):
    """
    Helper class that wraps values of certain types and produces a custom
    __repr__() that may be elided and is suitable for use during pretty
    printing
    """

    def __init__(self, v, threshold, indent):
        self.v = v
        self.indent = indent
        self.threshold = threshold

    @staticmethod
    def is_wrappable(v):
        numpy = get_module("numpy")
        if isinstance(v, (list, tuple)) and len(v) > 0 and not isinstance(v[0], dict):
            return True
        elif numpy and isinstance(v, numpy.ndarray):
            return True
        elif isinstance(v, str):
            return True
        else:
            return False

    def __repr__(self):
        numpy = get_module("numpy")
        if isinstance(self.v, (list, tuple)):
            # Handle lists/tuples
            res = _list_repr_elided(
                self.v, threshold=self.threshold, indent=self.indent
            )
            return res
        elif numpy and isinstance(self.v, numpy.ndarray):
            # Handle numpy arrays

            # Get original print opts
            orig_opts = numpy.get_printoptions()

            # Set threshold to self.max_list_elements
            numpy.set_printoptions(
                **dict(orig_opts, threshold=self.threshold, edgeitems=3, linewidth=80)
            )

            res = self.v.__repr__()

            # Add indent to all but the first line
            res_lines = res.split("\n")
            res = ("\n" + " " * self.indent).join(res_lines)

            # Restore print opts
            numpy.set_printoptions(**orig_opts)
            return res
        elif isinstance(self.v, str):
            # Handle strings
            if len(self.v) > 80:
                return "(" + repr(self.v[:30]) + " ... " + repr(self.v[-30:]) + ")"
            else:
                return self.v.__repr__()
        else:
            return self.v.__repr__()


class ElidedPrettyPrinter(PrettyPrinter):
    """
    PrettyPrinter subclass that elides long lists/arrays/strings
    """

    def __init__(self, *args, **kwargs):
        self.threshold = kwargs.pop("threshold", 200)
        PrettyPrinter.__init__(self, *args, **kwargs)

    def _format(self, val, stream, indent, allowance, context, level):
        if ElidedWrapper.is_wrappable(val):
            elided_val = ElidedWrapper(val, self.threshold, indent)

            return self._format(elided_val, stream, indent, allowance, context, level)
        else:
            return PrettyPrinter._format(
                self, val, stream, indent, allowance, context, level
            )


def node_generator(node, path=()):
    """
    General, node-yielding generator.

    Yields (node, path) tuples when it finds values that are dict
    instances.

    A path is a sequence of hashable values that can be used as either keys to
    a mapping (dict) or indices to a sequence (list). A path is always wrt to
    some object. Given an object, a path explains how to get from the top level
    of that object to a nested value in the object.

    :param (dict) node: Part of a dict to be traversed.
    :param (tuple[str]) path: Defines the path of the current node.
    :return: (Generator)

    Example:

        >>> for node, path in node_generator({'a': {'b': 5}}):
        ...     print(node, path)
        {'a': {'b': 5}} ()
        {'b': 5} ('a',)

    """
    if not isinstance(node, dict):
        return  # in case it's called with a non-dict node at top level
    yield node, path
    for key, val in node.items():
        if isinstance(val, dict):
            for item in node_generator(val, path + (key,)):
                yield item


def get_by_path(obj, path):
    """
    Iteratively get on obj for each key in path.

    :param (list|dict) obj: The top-level object.
    :param (tuple[str]|tuple[int]) path: Keys to access parts of obj.

    :return: (*)

    Example:

        >>> figure = {'data': [{'x': [5]}]}
        >>> path = ('data', 0, 'x')
        >>> get_by_path(figure, path)
        [5]
    """
    for key in path:
        obj = obj[key]
    return obj


def decode_unicode(coll):
    if isinstance(coll, list):
        for no, entry in enumerate(coll):
            if isinstance(entry, (dict, list)):
                coll[no] = decode_unicode(entry)
            else:
                if isinstance(entry, str):
                    try:
                        coll[no] = str(entry)
                    except UnicodeEncodeError:
                        pass
    elif isinstance(coll, dict):
        keys, vals = list(coll.keys()), list(coll.values())
        for key, val in zip(keys, vals):
            if isinstance(val, (dict, list)):
                coll[key] = decode_unicode(val)
            elif isinstance(val, str):
                try:
                    coll[key] = str(val)
                except UnicodeEncodeError:
                    pass
            coll[str(key)] = coll.pop(key)
    return coll
