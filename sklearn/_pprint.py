"""This module contains the _EstimatorPrettyPrinter class used in
BaseEstimator.__repr__ for pretty-printing estimators"""

from inspect import signature
import pprint

from .base import BaseEstimator


# Dummy classes for correctly rendering key-value tuples from dicts or
# parameters
class KeyValTuple(tuple):
    def __repr__(self):
        # needed for _dispatch[tuple.__repr__] not to be overridden
        return super().__repr__()


class KeyValTupleParam(KeyValTuple):
    pass


def _changed_params(estimator):
    """Return dict (param_name: value) of parameters that were given to
    estimator with non-default values."""

    params = estimator.get_params(deep=False)
    filtered_params = {}
    init_func = getattr(estimator.__init__, 'deprecated_original',
                        estimator.__init__)
    init_params = signature(init_func).parameters
    init_params = {name: param.default for (name, param) in
                   init_params.items()}
    for k, v in params.items():
        if v != init_params[k]:
            filtered_params[k] = v
    return filtered_params


class _EstimatorPrettyPrinter(pprint.PrettyPrinter):
    """Pretty Printer class for estimator objects.

    This extends the pprint.PrettyPrinter class, because:
    - we need estimators to be printed with their parameters, e.g.
      Estimator(param1=value1, ...) which is not supported by default.
    - the 'compact' parameter of PrettyPrinter is ignored for dicts, which
      may lead to very long representations that we want to avoid.

    Quick overview of pprint.PrettyPrinter (see also
    https://stackoverflow.com/questions/49565047/pprint-with-hex-numbers):

    - the entry point is the _format() method which calls format() (overridden
      here)
    - format() directly calls _safe_repr() for a first try at rendering the
      object
    - _safe_repr formats the whole object reccursively, only calling itself,
      not caring about line length or anything
    - back to _format(), if the output string is too long, _format() then calls
      the appropriate _pprint_TYPE() method (e.g. _pprint_list()) depending on
      the type of the object. This where the line length and the compact
      parameters are taken into account.
    - those _pprint_TYPE() methods will internally use the format() method for
      rendering the nested objects of an object (e.g. the elements of a list)

    In the end, everything has to be implemented twice: in _safe_repr and in
    the custom _pprint_TYPE methods. Unfortunately PrettyPrinter is really not
    straightforward to extend (especially when we want a compact output), so
    the code is a bit convoluted.

    This class overrides:
    - format() to support the changed_only parameter
    - _safe_repr to support printing of estimators (for when they fit on a
      single line)
    - _format_dict_items so that dict are correctly 'compacted'

    When estimators cannot be printed on a single line, the builtin _format()
    will call _pprint_estimator() because it was registered to do so (see
    _dispatch[BaseEstimator.__repr__] = _pprint_estimator).

    both _format_dict_items() and _pprint_estimator() use the
    _format_params_or_dict_items() method that will format parameters and
    key-value pairs respecting the compact parameter. This method needs another
    subroutine _pprint_key_val_tuple() used when a parameter or a key-value
    pair is too long to fit on a single line. This subroutine is called in
    _format() and is registered as well in the _dispatch dict (just like
    _pprint_estimator). We had to create the two classes KeyValTuple and
    KeyValTupleParam for this.
    """

    def __init__(self, indent=1, width=80, depth=None, stream=None, *,
                 compact=False, indent_at_name=True, changed_only=False):
        super().__init__(indent, width, depth, stream, compact=compact)
        self._indent_at_name = indent_at_name
        if self._indent_at_name:
            self._indent_per_level = 1  # ignore indent param
        self._changed_only = changed_only

    def format(self, object, context, maxlevels, level):
        return _safe_repr(object, context, maxlevels, level,
                          changed_only=self._changed_only)

    def _pprint_estimator(self, object, stream, indent, allowance, context,
                          level):
        stream.write(object.__class__.__name__ + '(')
        if self._indent_at_name:
            indent += len(object.__class__.__name__)

        if self._changed_only:
            params = _changed_params(object)
        else:
            params = object.get_params(deep=False)

        self._format_params(params.items(), stream, indent, allowance + 1,
                            context, level)
        stream.write(')')

    def _format_dict_items(self, items, stream, indent, allowance, context,
                           level):
        return self._format_params_or_dict_items(
            items, stream, indent, allowance, context, level, is_dict=True)

    def _format_params(self, items, stream, indent, allowance, context, level):
        return self._format_params_or_dict_items(
            items, stream, indent, allowance, context, level, is_dict=False)

    def _format_params_or_dict_items(self, object, stream, indent, allowance,
                                     context, level, is_dict):
        """Format dict items or parameters respecting the compact=True
        parameter. For some reason, the builtin rendering of dict items doesn't
        respect compact=True and will use one line per key-value if all cannot
        fit in a single line.
        Dict items will be rendered as <'key': value> while params will be
        rendered as <key=value>. The implementation is mostly copy/pasting from
        the builtin _format_items().
        """
        write = stream.write
        indent += self._indent_per_level
        delimnl = ',\n' + ' ' * indent
        delim = ''
        width = max_width = self._width - indent + 1
        it = iter(object)
        try:
            next_ent = next(it)
        except StopIteration:
            return
        last = False
        while not last:
            ent = next_ent
            try:
                next_ent = next(it)
            except StopIteration:
                last = True
                max_width -= allowance
                width -= allowance
            if self._compact:
                k, v = ent
                krepr = self._repr(k, context, level)
                vrepr = self._repr(v, context, level)
                if not is_dict:
                    krepr = krepr.strip("'")
                middle = ': ' if is_dict else '='
                rep = krepr + middle + vrepr
                w = len(rep) + 2
                if width < w:
                    width = max_width
                    if delim:
                        delim = delimnl
                if width >= w:
                    width -= w
                    write(delim)
                    delim = ', '
                    write(rep)
                    continue
            write(delim)
            delim = delimnl
            class_ = KeyValTuple if is_dict else KeyValTupleParam
            self._format(class_(ent), stream, indent,
                         allowance if last else 1, context, level)

    def _pprint_key_val_tuple(self, object, stream, indent, allowance, context,
                              level):
        """Pretty printing for key-value tuples from dict or parameters."""
        k, v = object
        rep = self._repr(k, context, level)
        if isinstance(object, KeyValTupleParam):
            rep = rep.strip("'")
            middle = '='
        else:
            middle = ': '
        stream.write(rep)
        stream.write(middle)
        self._format(v, stream, indent + len(rep) + len(middle), allowance,
                     context, level)

    _dispatch = pprint.PrettyPrinter._dispatch
    _dispatch[BaseEstimator.__repr__] = _pprint_estimator
    _dispatch[KeyValTuple.__repr__] = _pprint_key_val_tuple


def _safe_repr(object, context, maxlevels, level, changed_only=False):
    """Same as the builtin _safe_repr, with added support for Estimator
    objects."""
    typ = type(object)

    if typ in pprint._builtin_scalars:
        return repr(object), True, False

    r = getattr(typ, "__repr__", None)
    if issubclass(typ, dict) and r is dict.__repr__:
        if not object:
            return "{}", True, False
        objid = id(object)
        if maxlevels and level >= maxlevels:
            return "{...}", False, objid in context
        if objid in context:
            return pprint._recursion(object), False, True
        context[objid] = 1
        readable = True
        recursive = False
        components = []
        append = components.append
        level += 1
        saferepr = _safe_repr
        items = sorted(object.items(), key=pprint._safe_tuple)
        for k, v in items:
            krepr, kreadable, krecur = saferepr(
                k, context, maxlevels, level, changed_only=changed_only)
            vrepr, vreadable, vrecur = saferepr(
                v, context, maxlevels, level, changed_only=changed_only)
            append("%s: %s" % (krepr, vrepr))
            readable = readable and kreadable and vreadable
            if krecur or vrecur:
                recursive = True
        del context[objid]
        return "{%s}" % ", ".join(components), readable, recursive

    if (issubclass(typ, list) and r is list.__repr__) or \
       (issubclass(typ, tuple) and r is tuple.__repr__):
        if issubclass(typ, list):
            if not object:
                return "[]", True, False
            format = "[%s]"
        elif len(object) == 1:
            format = "(%s,)"
        else:
            if not object:
                return "()", True, False
            format = "(%s)"
        objid = id(object)
        if maxlevels and level >= maxlevels:
            return format % "...", False, objid in context
        if objid in context:
            return pprint._recursion(object), False, True
        context[objid] = 1
        readable = True
        recursive = False
        components = []
        append = components.append
        level += 1
        for o in object:
            orepr, oreadable, orecur = _safe_repr(
                o, context, maxlevels, level, changed_only=changed_only)
            append(orepr)
            if not oreadable:
                readable = False
            if orecur:
                recursive = True
        del context[objid]
        return format % ", ".join(components), readable, recursive

    if issubclass(typ, BaseEstimator):
        objid = id(object)
        if maxlevels and level >= maxlevels:
            return "{...}", False, objid in context
        if objid in context:
            return pprint._recursion(object), False, True
        context[objid] = 1
        readable = True
        recursive = False
        if changed_only:
            params = _changed_params(object)
        else:
            params = object.get_params(deep=False)
        components = []
        append = components.append
        level += 1
        saferepr = _safe_repr
        items = sorted(params.items(), key=pprint._safe_tuple)
        for k, v in items:
            krepr, kreadable, krecur = saferepr(
                k, context, maxlevels, level, changed_only=changed_only)
            vrepr, vreadable, vrecur = saferepr(
                v, context, maxlevels, level, changed_only=changed_only)
            append("%s=%s" % (krepr.strip("'"), vrepr))
            readable = readable and kreadable and vreadable
            if krecur or vrecur:
                recursive = True
        del context[objid]
        return ("%s(%s)" % (typ.__name__, ", ".join(components)), readable,
                recursive)

    rep = repr(object)
    return rep, (rep and not rep.startswith('<')), False
