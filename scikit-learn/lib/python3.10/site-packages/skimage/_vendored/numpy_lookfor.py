# Vendored subset of numpy/lib/utils.py in 1.26.3
# https://github.com/numpy/numpy/blob/b4bf93b936802618ebb49ee43e382b576b29a0a6/numpy/lib/utils.py
#
# Can be removed after deprecation of `skimage.lookfor` is completed.

import sys
import os
import re

from numpy import ufunc


# Cache for lookfor: {id(module): {name: (docstring, kind, index), ...}...}
# where kind: "func", "class", "module", "object"
# and index: index in breadth-first namespace traversal
_lookfor_caches = {}


# regexp whose match indicates that the string may contain a function
# signature
_function_signature_re = re.compile(r"[a-z0-9_]+\(.*[,=].*\)", re.I)


def _getmembers(item):
    import inspect

    try:
        members = inspect.getmembers(item)
    except Exception:
        members = [(x, getattr(item, x)) for x in dir(item) if hasattr(item, x)]
    return members


def _lookfor_generate_cache(module, import_modules, regenerate):
    """
    Generate docstring cache for given module.

    Parameters
    ----------
    module : str, None, module
        Module for which to generate docstring cache
    import_modules : bool
        Whether to import sub-modules in packages.
    regenerate : bool
        Re-generate the docstring cache

    Returns
    -------
    cache : dict {obj_full_name: (docstring, kind, index), ...}
        Docstring cache for the module, either cached one (regenerate=False)
        or newly generated.

    """
    # Local import to speed up numpy's import time.
    import inspect

    from io import StringIO

    if module is None:
        module = "numpy"

    if isinstance(module, str):
        try:
            __import__(module)
        except ImportError:
            return {}
        module = sys.modules[module]
    elif isinstance(module, list) or isinstance(module, tuple):
        cache = {}
        for mod in module:
            cache.update(_lookfor_generate_cache(mod, import_modules, regenerate))
        return cache

    if id(module) in _lookfor_caches and not regenerate:
        return _lookfor_caches[id(module)]

    # walk items and collect docstrings
    cache = {}
    _lookfor_caches[id(module)] = cache
    seen = {}
    index = 0
    stack = [(module.__name__, module)]
    while stack:
        name, item = stack.pop(0)
        if id(item) in seen:
            continue
        seen[id(item)] = True

        index += 1
        kind = "object"

        if inspect.ismodule(item):
            kind = "module"
            try:
                _all = item.__all__
            except AttributeError:
                _all = None

            # import sub-packages
            if import_modules and hasattr(item, '__path__'):
                for pth in item.__path__:
                    for mod_path in os.listdir(pth):
                        this_py = os.path.join(pth, mod_path)
                        init_py = os.path.join(pth, mod_path, '__init__.py')
                        if os.path.isfile(this_py) and mod_path.endswith('.py'):
                            to_import = mod_path[:-3]
                        elif os.path.isfile(init_py):
                            to_import = mod_path
                        else:
                            continue
                        if to_import == '__init__':
                            continue

                        try:
                            old_stdout = sys.stdout
                            old_stderr = sys.stderr
                            try:
                                sys.stdout = StringIO()
                                sys.stderr = StringIO()
                                __import__(f"{name}.{to_import}")
                            finally:
                                sys.stdout = old_stdout
                                sys.stderr = old_stderr
                        except KeyboardInterrupt:
                            # Assume keyboard interrupt came from a user
                            raise
                        except BaseException:
                            # Ignore also SystemExit and pytests.importorskip
                            # `Skipped` (these are BaseExceptions; gh-22345)
                            continue

            for n, v in _getmembers(item):
                try:
                    item_name = getattr(
                        v,
                        '__name__',
                        f"{name}.{n}",
                    )
                    mod_name = getattr(v, '__module__', None)
                except NameError:
                    # ref. SWIG's global cvars
                    #    NameError: Unknown C global variable
                    item_name = f"{name}.{n}"
                    mod_name = None
                if '.' not in item_name and mod_name:
                    item_name = f"{mod_name}.{item_name}"

                if not item_name.startswith(name + '.'):
                    # don't crawl "foreign" objects
                    if isinstance(v, ufunc):
                        # ... unless they are ufuncs
                        pass
                    else:
                        continue
                elif not (inspect.ismodule(v) or _all is None or n in _all):
                    continue
                stack.append((f"{name}.{n}", v))
        elif inspect.isclass(item):
            kind = "class"
            for n, v in _getmembers(item):
                stack.append((f"{name}.{n}", v))
        elif hasattr(item, "__call__"):
            kind = "func"

        try:
            doc = inspect.getdoc(item)
        except NameError:
            # ref SWIG's NameError: Unknown C global variable
            doc = None
        if doc is not None:
            cache[name] = (doc, kind, index)

    return cache


def lookfor(what, module=None, import_modules=True, regenerate=False, output=None):
    """
    Do a keyword search on docstrings.

    A list of objects that matched the search is displayed,
    sorted by relevance. All given keywords need to be found in the
    docstring for it to be returned as a result, but the order does
    not matter.

    Parameters
    ----------
    what : str
        String containing words to look for.
    module : str or list, optional
        Name of module(s) whose docstrings to go through.
    import_modules : bool, optional
        Whether to import sub-modules in packages. Default is True.
    regenerate : bool, optional
        Whether to re-generate the docstring cache. Default is False.
    output : file-like, optional
        File-like object to write the output to. If omitted, use a pager.

    See Also
    --------
    source, info

    Notes
    -----
    Relevance is determined only roughly, by checking if the keywords occur
    in the function name, at the start of a docstring, etc.

    Examples
    --------
    >>> np.lookfor('binary representation') # doctest: +SKIP
    Search results for 'binary representation'
    ------------------------------------------
    numpy.binary_repr
        Return the binary representation of the input number as a string.
    numpy.core.setup_common.long_double_representation
        Given a binary dump as given by GNU od -b, look for long double
    numpy.base_repr
        Return a string representation of a number in the given base system.
    ...

    """
    import pydoc

    # Cache
    cache = _lookfor_generate_cache(module, import_modules, regenerate)

    # Search
    # XXX: maybe using a real stemming search engine would be better?
    found = []
    whats = str(what).lower().split()
    if not whats:
        return

    for name, (docstring, kind, index) in cache.items():
        if kind in ('module', 'object'):
            # don't show modules or objects
            continue
        doc = docstring.lower()
        if all(w in doc for w in whats):
            found.append(name)

    # Relevance sort
    # XXX: this is full Harrison-Stetson heuristics now,
    # XXX: it probably could be improved

    kind_relevance = {'func': 1000, 'class': 1000, 'module': -1000, 'object': -1000}

    def relevance(name, docstr, kind, index):
        r = 0
        # do the keywords occur within the start of the docstring?
        first_doc = "\n".join(docstr.lower().strip().split("\n")[:3])
        r += sum([200 for w in whats if w in first_doc])
        # do the keywords occur in the function name?
        r += sum([30 for w in whats if w in name])
        # is the full name long?
        r += -len(name) * 5
        # is the object of bad type?
        r += kind_relevance.get(kind, -1000)
        # is the object deep in namespace hierarchy?
        r += -name.count('.') * 10
        r += max(-index / 100, -100)
        return r

    def relevance_value(a):
        return relevance(a, *cache[a])

    found.sort(key=relevance_value)

    # Pretty-print
    s = f"Search results for '{' '.join(whats)}'"
    help_text = [s, "-" * len(s)]
    for name in found[::-1]:
        doc, kind, ix = cache[name]

        doclines = [line.strip() for line in doc.strip().split("\n") if line.strip()]

        # find a suitable short description
        try:
            first_doc = doclines[0].strip()
            if _function_signature_re.search(first_doc):
                first_doc = doclines[1].strip()
        except IndexError:
            first_doc = ""
        help_text.append(f"{name}\n    {first_doc}")

    if not found:
        help_text.append("Nothing found.")

    # Output
    if output is not None:
        output.write("\n".join(help_text))
    elif len(help_text) > 10:
        pager = pydoc.getpager()
        pager("\n".join(help_text))
    else:
        print("\n".join(help_text))
