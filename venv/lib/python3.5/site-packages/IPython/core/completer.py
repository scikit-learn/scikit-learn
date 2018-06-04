"""Completion for IPython.

This module started as fork of the rlcompleter module in the Python standard
library.  The original enhancements made to rlcompleter have been sent
upstream and were accepted as of Python 2.3,

This module now support a wide variety of completion mechanism both available
for normal classic Python code, as well as completer for IPython specific
Syntax like magics.

Latex and Unicode completion
============================

IPython and compatible frontends not only can complete your code, but can help
you to input a wide range of characters. In particular we allow you to insert
a unicode character using the tab completion mechanism.

Forward latex/unicode completion
--------------------------------

Forward completion allows you to easily type a unicode character using its latex
name, or unicode long description. To do so type a backslash follow by the
relevant name and press tab:


Using latex completion:

.. code::

    \\alpha<tab>
    α

or using unicode completion:


.. code::

    \\greek small letter alpha<tab>
    α


Only valid Python identifiers will complete. Combining characters (like arrow or
dots) are also available, unlike latex they need to be put after the their
counterpart that is to say, `F\\\\vec<tab>` is correct, not `\\\\vec<tab>F`.

Some browsers are known to display combining characters incorrectly.

Backward latex completion
-------------------------

It is sometime challenging to know how to type a character, if you are using
IPython, or any compatible frontend you can prepend backslash to the character
and press `<tab>` to expand it to its latex form.

.. code::

    \\α<tab>
    \\alpha


Both forward and backward completions can be deactivated by setting the
``Completer.backslash_combining_completions`` option to ``False``.


Experimental
============

Starting with IPython 6.0, this module can make use of the Jedi library to
generate completions both using static analysis of the code, and dynamically
inspecting multiple namespaces. The APIs attached to this new mechanism is
unstable and will raise unless use in an :any:`provisionalcompleter` context
manager.

You will find that the following are experimental:

    - :any:`provisionalcompleter`
    - :any:`IPCompleter.completions`
    - :any:`Completion`
    - :any:`rectify_completions`

.. note::

    better name for :any:`rectify_completions` ?

We welcome any feedback on these new API, and we also encourage you to try this
module in debug mode (start IPython with ``--Completer.debug=True``) in order
to have extra logging information is :any:`jedi` is crashing, or if current
IPython completer pending deprecations are returning results not yet handled
by :any:`jedi`

Using Jedi for tab completion allow snippets like the following to work without
having to execute any code:

   >>> myvar = ['hello', 42]
   ... myvar[1].bi<tab>

Tab completion will be able to infer that ``myvar[1]`` is a real number without
executing any code unlike the previously available ``IPCompleter.greedy``
option.

Be sure to update :any:`jedi` to the latest stable version or to try the
current development version to get better completions.
"""


# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.
#
# Some of this code originated from rlcompleter in the Python standard library
# Copyright (C) 2001 Python Software Foundation, www.python.org


import __main__
import builtins as builtin_mod
import glob
import time
import inspect
import itertools
import keyword
import os
import re
import sys
import unicodedata
import string
import warnings

from contextlib import contextmanager
from importlib import import_module
from typing import Iterator, List, Tuple, Iterable, Union
from types import SimpleNamespace

from traitlets.config.configurable import Configurable
from IPython.core.error import TryNext
from IPython.core.inputsplitter import ESC_MAGIC
from IPython.core.latex_symbols import latex_symbols, reverse_latex_symbol
from IPython.core.oinspect import InspectColors
from IPython.utils import generics
from IPython.utils.dir2 import dir2, get_real_method
from IPython.utils.process import arg_split
from traitlets import Bool, Enum, observe, Int

# skip module docstests
skip_doctest = True

try:
    import jedi
    jedi.settings.case_insensitive_completion = False
    import jedi.api.helpers
    import jedi.api.classes
    JEDI_INSTALLED = True
except ImportError:
    JEDI_INSTALLED = False
#-----------------------------------------------------------------------------
# Globals
#-----------------------------------------------------------------------------

# Public API
__all__ = ['Completer','IPCompleter']

if sys.platform == 'win32':
    PROTECTABLES = ' '
else:
    PROTECTABLES = ' ()[]{}?=\\|;:\'#*"^&'

# Protect against returning an enormous number of completions which the frontend
# may have trouble processing.
MATCHES_LIMIT = 500

_deprecation_readline_sentinel = object()


class ProvisionalCompleterWarning(FutureWarning):
    """
    Exception raise by an experimental feature in this module.

    Wrap code in :any:`provisionalcompleter` context manager if you
    are certain you want to use an unstable feature.
    """
    pass

warnings.filterwarnings('error', category=ProvisionalCompleterWarning)

@contextmanager
def provisionalcompleter(action='ignore'):
    """


    This contest manager has to be used in any place where unstable completer
    behavior and API may be called.

    >>> with provisionalcompleter():
    ...     completer.do_experimetal_things() # works

    >>> completer.do_experimental_things() # raises.

    .. note:: Unstable

        By using this context manager you agree that the API in use may change
        without warning, and that you won't complain if they do so.

        You also understand that if the API is not to you liking you should report
        a bug to explain your use case upstream and improve the API and will loose
        credibility if you complain after the API is make stable.

        We'll be happy to get your feedback , feature request and improvement on
        any of the unstable APIs !
    """
    with warnings.catch_warnings():
        warnings.filterwarnings(action, category=ProvisionalCompleterWarning)
        yield


def has_open_quotes(s):
    """Return whether a string has open quotes.

    This simply counts whether the number of quote characters of either type in
    the string is odd.

    Returns
    -------
    If there is an open quote, the quote character is returned.  Else, return
    False.
    """
    # We check " first, then ', so complex cases with nested quotes will get
    # the " to take precedence.
    if s.count('"') % 2:
        return '"'
    elif s.count("'") % 2:
        return "'"
    else:
        return False


def protect_filename(s, protectables=PROTECTABLES):
    """Escape a string to protect certain characters."""
    if set(s) & set(protectables):
        if sys.platform == "win32":
            return '"' + s + '"'
        else:
            return "".join(("\\" + c if c in protectables else c) for c in s)
    else:
        return s


def expand_user(path:str) -> Tuple[str, bool, str]:
    """Expand ``~``-style usernames in strings.

    This is similar to :func:`os.path.expanduser`, but it computes and returns
    extra information that will be useful if the input was being used in
    computing completions, and you wish to return the completions with the
    original '~' instead of its expanded value.

    Parameters
    ----------
    path : str
      String to be expanded.  If no ~ is present, the output is the same as the
      input.

    Returns
    -------
    newpath : str
      Result of ~ expansion in the input path.
    tilde_expand : bool
      Whether any expansion was performed or not.
    tilde_val : str
      The value that ~ was replaced with.
    """
    # Default values
    tilde_expand = False
    tilde_val = ''
    newpath = path

    if path.startswith('~'):
        tilde_expand = True
        rest = len(path)-1
        newpath = os.path.expanduser(path)
        if rest:
            tilde_val = newpath[:-rest]
        else:
            tilde_val = newpath

    return newpath, tilde_expand, tilde_val


def compress_user(path:str, tilde_expand:bool, tilde_val:str) -> str:
    """Does the opposite of expand_user, with its outputs.
    """
    if tilde_expand:
        return path.replace(tilde_val, '~')
    else:
        return path


def completions_sorting_key(word):
    """key for sorting completions

    This does several things:

    - Demote any completions starting with underscores to the end
    - Insert any %magic and %%cellmagic completions in the alphabetical order
      by their name
    """
    prio1, prio2 = 0, 0

    if word.startswith('__'):
        prio1 = 2
    elif word.startswith('_'):
        prio1 = 1

    if word.endswith('='):
        prio1 = -1

    if word.startswith('%%'):
        # If there's another % in there, this is something else, so leave it alone
        if not "%" in word[2:]:
            word = word[2:]
            prio2 = 2
    elif word.startswith('%'):
        if not "%" in word[1:]:
            word = word[1:]
            prio2 = 1

    return prio1, word, prio2


class _FakeJediCompletion:
    """
    This is a workaround to communicate to the UI that Jedi has crashed and to
    report a bug. Will be used only id :any:`IPCompleter.debug` is set to true.

    Added in IPython 6.0 so should likely be removed for 7.0

    """

    def __init__(self, name):

        self.name = name
        self.complete = name
        self.type = 'crashed'
        self.name_with_symbols = name
        self.signature = ''
        self._origin = 'fake'

    def __repr__(self):
        return '<Fake completion object jedi has crashed>'


class Completion:
    """
    Completion object used and return by IPython completers.

    .. warning:: Unstable
        
        This function is unstable, API may change without warning.
        It will also raise unless use in proper context manager.

    This act as a middle ground :any:`Completion` object between the
    :any:`jedi.api.classes.Completion` object and the Prompt Toolkit completion
    object. While Jedi need a lot of information about evaluator and how the
    code should be ran/inspected, PromptToolkit (and other frontend) mostly
    need user facing information.

    - Which range should be replaced replaced by what.
    - Some metadata (like completion type), or meta information to displayed to
      the use user.

    For debugging purpose we can also store the origin of the completion (``jedi``,
    ``IPython.python_matches``, ``IPython.magics_matches``...).
    """

    __slots__ = ['start', 'end', 'text', 'type', 'signature', '_origin']

    def __init__(self, start: int, end: int, text: str, *, type: str=None, _origin='', signature='') -> None:
        warnings.warn("``Completion`` is a provisional API (as of IPython 6.0). "
                      "It may change without warnings. "
                      "Use in corresponding context manager.",
                      category=ProvisionalCompleterWarning, stacklevel=2)

        self.start = start
        self.end = end
        self.text = text
        self.type = type
        self.signature = signature
        self._origin = _origin

    def __repr__(self):
        return '<Completion start=%s end=%s text=%r type=%r, signature=%r,>' % \
                (self.start, self.end, self.text, self.type or '?', self.signature or '?')

    def __eq__(self, other)->Bool:
        """
        Equality and hash do not hash the type (as some completer may not be
        able to infer the type), but are use to (partially) de-duplicate
        completion.

        Completely de-duplicating completion is a bit tricker that just
        comparing as it depends on surrounding text, which Completions are not
        aware of.
        """
        return self.start == other.start and \
            self.end == other.end and \
            self.text == other.text

    def __hash__(self):
        return hash((self.start, self.end, self.text))


_IC = Iterable[Completion]


def _deduplicate_completions(text: str, completions: _IC)-> _IC:
    """
    Deduplicate a set of completions.

    .. warning:: Unstable

        This function is unstable, API may change without warning.

    Parameters
    ----------
    text: str
        text that should be completed.
    completions: Iterator[Completion]
        iterator over the completions to deduplicate

    Yields
    ------
    `Completions` objects


    Completions coming from multiple sources, may be different but end up having
    the same effect when applied to ``text``. If this is the case, this will
    consider completions as equal and only emit the first encountered.

    Not folded in `completions()` yet for debugging purpose, and to detect when
    the IPython completer does return things that Jedi does not, but should be
    at some point.
    """
    completions = list(completions)
    if not completions:
        return

    new_start = min(c.start for c in completions)
    new_end = max(c.end for c in completions)

    seen = set()
    for c in completions:
        new_text = text[new_start:c.start] + c.text + text[c.end:new_end]
        if new_text not in seen:
            yield c
            seen.add(new_text)


def rectify_completions(text: str, completions: _IC, *, _debug=False)->_IC:
    """
    Rectify a set of completions to all have the same ``start`` and ``end``

    .. warning:: Unstable

        This function is unstable, API may change without warning.
        It will also raise unless use in proper context manager.

    Parameters
    ----------
    text: str
        text that should be completed.
    completions: Iterator[Completion]
        iterator over the completions to rectify


    :any:`jedi.api.classes.Completion` s returned by Jedi may not have the same start and end, though
    the Jupyter Protocol requires them to behave like so. This will readjust
    the completion to have the same ``start`` and ``end`` by padding both
    extremities with surrounding text.

    During stabilisation should support a ``_debug`` option to log which
    completion are return by the IPython completer and not found in Jedi in
    order to make upstream bug report.
    """
    warnings.warn("`rectify_completions` is a provisional API (as of IPython 6.0). "
                 "It may change without warnings. "
                 "Use in corresponding context manager.",
                  category=ProvisionalCompleterWarning, stacklevel=2)

    completions = list(completions)
    if not completions:
        return
    starts = (c.start for c in completions)
    ends = (c.end for c in completions)

    new_start = min(starts)
    new_end = max(ends)

    seen_jedi = set()
    seen_python_matches = set()
    for c in completions:
        new_text = text[new_start:c.start] + c.text + text[c.end:new_end]
        if c._origin == 'jedi':
            seen_jedi.add(new_text)
        elif c._origin == 'IPCompleter.python_matches':
            seen_python_matches.add(new_text)
        yield Completion(new_start, new_end, new_text, type=c.type, _origin=c._origin, signature=c.signature)
    diff = seen_python_matches.difference(seen_jedi)
    if diff and _debug:
        print('IPython.python matches have extras:', diff)


if sys.platform == 'win32':
    DELIMS = ' \t\n`!@#$^&*()=+[{]}|;\'",<>?'
else:
    DELIMS = ' \t\n`!@#$^&*()=+[{]}\\|;:\'",<>?'

GREEDY_DELIMS = ' =\r\n'


class CompletionSplitter(object):
    """An object to split an input line in a manner similar to readline.

    By having our own implementation, we can expose readline-like completion in
    a uniform manner to all frontends.  This object only needs to be given the
    line of text to be split and the cursor position on said line, and it
    returns the 'word' to be completed on at the cursor after splitting the
    entire line.

    What characters are used as splitting delimiters can be controlled by
    setting the ``delims`` attribute (this is a property that internally
    automatically builds the necessary regular expression)"""

    # Private interface

    # A string of delimiter characters.  The default value makes sense for
    # IPython's most typical usage patterns.
    _delims = DELIMS

    # The expression (a normal string) to be compiled into a regular expression
    # for actual splitting.  We store it as an attribute mostly for ease of
    # debugging, since this type of code can be so tricky to debug.
    _delim_expr = None

    # The regular expression that does the actual splitting
    _delim_re = None

    def __init__(self, delims=None):
        delims = CompletionSplitter._delims if delims is None else delims
        self.delims = delims

    @property
    def delims(self):
        """Return the string of delimiter characters."""
        return self._delims

    @delims.setter
    def delims(self, delims):
        """Set the delimiters for line splitting."""
        expr = '[' + ''.join('\\'+ c for c in delims) + ']'
        self._delim_re = re.compile(expr)
        self._delims = delims
        self._delim_expr = expr

    def split_line(self, line, cursor_pos=None):
        """Split a line of text with a cursor at the given position.
        """
        l = line if cursor_pos is None else line[:cursor_pos]
        return self._delim_re.split(l)[-1]



class Completer(Configurable):

    greedy = Bool(False,
        help="""Activate greedy completion
        PENDING DEPRECTION. this is now mostly taken care of with Jedi.

        This will enable completion on elements of lists, results of function calls, etc.,
        but can be unsafe because the code is actually evaluated on TAB.
        """
    ).tag(config=True)

    use_jedi = Bool(default_value=False,
                    help="Experimental: Use Jedi to generate autocompletions. "
                    "Off by default.").tag(config=True)

    jedi_compute_type_timeout = Int(default_value=400,
        help="""Experimental: restrict time (in milliseconds) during which Jedi can compute types.
        Set to 0 to stop computing types. Non-zero value lower than 100ms may hurt
        performance by preventing jedi to build its cache.
        """).tag(config=True)

    debug = Bool(default_value=False,
                 help='Enable debug for the Completer. Mostly print extra '
                      'information for experimental jedi integration.')\
                      .tag(config=True)

    backslash_combining_completions = Bool(True, 
        help="Enable unicode completions, e.g. \\alpha<tab> . "
             "Includes completion of latex commands, unicode names, and expanding "
             "unicode characters back to latex commands.").tag(config=True)



    def __init__(self, namespace=None, global_namespace=None, **kwargs):
        """Create a new completer for the command line.

        Completer(namespace=ns, global_namespace=ns2) -> completer instance.

        If unspecified, the default namespace where completions are performed
        is __main__ (technically, __main__.__dict__). Namespaces should be
        given as dictionaries.

        An optional second namespace can be given.  This allows the completer
        to handle cases where both the local and global scopes need to be
        distinguished.
        """

        # Don't bind to namespace quite yet, but flag whether the user wants a
        # specific namespace or to use __main__.__dict__. This will allow us
        # to bind to __main__.__dict__ at completion time, not now.
        if namespace is None:
            self.use_main_ns = True
        else:
            self.use_main_ns = False
            self.namespace = namespace

        # The global namespace, if given, can be bound directly
        if global_namespace is None:
            self.global_namespace = {}
        else:
            self.global_namespace = global_namespace

        super(Completer, self).__init__(**kwargs)

    def complete(self, text, state):
        """Return the next possible completion for 'text'.

        This is called successively with state == 0, 1, 2, ... until it
        returns None.  The completion should begin with 'text'.

        """
        if self.use_main_ns:
            self.namespace = __main__.__dict__

        if state == 0:
            if "." in text:
                self.matches = self.attr_matches(text)
            else:
                self.matches = self.global_matches(text)
        try:
            return self.matches[state]
        except IndexError:
            return None

    def global_matches(self, text):
        """Compute matches when text is a simple name.

        Return a list of all keywords, built-in functions and names currently
        defined in self.namespace or self.global_namespace that match.

        """
        matches = []
        match_append = matches.append
        n = len(text)
        for lst in [keyword.kwlist,
                    builtin_mod.__dict__.keys(),
                    self.namespace.keys(),
                    self.global_namespace.keys()]:
            for word in lst:
                if word[:n] == text and word != "__builtins__":
                    match_append(word)

        snake_case_re = re.compile(r"[^_]+(_[^_]+)+?\Z")
        for lst in [self.namespace.keys(),
                    self.global_namespace.keys()]:
            shortened = {"_".join([sub[0] for sub in word.split('_')]) : word
                         for word in lst if snake_case_re.match(word)}
            for word in shortened.keys():
                if word[:n] == text and word != "__builtins__":
                    match_append(shortened[word])
        return matches

    def attr_matches(self, text):
        """Compute matches when text contains a dot.

        Assuming the text is of the form NAME.NAME....[NAME], and is
        evaluatable in self.namespace or self.global_namespace, it will be
        evaluated and its attributes (as revealed by dir()) are used as
        possible completions.  (For class instances, class members are
        also considered.)

        WARNING: this can still invoke arbitrary C code, if an object
        with a __getattr__ hook is evaluated.

        """

        # Another option, seems to work great. Catches things like ''.<tab>
        m = re.match(r"(\S+(\.\w+)*)\.(\w*)$", text)
    
        if m:
            expr, attr = m.group(1, 3)
        elif self.greedy:
            m2 = re.match(r"(.+)\.(\w*)$", self.line_buffer)
            if not m2:
                return []
            expr, attr = m2.group(1,2)
        else:
            return []
    
        try:
            obj = eval(expr, self.namespace)
        except:
            try:
                obj = eval(expr, self.global_namespace)
            except:
                return []

        if self.limit_to__all__ and hasattr(obj, '__all__'):
            words = get__all__entries(obj)
        else:
            words = dir2(obj)

        try:
            words = generics.complete_object(obj, words)
        except TryNext:
            pass
        except AssertionError:
            raise
        except Exception:
            # Silence errors from completion function
            #raise # dbg
            pass
        # Build match list to return
        n = len(attr)
        return [u"%s.%s" % (expr, w) for w in words if w[:n] == attr ]


def get__all__entries(obj):
    """returns the strings in the __all__ attribute"""
    try:
        words = getattr(obj, '__all__')
    except:
        return []
 
    return [w for w in words if isinstance(w, str)]


def match_dict_keys(keys: List[str], prefix: str, delims: str):
    """Used by dict_key_matches, matching the prefix to a list of keys

    Parameters
    ==========
    keys:
        list of keys in dictionary currently being completed.
    prefix:
        Part of the text already typed by the user. e.g. `mydict[b'fo`
    delims:
        String of delimiters to consider when finding the current key.

    Returns
    =======

    A tuple of three elements: ``quote``, ``token_start``, ``matched``, with
    ``quote`` being the quote that need to be used to close current string.
    ``token_start`` the position where the replacement should start occurring,
    ``matches`` a list of replacement/completion

    """
    if not prefix:
        return None, 0, [repr(k) for k in keys
                      if isinstance(k, (str, bytes))]
    quote_match = re.search('["\']', prefix)
    quote = quote_match.group()
    try:
        prefix_str = eval(prefix + quote, {})
    except Exception:
        return None, 0, []

    pattern = '[^' + ''.join('\\' + c for c in delims) + ']*$'
    token_match = re.search(pattern, prefix, re.UNICODE)
    token_start = token_match.start()
    token_prefix = token_match.group()

    matched = []
    for key in keys:
        try:
            if not key.startswith(prefix_str):
                continue
        except (AttributeError, TypeError, UnicodeError):
            # Python 3+ TypeError on b'a'.startswith('a') or vice-versa
            continue

        # reformat remainder of key to begin with prefix
        rem = key[len(prefix_str):]
        # force repr wrapped in '
        rem_repr = repr(rem + '"') if isinstance(rem, str) else repr(rem + b'"')
        if rem_repr.startswith('u') and prefix[0] not in 'uU':
            # Found key is unicode, but prefix is Py2 string.
            # Therefore attempt to interpret key as string.
            try:
                rem_repr = repr(rem.encode('ascii') + '"')
            except UnicodeEncodeError:
                continue

        rem_repr = rem_repr[1 + rem_repr.index("'"):-2]
        if quote == '"':
            # The entered prefix is quoted with ",
            # but the match is quoted with '.
            # A contained " hence needs escaping for comparison:
            rem_repr = rem_repr.replace('"', '\\"')

        # then reinsert prefix from start of token
        matched.append('%s%s' % (token_prefix, rem_repr))
    return quote, token_start, matched


def cursor_to_position(text:str, line:int, column:int)->int:
    """

    Convert the (line,column) position of the cursor in text to an offset in a
    string.

    Parameters
    ----------

    text : str
        The text in which to calculate the cursor offset
    line : int
        Line of the cursor; 0-indexed
    column : int
        Column of the cursor 0-indexed

    Return
    ------
        Position of the cursor in ``text``, 0-indexed.

    See Also
    --------
    position_to_cursor: reciprocal of this function

    """
    lines = text.split('\n')
    assert line <= len(lines), '{} <= {}'.format(str(line), str(len(lines)))

    return sum(len(l) + 1 for l in lines[:line]) + column

def position_to_cursor(text:str, offset:int)->Tuple[int, int]:
    """
    Convert the position of the cursor in text (0 indexed) to a line
    number(0-indexed) and a column number (0-indexed) pair

    Position should be a valid position in ``text``.

    Parameters
    ----------

    text : str
        The text in which to calculate the cursor offset
    offset : int
        Position of the cursor in ``text``, 0-indexed.

    Return
    ------
    (line, column) : (int, int)
        Line of the cursor; 0-indexed, column of the cursor 0-indexed


    See Also
    --------
    cursor_to_position : reciprocal of this function


    """

    assert 0 <= offset <= len(text) , "0 <= %s <= %s" % (offset , len(text))

    before = text[:offset]
    blines = before.split('\n')  # ! splitnes trim trailing \n
    line = before.count('\n')
    col = len(blines[-1])
    return line, col


def _safe_isinstance(obj, module, class_name):
    """Checks if obj is an instance of module.class_name if loaded
    """
    return (module in sys.modules and
            isinstance(obj, getattr(import_module(module), class_name)))


def back_unicode_name_matches(text):
    u"""Match unicode characters back to unicode name
    
    This does  ``☃`` -> ``\\snowman``

    Note that snowman is not a valid python3 combining character but will be expanded.
    Though it will not recombine back to the snowman character by the completion machinery.

    This will not either back-complete standard sequences like \\n, \\b ...
    
    Used on Python 3 only.
    """
    if len(text)<2:
        return u'', ()
    maybe_slash = text[-2]
    if maybe_slash != '\\':
        return u'', ()

    char = text[-1]
    # no expand on quote for completion in strings.
    # nor backcomplete standard ascii keys
    if char in string.ascii_letters or char in ['"',"'"]:
        return u'', ()
    try :
        unic = unicodedata.name(char)
        return '\\'+char,['\\'+unic]
    except KeyError:
        pass
    return u'', ()

def back_latex_name_matches(text:str):
    """Match latex characters back to unicode name
    
    This does ``\\ℵ`` -> ``\\aleph``

    Used on Python 3 only.
    """
    if len(text)<2:
        return u'', ()
    maybe_slash = text[-2]
    if maybe_slash != '\\':
        return u'', ()


    char = text[-1]
    # no expand on quote for completion in strings.
    # nor backcomplete standard ascii keys
    if char in string.ascii_letters or char in ['"',"'"]:
        return u'', ()
    try :
        latex = reverse_latex_symbol[char]
        # '\\' replace the \ as well
        return '\\'+char,[latex]
    except KeyError:
        pass
    return u'', ()


def _formatparamchildren(parameter) -> str:
    """
    Get parameter name and value from Jedi Private API

    Jedi does not expose a simple way to get `param=value` from its API.

    Parameter
    =========

    parameter:
        Jedi's function `Param`

    Returns
    =======

    A string like 'a', 'b=1', '*args', '**kwargs'


    """
    description = parameter.description
    if not description.startswith('param '):
        raise ValueError('Jedi function parameter description have change format.'
                         'Expected "param ...", found %r".' % description)
    return description[6:]

def _make_signature(completion)-> str:
    """
    Make the signature from a jedi completion

    Parameter
    =========

    completion: jedi.Completion
        object does not complete a function type

    Returns
    =======

    a string consisting of the function signature, with the parenthesis but
    without the function name. example:
    `(a, *args, b=1, **kwargs)`

    """

    return '(%s)'% ', '.join([f for f in (_formatparamchildren(p) for p in completion.params) if f])

class IPCompleter(Completer):
    """Extension of the completer class with IPython-specific features"""
    
    @observe('greedy')
    def _greedy_changed(self, change):
        """update the splitter and readline delims when greedy is changed"""
        if change['new']:
            self.splitter.delims = GREEDY_DELIMS
        else:
            self.splitter.delims = DELIMS
    
    merge_completions = Bool(True,
        help="""Whether to merge completion results into a single list
        
        If False, only the completion results from the first non-empty
        completer will be returned.
        """
    ).tag(config=True)
    omit__names = Enum((0,1,2), default_value=2,
        help="""Instruct the completer to omit private method names
        
        Specifically, when completing on ``object.<tab>``.
        
        When 2 [default]: all names that start with '_' will be excluded.
        
        When 1: all 'magic' names (``__foo__``) will be excluded.
        
        When 0: nothing will be excluded.
        """
    ).tag(config=True)
    limit_to__all__ = Bool(False,
        help="""
        DEPRECATED as of version 5.0.
        
        Instruct the completer to use __all__ for the completion
        
        Specifically, when completing on ``object.<tab>``.
        
        When True: only those names in obj.__all__ will be included.
        
        When False [default]: the __all__ attribute is ignored
        """,
    ).tag(config=True)

    @observe('limit_to__all__')
    def _limit_to_all_changed(self, change):
        warnings.warn('`IPython.core.IPCompleter.limit_to__all__` configuration '
            'value has been deprecated since IPython 5.0, will be made to have '
            'no effects and then removed in future version of IPython.',
            UserWarning)

    def __init__(self, shell=None, namespace=None, global_namespace=None,
                 use_readline=_deprecation_readline_sentinel, config=None, **kwargs):
        """IPCompleter() -> completer

        Return a completer object.

        Parameters
        ----------

        shell
            a pointer to the ipython shell itself.  This is needed
            because this completer knows about magic functions, and those can
            only be accessed via the ipython instance.

        namespace : dict, optional
            an optional dict where completions are performed.

        global_namespace : dict, optional
            secondary optional dict for completions, to
            handle cases (such as IPython embedded inside functions) where
            both Python scopes are visible.
       
        use_readline : bool, optional
            DEPRECATED, ignored since IPython 6.0, will have no effects
        """

        self.magic_escape = ESC_MAGIC
        self.splitter = CompletionSplitter()

        if use_readline is not _deprecation_readline_sentinel:
            warnings.warn('The `use_readline` parameter is deprecated and ignored since IPython 6.0.',
                          DeprecationWarning, stacklevel=2)

        # _greedy_changed() depends on splitter and readline being defined:
        Completer.__init__(self, namespace=namespace, global_namespace=global_namespace,
                            config=config, **kwargs)

        # List where completion matches will be stored
        self.matches = []
        self.shell = shell
        # Regexp to split filenames with spaces in them
        self.space_name_re = re.compile(r'([^\\] )')
        # Hold a local ref. to glob.glob for speed
        self.glob = glob.glob

        # Determine if we are running on 'dumb' terminals, like (X)Emacs
        # buffers, to avoid completion problems.
        term = os.environ.get('TERM','xterm')
        self.dumb_terminal = term in ['dumb','emacs']

        # Special handling of backslashes needed in win32 platforms
        if sys.platform == "win32":
            self.clean_glob = self._clean_glob_win32
        else:
            self.clean_glob = self._clean_glob

        #regexp to parse docstring for function signature
        self.docstring_sig_re = re.compile(r'^[\w|\s.]+\(([^)]*)\).*')
        self.docstring_kwd_re = re.compile(r'[\s|\[]*(\w+)(?:\s*=\s*.*)')
        #use this if positional argument name is also needed
        #= re.compile(r'[\s|\[]*(\w+)(?:\s*=?\s*.*)')

        self.magic_arg_matchers = [
            self.magic_config_matches,
            self.magic_color_matches,
        ]

        # This is set externally by InteractiveShell
        self.custom_completers = None

    @property
    def matchers(self):
        """All active matcher routines for completion"""
        if self.use_jedi:
            return [
                self.file_matches,
                self.magic_matches,
                self.dict_key_matches,
            ]
        else:
            return [
                self.python_matches,
                self.file_matches,
                self.magic_matches,
                self.python_func_kw_matches,
                self.dict_key_matches,
            ]

    def all_completions(self, text):
        """
        Wrapper around the complete method for the benefit of emacs.
        """
        return self.complete(text)[1]

    def _clean_glob(self, text):
        return self.glob("%s*" % text)

    def _clean_glob_win32(self,text):
        return [f.replace("\\","/")
                for f in self.glob("%s*" % text)]

    def file_matches(self, text):
        """Match filenames, expanding ~USER type strings.

        Most of the seemingly convoluted logic in this completer is an
        attempt to handle filenames with spaces in them.  And yet it's not
        quite perfect, because Python's readline doesn't expose all of the
        GNU readline details needed for this to be done correctly.

        For a filename with a space in it, the printed completions will be
        only the parts after what's already been typed (instead of the
        full completions, as is normally done).  I don't think with the
        current (as of Python 2.3) Python readline it's possible to do
        better."""

        # chars that require escaping with backslash - i.e. chars
        # that readline treats incorrectly as delimiters, but we
        # don't want to treat as delimiters in filename matching
        # when escaped with backslash
        if text.startswith('!'):
            text = text[1:]
            text_prefix = u'!'
        else:
            text_prefix = u''

        text_until_cursor = self.text_until_cursor
        # track strings with open quotes
        open_quotes = has_open_quotes(text_until_cursor)

        if '(' in text_until_cursor or '[' in text_until_cursor:
            lsplit = text
        else:
            try:
                # arg_split ~ shlex.split, but with unicode bugs fixed by us
                lsplit = arg_split(text_until_cursor)[-1]
            except ValueError:
                # typically an unmatched ", or backslash without escaped char.
                if open_quotes:
                    lsplit = text_until_cursor.split(open_quotes)[-1]
                else:
                    return []
            except IndexError:
                # tab pressed on empty line
                lsplit = ""

        if not open_quotes and lsplit != protect_filename(lsplit):
            # if protectables are found, do matching on the whole escaped name
            has_protectables = True
            text0,text = text,lsplit
        else:
            has_protectables = False
            text = os.path.expanduser(text)

        if text == "":
            return [text_prefix + protect_filename(f) for f in self.glob("*")]

        # Compute the matches from the filesystem
        if sys.platform == 'win32':
            m0 = self.clean_glob(text)
        else:
            m0 = self.clean_glob(text.replace('\\', ''))

        if has_protectables:
            # If we had protectables, we need to revert our changes to the
            # beginning of filename so that we don't double-write the part
            # of the filename we have so far
            len_lsplit = len(lsplit)
            matches = [text_prefix + text0 +
                       protect_filename(f[len_lsplit:]) for f in m0]
        else:
            if open_quotes:
                # if we have a string with an open quote, we don't need to
                # protect the names beyond the quote (and we _shouldn't_, as
                # it would cause bugs when the filesystem call is made).
                matches = m0 if sys.platform == "win32" else\
                    [protect_filename(f, open_quotes) for f in m0]
            else:
                matches = [text_prefix +
                           protect_filename(f) for f in m0]

        # Mark directories in input list by appending '/' to their names.
        return [x+'/' if os.path.isdir(x) else x for x in matches]

    def magic_matches(self, text):
        """Match magics"""
        # Get all shell magics now rather than statically, so magics loaded at
        # runtime show up too.
        lsm = self.shell.magics_manager.lsmagic()
        line_magics = lsm['line']
        cell_magics = lsm['cell']
        pre = self.magic_escape
        pre2 = pre+pre

        explicit_magic = text.startswith(pre)

        # Completion logic:
        # - user gives %%: only do cell magics
        # - user gives %: do both line and cell magics
        # - no prefix: do both
        # In other words, line magics are skipped if the user gives %% explicitly
        #
        # We also exclude magics that match any currently visible names:
        # https://github.com/ipython/ipython/issues/4877, unless the user has
        # typed a %:
        # https://github.com/ipython/ipython/issues/10754
        bare_text = text.lstrip(pre)
        global_matches = self.global_matches(bare_text)
        if not explicit_magic:
            def matches(magic):
                """
                Filter magics, in particular remove magics that match
                a name present in global namespace.
                """
                return ( magic.startswith(bare_text) and
                         magic not in global_matches )
        else:
            def matches(magic):
                return magic.startswith(bare_text)

        comp = [ pre2+m for m in cell_magics if matches(m)]
        if not text.startswith(pre2):
            comp += [ pre+m for m in line_magics if matches(m)]

        return comp

    def magic_config_matches(self, text:str) -> List[str]:
        """ Match class names and attributes for %config magic """
        texts = text.strip().split()

        if len(texts) > 0 and (texts[0] == 'config' or texts[0] == '%config'):
            # get all configuration classes
            classes = sorted(set([ c for c in self.shell.configurables
                                   if c.__class__.class_traits(config=True)
                                   ]), key=lambda x: x.__class__.__name__)
            classnames = [ c.__class__.__name__ for c in classes ]

            # return all classnames if config or %config is given
            if len(texts) == 1:
                return classnames

            # match classname
            classname_texts = texts[1].split('.')
            classname = classname_texts[0]
            classname_matches = [ c for c in classnames
                                  if c.startswith(classname) ]

            # return matched classes or the matched class with attributes
            if texts[1].find('.') < 0:
                return classname_matches
            elif len(classname_matches) == 1 and \
                            classname_matches[0] == classname:
                cls = classes[classnames.index(classname)].__class__
                help = cls.class_get_help()
                # strip leading '--' from cl-args:
                help = re.sub(re.compile(r'^--', re.MULTILINE), '', help)
                return [ attr.split('=')[0]
                         for attr in help.strip().splitlines()
                         if attr.startswith(texts[1]) ]
        return []

    def magic_color_matches(self, text:str) -> List[str] :
        """ Match color schemes for %colors magic"""
        texts = text.split()
        if text.endswith(' '):
            # .split() strips off the trailing whitespace. Add '' back
            # so that: '%colors ' -> ['%colors', '']
            texts.append('')

        if len(texts) == 2 and (texts[0] == 'colors' or texts[0] == '%colors'):
            prefix = texts[1]
            return [ color for color in InspectColors.keys()
                     if color.startswith(prefix) ]
        return []

    def _jedi_matches(self, cursor_column:int, cursor_line:int, text:str):
        """

        Return a list of :any:`jedi.api.Completions` object from a ``text`` and
        cursor position.

        Parameters
        ----------
        cursor_column : int
            column position of the cursor in ``text``, 0-indexed.
        cursor_line : int
            line position of the cursor in ``text``, 0-indexed
        text : str
            text to complete

        Debugging
        ---------

        If ``IPCompleter.debug`` is ``True`` may return a :any:`_FakeJediCompletion`
        object containing a string with the Jedi debug information attached.
        """
        namespaces = [self.namespace]
        if self.global_namespace is not None:
            namespaces.append(self.global_namespace)

        completion_filter = lambda x:x
        offset = cursor_to_position(text, cursor_line, cursor_column)
        # filter output if we are completing for object members
        if offset:
            pre = text[offset-1]
            if pre == '.':
                if self.omit__names == 2:
                    completion_filter = lambda c:not c.name.startswith('_')
                elif self.omit__names == 1:
                    completion_filter = lambda c:not (c.name.startswith('__') and c.name.endswith('__'))
                elif self.omit__names == 0:
                    completion_filter = lambda x:x
                else:
                    raise ValueError("Don't understand self.omit__names == {}".format(self.omit__names))

        interpreter = jedi.Interpreter(
            text[:offset], namespaces, column=cursor_column, line=cursor_line + 1)
        try_jedi = True

        try:
            # should we check the type of the node is Error ?
            try:
                # jedi < 0.11
                from jedi.parser.tree import ErrorLeaf
            except ImportError:
                # jedi >= 0.11
                from parso.tree import ErrorLeaf

            next_to_last_tree = interpreter._get_module().tree_node.children[-2]
            completing_string = False
            if isinstance(next_to_last_tree, ErrorLeaf):
                completing_string = next_to_last_tree.value.lstrip()[0] in {'"', "'"}
            # if we are in a string jedi is likely not the right candidate for
            # now. Skip it.
            try_jedi = not completing_string
        except Exception as e:
            # many of things can go wrong, we are using private API just don't crash.
            if self.debug:
                print("Error detecting if completing a non-finished string :", e, '|')

        if not try_jedi:
            return []
        try:
            return filter(completion_filter, interpreter.completions())
        except Exception as e:
            if self.debug:
                return [_FakeJediCompletion('Oops Jedi has crashed, please report a bug with the following:\n"""\n%s\ns"""' % (e))]
            else:
                return []

    def python_matches(self, text):
        """Match attributes or global python names"""
        if "." in text:
            try:
                matches = self.attr_matches(text)
                if text.endswith('.') and self.omit__names:
                    if self.omit__names == 1:
                        # true if txt is _not_ a __ name, false otherwise:
                        no__name = (lambda txt:
                                    re.match(r'.*\.__.*?__',txt) is None)
                    else:
                        # true if txt is _not_ a _ name, false otherwise:
                        no__name = (lambda txt:
                                    re.match(r'\._.*?',txt[txt.rindex('.'):]) is None)
                    matches = filter(no__name, matches)
            except NameError:
                # catches <undefined attributes>.<tab>
                matches = []
        else:
            matches = self.global_matches(text)
        return matches

    def _default_arguments_from_docstring(self, doc):
        """Parse the first line of docstring for call signature.

        Docstring should be of the form 'min(iterable[, key=func])\n'.
        It can also parse cython docstring of the form
        'Minuit.migrad(self, int ncall=10000, resume=True, int nsplit=1)'.
        """
        if doc is None:
            return []

        #care only the firstline
        line = doc.lstrip().splitlines()[0]

        #p = re.compile(r'^[\w|\s.]+\(([^)]*)\).*')
        #'min(iterable[, key=func])\n' -> 'iterable[, key=func]'
        sig = self.docstring_sig_re.search(line)
        if sig is None:
            return []
        # iterable[, key=func]' -> ['iterable[' ,' key=func]']
        sig = sig.groups()[0].split(',')
        ret = []
        for s in sig:
            #re.compile(r'[\s|\[]*(\w+)(?:\s*=\s*.*)')
            ret += self.docstring_kwd_re.findall(s)
        return ret

    def _default_arguments(self, obj):
        """Return the list of default arguments of obj if it is callable,
        or empty list otherwise."""
        call_obj = obj
        ret = []
        if inspect.isbuiltin(obj):
            pass
        elif not (inspect.isfunction(obj) or inspect.ismethod(obj)):
            if inspect.isclass(obj):
                #for cython embedsignature=True the constructor docstring
                #belongs to the object itself not __init__
                ret += self._default_arguments_from_docstring(
                            getattr(obj, '__doc__', ''))
                # for classes, check for __init__,__new__
                call_obj = (getattr(obj, '__init__', None) or
                       getattr(obj, '__new__', None))
            # for all others, check if they are __call__able
            elif hasattr(obj, '__call__'):
                call_obj = obj.__call__
        ret += self._default_arguments_from_docstring(
                 getattr(call_obj, '__doc__', ''))

        _keeps = (inspect.Parameter.KEYWORD_ONLY,
                  inspect.Parameter.POSITIONAL_OR_KEYWORD)

        try:
            sig = inspect.signature(call_obj)
            ret.extend(k for k, v in sig.parameters.items() if
                       v.kind in _keeps)
        except ValueError:
            pass

        return list(set(ret))

    def python_func_kw_matches(self,text):
        """Match named parameters (kwargs) of the last open function"""

        if "." in text: # a parameter cannot be dotted
            return []
        try: regexp = self.__funcParamsRegex
        except AttributeError:
            regexp = self.__funcParamsRegex = re.compile(r'''
                '.*?(?<!\\)' |    # single quoted strings or
                ".*?(?<!\\)" |    # double quoted strings or
                \w+          |    # identifier
                \S                # other characters
                ''', re.VERBOSE | re.DOTALL)
        # 1. find the nearest identifier that comes before an unclosed
        # parenthesis before the cursor
        # e.g. for "foo (1+bar(x), pa<cursor>,a=1)", the candidate is "foo"
        tokens = regexp.findall(self.text_until_cursor)
        iterTokens = reversed(tokens); openPar = 0

        for token in iterTokens:
            if token == ')':
                openPar -= 1
            elif token == '(':
                openPar += 1
                if openPar > 0:
                    # found the last unclosed parenthesis
                    break
        else:
            return []
        # 2. Concatenate dotted names ("foo.bar" for "foo.bar(x, pa" )
        ids = []
        isId = re.compile(r'\w+$').match

        while True:
            try:
                ids.append(next(iterTokens))
                if not isId(ids[-1]):
                    ids.pop(); break
                if not next(iterTokens) == '.':
                    break
            except StopIteration:
                break

        # Find all named arguments already assigned to, as to avoid suggesting
        # them again
        usedNamedArgs = set()
        par_level = -1
        for token, next_token in zip(tokens, tokens[1:]):
            if token == '(':
                par_level += 1
            elif token == ')':
                par_level -= 1

            if par_level != 0:
                continue

            if next_token != '=':
                continue

            usedNamedArgs.add(token)

        # lookup the candidate callable matches either using global_matches
        # or attr_matches for dotted names
        if len(ids) == 1:
            callableMatches = self.global_matches(ids[0])
        else:
            callableMatches = self.attr_matches('.'.join(ids[::-1]))
        argMatches = []
        for callableMatch in callableMatches:
            try:
                namedArgs = self._default_arguments(eval(callableMatch,
                                                        self.namespace))
            except:
                continue

            # Remove used named arguments from the list, no need to show twice
            for namedArg in set(namedArgs) - usedNamedArgs:
                if namedArg.startswith(text):
                    argMatches.append(u"%s=" %namedArg)
        return argMatches

    def dict_key_matches(self, text):
        "Match string keys in a dictionary, after e.g. 'foo[' "
        def get_keys(obj):
            # Objects can define their own completions by defining an
            # _ipy_key_completions_() method.
            method = get_real_method(obj, '_ipython_key_completions_')
            if method is not None:
                return method()

            # Special case some common in-memory dict-like types
            if isinstance(obj, dict) or\
               _safe_isinstance(obj, 'pandas', 'DataFrame'):
                try:
                    return list(obj.keys())
                except Exception:
                    return []
            elif _safe_isinstance(obj, 'numpy', 'ndarray') or\
                 _safe_isinstance(obj, 'numpy', 'void'):
                return obj.dtype.names or []
            return []

        try:
            regexps = self.__dict_key_regexps
        except AttributeError:
            dict_key_re_fmt = r'''(?x)
            (  # match dict-referring expression wrt greedy setting
                %s
            )
            \[   # open bracket
            \s*  # and optional whitespace
            ([uUbB]?  # string prefix (r not handled)
                (?:   # unclosed string
                    '(?:[^']|(?<!\\)\\')*
                |
                    "(?:[^"]|(?<!\\)\\")*
                )
            )?
            $
            '''
            regexps = self.__dict_key_regexps = {
                False: re.compile(dict_key_re_fmt % '''
                                  # identifiers separated by .
                                  (?!\d)\w+
                                  (?:\.(?!\d)\w+)*
                                  '''),
                True: re.compile(dict_key_re_fmt % '''
                                 .+
                                 ''')
            }

        match = regexps[self.greedy].search(self.text_until_cursor)
        if match is None:
            return []

        expr, prefix = match.groups()
        try:
            obj = eval(expr, self.namespace)
        except Exception:
            try:
                obj = eval(expr, self.global_namespace)
            except Exception:
                return []

        keys = get_keys(obj)
        if not keys:
            return keys
        closing_quote, token_offset, matches = match_dict_keys(keys, prefix, self.splitter.delims)
        if not matches:
            return matches
        
        # get the cursor position of
        # - the text being completed
        # - the start of the key text
        # - the start of the completion
        text_start = len(self.text_until_cursor) - len(text)
        if prefix:
            key_start = match.start(2)
            completion_start = key_start + token_offset
        else:
            key_start = completion_start = match.end()
        
        # grab the leading prefix, to make sure all completions start with `text`
        if text_start > key_start:
            leading = ''
        else:
            leading = text[text_start:completion_start]
        
        # the index of the `[` character
        bracket_idx = match.end(1)

        # append closing quote and bracket as appropriate
        # this is *not* appropriate if the opening quote or bracket is outside
        # the text given to this method
        suf = ''
        continuation = self.line_buffer[len(self.text_until_cursor):]
        if key_start > text_start and closing_quote:
            # quotes were opened inside text, maybe close them
            if continuation.startswith(closing_quote):
                continuation = continuation[len(closing_quote):]
            else:
                suf += closing_quote
        if bracket_idx > text_start:
            # brackets were opened inside text, maybe close them
            if not continuation.startswith(']'):
                suf += ']'
        
        return [leading + k + suf for k in matches]

    def unicode_name_matches(self, text):
        u"""Match Latex-like syntax for unicode characters base
        on the name of the character.
        
        This does  ``\\GREEK SMALL LETTER ETA`` -> ``η``

        Works only on valid python 3 identifier, or on combining characters that
        will combine to form a valid identifier.
        
        Used on Python 3 only.
        """
        slashpos = text.rfind('\\')
        if slashpos > -1:
            s = text[slashpos+1:]
            try :
                unic = unicodedata.lookup(s)
                # allow combining chars
                if ('a'+unic).isidentifier():
                    return '\\'+s,[unic]
            except KeyError:
                pass
        return u'', []


    def latex_matches(self, text):
        u"""Match Latex syntax for unicode characters.
        
        This does both ``\\alp`` -> ``\\alpha`` and ``\\alpha`` -> ``α``

        Used on Python 3 only.
        """
        slashpos = text.rfind('\\')
        if slashpos > -1:
            s = text[slashpos:]
            if s in latex_symbols:
                # Try to complete a full latex symbol to unicode
                # \\alpha -> α
                return s, [latex_symbols[s]]
            else:
                # If a user has partially typed a latex symbol, give them
                # a full list of options \al -> [\aleph, \alpha]
                matches = [k for k in latex_symbols if k.startswith(s)]
                return s, matches
        return u'', []

    def dispatch_custom_completer(self, text):
        if not self.custom_completers:
            return

        line = self.line_buffer
        if not line.strip():
            return None

        # Create a little structure to pass all the relevant information about
        # the current completion to any custom completer.
        event = SimpleNamespace()
        event.line = line
        event.symbol = text
        cmd = line.split(None,1)[0]
        event.command = cmd
        event.text_until_cursor = self.text_until_cursor

        # for foo etc, try also to find completer for %foo
        if not cmd.startswith(self.magic_escape):
            try_magic = self.custom_completers.s_matches(
                self.magic_escape + cmd)
        else:
            try_magic = []

        for c in itertools.chain(self.custom_completers.s_matches(cmd),
                 try_magic,
                 self.custom_completers.flat_matches(self.text_until_cursor)):
            try:
                res = c(event)
                if res:
                    # first, try case sensitive match
                    withcase = [r for r in res if r.startswith(text)]
                    if withcase:
                        return withcase
                    # if none, then case insensitive ones are ok too
                    text_low = text.lower()
                    return [r for r in res if r.lower().startswith(text_low)]
            except TryNext:
                pass
            except KeyboardInterrupt:
                """
                If custom completer take too long,
                let keyboard interrupt abort and return nothing.
                """
                break

        return None

    def completions(self, text: str, offset: int)->Iterator[Completion]:
        """
        Returns an iterator over the possible completions

        .. warning:: Unstable
            
            This function is unstable, API may change without warning.
            It will also raise unless use in proper context manager.

        Parameters
        ----------
        
        text:str
            Full text of the current input, multi line string.
        offset:int
            Integer representing the position of the cursor in ``text``. Offset
            is 0-based indexed.

        Yields
        ------
            :any:`Completion` object


        The cursor on a text can either be seen as being "in between"
        characters or "On" a character depending on the interface visible to
        the user. For consistency the cursor being on "in between" characters X
        and Y is equivalent to the cursor being "on" character Y, that is to say
        the character the cursor is on is considered as being after the cursor.

        Combining characters may span more that one position in the
        text.


        .. note::

            If ``IPCompleter.debug`` is :any:`True` will yield a ``--jedi/ipython--``
            fake Completion token to distinguish completion returned by Jedi
            and usual IPython completion.

        .. note::
            
            Completions are not completely deduplicated yet. If identical
            completions are coming from different sources this function does not
            ensure that each completion object will only be present once.
        """
        warnings.warn("_complete is a provisional API (as of IPython 6.0). "
                      "It may change without warnings. "
                      "Use in corresponding context manager.",
                      category=ProvisionalCompleterWarning, stacklevel=2)

        seen = set()
        try:
            for c in self._completions(text, offset, _timeout=self.jedi_compute_type_timeout/1000):
                if c and (c in seen):
                    continue
                yield c
                seen.add(c)
        except KeyboardInterrupt:
            """if completions take too long and users send keyboard interrupt,
            do not crash and return ASAP. """
            pass

    def _completions(self, full_text: str, offset: int, *, _timeout)->Iterator[Completion]:
        """
        Core completion module.Same signature as :any:`completions`, with the
        extra `timeout` parameter (in seconds).


        Computing jedi's completion ``.type`` can be quite expensive (it is a
        lazy property) and can require some warm-up, more warm up than just
        computing the ``name`` of a completion. The warm-up can be :

            - Long warm-up the first time a module is encountered after
            install/update: actually build parse/inference tree.

            - first time the module is encountered in a session: load tree from
            disk.

        We don't want to block completions for tens of seconds so we give the
        completer a "budget" of ``_timeout`` seconds per invocation to compute
        completions types, the completions that have not yet been computed will
        be marked as "unknown" an will have a chance to be computed next round
        are things get cached.

        Keep in mind that Jedi is not the only thing treating the completion so
        keep the timeout short-ish as if we take more than 0.3 second we still
        have lots of processing to do.

        """
        deadline = time.monotonic() + _timeout


        before = full_text[:offset]
        cursor_line, cursor_column = position_to_cursor(full_text, offset)

        matched_text, matches, matches_origin, jedi_matches = self._complete(
            full_text=full_text, cursor_line=cursor_line, cursor_pos=cursor_column)

        iter_jm = iter(jedi_matches)
        if _timeout:
            for jm in iter_jm:
                try:
                    type_ = jm.type
                except Exception:
                    if self.debug:
                        print("Error in Jedi getting type of ", jm)
                    type_ = None
                delta = len(jm.name_with_symbols) - len(jm.complete)
                if type_ == 'function':
                    signature = _make_signature(jm)
                else:
                    signature = ''
                yield Completion(start=offset - delta,
                                 end=offset,
                                 text=jm.name_with_symbols,
                                 type=type_,
                                 signature=signature,
                                 _origin='jedi')

                if time.monotonic() > deadline:
                    break

        for jm in iter_jm:
            delta = len(jm.name_with_symbols) - len(jm.complete)
            yield Completion(start=offset - delta,
                             end=offset,
                             text=jm.name_with_symbols,
                             type='<unknown>',  # don't compute type for speed
                             _origin='jedi',
                             signature='')


        start_offset = before.rfind(matched_text)

        # TODO:
        # Suppress this, right now just for debug.
        if jedi_matches and matches and self.debug:
            yield Completion(start=start_offset, end=offset, text='--jedi/ipython--',
                             _origin='debug', type='none', signature='')

        # I'm unsure if this is always true, so let's assert and see if it
        # crash
        assert before.endswith(matched_text)
        for m, t in zip(matches, matches_origin):
            yield Completion(start=start_offset, end=offset, text=m, _origin=t, signature='', type='<unknown>')


    def complete(self, text=None, line_buffer=None, cursor_pos=None):
        """Find completions for the given text and line context.

        Note that both the text and the line_buffer are optional, but at least
        one of them must be given.

        Parameters
        ----------
          text : string, optional
            Text to perform the completion on.  If not given, the line buffer
            is split using the instance's CompletionSplitter object.

          line_buffer : string, optional
            If not given, the completer attempts to obtain the current line
            buffer via readline.  This keyword allows clients which are
            requesting for text completions in non-readline contexts to inform
            the completer of the entire text.

          cursor_pos : int, optional
            Index of the cursor in the full line buffer.  Should be provided by
            remote frontends where kernel has no access to frontend state.

        Returns
        -------
        text : str
          Text that was actually used in the completion.

        matches : list
          A list of completion matches.


        .. note::

            This API is likely to be deprecated and replaced by
            :any:`IPCompleter.completions` in the future.


        """
        warnings.warn('`Completer.complete` is pending deprecation since '
                'IPython 6.0 and will be replaced by `Completer.completions`.',
                      PendingDeprecationWarning)
        # potential todo, FOLD the 3rd throw away argument of _complete
        # into the first 2 one.
        return self._complete(line_buffer=line_buffer, cursor_pos=cursor_pos, text=text, cursor_line=0)[:2]

    def _complete(self, *, cursor_line, cursor_pos, line_buffer=None, text=None,
                  full_text=None) -> Tuple[str, List[str], List[str], Iterable[_FakeJediCompletion]]:
        """

        Like complete but can also returns raw jedi completions as well as the
        origin of the completion text. This could (and should) be made much
        cleaner but that will be simpler once we drop the old (and stateful)
        :any:`complete` API.


        With current provisional API, cursor_pos act both (depending on the
        caller) as the offset in the ``text`` or ``line_buffer``, or as the
        ``column`` when passing multiline strings this could/should be renamed
        but would add extra noise.
        """

        # if the cursor position isn't given, the only sane assumption we can
        # make is that it's at the end of the line (the common case)
        if cursor_pos is None:
            cursor_pos = len(line_buffer) if text is None else len(text)

        if self.use_main_ns:
            self.namespace = __main__.__dict__

        # if text is either None or an empty string, rely on the line buffer
        if (not line_buffer) and full_text:
            line_buffer = full_text.split('\n')[cursor_line]
        if not text:
            text = self.splitter.split_line(line_buffer, cursor_pos)

        if self.backslash_combining_completions:
            # allow deactivation of these on windows.
            base_text = text if not line_buffer else line_buffer[:cursor_pos]
            latex_text, latex_matches = self.latex_matches(base_text)
            if latex_matches:
                return latex_text, latex_matches, ['latex_matches']*len(latex_matches), ()
            name_text = ''
            name_matches = []
            for meth in (self.unicode_name_matches, back_latex_name_matches, back_unicode_name_matches):
                name_text, name_matches = meth(base_text)
                if name_text:
                    return name_text, name_matches[:MATCHES_LIMIT], \
                           [meth.__qualname__]*min(len(name_matches), MATCHES_LIMIT), ()
        

        # If no line buffer is given, assume the input text is all there was
        if line_buffer is None:
            line_buffer = text

        self.line_buffer = line_buffer
        self.text_until_cursor = self.line_buffer[:cursor_pos]

        # Do magic arg matches
        for matcher in self.magic_arg_matchers:
            matches = list(matcher(line_buffer))[:MATCHES_LIMIT]
            if matches:
                origins = [matcher.__qualname__] * len(matches)
                return text, matches, origins, ()

        # Start with a clean slate of completions
        matches = []
        custom_res = self.dispatch_custom_completer(text)
        # FIXME: we should extend our api to return a dict with completions for
        # different types of objects.  The rlcomplete() method could then
        # simply collapse the dict into a list for readline, but we'd have
        # richer completion semantics in other environments.
        completions = ()
        if self.use_jedi:
            if not full_text:
                full_text = line_buffer
            completions = self._jedi_matches(
                cursor_pos, cursor_line, full_text)
        if custom_res is not None:
            # did custom completers produce something?
            matches = [(m, 'custom') for m in custom_res]
        else:
            # Extend the list of completions with the results of each
            # matcher, so we return results to the user from all
            # namespaces.
            if self.merge_completions:
                matches = []
                for matcher in self.matchers:
                    try:
                        matches.extend([(m, matcher.__qualname__)
                                        for m in matcher(text)])
                    except:
                        # Show the ugly traceback if the matcher causes an
                        # exception, but do NOT crash the kernel!
                        sys.excepthook(*sys.exc_info())
            else:
                for matcher in self.matchers:
                    matches = [(m, matcher.__qualname__)
                               for m in matcher(text)]
                    if matches:
                        break
        seen = set()
        filtered_matches = set()
        for m in matches:
            t, c = m
            if t not in seen:
                filtered_matches.add(m)
                seen.add(t)

        _filtered_matches = sorted(
            set(filtered_matches), key=lambda x: completions_sorting_key(x[0]))\
            [:MATCHES_LIMIT]

        _matches = [m[0] for m in _filtered_matches]
        origins = [m[1] for m in _filtered_matches]

        self.matches = _matches

        return text, _matches, origins, completions
