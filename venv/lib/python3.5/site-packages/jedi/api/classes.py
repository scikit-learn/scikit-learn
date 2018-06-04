"""
The :mod:`jedi.api.classes` module contains the return classes of the API.
These classes are the much bigger part of the whole API, because they contain
the interesting information about completion and goto operations.
"""
import re

from parso.python.tree import search_ancestor

from jedi import settings
from jedi.evaluate.utils import ignored, unite
from jedi.cache import memoize_method
from jedi.evaluate import imports
from jedi.evaluate import compiled
from jedi.evaluate.imports import ImportName
from jedi.evaluate.context import instance
from jedi.evaluate.context import ClassContext, FunctionContext, FunctionExecutionContext
from jedi.api.keywords import KeywordName


def _sort_names_by_start_pos(names):
    return sorted(names, key=lambda s: s.start_pos or (0, 0))


def defined_names(evaluator, context):
    """
    List sub-definitions (e.g., methods in class).

    :type scope: Scope
    :rtype: list of Definition
    """
    filter = next(context.get_filters(search_global=True))
    names = [name for name in filter.values()]
    return [Definition(evaluator, n) for n in _sort_names_by_start_pos(names)]


class BaseDefinition(object):
    _mapping = {
        'posixpath': 'os.path',
        'riscospath': 'os.path',
        'ntpath': 'os.path',
        'os2emxpath': 'os.path',
        'macpath': 'os.path',
        'genericpath': 'os.path',
        'posix': 'os',
        '_io': 'io',
        '_functools': 'functools',
        '_sqlite3': 'sqlite3',
        '__builtin__': '',
        'builtins': '',
    }

    _tuple_mapping = dict((tuple(k.split('.')), v) for (k, v) in {
        'argparse._ActionsContainer': 'argparse.ArgumentParser',
    }.items())

    def __init__(self, evaluator, name):
        self._evaluator = evaluator
        self._name = name
        """
        An instance of :class:`parso.reprsentation.Name` subclass.
        """
        self.is_keyword = isinstance(self._name, KeywordName)

        # generate a path to the definition
        self._module = name.get_root_context()
        if self.in_builtin_module():
            self.module_path = None
        else:
            self.module_path = self._module.py__file__()
            """Shows the file path of a module. e.g. ``/usr/lib/python2.7/os.py``"""

    @property
    def name(self):
        """
        Name of variable/function/class/module.

        For example, for ``x = None`` it returns ``'x'``.

        :rtype: str or None
        """
        return self._name.string_name

    @property
    def type(self):
        """
        The type of the definition.

        Here is an example of the value of this attribute.  Let's consider
        the following source.  As what is in ``variable`` is unambiguous
        to Jedi, :meth:`jedi.Script.goto_definitions` should return a list of
        definition for ``sys``, ``f``, ``C`` and ``x``.

        >>> from jedi import Script
        >>> source = '''
        ... import keyword
        ...
        ... class C:
        ...     pass
        ...
        ... class D:
        ...     pass
        ...
        ... x = D()
        ...
        ... def f():
        ...     pass
        ...
        ... for variable in [keyword, f, C, x]:
        ...     variable'''

        >>> script = Script(source)
        >>> defs = script.goto_definitions()

        Before showing what is in ``defs``, let's sort it by :attr:`line`
        so that it is easy to relate the result to the source code.

        >>> defs = sorted(defs, key=lambda d: d.line)
        >>> defs                           # doctest: +NORMALIZE_WHITESPACE
        [<Definition module keyword>, <Definition class C>,
         <Definition instance D>, <Definition def f>]

        Finally, here is what you can get from :attr:`type`:

        >>> defs = [str(d.type) for d in defs]  # It's unicode and in Py2 has u before it.
        >>> defs[0]
        'module'
        >>> defs[1]
        'class'
        >>> defs[2]
        'instance'
        >>> defs[3]
        'function'

        """
        tree_name = self._name.tree_name
        resolve = False
        if tree_name is not None:
            # TODO move this to their respective names.
            definition = tree_name.get_definition()
            if definition is not None and definition.type == 'import_from' and \
                    tree_name.is_definition():
                resolve = True

        if isinstance(self._name, imports.SubModuleName) or resolve:
            for context in self._name.infer():
                return context.api_type
        return self._name.api_type

    def _path(self):
        """The path to a module/class/function definition."""
        def to_reverse():
            name = self._name
            if name.api_type == 'module':
                try:
                    name = list(name.infer())[0].name
                except IndexError:
                    pass

            if name.api_type in 'module':
                module_contexts = name.infer()
                if module_contexts:
                    module_context, = module_contexts
                    for n in reversed(module_context.py__name__().split('.')):
                        yield n
                else:
                    # We don't really know anything about the path here. This
                    # module is just an import that would lead in an
                    # ImportError. So simply return the name.
                    yield name.string_name
                    return
            else:
                yield name.string_name

            parent_context = name.parent_context
            while parent_context is not None:
                try:
                    method = parent_context.py__name__
                except AttributeError:
                    try:
                        yield parent_context.name.string_name
                    except AttributeError:
                        pass
                else:
                    for name in reversed(method().split('.')):
                        yield name
                parent_context = parent_context.parent_context
        return reversed(list(to_reverse()))

    @property
    def module_name(self):
        """
        The module name.

        >>> from jedi import Script
        >>> source = 'import json'
        >>> script = Script(source, path='example.py')
        >>> d = script.goto_definitions()[0]
        >>> print(d.module_name)                       # doctest: +ELLIPSIS
        json
        """
        return self._module.name.string_name

    def in_builtin_module(self):
        """Whether this is a builtin module."""
        return isinstance(self._module, compiled.CompiledObject)

    @property
    def line(self):
        """The line where the definition occurs (starting with 1)."""
        start_pos = self._name.start_pos
        if start_pos is None:
            return None
        return start_pos[0]

    @property
    def column(self):
        """The column where the definition occurs (starting with 0)."""
        start_pos = self._name.start_pos
        if start_pos is None:
            return None
        return start_pos[1]

    def docstring(self, raw=False, fast=True):
        r"""
        Return a document string for this completion object.

        Example:

        >>> from jedi import Script
        >>> source = '''\
        ... def f(a, b=1):
        ...     "Document for function f."
        ... '''
        >>> script = Script(source, 1, len('def f'), 'example.py')
        >>> doc = script.goto_definitions()[0].docstring()
        >>> print(doc)
        f(a, b=1)
        <BLANKLINE>
        Document for function f.

        Notice that useful extra information is added to the actual
        docstring.  For function, it is call signature.  If you need
        actual docstring, use ``raw=True`` instead.

        >>> print(script.goto_definitions()[0].docstring(raw=True))
        Document for function f.

        :param fast: Don't follow imports that are only one level deep like
            ``import foo``, but follow ``from foo import bar``. This makes
            sense for speed reasons. Completing `import a` is slow if you use
            the ``foo.docstring(fast=False)`` on every object, because it
            parses all libraries starting with ``a``.
        """
        return _Help(self._name).docstring(fast=fast, raw=raw)

    @property
    def description(self):
        """A textual description of the object."""
        return self._name.string_name

    @property
    def full_name(self):
        """
        Dot-separated path of this object.

        It is in the form of ``<module>[.<submodule>[...]][.<object>]``.
        It is useful when you want to look up Python manual of the
        object at hand.

        Example:

        >>> from jedi import Script
        >>> source = '''
        ... import os
        ... os.path.join'''
        >>> script = Script(source, 3, len('os.path.join'), 'example.py')
        >>> print(script.goto_definitions()[0].full_name)
        os.path.join

        Notice that it returns ``'os.path.join'`` instead of (for example)
        ``'posixpath.join'``. This is not correct, since the modules name would
        be ``<module 'posixpath' ...>```. However most users find the latter
        more practical.
        """
        path = list(self._path())
        # TODO add further checks, the mapping should only occur on stdlib.
        if not path:
            return None  # for keywords the path is empty

        with ignored(KeyError):
            path[0] = self._mapping[path[0]]
        for key, repl in self._tuple_mapping.items():
            if tuple(path[:len(key)]) == key:
                path = [repl] + path[len(key):]

        return '.'.join(path if path[0] else path[1:])

    def goto_assignments(self):
        if self._name.tree_name is None:
            return self

        names = self._evaluator.goto(self._name.parent_context, self._name.tree_name)
        return [Definition(self._evaluator, n) for n in names]

    def _goto_definitions(self):
        # TODO make this function public.
        return [Definition(self._evaluator, d.name) for d in self._name.infer()]

    @property
    @memoize_method
    def params(self):
        """
        Raises an ``AttributeError``if the definition is not callable.
        Otherwise returns a list of `Definition` that represents the params.
        """
        def get_param_names(context):
            param_names = []
            if context.api_type == 'function':
                param_names = list(context.get_param_names())
                if isinstance(context, instance.BoundMethod):
                    param_names = param_names[1:]
            elif isinstance(context, (instance.AbstractInstanceContext, ClassContext)):
                if isinstance(context, ClassContext):
                    search = u'__init__'
                else:
                    search = u'__call__'
                names = context.get_function_slot_names(search)
                if not names:
                    return []

                # Just take the first one here, not optimal, but currently
                # there's no better solution.
                inferred = names[0].infer()
                param_names = get_param_names(next(iter(inferred)))
                if isinstance(context, ClassContext):
                    param_names = param_names[1:]
                return param_names
            elif isinstance(context, compiled.CompiledObject):
                return list(context.get_param_names())
            return param_names

        followed = list(self._name.infer())
        if not followed or not hasattr(followed[0], 'py__call__'):
            raise AttributeError()
        context = followed[0]  # only check the first one.

        return [Definition(self._evaluator, n) for n in get_param_names(context)]

    def parent(self):
        context = self._name.parent_context
        if context is None:
            return None

        if isinstance(context, FunctionExecutionContext):
            # TODO the function context should be a part of the function
            # execution context.
            context = FunctionContext(
                self._evaluator, context.parent_context, context.tree_node)
        return Definition(self._evaluator, context.name)

    def __repr__(self):
        return "<%s %s>" % (type(self).__name__, self.description)

    def get_line_code(self, before=0, after=0):
        """
        Returns the line of code where this object was defined.

        :param before: Add n lines before the current line to the output.
        :param after: Add n lines after the current line to the output.

        :return str: Returns the line(s) of code or an empty string if it's a
                     builtin.
        """
        if self.in_builtin_module():
            return ''

        lines = self._name.get_root_context().code_lines

        index = self._name.start_pos[0] - 1
        start_index = max(index - before, 0)
        return ''.join(lines[start_index:index + after + 1])


class Completion(BaseDefinition):
    """
    `Completion` objects are returned from :meth:`api.Script.completions`. They
    provide additional information about a completion.
    """
    def __init__(self, evaluator, name, stack, like_name_length):
        super(Completion, self).__init__(evaluator, name)

        self._like_name_length = like_name_length
        self._stack = stack

        # Completion objects with the same Completion name (which means
        # duplicate items in the completion)
        self._same_name_completions = []

    def _complete(self, like_name):
        append = ''
        if settings.add_bracket_after_function \
                and self.type == 'Function':
            append = '('

        if self._name.api_type == 'param' and self._stack is not None:
            node_names = list(self._stack.get_node_names(self._evaluator.grammar._pgen_grammar))
            if 'trailer' in node_names and 'argument' not in node_names:
                append += '='

        name = self._name.string_name
        if like_name:
            name = name[self._like_name_length:]
        return name + append

    @property
    def complete(self):
        """
        Return the rest of the word, e.g. completing ``isinstance``::

            isinstan# <-- Cursor is here

        would return the string 'ce'. It also adds additional stuff, depending
        on your `settings.py`.

        Assuming the following function definition::

            def foo(param=0):
                pass

        completing ``foo(par`` would give a ``Completion`` which `complete`
        would be `am=`


        """
        return self._complete(True)

    @property
    def name_with_symbols(self):
        """
        Similar to :attr:`name`, but like :attr:`name` returns also the
        symbols, for example assuming the following function definition::

            def foo(param=0):
                pass

        completing ``foo(`` would give a ``Completion`` which
        ``name_with_symbols`` would be "param=".

        """
        return self._complete(False)

    def docstring(self, raw=False, fast=True):
        if self._like_name_length >= 3:
            # In this case we can just resolve the like name, because we
            # wouldn't load like > 100 Python modules anymore.
            fast = False
        return super(Completion, self).docstring(raw=raw, fast=fast)

    @property
    def description(self):
        """Provide a description of the completion object."""
        # TODO improve the class structure.
        return Definition.description.__get__(self)

    def __repr__(self):
        return '<%s: %s>' % (type(self).__name__, self._name.string_name)

    @memoize_method
    def follow_definition(self):
        """
        Return the original definitions. I strongly recommend not using it for
        your completions, because it might slow down |jedi|. If you want to
        read only a few objects (<=20), it might be useful, especially to get
        the original docstrings. The basic problem of this function is that it
        follows all results. This means with 1000 completions (e.g.  numpy),
        it's just PITA-slow.
        """
        defs = self._name.infer()
        return [Definition(self._evaluator, d.name) for d in defs]


class Definition(BaseDefinition):
    """
    *Definition* objects are returned from :meth:`api.Script.goto_assignments`
    or :meth:`api.Script.goto_definitions`.
    """
    def __init__(self, evaluator, definition):
        super(Definition, self).__init__(evaluator, definition)

    @property
    def description(self):
        """
        A description of the :class:`.Definition` object, which is heavily used
        in testing. e.g. for ``isinstance`` it returns ``def isinstance``.

        Example:

        >>> from jedi import Script
        >>> source = '''
        ... def f():
        ...     pass
        ...
        ... class C:
        ...     pass
        ...
        ... variable = f if random.choice([0,1]) else C'''
        >>> script = Script(source, column=3)  # line is maximum by default
        >>> defs = script.goto_definitions()
        >>> defs = sorted(defs, key=lambda d: d.line)
        >>> defs
        [<Definition def f>, <Definition class C>]
        >>> str(defs[0].description)  # strip literals in python2
        'def f'
        >>> str(defs[1].description)
        'class C'

        """
        typ = self.type
        tree_name = self._name.tree_name
        if typ in ('function', 'class', 'module', 'instance') or tree_name is None:
            if typ == 'function':
                # For the description we want a short and a pythonic way.
                typ = 'def'
            return typ + ' ' + self._name.string_name
        elif typ == 'param':
            code = search_ancestor(tree_name, 'param').get_code(
                include_prefix=False,
                include_comma=False
            )
            return typ + ' ' + code

        definition = tree_name.get_definition() or tree_name
        # Remove the prefix, because that's not what we want for get_code
        # here.
        txt = definition.get_code(include_prefix=False)
        # Delete comments:
        txt = re.sub('#[^\n]+\n', ' ', txt)
        # Delete multi spaces/newlines
        txt = re.sub('\s+', ' ', txt).strip()
        return txt

    @property
    def desc_with_module(self):
        """
        In addition to the definition, also return the module.

        .. warning:: Don't use this function yet, its behaviour may change. If
            you really need it, talk to me.

        .. todo:: Add full path. This function is should return a
            `module.class.function` path.
        """
        position = '' if self.in_builtin_module else '@%s' % self.line
        return "%s:%s%s" % (self.module_name, self.description, position)

    @memoize_method
    def defined_names(self):
        """
        List sub-definitions (e.g., methods in class).

        :rtype: list of Definition
        """
        defs = self._name.infer()
        return sorted(
            unite(defined_names(self._evaluator, d) for d in defs),
            key=lambda s: s._name.start_pos or (0, 0)
        )

    def is_definition(self):
        """
        Returns True, if defined as a name in a statement, function or class.
        Returns False, if it's a reference to such a definition.
        """
        if self._name.tree_name is None:
            return True
        else:
            return self._name.tree_name.is_definition()

    def __eq__(self, other):
        return self._name.start_pos == other._name.start_pos \
            and self.module_path == other.module_path \
            and self.name == other.name \
            and self._evaluator == other._evaluator

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self._name.start_pos, self.module_path, self.name, self._evaluator))


class CallSignature(Definition):
    """
    `CallSignature` objects is the return value of `Script.function_definition`.
    It knows what functions you are currently in. e.g. `isinstance(` would
    return the `isinstance` function. without `(` it would return nothing.
    """
    def __init__(self, evaluator, executable_name, bracket_start_pos, index, key_name_str):
        super(CallSignature, self).__init__(evaluator, executable_name)
        self._index = index
        self._key_name_str = key_name_str
        self._bracket_start_pos = bracket_start_pos

    @property
    def index(self):
        """
        The Param index of the current call.
        Returns None if the index cannot be found in the curent call.
        """
        if self._key_name_str is not None:
            for i, param in enumerate(self.params):
                if self._key_name_str == param.name:
                    return i
            if self.params:
                param_name = self.params[-1]._name
                if param_name.tree_name is not None:
                    if param_name.tree_name.get_definition().star_count == 2:
                        return i
            return None

        if self._index >= len(self.params):
            for i, param in enumerate(self.params):
                tree_name = param._name.tree_name
                if tree_name is not None:
                    # *args case
                    if tree_name.get_definition().star_count == 1:
                        return i
            return None
        return self._index

    @property
    def bracket_start(self):
        """
        The indent of the bracket that is responsible for the last function
        call.
        """
        return self._bracket_start_pos

    def __repr__(self):
        return '<%s: %s index %s>' % \
            (type(self).__name__, self._name.string_name, self.index)


class _Help(object):
    """
    Temporary implementation, will be used as `Script.help() or something in
    the future.
    """
    def __init__(self, definition):
        self._name = definition

    @memoize_method
    def _get_contexts(self, fast):
        if isinstance(self._name, ImportName) and fast:
            return {}

        if self._name.api_type == 'statement':
            return {}

        return self._name.infer()

    def docstring(self, fast=True, raw=True):
        """
        The docstring ``__doc__`` for any object.

        See :attr:`doc` for example.
        """
        # TODO: Use all of the followed objects as output. Possibly divinding
        # them by a few dashes.
        for context in self._get_contexts(fast=fast):
            return context.py__doc__(include_call_signature=not raw)

        return ''
