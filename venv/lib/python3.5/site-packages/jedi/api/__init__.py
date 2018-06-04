"""
The API basically only provides one class. You can create a :class:`Script` and
use its methods.

Additionally you can add a debug function with :func:`set_debug_function`.
Alternatively, if you don't need a custom function and are happy with printing
debug messages to stdout, simply call :func:`set_debug_function` without
arguments.

.. warning:: Please, note that Jedi is **not thread safe**.
"""
import os
import sys

import parso
from parso.python import tree

from jedi._compatibility import force_unicode, is_py3
from jedi.parser_utils import get_executable_nodes
from jedi import debug
from jedi import settings
from jedi import cache
from jedi.api import classes
from jedi.api import interpreter
from jedi.api import helpers
from jedi.api.completion import Completion
from jedi.api.environment import InterpreterEnvironment
from jedi.api.project import get_default_project
from jedi.evaluate import Evaluator
from jedi.evaluate import imports
from jedi.evaluate import usages
from jedi.evaluate.arguments import try_iter_content
from jedi.evaluate.helpers import get_module_names, evaluate_call_of_leaf
from jedi.evaluate.sys_path import dotted_path_in_sys_path
from jedi.evaluate.filters import TreeNameDefinition, ParamName
from jedi.evaluate.syntax_tree import tree_name_to_contexts
from jedi.evaluate.context import ModuleContext
from jedi.evaluate.context.iterable import unpack_tuple_to_dict

# Jedi uses lots and lots of recursion. By setting this a little bit higher, we
# can remove some "maximum recursion depth" errors.
sys.setrecursionlimit(3000)


class Script(object):
    """
    A Script is the base for completions, goto or whatever you want to do with
    |jedi|.

    You can either use the ``source`` parameter or ``path`` to read a file.
    Usually you're going to want to use both of them (in an editor).

    The script might be analyzed in a different ``sys.path`` than |jedi|:

    - if `sys_path` parameter is not ``None``, it will be used as ``sys.path``
      for the script;

    - if `sys_path` parameter is ``None`` and ``VIRTUAL_ENV`` environment
      variable is defined, ``sys.path`` for the specified environment will be
      guessed (see :func:`jedi.evaluate.sys_path.get_venv_path`) and used for
      the script;

    - otherwise ``sys.path`` will match that of |jedi|.

    :param source: The source code of the current file, separated by newlines.
    :type source: str
    :param line: The line to perform actions on (starting with 1).
    :type line: int
    :param column: The column of the cursor (starting with 0).
    :type column: int
    :param path: The path of the file in the file system, or ``''`` if
        it hasn't been saved yet.
    :type path: str or None
    :param encoding: The encoding of ``source``, if it is not a
        ``unicode`` object (default ``'utf-8'``).
    :type encoding: str
    :param source_encoding: The encoding of ``source``, if it is not a
        ``unicode`` object (default ``'utf-8'``).
    :type encoding: str
    :param sys_path: ``sys.path`` to use during analysis of the script
    :type sys_path: list
    :param environment: TODO
    :type sys_path: Environment
    """
    def __init__(self, source=None, line=None, column=None, path=None,
                 encoding='utf-8', sys_path=None, environment=None):
        self._orig_path = path
        # An empty path (also empty string) should always result in no path.
        self.path = os.path.abspath(path) if path else None

        if source is None:
            # TODO add a better warning than the traceback!
            with open(path, 'rb') as f:
                source = f.read()

        # Load the Python grammar of the current interpreter.
        self._grammar = parso.load_grammar()

        if sys_path is not None and not is_py3:
            sys_path = list(map(force_unicode, sys_path))

        # Load the Python grammar of the current interpreter.
        project = get_default_project(
            os.path.dirname(self.path)if path else os.getcwd()
        )
        # TODO deprecate and remove sys_path from the Script API.
        if sys_path is not None:
            project._sys_path = sys_path
        self._evaluator = Evaluator(
            project, environment=environment, script_path=self.path
        )
        self._project = project
        debug.speed('init')
        self._module_node, source = self._evaluator.parse_and_get_code(
            code=source,
            path=self.path,
            cache=False,  # No disk cache, because the current script often changes.
            diff_cache=True,
            cache_path=settings.cache_directory
        )
        debug.speed('parsed')
        self._code_lines = parso.split_lines(source, keepends=True)
        self._code = source
        line = max(len(self._code_lines), 1) if line is None else line
        if not (0 < line <= len(self._code_lines)):
            raise ValueError('`line` parameter is not in a valid range.')

        line_string = self._code_lines[line - 1]
        line_len = len(line_string)
        if line_string.endswith('\r\n'):
            line_len -= 1
        if line_string.endswith('\n'):
            line_len -= 1

        column = line_len if column is None else column
        if not (0 <= column <= line_len):
            raise ValueError('`column` parameter is not in a valid range.')
        self._pos = line, column
        self._path = path

        cache.clear_time_caches()
        debug.reset_time()

    def _get_module(self):
        name = '__main__'
        if self.path is not None:
            n = dotted_path_in_sys_path(self._evaluator.get_sys_path(), self.path)
            if n is not None:
                name = n

        module = ModuleContext(
            self._evaluator, self._module_node, self.path,
            code_lines=self._code_lines
        )
        imports.add_module_to_cache(self._evaluator, name, module)
        return module

    def __repr__(self):
        return '<%s: %s>' % (self.__class__.__name__, repr(self._orig_path))

    def completions(self):
        """
        Return :class:`classes.Completion` objects. Those objects contain
        information about the completions, more than just names.

        :return: Completion objects, sorted by name and __ comes last.
        :rtype: list of :class:`classes.Completion`
        """
        debug.speed('completions start')
        completion = Completion(
            self._evaluator, self._get_module(), self._code_lines,
            self._pos, self.call_signatures
        )
        completions = completion.completions()
        debug.speed('completions end')
        return completions

    def goto_definitions(self):
        """
        Return the definitions of a the path under the cursor.  goto function!
        This follows complicated paths and returns the end, not the first
        definition. The big difference between :meth:`goto_assignments` and
        :meth:`goto_definitions` is that :meth:`goto_assignments` doesn't
        follow imports and statements. Multiple objects may be returned,
        because Python itself is a dynamic language, which means depending on
        an option you can have two different versions of a function.

        :rtype: list of :class:`classes.Definition`
        """
        leaf = self._module_node.get_name_of_position(self._pos)
        if leaf is None:
            leaf = self._module_node.get_leaf_for_position(self._pos)
            if leaf is None:
                return []

        context = self._evaluator.create_context(self._get_module(), leaf)
        definitions = helpers.evaluate_goto_definition(self._evaluator, context, leaf)

        names = [s.name for s in definitions]
        defs = [classes.Definition(self._evaluator, name) for name in names]
        # The additional set here allows the definitions to become unique in an
        # API sense. In the internals we want to separate more things than in
        # the API.
        return helpers.sorted_definitions(set(defs))

    def goto_assignments(self, follow_imports=False):
        """
        Return the first definition found, while optionally following imports.
        Multiple objects may be returned, because Python itself is a
        dynamic language, which means depending on an option you can have two
        different versions of a function.

        :rtype: list of :class:`classes.Definition`
        """
        def filter_follow_imports(names, check):
            for name in names:
                if check(name):
                    for result in filter_follow_imports(name.goto(), check):
                        yield result
                else:
                    yield name

        tree_name = self._module_node.get_name_of_position(self._pos)
        if tree_name is None:
            return []
        context = self._evaluator.create_context(self._get_module(), tree_name)
        names = list(self._evaluator.goto(context, tree_name))

        if follow_imports:
            def check(name):
                return name.is_import()
        else:
            def check(name):
                return isinstance(name, imports.SubModuleName)

        names = filter_follow_imports(names, check)

        defs = [classes.Definition(self._evaluator, d) for d in set(names)]
        return helpers.sorted_definitions(defs)

    def usages(self, additional_module_paths=()):
        """
        Return :class:`classes.Definition` objects, which contain all
        names that point to the definition of the name under the cursor. This
        is very useful for refactoring (renaming), or to show all usages of a
        variable.

        .. todo:: Implement additional_module_paths

        :rtype: list of :class:`classes.Definition`
        """
        tree_name = self._module_node.get_name_of_position(self._pos)
        if tree_name is None:
            # Must be syntax
            return []

        names = usages.usages(self._get_module(), tree_name)

        definitions = [classes.Definition(self._evaluator, n) for n in names]
        return helpers.sorted_definitions(definitions)

    def call_signatures(self):
        """
        Return the function object of the call you're currently in.

        E.g. if the cursor is here::

            abs(# <-- cursor is here

        This would return the ``abs`` function. On the other hand::

            abs()# <-- cursor is here

        This would return an empty list..

        :rtype: list of :class:`classes.CallSignature`
        """
        call_signature_details = \
            helpers.get_call_signature_details(self._module_node, self._pos)
        if call_signature_details is None:
            return []

        context = self._evaluator.create_context(
            self._get_module(),
            call_signature_details.bracket_leaf
        )
        definitions = helpers.cache_call_signatures(
            self._evaluator,
            context,
            call_signature_details.bracket_leaf,
            self._code_lines,
            self._pos
        )
        debug.speed('func_call followed')

        return [classes.CallSignature(self._evaluator, d.name,
                                      call_signature_details.bracket_leaf.start_pos,
                                      call_signature_details.call_index,
                                      call_signature_details.keyword_name_str)
                for d in definitions if hasattr(d, 'py__call__')]

    def _analysis(self):
        self._evaluator.is_analysis = True
        self._evaluator.analysis_modules = [self._module_node]
        module = self._get_module()
        try:
            for node in get_executable_nodes(self._module_node):
                context = module.create_context(node)
                if node.type in ('funcdef', 'classdef'):
                    # Resolve the decorators.
                    tree_name_to_contexts(self._evaluator, context, node.children[1])
                elif isinstance(node, tree.Import):
                    import_names = set(node.get_defined_names())
                    if node.is_nested():
                        import_names |= set(path[-1] for path in node.get_paths())
                    for n in import_names:
                        imports.infer_import(context, n)
                elif node.type == 'expr_stmt':
                    types = context.eval_node(node)
                    for testlist in node.children[:-1:2]:
                        # Iterate tuples.
                        unpack_tuple_to_dict(context, types, testlist)
                else:
                    if node.type == 'name':
                        defs = self._evaluator.goto_definitions(context, node)
                    else:
                        defs = evaluate_call_of_leaf(context, node)
                    try_iter_content(defs)
                self._evaluator.reset_recursion_limitations()

            ana = [a for a in self._evaluator.analysis if self.path == a.path]
            return sorted(set(ana), key=lambda x: x.line)
        finally:
            self._evaluator.is_analysis = False


class Interpreter(Script):
    """
    Jedi API for Python REPLs.

    In addition to completion of simple attribute access, Jedi
    supports code completion based on static code analysis.
    Jedi can complete attributes of object which is not initialized
    yet.

    >>> from os.path import join
    >>> namespace = locals()
    >>> script = Interpreter('join("").up', [namespace])
    >>> print(script.completions()[0].name)
    upper
    """

    def __init__(self, source, namespaces, **kwds):
        """
        Parse `source` and mixin interpreted Python objects from `namespaces`.

        :type source: str
        :arg  source: Code to parse.
        :type namespaces: list of dict
        :arg  namespaces: a list of namespace dictionaries such as the one
                          returned by :func:`locals`.

        Other optional arguments are same as the ones for :class:`Script`.
        If `line` and `column` are None, they are assumed be at the end of
        `source`.
        """
        try:
            namespaces = [dict(n) for n in namespaces]
        except Exception:
            raise TypeError("namespaces must be a non-empty list of dicts.")

        environment = kwds.get('environment', None)
        if environment is None:
            environment = InterpreterEnvironment()
        else:
            if not isinstance(environment, InterpreterEnvironment):
                raise TypeError("The environment needs to be an InterpreterEnvironment subclass.")

        super(Interpreter, self).__init__(source, environment=environment, **kwds)
        self.namespaces = namespaces

    def _get_module(self):
        return interpreter.MixedModuleContext(
            self._evaluator,
            self._module_node,
            self.namespaces,
            path=self.path,
            code_lines=self._code_lines,
        )


def names(source=None, path=None, encoding='utf-8', all_scopes=False,
          definitions=True, references=False, environment=None):
    """
    Returns a list of `Definition` objects, containing name parts.
    This means you can call ``Definition.goto_assignments()`` and get the
    reference of a name.
    The parameters are the same as in :py:class:`Script`, except or the
    following ones:

    :param all_scopes: If True lists the names of all scopes instead of only
        the module namespace.
    :param definitions: If True lists the names that have been defined by a
        class, function or a statement (``a = b`` returns ``a``).
    :param references: If True lists all the names that are not listed by
        ``definitions=True``. E.g. ``a = b`` returns ``b``.
    """
    def def_ref_filter(_def):
        is_def = _def._name.tree_name.is_definition()
        return definitions and is_def or references and not is_def

    def create_name(name):
        if name.parent.type == 'param':
            cls = ParamName
        else:
            cls = TreeNameDefinition
        is_module = name.parent.type == 'file_input'
        return cls(
            module_context.create_context(name if is_module else name.parent),
            name
        )

    # Set line/column to a random position, because they don't matter.
    script = Script(source, line=1, column=0, path=path, encoding=encoding, environment=environment)
    module_context = script._get_module()
    defs = [
        classes.Definition(
            script._evaluator,
            create_name(name)
        ) for name in get_module_names(script._module_node, all_scopes)
    ]
    return sorted(filter(def_ref_filter, defs), key=lambda x: (x.line, x.column))


def preload_module(*modules):
    """
    Preloading modules tells Jedi to load a module now, instead of lazy parsing
    of modules. Usful for IDEs, to control which modules to load on startup.

    :param modules: different module names, list of string.
    """
    for m in modules:
        s = "import %s as x; x." % m
        Script(s, 1, len(s), None).completions()


def set_debug_function(func_cb=debug.print_to_stdout, warnings=True,
                       notices=True, speed=True):
    """
    Define a callback debug function to get all the debug messages.

    If you don't specify any arguments, debug messages will be printed to stdout.

    :param func_cb: The callback function for debug messages, with n params.
    """
    debug.debug_function = func_cb
    debug.enable_warning = warnings
    debug.enable_notice = notices
    debug.enable_speed = speed
