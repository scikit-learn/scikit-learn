"""
Used only for REPL Completion.
"""

import inspect
import os

from jedi.parser_utils import get_cached_code_lines

from jedi import settings
from jedi.evaluate import compiled
from jedi.cache import underscore_memoization
from jedi.evaluate import imports
from jedi.evaluate.base_context import Context, ContextSet
from jedi.evaluate.context import ModuleContext
from jedi.evaluate.cache import evaluator_function_cache
from jedi.evaluate.compiled.getattr_static import getattr_static
from jedi.evaluate.compiled.access import compiled_objects_cache
from jedi.evaluate.compiled.context import create_cached_compiled_object


class MixedObject(object):
    """
    A ``MixedObject`` is used in two ways:

    1. It uses the default logic of ``parser.python.tree`` objects,
    2. except for getattr calls. The names dicts are generated in a fashion
       like ``CompiledObject``.

    This combined logic makes it possible to provide more powerful REPL
    completion. It allows side effects that are not noticable with the default
    parser structure to still be completeable.

    The biggest difference from CompiledObject to MixedObject is that we are
    generally dealing with Python code and not with C code. This will generate
    fewer special cases, because we in Python you don't have the same freedoms
    to modify the runtime.
    """
    def __init__(self, evaluator, parent_context, compiled_object, tree_context):
        self.evaluator = evaluator
        self.parent_context = parent_context
        self.compiled_object = compiled_object
        self._context = tree_context
        self.access_handle = compiled_object.access_handle

    # We have to overwrite everything that has to do with trailers, name
    # lookups and filters to make it possible to route name lookups towards
    # compiled objects and the rest towards tree node contexts.
    def py__getattribute__(*args, **kwargs):
        return Context.py__getattribute__(*args, **kwargs)

    def get_filters(self, *args, **kwargs):
        yield MixedObjectFilter(self.evaluator, self)

    def __repr__(self):
        return '<%s: %s>' % (type(self).__name__, self.access_handle.get_repr())

    def __getattr__(self, name):
        return getattr(self._context, name)


class MixedName(compiled.CompiledName):
    """
    The ``CompiledName._compiled_object`` is our MixedObject.
    """
    @property
    def start_pos(self):
        contexts = list(self.infer())
        if not contexts:
            # This means a start_pos that doesn't exist (compiled objects).
            return 0, 0
        return contexts[0].name.start_pos

    @start_pos.setter
    def start_pos(self, value):
        # Ignore the __init__'s start_pos setter call.
        pass

    @underscore_memoization
    def infer(self):
        access_handle = self.parent_context.access_handle
        # TODO use logic from compiled.CompiledObjectFilter
        access_handle = access_handle.getattr(self.string_name, default=None)
        return ContextSet(
            _create(self._evaluator, access_handle, parent_context=self.parent_context)
        )

    @property
    def api_type(self):
        return next(iter(self.infer())).api_type


class MixedObjectFilter(compiled.CompiledObjectFilter):
    name_class = MixedName

    def __init__(self, evaluator, mixed_object, is_instance=False):
        super(MixedObjectFilter, self).__init__(
            evaluator, mixed_object, is_instance)
        self._mixed_object = mixed_object

    #def _create(self, name):
        #return MixedName(self._evaluator, self._compiled_object, name)


@evaluator_function_cache()
def _load_module(evaluator, path):
    module_node = evaluator.grammar.parse(
        path=path,
        cache=True,
        diff_cache=True,
        cache_path=settings.cache_directory
    ).get_root_node()
    # python_module = inspect.getmodule(python_object)
    # TODO we should actually make something like this possible.
    #evaluator.modules[python_module.__name__] = module_node
    return module_node


def _get_object_to_check(python_object):
    """Check if inspect.getfile has a chance to find the source."""
    if (inspect.ismodule(python_object) or
            inspect.isclass(python_object) or
            inspect.ismethod(python_object) or
            inspect.isfunction(python_object) or
            inspect.istraceback(python_object) or
            inspect.isframe(python_object) or
            inspect.iscode(python_object)):
        return python_object

    try:
        return python_object.__class__
    except AttributeError:
        raise TypeError  # Prevents computation of `repr` within inspect.


def _find_syntax_node_name(evaluator, access_handle):
    # TODO accessing this is bad, but it probably doesn't matter that much,
    # because we're working with interpreteters only here.
    python_object = access_handle.access._obj
    try:
        python_object = _get_object_to_check(python_object)
        path = inspect.getsourcefile(python_object)
    except TypeError:
        # The type might not be known (e.g. class_with_dict.__weakref__)
        return None
    if path is None or not os.path.exists(path):
        # The path might not exist or be e.g. <stdin>.
        return None

    module_node = _load_module(evaluator, path)

    if inspect.ismodule(python_object):
        # We don't need to check names for modules, because there's not really
        # a way to write a module in a module in Python (and also __name__ can
        # be something like ``email.utils``).
        code_lines = get_cached_code_lines(evaluator.grammar, path)
        return module_node, module_node, path, code_lines

    try:
        name_str = python_object.__name__
    except AttributeError:
        # Stuff like python_function.__code__.
        return None

    if name_str == '<lambda>':
        return None  # It's too hard to find lambdas.

    # Doesn't always work (e.g. os.stat_result)
    try:
        names = module_node.get_used_names()[name_str]
    except KeyError:
        return None
    names = [n for n in names if n.is_definition()]

    try:
        code = python_object.__code__
        # By using the line number of a code object we make the lookup in a
        # file pretty easy. There's still a possibility of people defining
        # stuff like ``a = 3; foo(a); a = 4`` on the same line, but if people
        # do so we just don't care.
        line_nr = code.co_firstlineno
    except AttributeError:
        pass
    else:
        line_names = [name for name in names if name.start_pos[0] == line_nr]
        # There's a chance that the object is not available anymore, because
        # the code has changed in the background.
        if line_names:
            names = line_names

    code_lines = get_cached_code_lines(evaluator.grammar, path)
    # It's really hard to actually get the right definition, here as a last
    # resort we just return the last one. This chance might lead to odd
    # completions at some points but will lead to mostly correct type
    # inference, because people tend to define a public name in a module only
    # once.
    return module_node, names[-1].parent, path, code_lines


@compiled_objects_cache('mixed_cache')
def _create(evaluator, access_handle, parent_context, *args):
    compiled_object = create_cached_compiled_object(
        evaluator, access_handle, parent_context=parent_context.compiled_object)

    result = _find_syntax_node_name(evaluator, access_handle)
    if result is None:
        return compiled_object

    module_node, tree_node, path, code_lines = result

    if parent_context.tree_node.get_root_node() == module_node:
        module_context = parent_context.get_root_context()
    else:
        module_context = ModuleContext(
            evaluator, module_node,
            path=path,
            code_lines=code_lines,
        )
        # TODO this __name__ is probably wrong.
        name = compiled_object.get_root_context().py__name__()
        if name is not None:
            imports.add_module_to_cache(evaluator, name, module_context)

    tree_context = module_context.create_context(
        tree_node,
        node_is_context=True,
        node_is_object=True
    )
    if tree_node.type == 'classdef':
        if not access_handle.is_class():
            # Is an instance, not a class.
            tree_context, = tree_context.execute_evaluated()

    return MixedObject(
        evaluator,
        parent_context,
        compiled_object,
        tree_context=tree_context
    )
