"""
Recursions are the recipe of |jedi| to conquer Python code. However, someone
must stop recursions going mad. Some settings are here to make |jedi| stop at
the right time. You can read more about them :ref:`here <settings-recursion>`.

Next to :mod:`jedi.evaluate.cache` this module also makes |jedi| not
thread-safe. Why?  ``execution_recursion_decorator`` uses class variables to
count the function calls.

.. _settings-recursion:

Settings
~~~~~~~~~~

Recursion settings are important if you don't want extremly
recursive python code to go absolutely crazy.

The default values are based on experiments while completing the |jedi| library
itself (inception!). But I don't think there's any other Python library that
uses recursion in a similarly extreme way. Completion should also be fast and
therefore the quality might not always be maximal.

.. autodata:: recursion_limit
.. autodata:: total_function_execution_limit
.. autodata:: per_function_execution_limit
.. autodata:: per_function_recursion_limit
"""

from contextlib import contextmanager

from jedi import debug
from jedi.evaluate.base_context import NO_CONTEXTS


recursion_limit = 15
"""
Like ``sys.getrecursionlimit()``, just for |jedi|.
"""
total_function_execution_limit = 200
"""
This is a hard limit of how many non-builtin functions can be executed.
"""
per_function_execution_limit = 6
"""
The maximal amount of times a specific function may be executed.
"""
per_function_recursion_limit = 2
"""
A function may not be executed more than this number of times recursively.
"""


class RecursionDetector(object):
    def __init__(self):
        self.pushed_nodes = []


@contextmanager
def execution_allowed(evaluator, node):
    """
    A decorator to detect recursions in statements. In a recursion a statement
    at the same place, in the same module may not be executed two times.
    """
    pushed_nodes = evaluator.recursion_detector.pushed_nodes

    if node in pushed_nodes:
        debug.warning('catched stmt recursion: %s @%s', node,
                      node.start_pos)
        yield False
    else:
        try:
            pushed_nodes.append(node)
            yield True
        finally:
            pushed_nodes.pop()


def execution_recursion_decorator(default=NO_CONTEXTS):
    def decorator(func):
        def wrapper(execution, **kwargs):
            detector = execution.evaluator.execution_recursion_detector
            allowed = detector.push_execution(execution)
            try:
                if allowed:
                    result = default
                else:
                    result = func(execution, **kwargs)
            finally:
                detector.pop_execution()
            return result
        return wrapper
    return decorator


class ExecutionRecursionDetector(object):
    """
    Catches recursions of executions.
    """
    def __init__(self, evaluator):
        self._evaluator = evaluator

        self._recursion_level = 0
        self._parent_execution_funcs = []
        self._funcdef_execution_counts = {}
        self._execution_count = 0

    def pop_execution(self):
        self._parent_execution_funcs.pop()
        self._recursion_level -= 1

    def push_execution(self, execution):
        funcdef = execution.tree_node

        # These two will be undone in pop_execution.
        self._recursion_level += 1
        self._parent_execution_funcs.append(funcdef)

        module = execution.get_root_context()
        if module == self._evaluator.builtins_module:
            # We have control over builtins so we know they are not recursing
            # like crazy. Therefore we just let them execute always, because
            # they usually just help a lot with getting good results.
            return False

        if self._recursion_level > recursion_limit:
            return True

        if self._execution_count >= total_function_execution_limit:
            return True
        self._execution_count += 1

        if self._funcdef_execution_counts.setdefault(funcdef, 0) >= per_function_execution_limit:
            return True
        self._funcdef_execution_counts[funcdef] += 1

        if self._parent_execution_funcs.count(funcdef) > per_function_recursion_limit:
            return True
        return False
