"""
TODO Some parts of this module are still not well documented.
"""

from jedi.evaluate.context import ModuleContext
from jedi.evaluate import compiled
from jedi.evaluate.compiled import mixed
from jedi.evaluate.compiled.access import create_access_path
from jedi.evaluate.base_context import Context


def _create(evaluator, obj):
    return compiled.create_from_access_path(
        evaluator, create_access_path(evaluator, obj)
    )


class NamespaceObject(object):
    def __init__(self, dct):
        self.__dict__ = dct


class MixedModuleContext(Context):
    type = 'mixed_module'

    def __init__(self, evaluator, tree_module, namespaces, path, code_lines):
        self.evaluator = evaluator
        self._namespaces = namespaces

        self._namespace_objects = [NamespaceObject(n) for n in namespaces]
        self._module_context = ModuleContext(
            evaluator, tree_module,
            path=path,
            code_lines=code_lines
        )
        self.tree_node = tree_module

    def get_node(self):
        return self.tree_node

    def get_filters(self, *args, **kwargs):
        for filter in self._module_context.get_filters(*args, **kwargs):
            yield filter

        for namespace_obj in self._namespace_objects:
            compiled_object = _create(self.evaluator, namespace_obj)
            mixed_object = mixed.MixedObject(
                self.evaluator,
                parent_context=self,
                compiled_object=compiled_object,
                tree_context=self._module_context
            )
            for filter in mixed_object.get_filters(*args, **kwargs):
                yield filter

    @property
    def code_lines(self):
        return self._module_context.code_lines

    def __getattr__(self, name):
        return getattr(self._module_context, name)
