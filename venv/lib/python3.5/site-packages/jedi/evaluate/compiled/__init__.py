from jedi._compatibility import unicode
from jedi.evaluate.compiled.context import CompiledObject, CompiledName, \
    CompiledObjectFilter, CompiledContextName, create_from_access_path, \
    create_from_name


def builtin_from_name(evaluator, string):
    builtins = evaluator.builtins_module
    return create_from_name(evaluator, builtins, string)


def create_simple_object(evaluator, obj):
    """
    Only allows creations of objects that are easily picklable across Python
    versions.
    """
    assert isinstance(obj, (int, float, str, bytes, unicode, slice, complex))
    return create_from_access_path(
        evaluator,
        evaluator.compiled_subprocess.create_simple_object(obj)
    )


def get_special_object(evaluator, identifier):
    return create_from_access_path(
        evaluator,
        evaluator.compiled_subprocess.get_special_object(identifier)
    )


def get_string_context_set(evaluator):
    return builtin_from_name(evaluator, u'str').execute_evaluated()


def load_module(evaluator, **kwargs):
    access_path = evaluator.compiled_subprocess.load_module(**kwargs)
    if access_path is None:
        return None
    return create_from_access_path(evaluator, access_path)
