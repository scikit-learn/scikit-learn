import os
import inspect
from functools import partial

from .backend import LOKY_PICKLER

try:
    from cloudpickle import dumps, loads
    cloudpickle = True
except ImportError:
    cloudpickle = False

if not LOKY_PICKLER and cloudpickle:
    wrap_cache = dict()

    class CloudpickledObjectWrapper(object):
        def __init__(self, obj):
            self.pickled_obj = dumps(obj)

        def __reduce__(self):
            return loads, (self.pickled_obj,)

    def _wrap_non_picklable_objects(obj):
        need_wrap = "__main__" in getattr(obj, "__module__", "")
        if isinstance(obj, partial):
            return partial(
                _wrap_non_picklable_objects(obj.func),
                *[_wrap_non_picklable_objects(a) for a in obj.args],
                **{k: _wrap_non_picklable_objects(v)
                   for k, v in obj.keywords.items()}
            )
        if callable(obj):
            # Need wrap if the object is a function defined in a local scope of
            # another function.
            func_code = getattr(obj, "__code__", "")
            need_wrap |= getattr(func_code, "co_flags", 0) & inspect.CO_NESTED

            # Need wrap if the obj is a lambda expression
            func_name = getattr(obj, "__name__", "")
            need_wrap |= "<lambda>" in func_name

        if not need_wrap:
            return obj

        wrapped_obj = wrap_cache.get(obj)
        if wrapped_obj is None:
            wrapped_obj = CloudpickledObjectWrapper(obj)
            wrap_cache[obj] = wrapped_obj
        return wrapped_obj

else:
    def _wrap_non_picklable_objects(obj):
        return obj
