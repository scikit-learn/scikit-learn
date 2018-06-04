"""
- the popular ``_memoize_default`` works like a typical memoize and returns the
  default otherwise.
- ``CachedMetaClass`` uses ``_memoize_default`` to do the same with classes.
"""

_NO_DEFAULT = object()


def _memoize_default(default=_NO_DEFAULT, evaluator_is_first_arg=False, second_arg_is_evaluator=False):
    """ This is a typical memoization decorator, BUT there is one difference:
    To prevent recursion it sets defaults.

    Preventing recursion is in this case the much bigger use than speed. I
    don't think, that there is a big speed difference, but there are many cases
    where recursion could happen (think about a = b; b = a).
    """
    def func(function):
        def wrapper(obj, *args, **kwargs):
            # TODO These checks are kind of ugly and slow.
            if evaluator_is_first_arg:
                cache = obj.memoize_cache
            elif second_arg_is_evaluator:
                cache = args[0].memoize_cache  # needed for meta classes
            else:
                cache = obj.evaluator.memoize_cache

            try:
                memo = cache[function]
            except KeyError:
                memo = {}
                cache[function] = memo

            key = (obj, args, frozenset(kwargs.items()))
            if key in memo:
                return memo[key]
            else:
                if default is not _NO_DEFAULT:
                    memo[key] = default
                rv = function(obj, *args, **kwargs)
                memo[key] = rv
                return rv
        return wrapper

    return func


def evaluator_function_cache(default=_NO_DEFAULT):
    def decorator(func):
        return _memoize_default(default=default, evaluator_is_first_arg=True)(func)

    return decorator


def evaluator_method_cache(default=_NO_DEFAULT):
    def decorator(func):
        return _memoize_default(default=default)(func)

    return decorator


def evaluator_as_method_param_cache():
    def decorator(call):
        return _memoize_default(second_arg_is_evaluator=True)(call)

    return decorator


class CachedMetaClass(type):
    """
    This is basically almost the same than the decorator above, it just caches
    class initializations. Either you do it this way or with decorators, but
    with decorators you lose class access (isinstance, etc).
    """
    @evaluator_as_method_param_cache()
    def __call__(self, *args, **kwargs):
        return super(CachedMetaClass, self).__call__(*args, **kwargs)
