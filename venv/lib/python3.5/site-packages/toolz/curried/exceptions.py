import toolz


__all__ = ['merge_with', 'merge']


@toolz.curry
def merge_with(func, d, *dicts, **kwargs):
    return toolz.merge_with(func, d, *dicts, **kwargs)


@toolz.curry
def merge(d, *dicts, **kwargs):
    return toolz.merge(d, *dicts, **kwargs)


merge_with.__doc__ = toolz.merge_with.__doc__
merge.__doc__ = toolz.merge.__doc__
