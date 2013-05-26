from collections import defaultdict
from six import iteritems


def _default_iter_fits(est, fit, param_iter, *args, **kwargs):
    for params in param_iter:
        est.set_params(**params)
        yield params, fit(*args, **kwargs)


def iter_fits(est, param_iter, *args, **kwargs):
    if hasattr(est, 'iter_fits'):
        return est.iter_fits(param_iter, *args, **kwargs)
    return _default_iter_fits(est, est.fit, param_iter, *args, **kwargs)


def _fit_transform(est):
    if hasattr(est, 'fit_transform'):
        return est.fit_transform
    def fn(X, *args, **kwargs):
        est.fit(X, *args, **kwargs)
        return est.transform(X)
    return fn


def iter_fit_transforms(est, param_iter, X, *args, **kwargs):
    if hasattr(est, 'iter_fit_transforms'):
        return est.iter_fit_transforms(param_iter, X, *args, **kwargs)
    if hasattr(est, 'iter_fits'):
        return ((params, new_est.transform(X))
                for params, new_est
                in est.iter_fits(param_iter, X, *args, **kwargs))
    return _default_iter_fits(est, _fit_transform, param_iter, X, *args, **kwargs)


def group_params(items, key_name=lambda x: True):
    """bin by sub-dicts where values are compared by id if not hashable
    
    >>> a = ('boo', 6)
    >>> b = ('boo', 6)
    >>> id(a) == id(b)
    False
    >>> import numpy
    >>> c = numpy.array([1, 2, 3])
    >>> items = [{'p': a, 'q': 1}, {'p': b, 'q': 2}, {'p': c, 'q': 3}, {'p': c, 'q': 4}]
    >>> groups = list(group_params(items, lambda x: x == "p"))
    >>> len(groups)
    2
    >>> groups
    >>> sorted([[x['q'] for x in gitems] for g, gitems in groups if g['p'] == a][0])
    [1, 2]
    >>> sorted([[x['q'] for x in gitems] for g, gitems in groups if g['p'] is c][0])
    [3, 4]
    """
    # can reduce work if input is sorted tuples rather than dict

    groups = defaultdict(list)
    canonical = {}  # maps hashable x to a canonical instance of x
    id_to_obj = {}

    for params in items:
        group = []
        for k, v in iteritems(params):
            if key_name(k):
                try:
                    v_id = id(canonical.setdefault(v, v))
                except TypeError:
                    v_id = id(v)
                id_to_obj[v_id] = v
                group.append((k, v_id))
        groups[tuple(sorted(group))].append(params)

    return [({k: id_to_obj[v_id] for k, v_id in group}, group_items)
            for group, group_items in iteritems(groups)]


def find_group(haystack, needle):
    """Given the output of group_params, gets the entries for a single group value.
    
    Again designed to deal with non-hashable, non-comparable types..."""
    for group, entries in haystack:
        if len(group) != len(needle):
            continue
        for key, val in iteritems(group):
            try:
                other_val = needle[key]
            except KeyError:
                break
            if val is other_val:
                continue
            try:
                if val == other_val:
                    continue
            except (ValueError, TypeError):
                # can't convert numpy equality to boolean
                continue
            break
        else:
            return entries
    raise ValueError('Could not find {!r} in given groups:\n{!r}'.format(needle, haystack))
