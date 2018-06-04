import toolz
import toolz.curried
from toolz.curried import (take, first, second, sorted, merge_with, reduce,
                           merge, operator as cop)
from toolz.compatibility import import_module
from collections import defaultdict
from operator import add


def test_take():
    assert list(take(2)([1, 2, 3])) == [1, 2]


def test_first():
    assert first is toolz.itertoolz.first


def test_merge():
    assert merge(factory=lambda: defaultdict(int))({1: 1}) == {1: 1}
    assert merge({1: 1}) == {1: 1}
    assert merge({1: 1}, factory=lambda: defaultdict(int)) == {1: 1}


def test_merge_with():
    assert merge_with(sum)({1: 1}, {1: 2}) == {1: 3}


def test_merge_with_list():
    assert merge_with(sum, [{'a': 1}, {'a': 2}]) == {'a': 3}


def test_sorted():
    assert sorted(key=second)([(1, 2), (2, 1)]) == [(2, 1), (1, 2)]


def test_reduce():
    assert reduce(add)((1, 2, 3)) == 6


def test_module_name():
    assert toolz.curried.__name__ == 'toolz.curried'


def test_curried_operator():
    for k, v in vars(cop).items():
        if not callable(v):
            continue

        if not isinstance(v, toolz.curry):
            try:
                # Make sure it is unary
                v(1)
            except TypeError:
                try:
                    v('x')
                except TypeError:
                    pass
                else:
                    continue
                raise AssertionError(
                    'toolz.curried.operator.%s is not curried!' % k,
                )

    # Make sure this isn't totally empty.
    assert len(set(vars(cop)) & set(['add', 'sub', 'mul'])) == 3


def test_curried_namespace():
    exceptions = import_module('toolz.curried.exceptions')
    namespace = {}

    def should_curry(func):
        if not callable(func) or isinstance(func, toolz.curry):
            return False
        nargs = toolz.functoolz.num_required_args(func)
        if nargs is None or nargs > 1:
            return True
        return nargs == 1 and toolz.functoolz.has_keywords(func)


    def curry_namespace(ns):
        return dict(
            (name, toolz.curry(f) if should_curry(f) else f)
            for name, f in ns.items() if '__' not in name
        )

    from_toolz = curry_namespace(vars(toolz))
    from_exceptions = curry_namespace(vars(exceptions))
    namespace.update(toolz.merge(from_toolz, from_exceptions))

    namespace = toolz.valfilter(callable, namespace)
    curried_namespace = toolz.valfilter(callable, toolz.curried.__dict__)

    if namespace != curried_namespace:
        missing = set(namespace) - set(curried_namespace)
        if missing:
            raise AssertionError('There are missing functions in toolz.curried:\n    %s'
                                 % '    \n'.join(sorted(missing)))
        extra = set(curried_namespace) - set(namespace)
        if extra:
            raise AssertionError('There are extra functions in toolz.curried:\n    %s'
                                 % '    \n'.join(sorted(extra)))
        unequal = toolz.merge_with(list, namespace, curried_namespace)
        unequal = toolz.valfilter(lambda x: x[0] != x[1], unequal)
        messages = []
        for name, (orig_func, auto_func) in sorted(unequal.items()):
            if name in from_exceptions:
                messages.append('%s should come from toolz.curried.exceptions' % name)
            elif should_curry(getattr(toolz, name)):
                messages.append('%s should be curried from toolz' % name)
            else:
                messages.append('%s should come from toolz and NOT be curried' % name)
        raise AssertionError('\n'.join(messages))
