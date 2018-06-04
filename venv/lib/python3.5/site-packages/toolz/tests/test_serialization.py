from toolz import *
import toolz
import toolz.curried
import pickle
from toolz.compatibility import PY3, PY33, PY34
from toolz.utils import raises


def test_compose():
    f = compose(str, sum)
    g = pickle.loads(pickle.dumps(f))
    assert f((1, 2)) == g((1, 2))


def test_curry():
    f = curry(map)(str)
    g = pickle.loads(pickle.dumps(f))
    assert list(f((1, 2, 3))) == list(g((1, 2, 3)))


def test_juxt():
    f = juxt(str, int, bool)
    g = pickle.loads(pickle.dumps(f))
    assert f(1) == g(1)
    assert f.funcs == g.funcs


def test_complement():
    f = complement(bool)
    assert f(True) is False
    assert f(False) is True
    g = pickle.loads(pickle.dumps(f))
    assert f(True) == g(True)
    assert f(False) == g(False)


def test_instanceproperty():
    p = toolz.functoolz.InstanceProperty(bool)
    assert p.__get__(None) is None
    assert p.__get__(0) is False
    assert p.__get__(1) is True
    p2 = pickle.loads(pickle.dumps(p))
    assert p2.__get__(None) is None
    assert p2.__get__(0) is False
    assert p2.__get__(1) is True


def f(x, y):
    return x, y


def test_flip():
    flip = pickle.loads(pickle.dumps(toolz.functoolz.flip))
    assert flip is toolz.functoolz.flip
    g1 = flip(f)
    g2 = pickle.loads(pickle.dumps(g1))
    assert g1(1, 2) == g2(1, 2) == f(2, 1)
    g1 = flip(f)(1)
    g2 = pickle.loads(pickle.dumps(g1))
    assert g1(2) == g2(2) == f(2, 1)


def test_curried_exceptions():
    # This tests a global curried object that isn't defined in toolz.functoolz
    merge = pickle.loads(pickle.dumps(toolz.curried.merge))
    assert merge is toolz.curried.merge


@toolz.curry
class GlobalCurried(object):
    def __init__(self, x, y):
        self.x = x
        self.y = y

    @toolz.curry
    def f1(self, a, b):
        return self.x + self.y + a + b

    def g1(self):
        pass

    def __reduce__(self):
        """Allow us to serialize instances of GlobalCurried"""
        return (GlobalCurried, (self.x, self.y))

    @toolz.curry
    class NestedCurried(object):
        def __init__(self, x, y):
            self.x = x
            self.y = y

        @toolz.curry
        def f2(self, a, b):
            return self.x + self.y + a + b

        def g2(self):
            pass

        def __reduce__(self):
            """Allow us to serialize instances of NestedCurried"""
            return (GlobalCurried.NestedCurried, (self.x, self.y))

    class Nested(object):
        def __init__(self, x, y):
            self.x = x
            self.y = y

        @toolz.curry
        def f3(self, a, b):
            return self.x + self.y + a + b

        def g3(self):
            pass


def test_curried_qualname():
    if not PY3:
        return

    def preserves_identity(obj):
        return pickle.loads(pickle.dumps(obj)) is obj

    assert preserves_identity(GlobalCurried)
    assert preserves_identity(GlobalCurried.func.f1)
    assert preserves_identity(GlobalCurried.func.NestedCurried)
    assert preserves_identity(GlobalCurried.func.NestedCurried.func.f2)
    assert preserves_identity(GlobalCurried.func.Nested.f3)

    global_curried1 = GlobalCurried(1)
    global_curried2 = pickle.loads(pickle.dumps(global_curried1))
    assert global_curried1 is not global_curried2
    assert global_curried1(2).f1(3, 4) == global_curried2(2).f1(3, 4) == 10

    global_curried3 = global_curried1(2)
    global_curried4 = pickle.loads(pickle.dumps(global_curried3))
    assert global_curried3 is not global_curried4
    assert global_curried3.f1(3, 4) == global_curried4.f1(3, 4) == 10

    func1 = global_curried1(2).f1(3)
    func2 = pickle.loads(pickle.dumps(func1))
    assert func1 is not func2
    assert func1(4) == func2(4) == 10

    nested_curried1 = GlobalCurried.func.NestedCurried(1)
    nested_curried2 = pickle.loads(pickle.dumps(nested_curried1))
    assert nested_curried1 is not nested_curried2
    assert nested_curried1(2).f2(3, 4) == nested_curried2(2).f2(3, 4) == 10

    # If we add `curry.__getattr__` forwarding, the following tests will pass

    # if not PY33 and not PY34:
    #     assert preserves_identity(GlobalCurried.func.g1)
    #     assert preserves_identity(GlobalCurried.func.NestedCurried.func.g2)
    #     assert preserves_identity(GlobalCurried.func.Nested)
    #     assert preserves_identity(GlobalCurried.func.Nested.g3)
    #
    # # Rely on curry.__getattr__
    # assert preserves_identity(GlobalCurried.f1)
    # assert preserves_identity(GlobalCurried.NestedCurried)
    # assert preserves_identity(GlobalCurried.NestedCurried.f2)
    # assert preserves_identity(GlobalCurried.Nested.f3)
    # if not PY33 and not PY34:
    #     assert preserves_identity(GlobalCurried.g1)
    #     assert preserves_identity(GlobalCurried.NestedCurried.g2)
    #     assert preserves_identity(GlobalCurried.Nested)
    #     assert preserves_identity(GlobalCurried.Nested.g3)
    #
    # nested_curried3 = nested_curried1(2)
    # nested_curried4 = pickle.loads(pickle.dumps(nested_curried3))
    # assert nested_curried3 is not nested_curried4
    # assert nested_curried3.f2(3, 4) == nested_curried4.f2(3, 4) == 10
    #
    # func1 = nested_curried1(2).f2(3)
    # func2 = pickle.loads(pickle.dumps(func1))
    # assert func1 is not func2
    # assert func1(4) == func2(4) == 10
    #
    # if not PY33 and not PY34:
    #     nested3 = GlobalCurried.func.Nested(1, 2)
    #     nested4 = pickle.loads(pickle.dumps(nested3))
    #     assert nested3 is not nested4
    #     assert nested3.f3(3, 4) == nested4.f3(3, 4) == 10
    #
    #     func1 = nested3.f3(3)
    #     func2 = pickle.loads(pickle.dumps(func1))
    #     assert func1 is not func2
    #     assert func1(4) == func2(4) == 10


def test_curried_bad_qualname():
    @toolz.curry
    class Bad(object):
        __qualname__ = 'toolz.functoolz.not.a.valid.path'

    assert raises(pickle.PicklingError, lambda: pickle.dumps(Bad))
