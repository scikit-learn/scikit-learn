import nose.tools as nt
from IPython.utils.dir2 import dir2


class Base(object):
    x = 1
    z = 23


def test_base():
    res = dir2(Base())
    assert ('x' in res)
    assert ('z' in res)
    assert ('y' not in res)
    assert ('__class__' in res)
    nt.assert_equal(res.count('x'), 1)
    nt.assert_equal(res.count('__class__'), 1)

def test_SubClass():

    class SubClass(Base):
        y = 2

    res = dir2(SubClass())
    assert ('y' in res)
    nt.assert_equal(res.count('y'), 1)
    nt.assert_equal(res.count('x'), 1)


def test_SubClass_with_trait_names_attr():
    # usecase: trait_names is used in a class describing psychological classification

    class SubClass(Base):
        y = 2
        trait_names = 44

    res = dir2(SubClass())
    assert('trait_names' in res)


def test_misbehaving_object_without_trait_names():
    # dir2 shouldn't raise even when objects are dumb and raise
    # something other than AttribteErrors on bad getattr.

    class MisbehavingGetattr(object):
        def __getattr__(self):
            raise KeyError("I should be caught")

        def some_method(self):
            pass

    class SillierWithDir(MisbehavingGetattr):
        def __dir__(self):
            return ['some_method']

    for bad_klass in (MisbehavingGetattr, SillierWithDir):
        res = dir2(bad_klass())
        assert('some_method' in res)
