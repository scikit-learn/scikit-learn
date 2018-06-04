import pytest

import dask
from dask.order import ndependencies, order
from dask.core import get_deps
from dask.utils_test import add, inc


@pytest.fixture(params=['abcde', 'edcba'])
def abcde(request):
    return request.param


def issorted(L, reverse=False):
    return sorted(L, reverse=reverse) == L


def f(*args):
    pass


def test_ordering_keeps_groups_together(abcde):
    a, b, c, d, e = abcde
    d = dict(((a, i), (f,)) for i in range(4))
    d.update({(b, 0): (f, (a, 0), (a, 1)),
              (b, 1): (f, (a, 2), (a, 3))})
    o = order(d)

    assert abs(o[(a, 0)] - o[(a, 1)]) == 1
    assert abs(o[(a, 2)] - o[(a, 3)]) == 1

    d = dict(((a, i), (f,)) for i in range(4))
    d.update({(b, 0): (f, (a, 0), (a, 2)),
              (b, 1): (f, (a, 1), (a, 3))})
    o = order(d)

    assert abs(o[(a, 0)] - o[(a, 2)]) == 1
    assert abs(o[(a, 1)] - o[(a, 3)]) == 1


@pytest.mark.xfail(reason="Can't please 'em all")
def test_avoid_broker_nodes(abcde):
    """

    b0    b1  b2
    |      \  /
    a0      a1

    a0 should be run before a1
    """
    a, b, c, d, e = abcde
    dsk = {(a, 0): (f,), (a, 1): (f,),
           (b, 0): (f, (a, 0)), (b, 1): (f, (a, 1)), (b, 2): (f, (a, 1))}

    o = order(dsk)

    assert o[(a, 0)] < o[(a, 1)]

    # Switch name of 0, 1 to ensure that this isn't due to string comparison
    dsk = {(a, 1): (f,), (a, 0): (f,),
           (b, 0): (f, (a, 1)), (b, 1): (f, (a, 0)), (b, 2): (f, (a, 0))}

    o = order(dsk)

    assert o[(a, 0)] > o[(a, 1)]


def test_base_of_reduce_preferred(abcde):
    """
               a3
              /|
            a2 |
           /|  |
         a1 |  |
        /|  |  |
      a0 |  |  |
      |  |  |  |
      b0 b1 b2 b3
        \ \ / /
           c

    We really want to run b0 quickly
    """
    a, b, c, d, e = abcde
    dsk = {(a, i): (f, (a, i - 1), (b, i)) for i in [1, 2, 3]}
    dsk[(a, 0)] = (f, (b, 0))
    dsk.update({(b, i): (f, c, 1) for i in [0, 1, 2, 3]})
    dsk[c] = 1

    o = order(dsk)

    assert o[(b, 0)] <= 4
    assert o[(b, 1)] <= 6


@pytest.mark.xfail(reason="Can't please 'em all")
def test_avoid_upwards_branching(abcde):
    """
         a1
         |
         a2
         |
         a3    d1
        /  \  /
      b1    c1
      |     |
      b2    c2
            |
            c3

    Prefer b1 over c1 because it won't stick around waiting for d1 to complete
    """
    a, b, c, d, e = abcde
    dsk = {(a, 1): (f, (a, 2)),
           (a, 2): (f, (a, 3)),
           (a, 3): (f, (b, 1), (c, 1)),
           (b, 1): (f, (b, 2)),
           (c, 1): (f, (c, 2)),
           (c, 2): (f, (c, 3)),
           (d, 1): (f, (c, 1))}

    o = order(dsk)

    assert o[(b, 1)] < o[(c, 1)]


def test_avoid_upwards_branching_complex(abcde):
    """
         a1
         |
    e2   a2  d2  d3
    |    |    \  /
    e1   a3    d1
     \  /  \  /
      b1    c1
      |     |
      b2    c2
            |
            c3

    Prefer c1 over b1 because c1 will stay in memory less long while b1
    computes
    """
    a, b, c, d, e = abcde
    dsk = {(a, 1): (f, (a, 2)),
           (a, 2): (f, (a, 3)),
           (a, 3): (f, (b, 1), (c, 1)),
           (b, 1): (f, (b, 2)),
           (b, 2): (f,),
           (c, 1): (f, (c, 2)),
           (c, 2): (f, (c, 3)),
           (c, 3): (f,),
           (d, 1): (f, (c, 1)),
           (d, 2): (f, (d, 1)),
           (d, 3): (f, (d, 1)),
           (e, 1): (f, (b, 1)),
           (e, 2): (f, (e, 1))}

    o = order(dsk)

    assert o[(c, 1)] < o[(b, 1)]


@pytest.mark.xfail(reason="this case is ambiguous")
def test_deep_bases_win_over_dependents(abcde):
    """
    It's not clear who should run first, e or d

    1.  d is nicer because it exposes parallelism
    2.  e is nicer (hypothetically) because it will be sooner released
        (though in this case we need d to run first regardless)

            a
          / | \   .
         b  c |
        / \ | /
       e    d
    """
    a, b, c, d, e = abcde
    dsk = {a: (f, b, c, d), b: (f, d, e), c: (f, d), d: 1, e: 2}

    o = order(dsk)
    assert o[e] < o[d]
    assert o[d] < o[b] or o[d] < o[c]


def test_prefer_deep(abcde):
    """
        c
        |
    e   b
    |   |
    d   a

    Prefer longer chains first so we should start with c
    """
    a, b, c, d, e = abcde
    dsk = {a: 1, b: (f, a), c: (f, b),
           d: 1, e: (f, d)}

    o = order(dsk)
    assert o[a] < o[d]
    assert o[b] < o[d]


def test_stacklimit(abcde):
    dsk = dict(('x%s' % (i + 1), (inc, 'x%s' % i)) for i in range(10000))
    dependencies, dependents = get_deps(dsk)
    ndependencies(dependencies, dependents)


@pytest.mark.xfail(reason="Can't please 'em all")
def test_break_ties_by_str(abcde):
    a, b, c, d, e = abcde
    dsk = {('x', i): (inc, i) for i in range(10)}
    x_keys = sorted(dsk)
    dsk['y'] = list(x_keys)

    o = order(dsk)
    expected = {'y': 0}
    expected.update({k: i + 1 for i, k in enumerate(x_keys)})

    assert o == expected


def test_order_doesnt_fail_on_mixed_type_keys(abcde):
    order({'x': (inc, 1),
           ('y', 0): (inc, 2),
           'z': (add, 'x', ('y', 0))})


def test_gh_3055():
    da = pytest.importorskip('dask.array')
    A, B = 20, 99
    x = da.random.normal(size=(A, B), chunks=(1, None))
    for _ in range(2):
        y = (x[:, None, :] * x[:, :, None]).cumsum(axis=0)
        x = x.cumsum(axis=0)
    w = (y * x[:, None]).sum(axis=(1,2))

    dsk = dict(w.__dask_graph__())
    o = order(dsk)
    L = [o[k] for k in w.__dask_keys__()]
    assert sum(x < len(o) / 2 for x in L) > len(L) / 3  # some complete quickly
    assert sorted(L) == L  # operate in order


def test_type_comparisions_ok(abcde):
    a, b, c, d, e = abcde
    dsk = {a: 1, (a, 1): 2, (a, b, 1): 3}
    order(dsk)  # this doesn't err


def test_prefer_short_dependents(abcde):
    """

         a
         |
      d  b  e
       \ | /
         c

    Prefer to finish d and e before starting b.  That way c can be released
    during the long computations.
    """
    a, b, c, d, e = abcde
    dsk = {c: (f,), d: (f, c), e: (f, c), b: (f, c), a: (f, b)}

    o = order(dsk)
    assert o[d] < o[b]
    assert o[e] < o[b]


@pytest.mark.xfail(reason="This is challenging to do precisely")
def test_run_smaller_sections(abcde):
    """
            aa
           / |
      b   d  bb dd
     / \ /|  | /
    a   c e  cc

    Prefer to run acb first because then we can get that out of the way
    """
    a, b, c, d, e = abcde
    aa, bb, cc, dd = [x * 2 for x in [a, b, c, d]]

    expected = [a, c, b, e, d, cc, bb, aa, dd]

    log = []

    def f(x):
        def _(*args):
            log.append(x)
        return _

    dsk = {a: (f(a),),
           c: (f(c),),
           e: (f(e),),
           cc: (f(cc),),
           b: (f(b), a, c),
           d: (f(d), c, e),
           bb: (f(bb), cc),
           aa: (f(aa), d, bb),
           dd: (f(dd), cc)}

    dask.get(dsk, [aa, b, dd])  # trigger computation

    assert log == expected


def test_local_parents_of_reduction(abcde):
    """

            c1
            |
        b1  c2
        |  /|
    a1  b2  c3
    |  /|
    a2  b3
    |
    a3

    Prefer to finish a1 stack before proceding to b2
    """
    a, b, c, d, e = abcde
    a1, a2, a3 = [a + i for i in '123']
    b1, b2, b3 = [b + i for i in '123']
    c1, c2, c3 = [c + i for i in '123']

    expected = [a3, a2, a1,
                b3, b2, b1,
                c3, c2, c1]

    log = []

    def f(x):
        def _(*args):
            log.append(x)
        return _

    dsk = {a3: (f(a3),),
           a2: (f(a2), a3),
           a1: (f(a1), a2),
           b3: (f(b3),),
           b2: (f(b2), b3, a2),
           b1: (f(b1), b2),
           c3: (f(c3),),
           c2: (f(c2), c3, b2),
           c1: (f(c1), c2)}

    order(dsk)
    dask.get(dsk, [a1, b1, c1])  # trigger computation

    assert log == expected


def test_nearest_neighbor(abcde):
    """

    a1  a2  a3  a4  a5  a6  a7 a8  a9
     \  |  /  \ |  /  \ |  / \ |  /
        b1      b2      b3     b4

    Want to finish off a local group before moving on.
    This is difficult because all groups are connected.
    """
    a, b, c, _, _ = abcde
    a1, a2, a3, a4, a5, a6, a7, a8, a9 = [a + i for i in '123456789']
    b1, b2, b3, b4 = [b + i for i in '1234']

    dsk = {b1: (f,),
           b2: (f,),
           b3: (f,),
           b4: (f,),
           a1: (f, b1),
           a2: (f, b1),
           a3: (f, b1, b2),
           a4: (f, b2),
           a5: (f, b2, b3),
           a6: (f, b3),
           a7: (f, b3, b4),
           a8: (f, b4),
           a9: (f, b4)}

    o = order(dsk)

    assert 3 < sum(o[a + i] < len(o) / 2 for i in '123456789') < 7
    assert 1 < sum(o[b + i] < len(o) / 2 for i in '1234') < 4
    assert o[min([b1, b2, b3, b4])] == 0


def test_string_ordering():
    """ Prefer ordering tasks by name first """
    dsk = {('a', 1): (f,), ('a', 2): (f,), ('a', 3): (f,)}
    o = order(dsk)
    assert o == {('a', 1): 0,
                 ('a', 2): 1,
                 ('a', 3): 2}


def test_string_ordering_dependents():
    """ Prefer ordering tasks by name first even when in dependencies """
    dsk = {('a', 1): (f, 'b'), ('a', 2): (f, 'b'), ('a', 3): (f, 'b'),
           'b': (f,)}
    o = order(dsk)
    assert o == {'b': 0,
                 ('a', 1): 1,
                 ('a', 2): 2,
                 ('a', 3): 3}
