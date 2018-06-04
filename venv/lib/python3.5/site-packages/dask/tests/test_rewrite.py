from dask.rewrite import RewriteRule, RuleSet, head, args, VAR, Traverser
from dask.utils_test import inc, add


def double(x):
    return x * 2


def test_head():
    assert head((inc, 1)) == inc
    assert head((add, 1, 2)) == add
    assert head((add, (inc, 1), (inc, 1))) == add
    assert head([1, 2, 3]) == list


def test_args():
    assert args((inc, 1)) == (1,)
    assert args((add, 1, 2)) == (1, 2)
    assert args(1) == ()
    assert args([1, 2, 3]) == [1, 2, 3]


def test_traverser():
    term = (add, (inc, 1), (double, (inc, 1), 2))
    t = Traverser(term)
    t2 = t.copy()
    assert t.current == add
    t.next()
    assert t.current == inc
    # Ensure copies aren't advanced when the original advances
    assert t2.current == add
    t.skip()
    assert t.current == double
    t.next()
    assert t.current == inc
    assert list(t2) == [add, inc, 1, double, inc, 1, 2]


vars = ("a", "b", "c")
# add(a, 1) -> inc(a)
rule1 = RewriteRule((add, "a", 1), (inc, "a"), vars)
# add(a, a) -> double(a)
rule2 = RewriteRule((add, "a", "a"), (double, "a"), vars)
# add(inc(a), inc(a)) -> add(double(a), 2)
rule3 = RewriteRule((add, (inc, "a"), (inc, "a")), (add, (double, "a"), 2), vars)
# add(inc(b), inc(a)) -> add(add(a, b), 2)
rule4 = RewriteRule((add, (inc, "b"), (inc, "a")), (add, (add, "a", "b"), 2), vars)
# sum([c, b, a]) -> add(add(a, b), c)
rule5 = RewriteRule((sum, ["c", "b", "a"]), (add, (add, "a", "b"), "c"), vars)
# list(x) -> x if x is a list


def repl_list(sd):
    x = sd['x']
    if isinstance(x, list):
        return x
    else:
        return (list, x)


rule6 = RewriteRule((list, 'x'), repl_list, ('x',))


def test_RewriteRule():
    # Test extraneous vars are removed, varlist is correct
    assert rule1.vars == ("a",)
    assert rule1._varlist == ["a"]
    assert rule2.vars == ("a",)
    assert rule2._varlist == ["a", "a"]
    assert rule3.vars == ("a",)
    assert rule3._varlist == ["a", "a"]
    assert rule4.vars == ("a", "b")
    assert rule4._varlist == ["b", "a"]
    assert rule5.vars == ("a", "b", "c")
    assert rule5._varlist == ["c", "b", "a"]


def test_RewriteRuleSubs():
    # Test both rhs substitution and callable rhs
    assert rule1.subs({'a': 1}) == (inc, 1)
    assert rule6.subs({'x': [1, 2, 3]}) == [1, 2, 3]


rules = [rule1, rule2, rule3, rule4, rule5, rule6]
rs = RuleSet(*rules)


def test_RuleSet():
    net = ({add: ({VAR: ({VAR: ({}, [1]), 1: ({}, [0])}, []),
            inc: ({VAR: ({inc: ({VAR: ({}, [2, 3])}, [])}, [])}, [])}, []),
            list: ({VAR: ({}, [5])}, []),
            sum: ({list: ({VAR: ({VAR: ({VAR: ({}, [4])}, [])}, [])}, [])}, [])}, [])
    assert rs._net == net
    assert rs.rules == rules


def test_matches():
    term = (add, 2, 1)
    matches = list(rs.iter_matches(term))
    assert len(matches) == 1
    assert matches[0] == (rule1, {'a': 2})
    # Test matches specific before general
    term = (add, 1, 1)
    matches = list(rs.iter_matches(term))
    assert len(matches) == 2
    assert matches[0] == (rule1, {'a': 1})
    assert matches[1] == (rule2, {'a': 1})
    # Test matches unhashable. What it's getting rewritten to doesn't make
    # sense, this is just to test that it works. :)
    term = (add, [1], [1])
    matches = list(rs.iter_matches(term))
    assert len(matches) == 1
    assert matches[0] == (rule2, {'a': [1]})
    # Test match at depth
    term = (add, (inc, 1), (inc, 1))
    matches = list(rs.iter_matches(term))
    assert len(matches) == 3
    assert matches[0] == (rule3, {'a': 1})
    assert matches[1] == (rule4, {'a': 1, 'b': 1})
    assert matches[2] == (rule2, {'a': (inc, 1)})
    # Test non-linear pattern checking
    term = (add, 2, 3)
    matches = list(rs.iter_matches(term))
    assert len(matches) == 0


def test_rewrite():
    # Rewrite inside list
    term = (sum, [(add, 1, 1), (add, 1, 1), (add, 1, 1)])
    new_term = rs.rewrite(term)
    assert new_term == (add, (add, (inc, 1), (inc, 1)), (inc, 1))
    # Rules aren't applied to exhaustion, this can be further simplified
    new_term = rs.rewrite(new_term)
    assert new_term == (add, (add, (double, 1), 2), (inc, 1))
    term = (add, (add, (add, (add, 1, 2), (add, 1, 2)), (add, (add, 1, 2), (add, 1, 2))), 1)
    assert rs.rewrite(term) == (inc, (double, (double, (add, 1, 2))))
    # Callable RewriteRule rhs
    term = (list, [1, 2, 3])
    assert rs.rewrite(term) == [1, 2, 3]
    term = (list, (map, inc, [1, 2, 3]))
    assert rs.rewrite(term) == term
