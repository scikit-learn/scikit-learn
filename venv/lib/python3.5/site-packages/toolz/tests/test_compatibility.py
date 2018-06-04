from toolz.compatibility import map, filter, iteritems, iterkeys, itervalues


def test_map_filter_are_lazy():
    def bad(x):
        raise Exception()
    map(bad, [1, 2, 3])
    filter(bad, [1, 2, 3])


def test_dict_iteration():
    d = {'a': 1, 'b': 2, 'c': 3}
    assert not isinstance(iteritems(d), list)
    assert not isinstance(iterkeys(d), list)
    assert not isinstance(itervalues(d), list)
    assert set(iteritems(d)) == set(d.items())
    assert set(iterkeys(d)) == set(d.keys())
    assert set(itervalues(d)) == set(d.values())
