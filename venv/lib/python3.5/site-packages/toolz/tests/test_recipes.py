from toolz import first, identity, countby, partitionby


def iseven(x):
    return x % 2 == 0


def test_countby():
    assert countby(iseven, [1, 2, 3]) == {True: 1, False: 2}
    assert countby(len, ['cat', 'dog', 'mouse']) == {3: 2, 5: 1}
    assert countby(0, ('ab', 'ac', 'bc')) == {'a': 2, 'b': 1}


def test_partitionby():
    assert list(partitionby(identity, [])) == []

    vowels = "aeiou"
    assert (list(partitionby(vowels.__contains__, "abcdefghi")) ==
            [("a",), ("b", "c", "d"), ("e",), ("f", "g", "h"), ("i",)])

    assert (list(map(first,
                     partitionby(identity,
                                 [1, 1, 1, 2, 3, 3, 2, 2, 3]))) ==
            [1, 2, 3, 2, 3])

    assert ''.join(map(first,
                       partitionby(identity, "Khhhaaaaannnnn!!!!"))) == 'Khan!'
