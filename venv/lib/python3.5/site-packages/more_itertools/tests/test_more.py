from __future__ import division, print_function, unicode_literals

from decimal import Decimal
from doctest import DocTestSuite
from fractions import Fraction
from functools import partial, reduce
from heapq import merge
from io import StringIO
from itertools import (
    chain,
    count,
    groupby,
    islice,
    permutations,
    product,
    repeat,
)
from operator import add, itemgetter
from unittest import TestCase

from six.moves import filter, map, range, zip

import more_itertools as mi


def load_tests(loader, tests, ignore):
    # Add the doctests
    tests.addTests(DocTestSuite('more_itertools.more'))
    return tests


class CollateTests(TestCase):
    """Unit tests for ``collate()``"""
    # Also accidentally tests peekable, though that could use its own tests

    def test_default(self):
        """Test with the default `key` function."""
        iterables = [range(4), range(7), range(3, 6)]
        self.assertEqual(
            sorted(reduce(list.__add__, [list(it) for it in iterables])),
            list(mi.collate(*iterables))
        )

    def test_key(self):
        """Test using a custom `key` function."""
        iterables = [range(5, 0, -1), range(4, 0, -1)]
        actual = sorted(
            reduce(list.__add__, [list(it) for it in iterables]), reverse=True
        )
        expected = list(mi.collate(*iterables, key=lambda x: -x))
        self.assertEqual(actual, expected)

    def test_empty(self):
        """Be nice if passed an empty list of iterables."""
        self.assertEqual([], list(mi.collate()))

    def test_one(self):
        """Work when only 1 iterable is passed."""
        self.assertEqual([0, 1], list(mi.collate(range(2))))

    def test_reverse(self):
        """Test the `reverse` kwarg."""
        iterables = [range(4, 0, -1), range(7, 0, -1), range(3, 6, -1)]

        actual = sorted(
            reduce(list.__add__, [list(it) for it in iterables]), reverse=True
        )
        expected = list(mi.collate(*iterables, reverse=True))
        self.assertEqual(actual, expected)

    def test_alias(self):
        self.assertNotEqual(merge.__doc__, mi.collate.__doc__)
        self.assertNotEqual(partial.__doc__, mi.collate.__doc__)


class ChunkedTests(TestCase):
    """Tests for ``chunked()``"""

    def test_even(self):
        """Test when ``n`` divides evenly into the length of the iterable."""
        self.assertEqual(
            list(mi.chunked('ABCDEF', 3)), [['A', 'B', 'C'], ['D', 'E', 'F']]
        )

    def test_odd(self):
        """Test when ``n`` does not divide evenly into the length of the
        iterable.

        """
        self.assertEqual(
            list(mi.chunked('ABCDE', 3)), [['A', 'B', 'C'], ['D', 'E']]
        )


class FirstTests(TestCase):
    """Tests for ``first()``"""

    def test_many(self):
        """Test that it works on many-item iterables."""
        # Also try it on a generator expression to make sure it works on
        # whatever those return, across Python versions.
        self.assertEqual(mi.first(x for x in range(4)), 0)

    def test_one(self):
        """Test that it doesn't raise StopIteration prematurely."""
        self.assertEqual(mi.first([3]), 3)

    def test_empty_stop_iteration(self):
        """It should raise StopIteration for empty iterables."""
        self.assertRaises(ValueError, lambda: mi.first([]))

    def test_default(self):
        """It should return the provided default arg for empty iterables."""
        self.assertEqual(mi.first([], 'boo'), 'boo')


class PeekableTests(TestCase):
    """Tests for ``peekable()`` behavor not incidentally covered by testing
    ``collate()``

    """
    def test_peek_default(self):
        """Make sure passing a default into ``peek()`` works."""
        p = mi.peekable([])
        self.assertEqual(p.peek(7), 7)

    def test_truthiness(self):
        """Make sure a ``peekable`` tests true iff there are items remaining in
        the iterable.

        """
        p = mi.peekable([])
        self.assertFalse(p)

        p = mi.peekable(range(3))
        self.assertTrue(p)

    def test_simple_peeking(self):
        """Make sure ``next`` and ``peek`` advance and don't advance the
        iterator, respectively.

        """
        p = mi.peekable(range(10))
        self.assertEqual(next(p), 0)
        self.assertEqual(p.peek(), 1)
        self.assertEqual(next(p), 1)

    def test_indexing(self):
        """
        Indexing into the peekable shouldn't advance the iterator.
        """
        p = mi.peekable('abcdefghijkl')

        # The 0th index is what ``next()`` will return
        self.assertEqual(p[0], 'a')
        self.assertEqual(next(p), 'a')

        # Indexing further into the peekable shouldn't advance the itertor
        self.assertEqual(p[2], 'd')
        self.assertEqual(next(p), 'b')

        # The 0th index moves up with the iterator; the last index follows
        self.assertEqual(p[0], 'c')
        self.assertEqual(p[9], 'l')

        self.assertEqual(next(p), 'c')
        self.assertEqual(p[8], 'l')

        # Negative indexing should work too
        self.assertEqual(p[-2], 'k')
        self.assertEqual(p[-9], 'd')
        self.assertRaises(IndexError, lambda: p[-10])

    def test_slicing(self):
        """Slicing the peekable shouldn't advance the iterator."""
        seq = list('abcdefghijkl')
        p = mi.peekable(seq)

        # Slicing the peekable should just be like slicing a re-iterable
        self.assertEqual(p[1:4], seq[1:4])

        # Advancing the iterator moves the slices up also
        self.assertEqual(next(p), 'a')
        self.assertEqual(p[1:4], seq[1:][1:4])

        # Implicit starts and stop should work
        self.assertEqual(p[:5], seq[1:][:5])
        self.assertEqual(p[:], seq[1:][:])

        # Indexing past the end should work
        self.assertEqual(p[:100], seq[1:][:100])

        # Steps should work, including negative
        self.assertEqual(p[::2], seq[1:][::2])
        self.assertEqual(p[::-1], seq[1:][::-1])

    def test_slicing_reset(self):
        """Test slicing on a fresh iterable each time"""
        iterable = ['0', '1', '2', '3', '4', '5']
        indexes = list(range(-4, len(iterable) + 4)) + [None]
        steps = [1, 2, 3, 4, -1, -2, -3, 4]
        for slice_args in product(indexes, indexes, steps):
            it = iter(iterable)
            p = mi.peekable(it)
            next(p)
            index = slice(*slice_args)
            actual = p[index]
            expected = iterable[1:][index]
            self.assertEqual(actual, expected, slice_args)

    def test_slicing_error(self):
        iterable = '01234567'
        p = mi.peekable(iter(iterable))

        # Prime the cache
        p.peek()
        old_cache = list(p._cache)

        # Illegal slice
        with self.assertRaises(ValueError):
            p[1:-1:0]

        # Neither the cache nor the iteration should be affected
        self.assertEqual(old_cache, list(p._cache))
        self.assertEqual(list(p), list(iterable))

    def test_passthrough(self):
        """Iterating a peekable without using ``peek()`` or ``prepend()``
        should just give the underlying iterable's elements (a trivial test but
        useful to set a baseline in case something goes wrong)"""
        expected = [1, 2, 3, 4, 5]
        actual = list(mi.peekable(expected))
        self.assertEqual(actual, expected)

    # prepend() behavior tests

    def test_prepend(self):
        """Tests intersperesed ``prepend()`` and ``next()`` calls"""
        it = mi.peekable(range(2))
        actual = []

        # Test prepend() before next()
        it.prepend(10)
        actual += [next(it), next(it)]

        # Test prepend() between next()s
        it.prepend(11)
        actual += [next(it), next(it)]

        # Test prepend() after source iterable is consumed
        it.prepend(12)
        actual += [next(it)]

        expected = [10, 0, 11, 1, 12]
        self.assertEqual(actual, expected)

    def test_multi_prepend(self):
        """Tests prepending multiple items and getting them in proper order"""
        it = mi.peekable(range(5))
        actual = [next(it), next(it)]
        it.prepend(10, 11, 12)
        it.prepend(20, 21)
        actual += list(it)
        expected = [0, 1, 20, 21, 10, 11, 12, 2, 3, 4]
        self.assertEqual(actual, expected)

    def test_empty(self):
        """Tests prepending in front of an empty iterable"""
        it = mi.peekable([])
        it.prepend(10)
        actual = list(it)
        expected = [10]
        self.assertEqual(actual, expected)

    def test_prepend_truthiness(self):
        """Tests that ``__bool__()`` or ``__nonzero__()`` works properly
        with ``prepend()``"""
        it = mi.peekable(range(5))
        self.assertTrue(it)
        actual = list(it)
        self.assertFalse(it)
        it.prepend(10)
        self.assertTrue(it)
        actual += [next(it)]
        self.assertFalse(it)
        expected = [0, 1, 2, 3, 4, 10]
        self.assertEqual(actual, expected)

    def test_multi_prepend_peek(self):
        """Tests prepending multiple elements and getting them in reverse order
        while peeking"""
        it = mi.peekable(range(5))
        actual = [next(it), next(it)]
        self.assertEqual(it.peek(), 2)
        it.prepend(10, 11, 12)
        self.assertEqual(it.peek(), 10)
        it.prepend(20, 21)
        self.assertEqual(it.peek(), 20)
        actual += list(it)
        self.assertFalse(it)
        expected = [0, 1, 20, 21, 10, 11, 12, 2, 3, 4]
        self.assertEqual(actual, expected)

    def test_prepend_after_stop(self):
        """Test resuming iteration after a previous exhaustion"""
        it = mi.peekable(range(3))
        self.assertEqual(list(it), [0, 1, 2])
        self.assertRaises(StopIteration, lambda: next(it))
        it.prepend(10)
        self.assertEqual(next(it), 10)
        self.assertRaises(StopIteration, lambda: next(it))

    def test_prepend_slicing(self):
        """Tests interaction between prepending and slicing"""
        seq = list(range(20))
        p = mi.peekable(seq)

        p.prepend(30, 40, 50)
        pseq = [30, 40, 50] + seq  # pseq for prepended_seq

        # adapt the specific tests from test_slicing
        self.assertEqual(p[0], 30)
        self.assertEqual(p[1:8], pseq[1:8])
        self.assertEqual(p[1:], pseq[1:])
        self.assertEqual(p[:5], pseq[:5])
        self.assertEqual(p[:], pseq[:])
        self.assertEqual(p[:100], pseq[:100])
        self.assertEqual(p[::2], pseq[::2])
        self.assertEqual(p[::-1], pseq[::-1])

    def test_prepend_indexing(self):
        """Tests interaction between prepending and indexing"""
        seq = list(range(20))
        p = mi.peekable(seq)

        p.prepend(30, 40, 50)

        self.assertEqual(p[0], 30)
        self.assertEqual(next(p), 30)
        self.assertEqual(p[2], 0)
        self.assertEqual(next(p), 40)
        self.assertEqual(p[0], 50)
        self.assertEqual(p[9], 8)
        self.assertEqual(next(p), 50)
        self.assertEqual(p[8], 8)
        self.assertEqual(p[-2], 18)
        self.assertEqual(p[-9], 11)
        self.assertRaises(IndexError, lambda: p[-21])

    def test_prepend_iterable(self):
        """Tests prepending from an iterable"""
        it = mi.peekable(range(5))
        # Don't directly use the range() object to avoid any range-specific
        # optimizations
        it.prepend(*(x for x in range(5)))
        actual = list(it)
        expected = list(chain(range(5), range(5)))
        self.assertEqual(actual, expected)

    def test_prepend_many(self):
        """Tests that prepending a huge number of elements works"""
        it = mi.peekable(range(5))
        # Don't directly use the range() object to avoid any range-specific
        # optimizations
        it.prepend(*(x for x in range(20000)))
        actual = list(it)
        expected = list(chain(range(20000), range(5)))
        self.assertEqual(actual, expected)

    def test_prepend_reversed(self):
        """Tests prepending from a reversed iterable"""
        it = mi.peekable(range(3))
        it.prepend(*reversed((10, 11, 12)))
        actual = list(it)
        expected = [12, 11, 10, 0, 1, 2]
        self.assertEqual(actual, expected)


class ConsumerTests(TestCase):
    """Tests for ``consumer()``"""

    def test_consumer(self):
        @mi.consumer
        def eater():
            while True:
                x = yield  # noqa

        e = eater()
        e.send('hi')  # without @consumer, would raise TypeError


class DistinctPermutationsTests(TestCase):
    def test_distinct_permutations(self):
        """Make sure the output for ``distinct_permutations()`` is the same as
        set(permutations(it)).

        """
        iterable = ['z', 'a', 'a', 'q', 'q', 'q', 'y']
        test_output = sorted(mi.distinct_permutations(iterable))
        ref_output = sorted(set(permutations(iterable)))
        self.assertEqual(test_output, ref_output)


class IlenTests(TestCase):
    def test_ilen(self):
        """Sanity-checks for ``ilen()``."""
        # Non-empty
        self.assertEqual(
            mi.ilen(filter(lambda x: x % 10 == 0, range(101))), 11
        )

        # Empty
        self.assertEqual(mi.ilen((x for x in range(0))), 0)

        # Iterable with __len__
        self.assertEqual(mi.ilen(list(range(6))), 6)


class WithIterTests(TestCase):
    def test_with_iter(self):
        s = StringIO('One fish\nTwo fish')
        initial_words = [line.split()[0] for line in mi.with_iter(s)]

        # Iterable's items should be faithfully represented
        self.assertEqual(initial_words, ['One', 'Two'])
        # The file object should be closed
        self.assertEqual(s.closed, True)


class OneTests(TestCase):
    def test_basic(self):
        it = iter(['item'])
        self.assertEqual(mi.one(it), 'item')

    def test_too_short(self):
        it = iter([])
        self.assertRaises(ValueError, lambda: mi.one(it))
        self.assertRaises(IndexError, lambda: mi.one(it, too_short=IndexError))

    def test_too_long(self):
        it = count()
        self.assertRaises(ValueError, lambda: mi.one(it))  # burn 0 and 1
        self.assertEqual(next(it), 2)
        self.assertRaises(
            OverflowError, lambda: mi.one(it, too_long=OverflowError)
        )


class IntersperseTest(TestCase):
    """ Tests for intersperse() """

    def test_even(self):
        iterable = (x for x in '01')
        self.assertEqual(
            list(mi.intersperse(None, iterable)), ['0', None, '1']
        )

    def test_odd(self):
        iterable = (x for x in '012')
        self.assertEqual(
            list(mi.intersperse(None, iterable)), ['0', None, '1', None, '2']
        )

    def test_nested(self):
        element = ('a', 'b')
        iterable = (x for x in '012')
        actual = list(mi.intersperse(element, iterable))
        expected = ['0', ('a', 'b'), '1', ('a', 'b'), '2']
        self.assertEqual(actual, expected)

    def test_not_iterable(self):
        self.assertRaises(TypeError, lambda: mi.intersperse('x', 1))

    def test_n(self):
        for n, element, expected in [
            (1, '_', ['0', '_', '1', '_', '2', '_', '3', '_', '4', '_', '5']),
            (2, '_', ['0', '1', '_', '2', '3', '_', '4', '5']),
            (3, '_', ['0', '1', '2', '_', '3', '4', '5']),
            (4, '_', ['0', '1', '2', '3', '_', '4', '5']),
            (5, '_', ['0', '1', '2', '3', '4', '_', '5']),
            (6, '_', ['0', '1', '2', '3', '4', '5']),
            (7, '_', ['0', '1', '2', '3', '4', '5']),
            (3, ['a', 'b'], ['0', '1', '2', ['a', 'b'], '3', '4', '5']),
        ]:
            iterable = (x for x in '012345')
            actual = list(mi.intersperse(element, iterable, n=n))
            self.assertEqual(actual, expected)

    def test_n_zero(self):
        self.assertRaises(
            ValueError, lambda: list(mi.intersperse('x', '012', n=0))
        )


class UniqueToEachTests(TestCase):
    """Tests for ``unique_to_each()``"""

    def test_all_unique(self):
        """When all the input iterables are unique the output should match
        the input."""
        iterables = [[1, 2], [3, 4, 5], [6, 7, 8]]
        self.assertEqual(mi.unique_to_each(*iterables), iterables)

    def test_duplicates(self):
        """When there are duplicates in any of the input iterables that aren't
        in the rest, those duplicates should be emitted."""
        iterables = ["mississippi", "missouri"]
        self.assertEqual(
            mi.unique_to_each(*iterables), [['p', 'p'], ['o', 'u', 'r']]
        )

    def test_mixed(self):
        """When the input iterables contain different types the function should
        still behave properly"""
        iterables = ['x', (i for i in range(3)), [1, 2, 3], tuple()]
        self.assertEqual(mi.unique_to_each(*iterables), [['x'], [0], [3], []])


class WindowedTests(TestCase):
    """Tests for ``windowed()``"""

    def test_basic(self):
        actual = list(mi.windowed([1, 2, 3, 4, 5], 3))
        expected = [(1, 2, 3), (2, 3, 4), (3, 4, 5)]
        self.assertEqual(actual, expected)

    def test_large_size(self):
        """
        When the window size is larger than the iterable, and no fill value is
        given,``None`` should be filled in.
        """
        actual = list(mi.windowed([1, 2, 3, 4, 5], 6))
        expected = [(1, 2, 3, 4, 5, None)]
        self.assertEqual(actual, expected)

    def test_fillvalue(self):
        """
        When sizes don't match evenly, the given fill value should be used.
        """
        iterable = [1, 2, 3, 4, 5]

        for n, kwargs, expected in [
            (6, {}, [(1, 2, 3, 4, 5, '!')]),  # n > len(iterable)
            (3, {'step': 3}, [(1, 2, 3), (4, 5, '!')]),  # using ``step``
        ]:
            actual = list(mi.windowed(iterable, n, fillvalue='!', **kwargs))
            self.assertEqual(actual, expected)

    def test_zero(self):
        """When the window size is zero, an empty tuple should be emitted."""
        actual = list(mi.windowed([1, 2, 3, 4, 5], 0))
        expected = [tuple()]
        self.assertEqual(actual, expected)

    def test_negative(self):
        """When the window size is negative, ValueError should be raised."""
        with self.assertRaises(ValueError):
            list(mi.windowed([1, 2, 3, 4, 5], -1))

    def test_step(self):
        """The window should advance by the number of steps provided"""
        iterable = [1, 2, 3, 4, 5, 6, 7]
        for n, step, expected in [
            (3, 2, [(1, 2, 3), (3, 4, 5), (5, 6, 7)]),  # n > step
            (3, 3, [(1, 2, 3), (4, 5, 6), (7, None, None)]),  # n == step
            (3, 4, [(1, 2, 3), (5, 6, 7)]),  # line up nicely
            (3, 5, [(1, 2, 3), (6, 7, None)]),  # off by one
            (3, 6, [(1, 2, 3), (7, None, None)]),  # off by two
            (3, 7, [(1, 2, 3)]),  # step past the end
            (7, 8, [(1, 2, 3, 4, 5, 6, 7)]),  # step > len(iterable)
        ]:
            actual = list(mi.windowed(iterable, n, step=step))
            self.assertEqual(actual, expected)

        # Step must be greater than or equal to 1
        with self.assertRaises(ValueError):
            list(mi.windowed(iterable, 3, step=0))


class BucketTests(TestCase):
    """Tests for ``bucket()``"""

    def test_basic(self):
        iterable = [10, 20, 30, 11, 21, 31, 12, 22, 23, 33]
        D = mi.bucket(iterable, key=lambda x: 10 * (x // 10))

        # In-order access
        self.assertEqual(list(D[10]), [10, 11, 12])

        # Out of order access
        self.assertEqual(list(D[30]), [30, 31, 33])
        self.assertEqual(list(D[20]), [20, 21, 22, 23])

        self.assertEqual(list(D[40]), [])  # Nothing in here!

    def test_in(self):
        iterable = [10, 20, 30, 11, 21, 31, 12, 22, 23, 33]
        D = mi.bucket(iterable, key=lambda x: 10 * (x // 10))

        self.assertTrue(10 in D)
        self.assertFalse(40 in D)
        self.assertTrue(20 in D)
        self.assertFalse(21 in D)

        # Checking in-ness shouldn't advance the iterator
        self.assertEqual(next(D[10]), 10)

    def test_validator(self):
        iterable = count(0)
        key = lambda x: int(str(x)[0])  # First digit of each number
        validator = lambda x: 0 < x < 10  # No leading zeros
        D = mi.bucket(iterable, key, validator=validator)
        self.assertEqual(mi.take(3, D[1]), [1, 10, 11])
        self.assertNotIn(0, D)  # Non-valid entries don't return True
        self.assertNotIn(0, D._cache)  # Don't store non-valid entries
        self.assertEqual(list(D[0]), [])


class SpyTests(TestCase):
    """Tests for ``spy()``"""

    def test_basic(self):
        original_iterable = iter('abcdefg')
        head, new_iterable = mi.spy(original_iterable)
        self.assertEqual(head, ['a'])
        self.assertEqual(
            list(new_iterable), ['a', 'b', 'c', 'd', 'e', 'f', 'g']
        )

    def test_unpacking(self):
        original_iterable = iter('abcdefg')
        (first, second, third), new_iterable = mi.spy(original_iterable, 3)
        self.assertEqual(first, 'a')
        self.assertEqual(second, 'b')
        self.assertEqual(third, 'c')
        self.assertEqual(
            list(new_iterable), ['a', 'b', 'c', 'd', 'e', 'f', 'g']
        )

    def test_too_many(self):
        original_iterable = iter('abc')
        head, new_iterable = mi.spy(original_iterable, 4)
        self.assertEqual(head, ['a', 'b', 'c'])
        self.assertEqual(list(new_iterable), ['a', 'b', 'c'])

    def test_zero(self):
        original_iterable = iter('abc')
        head, new_iterable = mi.spy(original_iterable, 0)
        self.assertEqual(head, [])
        self.assertEqual(list(new_iterable), ['a', 'b', 'c'])


class InterleaveTests(TestCase):
    def test_even(self):
        actual = list(mi.interleave([1, 4, 7], [2, 5, 8], [3, 6, 9]))
        expected = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        self.assertEqual(actual, expected)

    def test_short(self):
        actual = list(mi.interleave([1, 4], [2, 5, 7], [3, 6, 8]))
        expected = [1, 2, 3, 4, 5, 6]
        self.assertEqual(actual, expected)

    def test_mixed_types(self):
        it_list = ['a', 'b', 'c', 'd']
        it_str = '12345'
        it_inf = count()
        actual = list(mi.interleave(it_list, it_str, it_inf))
        expected = ['a', '1', 0, 'b', '2', 1, 'c', '3', 2, 'd', '4', 3]
        self.assertEqual(actual, expected)


class InterleaveLongestTests(TestCase):
    def test_even(self):
        actual = list(mi.interleave_longest([1, 4, 7], [2, 5, 8], [3, 6, 9]))
        expected = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        self.assertEqual(actual, expected)

    def test_short(self):
        actual = list(mi.interleave_longest([1, 4], [2, 5, 7], [3, 6, 8]))
        expected = [1, 2, 3, 4, 5, 6, 7, 8]
        self.assertEqual(actual, expected)

    def test_mixed_types(self):
        it_list = ['a', 'b', 'c', 'd']
        it_str = '12345'
        it_gen = (x for x in range(3))
        actual = list(mi.interleave_longest(it_list, it_str, it_gen))
        expected = ['a', '1', 0, 'b', '2', 1, 'c', '3', 2, 'd', '4', '5']
        self.assertEqual(actual, expected)


class TestCollapse(TestCase):
    """Tests for ``collapse()``"""

    def test_collapse(self):
        l = [[1], 2, [[3], 4], [[[5]]]]
        self.assertEqual(list(mi.collapse(l)), [1, 2, 3, 4, 5])

    def test_collapse_to_string(self):
        l = [["s1"], "s2", [["s3"], "s4"], [[["s5"]]]]
        self.assertEqual(list(mi.collapse(l)), ["s1", "s2", "s3", "s4", "s5"])

    def test_collapse_flatten(self):
        l = [[1], [2], [[3], 4], [[[5]]]]
        self.assertEqual(list(mi.collapse(l, levels=1)), list(mi.flatten(l)))

    def test_collapse_to_level(self):
        l = [[1], 2, [[3], 4], [[[5]]]]
        self.assertEqual(list(mi.collapse(l, levels=2)), [1, 2, 3, 4, [5]])
        self.assertEqual(
            list(mi.collapse(mi.collapse(l, levels=1), levels=1)),
            list(mi.collapse(l, levels=2))
        )

    def test_collapse_to_list(self):
        l = (1, [2], (3, [4, (5,)], 'ab'))
        actual = list(mi.collapse(l, base_type=list))
        expected = [1, [2], 3, [4, (5,)], 'ab']
        self.assertEqual(actual, expected)


class SideEffectTests(TestCase):
    """Tests for ``side_effect()``"""

    def test_individual(self):
        # The function increments the counter for each call
        counter = [0]

        def func(arg):
            counter[0] += 1

        result = list(mi.side_effect(func, range(10)))
        self.assertEqual(result, list(range(10)))
        self.assertEqual(counter[0], 10)

    def test_chunked(self):
        # The function increments the counter for each call
        counter = [0]

        def func(arg):
            counter[0] += 1

        result = list(mi.side_effect(func, range(10), 2))
        self.assertEqual(result, list(range(10)))
        self.assertEqual(counter[0], 5)

    def test_before_after(self):
        f = StringIO()
        collector = []

        def func(item):
            print(item, file=f)
            collector.append(f.getvalue())

        def it():
            yield u'a'
            yield u'b'
            raise RuntimeError('kaboom')

        before = lambda: print('HEADER', file=f)
        after = f.close

        try:
            mi.consume(mi.side_effect(func, it(), before=before, after=after))
        except RuntimeError:
            pass

        # The iterable should have been written to the file
        self.assertEqual(collector, [u'HEADER\na\n', u'HEADER\na\nb\n'])

        # The file should be closed even though something bad happened
        self.assertTrue(f.closed)

    def test_before_fails(self):
        f = StringIO()
        func = lambda x: print(x, file=f)

        def before():
            raise RuntimeError('ouch')

        try:
            mi.consume(
                mi.side_effect(func, u'abc', before=before, after=f.close)
            )
        except RuntimeError:
            pass

        # The file should be closed even though something bad happened in the
        # before function
        self.assertTrue(f.closed)


class SlicedTests(TestCase):
    """Tests for ``sliced()``"""

    def test_even(self):
        """Test when the length of the sequence is divisible by *n*"""
        seq = 'ABCDEFGHI'
        self.assertEqual(list(mi.sliced(seq, 3)), ['ABC', 'DEF', 'GHI'])

    def test_odd(self):
        """Test when the length of the sequence is not divisible by *n*"""
        seq = 'ABCDEFGHI'
        self.assertEqual(list(mi.sliced(seq, 4)), ['ABCD', 'EFGH', 'I'])

    def test_not_sliceable(self):
        seq = (x for x in 'ABCDEFGHI')

        with self.assertRaises(TypeError):
            list(mi.sliced(seq, 3))


class SplitAtTests(TestCase):
    """Tests for ``split()``"""

    def comp_with_str_split(self, str_to_split, delim):
        pred = lambda c: c == delim
        actual = list(map(''.join, mi.split_at(str_to_split, pred)))
        expected = str_to_split.split(delim)
        self.assertEqual(actual, expected)

    def test_seperators(self):
        test_strs = ['', 'abcba', 'aaabbbcccddd', 'e']
        for s, delim in product(test_strs, 'abcd'):
            self.comp_with_str_split(s, delim)


class SplitBeforeTest(TestCase):
    """Tests for ``split_before()``"""

    def test_starts_with_sep(self):
        actual = list(mi.split_before('xooxoo', lambda c: c == 'x'))
        expected = [['x', 'o', 'o'], ['x', 'o', 'o']]
        self.assertEqual(actual, expected)

    def test_ends_with_sep(self):
        actual = list(mi.split_before('ooxoox', lambda c: c == 'x'))
        expected = [['o', 'o'], ['x', 'o', 'o'], ['x']]
        self.assertEqual(actual, expected)

    def test_no_sep(self):
        actual = list(mi.split_before('ooo', lambda c: c == 'x'))
        expected = [['o', 'o', 'o']]
        self.assertEqual(actual, expected)


class SplitAfterTest(TestCase):
    """Tests for ``split_after()``"""

    def test_starts_with_sep(self):
        actual = list(mi.split_after('xooxoo', lambda c: c == 'x'))
        expected = [['x'], ['o', 'o', 'x'], ['o', 'o']]
        self.assertEqual(actual, expected)

    def test_ends_with_sep(self):
        actual = list(mi.split_after('ooxoox', lambda c: c == 'x'))
        expected = [['o', 'o', 'x'], ['o', 'o', 'x']]
        self.assertEqual(actual, expected)

    def test_no_sep(self):
        actual = list(mi.split_after('ooo', lambda c: c == 'x'))
        expected = [['o', 'o', 'o']]
        self.assertEqual(actual, expected)


class PaddedTest(TestCase):
    """Tests for ``padded()``"""

    def test_no_n(self):
        seq = [1, 2, 3]

        # No fillvalue
        self.assertEqual(mi.take(5, mi.padded(seq)), [1, 2, 3, None, None])

        # With fillvalue
        self.assertEqual(
            mi.take(5, mi.padded(seq, fillvalue='')), [1, 2, 3, '', '']
        )

    def test_invalid_n(self):
        self.assertRaises(ValueError, lambda: list(mi.padded([1, 2, 3], n=-1)))
        self.assertRaises(ValueError, lambda: list(mi.padded([1, 2, 3], n=0)))

    def test_valid_n(self):
        seq = [1, 2, 3, 4, 5]

        # No need for padding: len(seq) <= n
        self.assertEqual(list(mi.padded(seq, n=4)), [1, 2, 3, 4, 5])
        self.assertEqual(list(mi.padded(seq, n=5)), [1, 2, 3, 4, 5])

        # No fillvalue
        self.assertEqual(
            list(mi.padded(seq, n=7)), [1, 2, 3, 4, 5, None, None]
        )

        # With fillvalue
        self.assertEqual(
            list(mi.padded(seq, fillvalue='', n=7)), [1, 2, 3, 4, 5, '', '']
        )

    def test_next_multiple(self):
        seq = [1, 2, 3, 4, 5, 6]

        # No need for padding: len(seq) % n == 0
        self.assertEqual(
            list(mi.padded(seq, n=3, next_multiple=True)), [1, 2, 3, 4, 5, 6]
        )

        # Padding needed: len(seq) < n
        self.assertEqual(
            list(mi.padded(seq, n=8, next_multiple=True)),
            [1, 2, 3, 4, 5, 6, None, None]
        )

        # No padding needed: len(seq) == n
        self.assertEqual(
            list(mi.padded(seq, n=6, next_multiple=True)), [1, 2, 3, 4, 5, 6]
        )

        # Padding needed: len(seq) > n
        self.assertEqual(
            list(mi.padded(seq, n=4, next_multiple=True)),
            [1, 2, 3, 4, 5, 6, None, None]
        )

        # With fillvalue
        self.assertEqual(
            list(mi.padded(seq, fillvalue='', n=4, next_multiple=True)),
            [1, 2, 3, 4, 5, 6, '', '']
        )


class DistributeTest(TestCase):
    """Tests for distribute()"""

    def test_invalid_n(self):
        self.assertRaises(ValueError, lambda: mi.distribute(-1, [1, 2, 3]))
        self.assertRaises(ValueError, lambda: mi.distribute(0, [1, 2, 3]))

    def test_basic(self):
        iterable = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

        for n, expected in [
            (1, [iterable]),
            (2, [[1, 3, 5, 7, 9], [2, 4, 6, 8, 10]]),
            (3, [[1, 4, 7, 10], [2, 5, 8], [3, 6, 9]]),
            (10, [[n] for n in range(1, 10 + 1)]),
        ]:
            self.assertEqual(
                [list(x) for x in mi.distribute(n, iterable)], expected
            )

    def test_large_n(self):
        iterable = [1, 2, 3, 4]
        self.assertEqual(
            [list(x) for x in mi.distribute(6, iterable)],
            [[1], [2], [3], [4], [], []]
        )


class StaggerTest(TestCase):
    """Tests for ``stagger()``"""

    def test_default(self):
        iterable = [0, 1, 2, 3]
        actual = list(mi.stagger(iterable))
        expected = [(None, 0, 1), (0, 1, 2), (1, 2, 3)]
        self.assertEqual(actual, expected)

    def test_offsets(self):
        iterable = [0, 1, 2, 3]
        for offsets, expected in [
            ((-2, 0, 2), [('', 0, 2), ('', 1, 3)]),
            ((-2, -1), [('', ''), ('', 0), (0, 1), (1, 2), (2, 3)]),
            ((1, 2), [(1, 2), (2, 3)]),
        ]:
            all_groups = mi.stagger(iterable, offsets=offsets, fillvalue='')
            self.assertEqual(list(all_groups), expected)

    def test_longest(self):
        iterable = [0, 1, 2, 3]
        for offsets, expected in [
            (
                (-1, 0, 1),
                [('', 0, 1), (0, 1, 2), (1, 2, 3), (2, 3, ''), (3, '', '')]
            ),
            ((-2, -1), [('', ''), ('', 0), (0, 1), (1, 2), (2, 3), (3, '')]),
            ((1, 2), [(1, 2), (2, 3), (3, '')]),
        ]:
            all_groups = mi.stagger(
                iterable, offsets=offsets, fillvalue='', longest=True
            )
            self.assertEqual(list(all_groups), expected)


class ZipOffsetTest(TestCase):
    """Tests for ``zip_offset()``"""

    def test_shortest(self):
        a_1 = [0, 1, 2, 3]
        a_2 = [0, 1, 2, 3, 4, 5]
        a_3 = [0, 1, 2, 3, 4, 5, 6, 7]
        actual = list(
            mi.zip_offset(a_1, a_2, a_3, offsets=(-1, 0, 1), fillvalue='')
        )
        expected = [('', 0, 1), (0, 1, 2), (1, 2, 3), (2, 3, 4), (3, 4, 5)]
        self.assertEqual(actual, expected)

    def test_longest(self):
        a_1 = [0, 1, 2, 3]
        a_2 = [0, 1, 2, 3, 4, 5]
        a_3 = [0, 1, 2, 3, 4, 5, 6, 7]
        actual = list(
            mi.zip_offset(a_1, a_2, a_3, offsets=(-1, 0, 1), longest=True)
        )
        expected = [
            (None, 0, 1),
            (0, 1, 2),
            (1, 2, 3),
            (2, 3, 4),
            (3, 4, 5),
            (None, 5, 6),
            (None, None, 7),
        ]
        self.assertEqual(actual, expected)

    def test_mismatch(self):
        iterables = [0, 1, 2], [2, 3, 4]
        offsets = (-1, 0, 1)
        self.assertRaises(
            ValueError,
            lambda: list(mi.zip_offset(*iterables, offsets=offsets))
        )


class SortTogetherTest(TestCase):
    """Tests for sort_together()"""

    def test_key_list(self):
        """tests `key_list` including default, iterables include duplicates"""
        iterables = [
            ['GA', 'GA', 'GA', 'CT', 'CT', 'CT'],
            ['May', 'Aug.', 'May', 'June', 'July', 'July'],
            [97, 20, 100, 70, 100, 20]
        ]

        self.assertEqual(
            mi.sort_together(iterables),
            [
                ('CT', 'CT', 'CT', 'GA', 'GA', 'GA'),
                ('June', 'July', 'July', 'May', 'Aug.', 'May'),
                (70, 100, 20, 97, 20, 100)
            ]
        )

        self.assertEqual(
            mi.sort_together(iterables, key_list=(0, 1)),
            [
                ('CT', 'CT', 'CT', 'GA', 'GA', 'GA'),
                ('July', 'July', 'June', 'Aug.', 'May', 'May'),
                (100, 20, 70, 20, 97, 100)
            ]
        )

        self.assertEqual(
            mi.sort_together(iterables, key_list=(0, 1, 2)),
            [
                ('CT', 'CT', 'CT', 'GA', 'GA', 'GA'),
                ('July', 'July', 'June', 'Aug.', 'May', 'May'),
                (20, 100, 70, 20, 97, 100)
            ]
        )

        self.assertEqual(
            mi.sort_together(iterables, key_list=(2,)),
            [
                ('GA', 'CT', 'CT', 'GA', 'GA', 'CT'),
                ('Aug.', 'July', 'June', 'May', 'May', 'July'),
                (20, 20, 70, 97, 100, 100)
            ]
        )

    def test_invalid_key_list(self):
        """tests `key_list` for indexes not available in `iterables`"""
        iterables = [
            ['GA', 'GA', 'GA', 'CT', 'CT', 'CT'],
            ['May', 'Aug.', 'May', 'June', 'July', 'July'],
            [97, 20, 100, 70, 100, 20]
        ]

        self.assertRaises(
            IndexError, lambda: mi.sort_together(iterables, key_list=(5,))
        )

    def test_reverse(self):
        """tests `reverse` to ensure a reverse sort for `key_list` iterables"""
        iterables = [
            ['GA', 'GA', 'GA', 'CT', 'CT', 'CT'],
            ['May', 'Aug.', 'May', 'June', 'July', 'July'],
            [97, 20, 100, 70, 100, 20]
        ]

        self.assertEqual(
            mi.sort_together(iterables, key_list=(0, 1, 2), reverse=True),
            [('GA', 'GA', 'GA', 'CT', 'CT', 'CT'),
             ('May', 'May', 'Aug.', 'June', 'July', 'July'),
             (100, 97, 20, 70, 100, 20)]
        )

    def test_uneven_iterables(self):
        """tests trimming of iterables to the shortest length before sorting"""
        iterables = [['GA', 'GA', 'GA', 'CT', 'CT', 'CT', 'MA'],
                     ['May', 'Aug.', 'May', 'June', 'July', 'July'],
                     [97, 20, 100, 70, 100, 20, 0]]

        self.assertEqual(
            mi.sort_together(iterables),
            [
                ('CT', 'CT', 'CT', 'GA', 'GA', 'GA'),
                ('June', 'July', 'July', 'May', 'Aug.', 'May'),
                (70, 100, 20, 97, 20, 100)
            ]
        )


class DivideTest(TestCase):
    """Tests for divide()"""

    def test_invalid_n(self):
        self.assertRaises(ValueError, lambda: mi.divide(-1, [1, 2, 3]))
        self.assertRaises(ValueError, lambda: mi.divide(0, [1, 2, 3]))

    def test_basic(self):
        iterable = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

        for n, expected in [
            (1, [iterable]),
            (2, [[1, 2, 3, 4, 5], [6, 7, 8, 9, 10]]),
            (3, [[1, 2, 3, 4], [5, 6, 7], [8, 9, 10]]),
            (10, [[n] for n in range(1, 10 + 1)]),
        ]:
            self.assertEqual(
                [list(x) for x in mi.divide(n, iterable)], expected
            )

    def test_large_n(self):
        iterable = [1, 2, 3, 4]
        self.assertEqual(
            [list(x) for x in mi.divide(6, iterable)],
            [[1], [2], [3], [4], [], []]
        )


class TestAlwaysIterable(TestCase):
    """Tests for always_iterable()"""
    def test_single(self):
        self.assertEqual(list(mi.always_iterable(1)), [1])

    def test_strings(self):
        for obj in ['foo', b'bar', u'baz']:
            actual = list(mi.always_iterable(obj))
            expected = [obj]
            self.assertEqual(actual, expected)

    def test_base_type(self):
        dict_obj = {'a': 1, 'b': 2}
        str_obj = '123'

        # Default: dicts are iterable like they normally are
        default_actual = list(mi.always_iterable(dict_obj))
        default_expected = list(dict_obj)
        self.assertEqual(default_actual, default_expected)

        # Unitary types set: dicts are not iterable
        custom_actual = list(mi.always_iterable(dict_obj, base_type=dict))
        custom_expected = [dict_obj]
        self.assertEqual(custom_actual, custom_expected)

        # With unitary types set, strings are iterable
        str_actual = list(mi.always_iterable(str_obj, base_type=None))
        str_expected = list(str_obj)
        self.assertEqual(str_actual, str_expected)

    def test_iterables(self):
        self.assertEqual(list(mi.always_iterable([0, 1])), [0, 1])
        self.assertEqual(
            list(mi.always_iterable([0, 1], base_type=list)), [[0, 1]]
        )
        self.assertEqual(
            list(mi.always_iterable(iter('foo'))), ['f', 'o', 'o']
        )
        self.assertEqual(list(mi.always_iterable([])), [])

    def test_none(self):
        self.assertEqual(list(mi.always_iterable(None)), [])

    def test_generator(self):
        def _gen():
            yield 0
            yield 1

        self.assertEqual(list(mi.always_iterable(_gen())), [0, 1])


class AdjacentTests(TestCase):
    def test_typical(self):
        actual = list(mi.adjacent(lambda x: x % 5 == 0, range(10)))
        expected = [(True, 0), (True, 1), (False, 2), (False, 3), (True, 4),
                    (True, 5), (True, 6), (False, 7), (False, 8), (False, 9)]
        self.assertEqual(actual, expected)

    def test_empty_iterable(self):
        actual = list(mi.adjacent(lambda x: x % 5 == 0, []))
        expected = []
        self.assertEqual(actual, expected)

    def test_length_one(self):
        actual = list(mi.adjacent(lambda x: x % 5 == 0, [0]))
        expected = [(True, 0)]
        self.assertEqual(actual, expected)

        actual = list(mi.adjacent(lambda x: x % 5 == 0, [1]))
        expected = [(False, 1)]
        self.assertEqual(actual, expected)

    def test_consecutive_true(self):
        """Test that when the predicate matches multiple consecutive elements
        it doesn't repeat elements in the output"""
        actual = list(mi.adjacent(lambda x: x % 5 < 2, range(10)))
        expected = [(True, 0), (True, 1), (True, 2), (False, 3), (True, 4),
                    (True, 5), (True, 6), (True, 7), (False, 8), (False, 9)]
        self.assertEqual(actual, expected)

    def test_distance(self):
        actual = list(mi.adjacent(lambda x: x % 5 == 0, range(10), distance=2))
        expected = [(True, 0), (True, 1), (True, 2), (True, 3), (True, 4),
                    (True, 5), (True, 6), (True, 7), (False, 8), (False, 9)]
        self.assertEqual(actual, expected)

        actual = list(mi.adjacent(lambda x: x % 5 == 0, range(10), distance=3))
        expected = [(True, 0), (True, 1), (True, 2), (True, 3), (True, 4),
                    (True, 5), (True, 6), (True, 7), (True, 8), (False, 9)]
        self.assertEqual(actual, expected)

    def test_large_distance(self):
        """Test distance larger than the length of the iterable"""
        iterable = range(10)
        actual = list(mi.adjacent(lambda x: x % 5 == 4, iterable, distance=20))
        expected = list(zip(repeat(True), iterable))
        self.assertEqual(actual, expected)

        actual = list(mi.adjacent(lambda x: False, iterable, distance=20))
        expected = list(zip(repeat(False), iterable))
        self.assertEqual(actual, expected)

    def test_zero_distance(self):
        """Test that adjacent() reduces to zip+map when distance is 0"""
        iterable = range(1000)
        predicate = lambda x: x % 4 == 2
        actual = mi.adjacent(predicate, iterable, 0)
        expected = zip(map(predicate, iterable), iterable)
        self.assertTrue(all(a == e for a, e in zip(actual, expected)))

    def test_negative_distance(self):
        """Test that adjacent() raises an error with negative distance"""
        pred = lambda x: x
        self.assertRaises(
            ValueError, lambda: mi.adjacent(pred, range(1000), -1)
        )
        self.assertRaises(
            ValueError, lambda: mi.adjacent(pred, range(10), -10)
        )

    def test_grouping(self):
        """Test interaction of adjacent() with groupby_transform()"""
        iterable = mi.adjacent(lambda x: x % 5 == 0, range(10))
        grouper = mi.groupby_transform(iterable, itemgetter(0), itemgetter(1))
        actual = [(k, list(g)) for k, g in grouper]
        expected = [
            (True, [0, 1]),
            (False, [2, 3]),
            (True, [4, 5, 6]),
            (False, [7, 8, 9]),
        ]
        self.assertEqual(actual, expected)

    def test_call_once(self):
        """Test that the predicate is only called once per item."""
        already_seen = set()
        iterable = range(10)

        def predicate(item):
            self.assertNotIn(item, already_seen)
            already_seen.add(item)
            return True

        actual = list(mi.adjacent(predicate, iterable))
        expected = [(True, x) for x in iterable]
        self.assertEqual(actual, expected)


class GroupByTransformTests(TestCase):
    def assertAllGroupsEqual(self, groupby1, groupby2):
        """Compare two groupby objects for equality, both keys and groups."""
        for a, b in zip(groupby1, groupby2):
            key1, group1 = a
            key2, group2 = b
            self.assertEqual(key1, key2)
            self.assertListEqual(list(group1), list(group2))
        self.assertRaises(StopIteration, lambda: next(groupby1))
        self.assertRaises(StopIteration, lambda: next(groupby2))

    def test_default_funcs(self):
        """Test that groupby_transform() with default args mimics groupby()"""
        iterable = [(x // 5, x) for x in range(1000)]
        actual = mi.groupby_transform(iterable)
        expected = groupby(iterable)
        self.assertAllGroupsEqual(actual, expected)

    def test_valuefunc(self):
        iterable = [(int(x / 5), int(x / 3), x) for x in range(10)]

        # Test the standard usage of grouping one iterable using another's keys
        grouper = mi.groupby_transform(
            iterable, keyfunc=itemgetter(0), valuefunc=itemgetter(-1)
        )
        actual = [(k, list(g)) for k, g in grouper]
        expected = [(0, [0, 1, 2, 3, 4]), (1, [5, 6, 7, 8, 9])]
        self.assertEqual(actual, expected)

        grouper = mi.groupby_transform(
            iterable, keyfunc=itemgetter(1), valuefunc=itemgetter(-1)
        )
        actual = [(k, list(g)) for k, g in grouper]
        expected = [(0, [0, 1, 2]), (1, [3, 4, 5]), (2, [6, 7, 8]), (3, [9])]
        self.assertEqual(actual, expected)

        # and now for something a little different
        d = dict(zip(range(10), 'abcdefghij'))
        grouper = mi.groupby_transform(
            range(10), keyfunc=lambda x: x // 5, valuefunc=d.get
        )
        actual = [(k, ''.join(g)) for k, g in grouper]
        expected = [(0, 'abcde'), (1, 'fghij')]
        self.assertEqual(actual, expected)

    def test_no_valuefunc(self):
        iterable = range(1000)

        def key(x):
            return x // 5

        actual = mi.groupby_transform(iterable, key, valuefunc=None)
        expected = groupby(iterable, key)
        self.assertAllGroupsEqual(actual, expected)

        actual = mi.groupby_transform(iterable, key)  # default valuefunc
        expected = groupby(iterable, key)
        self.assertAllGroupsEqual(actual, expected)


class NumericRangeTests(TestCase):
    def test_basic(self):
        for args, expected in [
            ((4,), [0, 1, 2, 3]),
            ((4.0,), [0.0, 1.0, 2.0, 3.0]),
            ((1.0, 4), [1.0, 2.0, 3.0]),
            ((1, 4.0), [1, 2, 3]),
            ((1.0, 5), [1.0, 2.0, 3.0, 4.0]),
            ((0, 20, 5), [0, 5, 10, 15]),
            ((0, 20, 5.0), [0.0, 5.0, 10.0, 15.0]),
            ((0, 10, 3), [0, 3, 6, 9]),
            ((0, 10, 3.0), [0.0, 3.0, 6.0, 9.0]),
            ((0, -5, -1), [0, -1, -2, -3, -4]),
            ((0.0, -5, -1), [0.0, -1.0, -2.0, -3.0, -4.0]),
            ((1, 2, Fraction(1, 2)), [Fraction(1, 1), Fraction(3, 2)]),
            ((0,), []),
            ((0.0,), []),
            ((1, 0), []),
            ((1.0, 0.0), []),
            ((Fraction(2, 1),), [Fraction(0, 1), Fraction(1, 1)]),
            ((Decimal('2.0'),), [Decimal('0.0'), Decimal('1.0')]),
        ]:
            actual = list(mi.numeric_range(*args))
            self.assertEqual(actual, expected)
            self.assertTrue(
                all(type(a) == type(e) for a, e in zip(actual, expected))
            )

    def test_arg_count(self):
        self.assertRaises(TypeError, lambda: list(mi.numeric_range()))
        self.assertRaises(
            TypeError, lambda: list(mi.numeric_range(0, 1, 2, 3))
        )

    def test_zero_step(self):
        self.assertRaises(
            ValueError, lambda: list(mi.numeric_range(1, 2, 0))
        )


class CountCycleTests(TestCase):
    def test_basic(self):
        expected = [
            (0, 'a'), (0, 'b'), (0, 'c'),
            (1, 'a'), (1, 'b'), (1, 'c'),
            (2, 'a'), (2, 'b'), (2, 'c'),
        ]
        for actual in [
            mi.take(9, mi.count_cycle('abc')),  # n=None
            list(mi.count_cycle('abc', 3)),  # n=3
        ]:
            self.assertEqual(actual, expected)

    def test_empty(self):
        self.assertEqual(list(mi.count_cycle('')), [])
        self.assertEqual(list(mi.count_cycle('', 2)), [])

    def test_negative(self):
        self.assertEqual(list(mi.count_cycle('abc', -3)), [])


class LocateTests(TestCase):
    def test_default_pred(self):
        iterable = [0, 1, 1, 0, 1, 0, 0]
        actual = list(mi.locate(iterable))
        expected = [1, 2, 4]
        self.assertEqual(actual, expected)

    def test_no_matches(self):
        iterable = [0, 0, 0]
        actual = list(mi.locate(iterable))
        expected = []
        self.assertEqual(actual, expected)

    def test_custom_pred(self):
        iterable = ['0', 1, 1, '0', 1, '0', '0']
        pred = lambda x: x == '0'
        actual = list(mi.locate(iterable, pred))
        expected = [0, 3, 5, 6]
        self.assertEqual(actual, expected)


class StripFunctionTests(TestCase):
    def test_hashable(self):
        iterable = list('www.example.com')
        pred = lambda x: x in set('cmowz.')

        self.assertEqual(list(mi.lstrip(iterable, pred)), list('example.com'))
        self.assertEqual(list(mi.rstrip(iterable, pred)), list('www.example'))
        self.assertEqual(list(mi.strip(iterable, pred)), list('example'))

    def test_not_hashable(self):
        iterable = [
            list('http://'), list('www'), list('.example'), list('.com')
        ]
        pred = lambda x: x in [list('http://'), list('www'), list('.com')]

        self.assertEqual(list(mi.lstrip(iterable, pred)), iterable[2:])
        self.assertEqual(list(mi.rstrip(iterable, pred)), iterable[:3])
        self.assertEqual(list(mi.strip(iterable, pred)), iterable[2: 3])

    def test_math(self):
        iterable = [0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2]
        pred = lambda x: x <= 2

        self.assertEqual(list(mi.lstrip(iterable, pred)), iterable[3:])
        self.assertEqual(list(mi.rstrip(iterable, pred)), iterable[:-3])
        self.assertEqual(list(mi.strip(iterable, pred)), iterable[3:-3])


class IsliceExtendedTests(TestCase):
    def test_all(self):
        iterable = ['0', '1', '2', '3', '4', '5']
        indexes = list(range(-4, len(iterable) + 4)) + [None]
        steps = [1, 2, 3, 4, -1, -2, -3, 4]
        for slice_args in product(indexes, indexes, steps):
            try:
                actual = list(mi.islice_extended(iterable, *slice_args))
            except Exception as e:
                self.fail((slice_args, e))

            expected = iterable[slice(*slice_args)]
            self.assertEqual(actual, expected, slice_args)

    def test_zero_step(self):
        with self.assertRaises(ValueError):
            list(mi.islice_extended([1, 2, 3], 0, 1, 0))


class ConsecutiveGroupsTest(TestCase):
    def test_numbers(self):
        iterable = [-10, -8, -7, -6, 1, 2, 4, 5, -1, 7]
        actual = [list(g) for g in mi.consecutive_groups(iterable)]
        expected = [[-10], [-8, -7, -6], [1, 2], [4, 5], [-1], [7]]
        self.assertEqual(actual, expected)

    def test_custom_ordering(self):
        iterable = ['1', '10', '11', '20', '21', '22', '30', '31']
        ordering = lambda x: int(x)
        actual = [list(g) for g in mi.consecutive_groups(iterable, ordering)]
        expected = [['1'], ['10', '11'], ['20', '21', '22'], ['30', '31']]
        self.assertEqual(actual, expected)

    def test_exotic_ordering(self):
        iterable = [
            ('a', 'b', 'c', 'd'),
            ('a', 'c', 'b', 'd'),
            ('a', 'c', 'd', 'b'),
            ('a', 'd', 'b', 'c'),
            ('d', 'b', 'c', 'a'),
            ('d', 'c', 'a', 'b'),
        ]
        ordering = list(permutations('abcd')).index
        actual = [list(g) for g in mi.consecutive_groups(iterable, ordering)]
        expected = [
            [('a', 'b', 'c', 'd')],
            [('a', 'c', 'b', 'd'), ('a', 'c', 'd', 'b'), ('a', 'd', 'b', 'c')],
            [('d', 'b', 'c', 'a'), ('d', 'c', 'a', 'b')],
        ]
        self.assertEqual(actual, expected)


class DifferenceTest(TestCase):
    def test_normal(self):
        iterable = [10, 20, 30, 40, 50]
        actual = list(mi.difference(iterable))
        expected = [10, 10, 10, 10, 10]
        self.assertEqual(actual, expected)

    def test_custom(self):
        iterable = [10, 20, 30, 40, 50]
        actual = list(mi.difference(iterable, add))
        expected = [10, 30, 50, 70, 90]
        self.assertEqual(actual, expected)

    def test_roundtrip(self):
        original = list(range(100))
        accumulated = mi.accumulate(original)
        actual = list(mi.difference(accumulated))
        self.assertEqual(actual, original)

    def test_one(self):
        self.assertEqual(list(mi.difference([0])), [0])

    def test_empty(self):
        self.assertEqual(list(mi.difference([])), [])


class SeekableTest(TestCase):
    def test_exhaustion_reset(self):
        iterable = [str(n) for n in range(10)]

        s = mi.seekable(iterable)
        self.assertEqual(list(s), iterable)  # Normal iteration
        self.assertEqual(list(s), [])  # Iterable is exhausted

        s.seek(0)
        self.assertEqual(list(s), iterable)  # Back in action

    def test_partial_reset(self):
        iterable = [str(n) for n in range(10)]

        s = mi.seekable(iterable)
        self.assertEqual(mi.take(5, s), iterable[:5])  # Normal iteration

        s.seek(1)
        self.assertEqual(list(s), iterable[1:])  # Get the rest of the iterable

    def test_forward(self):
        iterable = [str(n) for n in range(10)]

        s = mi.seekable(iterable)
        self.assertEqual(mi.take(1, s), iterable[:1])  # Normal iteration

        s.seek(3)  # Skip over index 2
        self.assertEqual(list(s), iterable[3:])  # Result is similar to slicing

        s.seek(0)  # Back to 0
        self.assertEqual(list(s), iterable)  # No difference in result

    def test_past_end(self):
        iterable = [str(n) for n in range(10)]

        s = mi.seekable(iterable)
        self.assertEqual(mi.take(1, s), iterable[:1])  # Normal iteration

        s.seek(20)
        self.assertEqual(list(s), [])  # Iterable is exhausted

        s.seek(0)  # Back to 0
        self.assertEqual(list(s), iterable)  # No difference in result

    def test_elements(self):
        iterable = map(str, count())

        s = mi.seekable(iterable)
        mi.take(10, s)

        elements = s.elements()
        self.assertEqual(
            [elements[i] for i in range(10)], [str(n) for n in range(10)]
        )
        self.assertEqual(len(elements), 10)

        mi.take(10, s)
        self.assertEqual(list(elements), [str(n) for n in range(20)])


class SequenceViewTests(TestCase):
    def test_init(self):
        view = mi.SequenceView((1, 2, 3))
        self.assertEqual(repr(view), "SequenceView((1, 2, 3))")
        self.assertRaises(TypeError, lambda: mi.SequenceView({}))

    def test_update(self):
        seq = [1, 2, 3]
        view = mi.SequenceView(seq)
        self.assertEqual(len(view), 3)
        self.assertEqual(repr(view), "SequenceView([1, 2, 3])")

        seq.pop()
        self.assertEqual(len(view), 2)
        self.assertEqual(repr(view), "SequenceView([1, 2])")

    def test_indexing(self):
        seq = ('a', 'b', 'c', 'd', 'e', 'f')
        view = mi.SequenceView(seq)
        for i in range(-len(seq), len(seq)):
            self.assertEqual(view[i], seq[i])

    def test_slicing(self):
        seq = ('a', 'b', 'c', 'd', 'e', 'f')
        view = mi.SequenceView(seq)
        n = len(seq)
        indexes = list(range(-n - 1, n + 1)) + [None]
        steps = list(range(-n, n + 1))
        steps.remove(0)
        for slice_args in product(indexes, indexes, steps):
            i = slice(*slice_args)
            self.assertEqual(view[i], seq[i])

    def test_abc_methods(self):
        # collections.Sequence should provide all of this functionality
        seq = ('a', 'b', 'c', 'd', 'e', 'f', 'f')
        view = mi.SequenceView(seq)

        # __contains__
        self.assertIn('b', view)
        self.assertNotIn('g', view)

        # __iter__
        self.assertEqual(list(iter(view)), list(seq))

        # __reversed__
        self.assertEqual(list(reversed(view)), list(reversed(seq)))

        # index
        self.assertEqual(view.index('b'), 1)

        # count
        self.assertEqual(seq.count('f'), 2)


class RunLengthTest(TestCase):
    def test_encode(self):
        iterable = (int(str(n)[0]) for n in count(800))
        actual = mi.take(4, mi.run_length.encode(iterable))
        expected = [(8, 100), (9, 100), (1, 1000), (2, 1000)]
        self.assertEqual(actual, expected)

    def test_decode(self):
        iterable = [('d', 4), ('c', 3), ('b', 2), ('a', 1)]
        actual = ''.join(mi.run_length.decode(iterable))
        expected = 'ddddcccbba'
        self.assertEqual(actual, expected)


class ExactlyNTests(TestCase):
    """Tests for ``exactly_n()``"""

    def test_true(self):
        """Iterable has ``n`` ``True`` elements"""
        self.assertTrue(mi.exactly_n([True, False, True], 2))
        self.assertTrue(mi.exactly_n([1, 1, 1, 0], 3))
        self.assertTrue(mi.exactly_n([False, False], 0))
        self.assertTrue(mi.exactly_n(range(100), 10, lambda x: x < 10))

    def test_false(self):
        """Iterable does not have ``n`` ``True`` elements"""
        self.assertFalse(mi.exactly_n([True, False, False], 2))
        self.assertFalse(mi.exactly_n([True, True, False], 1))
        self.assertFalse(mi.exactly_n([False], 1))
        self.assertFalse(mi.exactly_n([True], -1))
        self.assertFalse(mi.exactly_n(repeat(True), 100))

    def test_empty(self):
        """Return ``True`` if the iterable is empty and ``n`` is 0"""
        self.assertTrue(mi.exactly_n([], 0))
        self.assertFalse(mi.exactly_n([], 1))


class AlwaysReversibleTests(TestCase):
    """Tests for ``always_reversible()``"""

    def test_regular_reversed(self):
        self.assertEqual(list(reversed(range(10))),
                         list(mi.always_reversible(range(10))))
        self.assertEqual(list(reversed([1, 2, 3])),
                         list(mi.always_reversible([1, 2, 3])))
        self.assertEqual(reversed([1, 2, 3]).__class__,
                         mi.always_reversible([1, 2, 3]).__class__)

    def test_nonseq_reversed(self):
        # Create a non-reversible generator from a sequence
        with self.assertRaises(TypeError):
            reversed(x for x in range(10))

        self.assertEqual(list(reversed(range(10))),
                         list(mi.always_reversible(x for x in range(10))))
        self.assertEqual(list(reversed([1, 2, 3])),
                         list(mi.always_reversible(x for x in [1, 2, 3])))
        self.assertNotEqual(reversed((1, 2)).__class__,
                            mi.always_reversible(x for x in (1, 2)).__class__)


class CircularShiftsTests(TestCase):
    def test_empty(self):
        # empty iterable -> empty list
        self.assertEqual(list(mi.circular_shifts([])), [])

    def test_simple_circular_shifts(self):
        # test the a simple iterator case
        self.assertEqual(
            mi.circular_shifts(range(4)),
            [(0, 1, 2, 3), (1, 2, 3, 0), (2, 3, 0, 1), (3, 0, 1, 2)]
        )

    def test_duplicates(self):
        # test non-distinct entries
        self.assertEqual(
            mi.circular_shifts([0, 1, 0, 1]),
            [(0, 1, 0, 1), (1, 0, 1, 0), (0, 1, 0, 1), (1, 0, 1, 0)]
        )


class MakeDecoratorTests(TestCase):
    def test_basic(self):
        slicer = mi.make_decorator(islice)

        @slicer(1, 10, 2)
        def user_function(arg_1, arg_2, kwarg_1=None):
            self.assertEqual(arg_1, 'arg_1')
            self.assertEqual(arg_2, 'arg_2')
            self.assertEqual(kwarg_1, 'kwarg_1')
            return map(str, count())

        it = user_function('arg_1', 'arg_2', kwarg_1='kwarg_1')
        actual = list(it)
        expected = ['1', '3', '5', '7', '9']
        self.assertEqual(actual, expected)

    def test_result_index(self):
        def stringify(*args, **kwargs):
            self.assertEqual(args[0], 'arg_0')
            iterable = args[1]
            self.assertEqual(args[2], 'arg_2')
            self.assertEqual(kwargs['kwarg_1'], 'kwarg_1')
            return map(str, iterable)

        stringifier = mi.make_decorator(stringify, result_index=1)

        @stringifier('arg_0', 'arg_2', kwarg_1='kwarg_1')
        def user_function(n):
            return count(n)

        it = user_function(1)
        actual = mi.take(5, it)
        expected = ['1', '2', '3', '4', '5']
        self.assertEqual(actual, expected)

    def test_wrap_class(self):
        seeker = mi.make_decorator(mi.seekable)

        @seeker()
        def user_function(n):
            return map(str, range(n))

        it = user_function(5)
        self.assertEqual(list(it), ['0', '1', '2', '3', '4'])

        it.seek(0)
        self.assertEqual(list(it), ['0', '1', '2', '3', '4'])
