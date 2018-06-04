from doctest import DocTestSuite
from unittest import TestCase

from itertools import combinations
from six.moves import range

import more_itertools as mi


def load_tests(loader, tests, ignore):
    # Add the doctests
    tests.addTests(DocTestSuite('more_itertools.recipes'))
    return tests


class AccumulateTests(TestCase):
    """Tests for ``accumulate()``"""

    def test_empty(self):
        """Test that an empty input returns an empty output"""
        self.assertEqual(list(mi.accumulate([])), [])

    def test_default(self):
        """Test accumulate with the default function (addition)"""
        self.assertEqual(list(mi.accumulate([1, 2, 3])), [1, 3, 6])

    def test_bogus_function(self):
        """Test accumulate with an invalid function"""
        with self.assertRaises(TypeError):
            list(mi.accumulate([1, 2, 3], func=lambda x: x))

    def test_custom_function(self):
        """Test accumulate with a custom function"""
        self.assertEqual(
            list(mi.accumulate((1, 2, 3, 2, 1), func=max)), [1, 2, 3, 3, 3]
        )


class TakeTests(TestCase):
    """Tests for ``take()``"""

    def test_simple_take(self):
        """Test basic usage"""
        t = mi.take(5, range(10))
        self.assertEqual(t, [0, 1, 2, 3, 4])

    def test_null_take(self):
        """Check the null case"""
        t = mi.take(0, range(10))
        self.assertEqual(t, [])

    def test_negative_take(self):
        """Make sure taking negative items results in a ValueError"""
        self.assertRaises(ValueError, lambda: mi.take(-3, range(10)))

    def test_take_too_much(self):
        """Taking more than an iterator has remaining should return what the
        iterator has remaining.

        """
        t = mi.take(10, range(5))
        self.assertEqual(t, [0, 1, 2, 3, 4])


class TabulateTests(TestCase):
    """Tests for ``tabulate()``"""

    def test_simple_tabulate(self):
        """Test the happy path"""
        t = mi.tabulate(lambda x: x)
        f = tuple([next(t) for _ in range(3)])
        self.assertEqual(f, (0, 1, 2))

    def test_count(self):
        """Ensure tabulate accepts specific count"""
        t = mi.tabulate(lambda x: 2 * x, -1)
        f = (next(t), next(t), next(t))
        self.assertEqual(f, (-2, 0, 2))


class TailTests(TestCase):
    """Tests for ``tail()``"""

    def test_greater(self):
        """Length of iterable is greather than requested tail"""
        self.assertEqual(list(mi.tail(3, 'ABCDEFG')), ['E', 'F', 'G'])

    def test_equal(self):
        """Length of iterable is equal to the requested tail"""
        self.assertEqual(
            list(mi.tail(7, 'ABCDEFG')), ['A', 'B', 'C', 'D', 'E', 'F', 'G']
        )

    def test_less(self):
        """Length of iterable is less than requested tail"""
        self.assertEqual(
            list(mi.tail(8, 'ABCDEFG')), ['A', 'B', 'C', 'D', 'E', 'F', 'G']
        )


class ConsumeTests(TestCase):
    """Tests for ``consume()``"""

    def test_sanity(self):
        """Test basic functionality"""
        r = (x for x in range(10))
        mi.consume(r, 3)
        self.assertEqual(3, next(r))

    def test_null_consume(self):
        """Check the null case"""
        r = (x for x in range(10))
        mi.consume(r, 0)
        self.assertEqual(0, next(r))

    def test_negative_consume(self):
        """Check that negative consumsion throws an error"""
        r = (x for x in range(10))
        self.assertRaises(ValueError, lambda: mi.consume(r, -1))

    def test_total_consume(self):
        """Check that iterator is totally consumed by default"""
        r = (x for x in range(10))
        mi.consume(r)
        self.assertRaises(StopIteration, lambda: next(r))


class NthTests(TestCase):
    """Tests for ``nth()``"""

    def test_basic(self):
        """Make sure the nth item is returned"""
        l = range(10)
        for i, v in enumerate(l):
            self.assertEqual(mi.nth(l, i), v)

    def test_default(self):
        """Ensure a default value is returned when nth item not found"""
        l = range(3)
        self.assertEqual(mi.nth(l, 100, "zebra"), "zebra")

    def test_negative_item_raises(self):
        """Ensure asking for a negative item raises an exception"""
        self.assertRaises(ValueError, lambda: mi.nth(range(10), -3))


class AllEqualTests(TestCase):
    """Tests for ``all_equal()``"""

    def test_true(self):
        """Everything is equal"""
        self.assertTrue(mi.all_equal('aaaaaa'))
        self.assertTrue(mi.all_equal([0, 0, 0, 0]))

    def test_false(self):
        """Not everything is equal"""
        self.assertFalse(mi.all_equal('aaaaab'))
        self.assertFalse(mi.all_equal([0, 0, 0, 1]))

    def test_tricky(self):
        """Not everything is identical, but everything is equal"""
        items = [1, complex(1, 0), 1.0]
        self.assertTrue(mi.all_equal(items))

    def test_empty(self):
        """Return True if the iterable is empty"""
        self.assertTrue(mi.all_equal(''))
        self.assertTrue(mi.all_equal([]))

    def test_one(self):
        """Return True if the iterable is singular"""
        self.assertTrue(mi.all_equal('0'))
        self.assertTrue(mi.all_equal([0]))


class QuantifyTests(TestCase):
    """Tests for ``quantify()``"""

    def test_happy_path(self):
        """Make sure True count is returned"""
        q = [True, False, True]
        self.assertEqual(mi.quantify(q), 2)

    def test_custom_predicate(self):
        """Ensure non-default predicates return as expected"""
        q = range(10)
        self.assertEqual(mi.quantify(q, lambda x: x % 2 == 0), 5)


class PadnoneTests(TestCase):
    """Tests for ``padnone()``"""

    def test_happy_path(self):
        """wrapper iterator should return None indefinitely"""
        r = range(2)
        p = mi.padnone(r)
        self.assertEqual([0, 1, None, None], [next(p) for _ in range(4)])


class NcyclesTests(TestCase):
    """Tests for ``nyclces()``"""

    def test_happy_path(self):
        """cycle a sequence three times"""
        r = ["a", "b", "c"]
        n = mi.ncycles(r, 3)
        self.assertEqual(
            ["a", "b", "c", "a", "b", "c", "a", "b", "c"],
            list(n)
        )

    def test_null_case(self):
        """asking for 0 cycles should return an empty iterator"""
        n = mi.ncycles(range(100), 0)
        self.assertRaises(StopIteration, lambda: next(n))

    def test_pathalogical_case(self):
        """asking for negative cycles should return an empty iterator"""
        n = mi.ncycles(range(100), -10)
        self.assertRaises(StopIteration, lambda: next(n))


class DotproductTests(TestCase):
    """Tests for ``dotproduct()``'"""

    def test_happy_path(self):
        """simple dotproduct example"""
        self.assertEqual(400, mi.dotproduct([10, 10], [20, 20]))


class FlattenTests(TestCase):
    """Tests for ``flatten()``"""

    def test_basic_usage(self):
        """ensure list of lists is flattened one level"""
        f = [[0, 1, 2], [3, 4, 5]]
        self.assertEqual(list(range(6)), list(mi.flatten(f)))

    def test_single_level(self):
        """ensure list of lists is flattened only one level"""
        f = [[0, [1, 2]], [[3, 4], 5]]
        self.assertEqual([0, [1, 2], [3, 4], 5], list(mi.flatten(f)))


class RepeatfuncTests(TestCase):
    """Tests for ``repeatfunc()``"""

    def test_simple_repeat(self):
        """test simple repeated functions"""
        r = mi.repeatfunc(lambda: 5)
        self.assertEqual([5, 5, 5, 5, 5], [next(r) for _ in range(5)])

    def test_finite_repeat(self):
        """ensure limited repeat when times is provided"""
        r = mi.repeatfunc(lambda: 5, times=5)
        self.assertEqual([5, 5, 5, 5, 5], list(r))

    def test_added_arguments(self):
        """ensure arguments are applied to the function"""
        r = mi.repeatfunc(lambda x: x, 2, 3)
        self.assertEqual([3, 3], list(r))

    def test_null_times(self):
        """repeat 0 should return an empty iterator"""
        r = mi.repeatfunc(range, 0, 3)
        self.assertRaises(StopIteration, lambda: next(r))


class PairwiseTests(TestCase):
    """Tests for ``pairwise()``"""

    def test_base_case(self):
        """ensure an iterable will return pairwise"""
        p = mi.pairwise([1, 2, 3])
        self.assertEqual([(1, 2), (2, 3)], list(p))

    def test_short_case(self):
        """ensure an empty iterator if there's not enough values to pair"""
        p = mi.pairwise("a")
        self.assertRaises(StopIteration, lambda: next(p))


class GrouperTests(TestCase):
    """Tests for ``grouper()``"""

    def test_even(self):
        """Test when group size divides evenly into the length of
        the iterable.

        """
        self.assertEqual(
            list(mi.grouper(3, 'ABCDEF')), [('A', 'B', 'C'), ('D', 'E', 'F')]
        )

    def test_odd(self):
        """Test when group size does not divide evenly into the length of the
        iterable.

        """
        self.assertEqual(
            list(mi.grouper(3, 'ABCDE')), [('A', 'B', 'C'), ('D', 'E', None)]
        )

    def test_fill_value(self):
        """Test that the fill value is used to pad the final group"""
        self.assertEqual(
            list(mi.grouper(3, 'ABCDE', 'x')),
            [('A', 'B', 'C'), ('D', 'E', 'x')]
        )


class RoundrobinTests(TestCase):
    """Tests for ``roundrobin()``"""

    def test_even_groups(self):
        """Ensure ordered output from evenly populated iterables"""
        self.assertEqual(
            list(mi.roundrobin('ABC', [1, 2, 3], range(3))),
            ['A', 1, 0, 'B', 2, 1, 'C', 3, 2]
        )

    def test_uneven_groups(self):
        """Ensure ordered output from unevenly populated iterables"""
        self.assertEqual(
            list(mi.roundrobin('ABCD', [1, 2], range(0))),
            ['A', 1, 'B', 2, 'C', 'D']
        )


class PartitionTests(TestCase):
    """Tests for ``partition()``"""

    def test_bool(self):
        """Test when pred() returns a boolean"""
        lesser, greater = mi.partition(lambda x: x > 5, range(10))
        self.assertEqual(list(lesser), [0, 1, 2, 3, 4, 5])
        self.assertEqual(list(greater), [6, 7, 8, 9])

    def test_arbitrary(self):
        """Test when pred() returns an integer"""
        divisibles, remainders = mi.partition(lambda x: x % 3, range(10))
        self.assertEqual(list(divisibles), [0, 3, 6, 9])
        self.assertEqual(list(remainders), [1, 2, 4, 5, 7, 8])


class PowersetTests(TestCase):
    """Tests for ``powerset()``"""

    def test_combinatorics(self):
        """Ensure a proper enumeration"""
        p = mi.powerset([1, 2, 3])
        self.assertEqual(
            list(p),
            [(), (1,), (2,), (3,), (1, 2), (1, 3), (2, 3), (1, 2, 3)]
        )


class UniqueEverseenTests(TestCase):
    """Tests for ``unique_everseen()``"""

    def test_everseen(self):
        """ensure duplicate elements are ignored"""
        u = mi.unique_everseen('AAAABBBBCCDAABBB')
        self.assertEqual(
            ['A', 'B', 'C', 'D'],
            list(u)
        )

    def test_custom_key(self):
        """ensure the custom key comparison works"""
        u = mi.unique_everseen('aAbACCc', key=str.lower)
        self.assertEqual(list('abC'), list(u))

    def test_unhashable(self):
        """ensure things work for unhashable items"""
        iterable = ['a', [1, 2, 3], [1, 2, 3], 'a']
        u = mi.unique_everseen(iterable)
        self.assertEqual(list(u), ['a', [1, 2, 3]])

    def test_unhashable_key(self):
        """ensure things work for unhashable items with a custom key"""
        iterable = ['a', [1, 2, 3], [1, 2, 3], 'a']
        u = mi.unique_everseen(iterable, key=lambda x: x)
        self.assertEqual(list(u), ['a', [1, 2, 3]])


class UniqueJustseenTests(TestCase):
    """Tests for ``unique_justseen()``"""

    def test_justseen(self):
        """ensure only last item is remembered"""
        u = mi.unique_justseen('AAAABBBCCDABB')
        self.assertEqual(list('ABCDAB'), list(u))

    def test_custom_key(self):
        """ensure the custom key comparison works"""
        u = mi.unique_justseen('AABCcAD', str.lower)
        self.assertEqual(list('ABCAD'), list(u))


class IterExceptTests(TestCase):
    """Tests for ``iter_except()``"""

    def test_exact_exception(self):
        """ensure the exact specified exception is caught"""
        l = [1, 2, 3]
        i = mi.iter_except(l.pop, IndexError)
        self.assertEqual(list(i), [3, 2, 1])

    def test_generic_exception(self):
        """ensure the generic exception can be caught"""
        l = [1, 2]
        i = mi.iter_except(l.pop, Exception)
        self.assertEqual(list(i), [2, 1])

    def test_uncaught_exception_is_raised(self):
        """ensure a non-specified exception is raised"""
        l = [1, 2, 3]
        i = mi.iter_except(l.pop, KeyError)
        self.assertRaises(IndexError, lambda: list(i))

    def test_first(self):
        """ensure first is run before the function"""
        l = [1, 2, 3]
        f = lambda: 25
        i = mi.iter_except(l.pop, IndexError, f)
        self.assertEqual(list(i), [25, 3, 2, 1])


class FirstTrueTests(TestCase):
    """Tests for ``first_true()``"""

    def test_something_true(self):
        """Test with no keywords"""
        self.assertEqual(mi.first_true(range(10)), 1)

    def test_nothing_true(self):
        """Test default return value."""
        self.assertEqual(mi.first_true([0, 0, 0]), False)

    def test_default(self):
        """Test with a default keyword"""
        self.assertEqual(mi.first_true([0, 0, 0], default='!'), '!')

    def test_pred(self):
        """Test with a custom predicate"""
        self.assertEqual(
            mi.first_true([2, 4, 6], pred=lambda x: x % 3 == 0), 6
        )


class RandomProductTests(TestCase):
    """Tests for ``random_product()``

    Since random.choice() has different results with the same seed across
    python versions 2.x and 3.x, these tests use highly probably events to
    create predictable outcomes across platforms.
    """

    def test_simple_lists(self):
        """Ensure that one item is chosen from each list in each pair.
        Also ensure that each item from each list eventually appears in
        the chosen combinations.

        Odds are roughly 1 in 7.1 * 10e16 that one item from either list will
        not be chosen after 100 samplings of one item from each list. Just to
        be safe, better use a known random seed, too.

        """
        nums = [1, 2, 3]
        lets = ['a', 'b', 'c']
        n, m = zip(*[mi.random_product(nums, lets) for _ in range(100)])
        n, m = set(n), set(m)
        self.assertEqual(n, set(nums))
        self.assertEqual(m, set(lets))
        self.assertEqual(len(n), len(nums))
        self.assertEqual(len(m), len(lets))

    def test_list_with_repeat(self):
        """ensure multiple items are chosen, and that they appear to be chosen
        from one list then the next, in proper order.

        """
        nums = [1, 2, 3]
        lets = ['a', 'b', 'c']
        r = list(mi.random_product(nums, lets, repeat=100))
        self.assertEqual(2 * 100, len(r))
        n, m = set(r[::2]), set(r[1::2])
        self.assertEqual(n, set(nums))
        self.assertEqual(m, set(lets))
        self.assertEqual(len(n), len(nums))
        self.assertEqual(len(m), len(lets))


class RandomPermutationTests(TestCase):
    """Tests for ``random_permutation()``"""

    def test_full_permutation(self):
        """ensure every item from the iterable is returned in a new ordering

        15 elements have a 1 in 1.3 * 10e12 of appearing in sorted order, so
        we fix a seed value just to be sure.

        """
        i = range(15)
        r = mi.random_permutation(i)
        self.assertEqual(set(i), set(r))
        if i == r:
            raise AssertionError("Values were not permuted")

    def test_partial_permutation(self):
        """ensure all returned items are from the iterable, that the returned
        permutation is of the desired length, and that all items eventually
        get returned.

        Sampling 100 permutations of length 5 from a set of 15 leaves a
        (2/3)^100 chance that an item will not be chosen. Multiplied by 15
        items, there is a 1 in 2.6e16 chance that at least 1 item will not
        show up in the resulting output. Using a random seed will fix that.

        """
        items = range(15)
        item_set = set(items)
        all_items = set()
        for _ in range(100):
            permutation = mi.random_permutation(items, 5)
            self.assertEqual(len(permutation), 5)
            permutation_set = set(permutation)
            self.assertLessEqual(permutation_set, item_set)
            all_items |= permutation_set
        self.assertEqual(all_items, item_set)


class RandomCombinationTests(TestCase):
    """Tests for ``random_combination()``"""

    def test_psuedorandomness(self):
        """ensure different subsets of the iterable get returned over many
        samplings of random combinations"""
        items = range(15)
        all_items = set()
        for _ in range(50):
            combination = mi.random_combination(items, 5)
            all_items |= set(combination)
        self.assertEqual(all_items, set(items))

    def test_no_replacement(self):
        """ensure that elements are sampled without replacement"""
        items = range(15)
        for _ in range(50):
            combination = mi.random_combination(items, len(items))
            self.assertEqual(len(combination), len(set(combination)))
        self.assertRaises(
            ValueError, lambda: mi.random_combination(items, len(items) + 1)
        )


class RandomCombinationWithReplacementTests(TestCase):
    """Tests for ``random_combination_with_replacement()``"""

    def test_replacement(self):
        """ensure that elements are sampled with replacement"""
        items = range(5)
        combo = mi.random_combination_with_replacement(items, len(items) * 2)
        self.assertEqual(2 * len(items), len(combo))
        if len(set(combo)) == len(combo):
            raise AssertionError("Combination contained no duplicates")

    def test_pseudorandomness(self):
        """ensure different subsets of the iterable get returned over many
        samplings of random combinations"""
        items = range(15)
        all_items = set()
        for _ in range(50):
            combination = mi.random_combination_with_replacement(items, 5)
            all_items |= set(combination)
        self.assertEqual(all_items, set(items))


class NthCombinationTests(TestCase):
    def test_basic(self):
        iterable = 'abcdefg'
        r = 4
        for index, expected in enumerate(combinations(iterable, r)):
            actual = mi.nth_combination(iterable, r, index)
            self.assertEqual(actual, expected)

    def test_long(self):
        actual = mi.nth_combination(range(180), 4, 2000000)
        expected = (2, 12, 35, 126)
        self.assertEqual(actual, expected)
