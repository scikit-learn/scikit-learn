import pytest
from numpy.testing import assert_allclose

from sklearn.association._apriori import Apriori


def _as_sorted_tuple(s):
    return tuple(sorted(s))


def test_mines_frequent_itemsets_from_transactions():
    X = [
        ["a", "b"],
        ["b", "c"],
        ["a", "b", "c"],
        ["a"],
    ]
    apr = Apriori(min_support=0.5)
    apr.fit(X)

    # expected supports (fractions over 4 transactions)
    expected = {
        _as_sorted_tuple(["a"]): 3 / 4,
        _as_sorted_tuple(["b"]): 3 / 4,
        _as_sorted_tuple(["c"]): 2 / 4,
        _as_sorted_tuple(["a", "b"]): 2 / 4,
        _as_sorted_tuple(["b", "c"]): 2 / 4,
    }

    # convert keys to sorted tuples for comparison
    supports = {tuple(sorted(k)): v for k, v in apr.frequent_itemsets_.items()}

    for key, val in expected.items():
        assert key in supports
        assert_allclose(supports[key], val)


def test_generate_rules_confidence_lift_leverage():
    X = [
        ["a", "b"],
        ["b", "c"],
        ["a", "b", "c"],
        ["a"],
    ]
    apr = Apriori(min_support=0.5)
    apr.fit(X)

    rules = apr.generate_rules(metric="confidence", min_threshold=0.7)

    # normalize rules to list of dicts whether pandas present or not
    try:
        import pandas as pd

        if isinstance(rules, pd.DataFrame):
            records = rules.to_dict("records")
        else:
            records = rules
    except Exception:
        records = rules

    # find rule antecedent ('c',) -> consequent ('b',)
    match = None
    for r in records:
        if tuple(sorted(r["antecedent"])) == ("c",) and tuple(
            sorted(r["consequent"])
        ) == ("b",):
            match = r
            break

    assert match is not None, "expected rule (c) -> (b) not found"

    # check metrics: support itemset {'b','c'} = 0.5, antecedent 'c' support 0.5
    assert_allclose(match["support"], 0.5, rtol=1e-12)
    assert_allclose(match["confidence"], 1.0, rtol=1e-12)
    # lift = confidence / support(consequent) -> 1.0 / 0.75 = 4/3
    assert_allclose(match["lift"], pytest.approx(4 / 3), rtol=1e-12)
    # leverage = sup({b,c}) - sup(b)*sup(c) = 0.5 - 0.75*0.5 = 0.125
    assert_allclose(match["leverage"], pytest.approx(0.125), rtol=1e-12)


def test_onehot_dataframe_input_uses_column_names():
    pd = pytest.importorskip("pandas")

    df = pd.DataFrame(
        [
            {"a": 1, "b": 1, "c": 0},
            {"a": 0, "b": 1, "c": 1},
            {"a": 1, "b": 1, "c": 1},
            {"a": 1, "b": 0, "c": 0},
        ]
    )

    apr = Apriori(min_support=0.5, use_colnames=True)
    apr.fit(df)

    supports = {tuple(sorted(k)): v for k, v in apr.frequent_itemsets_.items()}

    # column names should appear as items
    assert ("a",) in supports
    assert ("b",) in supports
    assert ("c",) in supports
    # check support for 'a' is 3/4
    assert_allclose(supports[("a",)], 3 / 4)


def test_fit_empty_transactions_raises():
    apr = Apriori(min_support=0.1)
    with pytest.raises(ValueError, match="No transactions provided"):
        apr.fit([])
