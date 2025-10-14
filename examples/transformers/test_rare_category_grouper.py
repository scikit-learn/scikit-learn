
import pandas as pd
import pytest

from .rare_category_grouper import RareCategoryGrouper


def test_basic_min_freq():
    df = pd.DataFrame({"c": ["a", "a", "b", "b", "c", "d"]})
    rcg = RareCategoryGrouper(min_freq=2, columns=["c"], other_label="Other").fit(df)
    out = rcg.transform(df)
    assert out["c"].tolist() == ["a", "a", "b", "b", "Other", "Other"]
    assert set(rcg.get_kept_categories()["c"]) == {"a", "b"}


def test_min_prop_or_min_freq_logic():

    df = pd.DataFrame({"c": ["x", "y", "y", "z", "z", "z"]})

    rcg = RareCategoryGrouper(min_freq=3, min_prop=0.4, columns=["c"]).fit(df)
    kept = set(rcg.get_kept_categories()["c"])
    assert kept == {"z"}  

    out = rcg.transform(df)
    assert out["c"].tolist() == ["Other", "Other", "Other", "z", "z", "z"]


def test_unseen_categories_mapped_to_other():
    train = pd.DataFrame({"city": ["LA", "LA", "SF", "SF", "NY"]})
    test = pd.DataFrame({"city": ["LA", "SEA", "PDX", "SF"]})
    rcg = RareCategoryGrouper(min_freq=2, columns=["city"]).fit(train)
    out = rcg.transform(test)

    assert out["city"].tolist() == ["LA", "Other", "Other", "SF"]


def test_numeric_columns_are_untouched():
    df = pd.DataFrame({"city": ["A", "A", "B", "C"], "age": [20, 22, 21, 23]})
    rcg = RareCategoryGrouper(min_freq=2, columns=["city"]).fit(df)
    out = rcg.transform(df)
    assert out["age"].equals(df["age"])  


def test_default_column_selection_object_and_categorical():
    df = pd.DataFrame(
        {
            "city": pd.Series(["A", "A", "B", "C"], dtype="category"),
            "brand": ["x", "y", "x", "z"],
            "age": [1, 2, 3, 4], 
        }
    )
    rcg = RareCategoryGrouper(min_freq=2)  
    rcg.fit(df)
    out = rcg.transform(df)

    assert out["city"].tolist() == ["A", "A", "Other", "Other"]
    assert out["brand"].tolist() == ["x", "Other", "x", "Other"]
    assert out["age"].equals(df["age"])


def test_custom_other_label():
    df = pd.DataFrame({"c": ["a", "a", "b", "c"]})
    rcg = RareCategoryGrouper(min_freq=2, columns=["c"], other_label="<RARE>").fit(df)
    out = rcg.transform(df)
    assert out["c"].tolist() == ["a", "a", "<RARE>", "<RARE>"]


def test_input_validation_dataframe_only():
    rcg = RareCategoryGrouper(min_freq=1)
    with pytest.raises(TypeError):
        rcg.fit([["a"], ["b"]])  

    with pytest.raises(TypeError):
        rcg.transform([["a"], ["b"]]) 


def test_not_fitted_error():
    df = pd.DataFrame({"c": ["a", "b"]})
    rcg = RareCategoryGrouper(min_freq=1)
    with pytest.raises(AttributeError):
        rcg.transform(df)  
