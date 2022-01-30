import numpy as np
import pandas as pd
import datetime
import matplotlib as mpl
import matplotlib.pyplot as plt

import pytest


def has_verdana():
    """Helper to verify if Verdana font is present"""
    # This import is relatively lengthy, so to prevent its import for
    # testing other tests in this module not requiring this knowledge,
    # import font_manager here
    import matplotlib.font_manager as mplfm
    try:
        verdana_font = mplfm.findfont('Verdana', fallback_to_default=False)
    except:  # noqa
        # if https://github.com/matplotlib/matplotlib/pull/3435
        # gets accepted
        return False
    # otherwise check if not matching the logic for a 'default' one
    try:
        unlikely_font = mplfm.findfont("very_unlikely_to_exist1234",
                                       fallback_to_default=False)
    except:  # noqa
        # if matched verdana but not unlikely, Verdana must exist
        return True
    # otherwise -- if they match, must be the same default
    return verdana_font != unlikely_font


@pytest.fixture(scope="session", autouse=True)
def remove_pandas_unit_conversion():
    # Prior to pandas 1.0, it registered its own datetime converters,
    # but they are less powerful than what matplotlib added in 2.2,
    # and we rely on that functionality in seaborn.
    # https://github.com/matplotlib/matplotlib/pull/9779
    # https://github.com/pandas-dev/pandas/issues/27036
    mpl.units.registry[np.datetime64] = mpl.dates.DateConverter()
    mpl.units.registry[datetime.date] = mpl.dates.DateConverter()
    mpl.units.registry[datetime.datetime] = mpl.dates.DateConverter()


@pytest.fixture(autouse=True)
def close_figs():
    yield
    plt.close("all")


@pytest.fixture(autouse=True)
def random_seed():
    seed = sum(map(ord, "seaborn random global"))
    np.random.seed(seed)


@pytest.fixture()
def rng():
    seed = sum(map(ord, "seaborn random object"))
    return np.random.RandomState(seed)


@pytest.fixture
def wide_df(rng):

    columns = list("abc")
    index = pd.Int64Index(np.arange(10, 50, 2), name="wide_index")
    values = rng.normal(size=(len(index), len(columns)))
    return pd.DataFrame(values, index=index, columns=columns)


@pytest.fixture
def wide_array(wide_df):

    # Requires panads >= 0.24
    # return wide_df.to_numpy()
    return np.asarray(wide_df)


@pytest.fixture
def flat_series(rng):

    index = pd.Int64Index(np.arange(10, 30), name="t")
    return pd.Series(rng.normal(size=20), index, name="s")


@pytest.fixture
def flat_array(flat_series):

    # Requires panads >= 0.24
    # return flat_series.to_numpy()
    return np.asarray(flat_series)


@pytest.fixture
def flat_list(flat_series):

    # Requires panads >= 0.24
    # return flat_series.to_list()
    return flat_series.tolist()


@pytest.fixture(params=["series", "array", "list"])
def flat_data(rng, request):

    index = pd.Int64Index(np.arange(10, 30), name="t")
    series = pd.Series(rng.normal(size=20), index, name="s")
    if request.param == "series":
        data = series
    elif request.param == "array":
        try:
            data = series.to_numpy()  # Requires pandas >= 0.24
        except AttributeError:
            data = np.asarray(series)
    elif request.param == "list":
        try:
            data = series.to_list()  # Requires pandas >= 0.24
        except AttributeError:
            data = series.tolist()
    return data


@pytest.fixture
def wide_list_of_series(rng):

    return [pd.Series(rng.normal(size=20), np.arange(20), name="a"),
            pd.Series(rng.normal(size=10), np.arange(5, 15), name="b")]


@pytest.fixture
def wide_list_of_arrays(wide_list_of_series):

    # Requires pandas >= 0.24
    # return [s.to_numpy() for s in wide_list_of_series]
    return [np.asarray(s) for s in wide_list_of_series]


@pytest.fixture
def wide_list_of_lists(wide_list_of_series):

    # Requires pandas >= 0.24
    # return [s.to_list() for s in wide_list_of_series]
    return [s.tolist() for s in wide_list_of_series]


@pytest.fixture
def wide_dict_of_series(wide_list_of_series):

    return {s.name: s for s in wide_list_of_series}


@pytest.fixture
def wide_dict_of_arrays(wide_list_of_series):

    # Requires pandas >= 0.24
    # return {s.name: s.to_numpy() for s in wide_list_of_series}
    return {s.name: np.asarray(s) for s in wide_list_of_series}


@pytest.fixture
def wide_dict_of_lists(wide_list_of_series):

    # Requires pandas >= 0.24
    # return {s.name: s.to_list() for s in wide_list_of_series}
    return {s.name: s.tolist() for s in wide_list_of_series}


@pytest.fixture
def long_df(rng):

    n = 100
    df = pd.DataFrame(dict(
        x=rng.uniform(0, 20, n).round().astype("int"),
        y=rng.normal(size=n),
        z=rng.lognormal(size=n),
        a=rng.choice(list("abc"), n),
        b=rng.choice(list("mnop"), n),
        c=rng.choice([0, 1], n, [.3, .7]),
        t=rng.choice(np.arange("2004-07-30", "2007-07-30", dtype="datetime64[Y]"), n),
        s=rng.choice([2, 4, 8], n),
        f=rng.choice([0.2, 0.3], n),
    ))

    a_cat = df["a"].astype("category")
    new_categories = np.roll(a_cat.cat.categories, 1)
    df["a_cat"] = a_cat.cat.reorder_categories(new_categories)

    df["s_cat"] = df["s"].astype("category")
    df["s_str"] = df["s"].astype(str)

    return df


@pytest.fixture
def long_dict(long_df):

    return long_df.to_dict()


@pytest.fixture
def repeated_df(rng):

    n = 100
    return pd.DataFrame(dict(
        x=np.tile(np.arange(n // 2), 2),
        y=rng.normal(size=n),
        a=rng.choice(list("abc"), n),
        u=np.repeat(np.arange(2), n // 2),
    ))


@pytest.fixture
def missing_df(rng, long_df):

    df = long_df.copy()
    for col in df:
        idx = rng.permutation(df.index)[:10]
        df.loc[idx, col] = np.nan
    return df


@pytest.fixture
def object_df(rng, long_df):

    df = long_df.copy()
    # objectify numeric columns
    for col in ["c", "s", "f"]:
        df[col] = df[col].astype(object)
    return df


@pytest.fixture
def null_series(flat_series):

    return pd.Series(index=flat_series.index, dtype='float64')
