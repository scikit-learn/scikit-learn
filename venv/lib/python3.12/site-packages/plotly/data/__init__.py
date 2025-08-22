"""
Built-in datasets for demonstration, educational and test purposes.
"""

import os
from importlib import import_module

import narwhals.stable.v1 as nw

AVAILABLE_BACKENDS = {"pandas", "polars", "pyarrow", "modin", "cudf"}
BACKENDS_WITH_INDEX_SUPPORT = {"pandas", "modin", "cudf"}


def gapminder(
    datetimes=False,
    centroids=False,
    year=None,
    pretty_names=False,
    return_type="pandas",
):
    """
    Each row represents a country on a given year.

    https://www.gapminder.org/data/

    Parameters
    ----------
    datetimes: bool
        Whether or not 'year' column will converted to datetime type

    centroids: bool
        If True, ['centroid_lat', 'centroid_lon'] columns are added

    year: int | None
        If provided, the dataset will be filtered for that year

    pretty_names: bool
        If True, prettifies the column names

    return_type: {'pandas', 'polars', 'pyarrow', 'modin', 'cudf'}
        Type of the resulting dataframe

    Returns
    -------
    Dataframe of `return_type` type
        Dataframe with 1704 rows and the following columns:
        `['country', 'continent', 'year', 'lifeExp', 'pop', 'gdpPercap',
        'iso_alpha', 'iso_num']`.

        If `datetimes` is True, the 'year' column will be a datetime column
        If `centroids` is True, two new columns are added: ['centroid_lat', 'centroid_lon']
        If `year` is an integer, the dataset will be filtered for that year
    """
    df = nw.from_native(
        _get_dataset("gapminder", return_type=return_type), eager_only=True
    )
    if year:
        df = df.filter(nw.col("year") == year)
    if datetimes:
        df = df.with_columns(
            # Concatenate the year value with the literal "-01-01" so that it can be
            # casted to datetime from "%Y-%m-%d" format
            nw.concat_str(
                [nw.col("year").cast(nw.String()), nw.lit("-01-01")]
            ).str.to_datetime(format="%Y-%m-%d")
        )
    if not centroids:
        df = df.drop("centroid_lat", "centroid_lon")
    if pretty_names:
        df = df.rename(
            dict(
                country="Country",
                continent="Continent",
                year="Year",
                lifeExp="Life Expectancy",
                gdpPercap="GDP per Capita",
                pop="Population",
                iso_alpha="ISO Alpha Country Code",
                iso_num="ISO Numeric Country Code",
                centroid_lat="Centroid Latitude",
                centroid_lon="Centroid Longitude",
            )
        )
    return df.to_native()


def tips(pretty_names=False, return_type="pandas"):
    """
    Each row represents a restaurant bill.

    https://vincentarelbundock.github.io/Rdatasets/doc/reshape2/tips.html

    Parameters
    ----------
    pretty_names: bool
        If True, prettifies the column names

    return_type: {'pandas', 'polars', 'pyarrow', 'modin', 'cudf'}
        Type of the resulting dataframe

    Returns
    -------
    Dataframe of `return_type` type
        Dataframe with 244 rows and the following columns:
        `['total_bill', 'tip', 'sex', 'smoker', 'day', 'time', 'size']`.
    """

    df = nw.from_native(_get_dataset("tips", return_type=return_type), eager_only=True)
    if pretty_names:
        df = df.rename(
            dict(
                total_bill="Total Bill",
                tip="Tip",
                sex="Payer Gender",
                smoker="Smokers at Table",
                day="Day of Week",
                time="Meal",
                size="Party Size",
            )
        )
    return df.to_native()


def iris(return_type="pandas"):
    """
    Each row represents a flower.

    https://en.wikipedia.org/wiki/Iris_flower_data_set

    Parameters
    ----------
    return_type: {'pandas', 'polars', 'pyarrow', 'modin', 'cudf'}
        Type of the resulting dataframe

    Returns
    -------
    Dataframe of `return_type` type
        Dataframe with 150 rows and the following columns:
        `['sepal_length', 'sepal_width', 'petal_length', 'petal_width', 'species', 'species_id']`.
    """
    return _get_dataset("iris", return_type=return_type)


def wind(return_type="pandas"):
    """
    Each row represents a level of wind intensity in a cardinal direction, and its frequency.

    Parameters
    ----------
    return_type: {'pandas', 'polars', 'pyarrow', 'modin', 'cudf'}
        Type of the resulting dataframe

    Returns
    -------
    Dataframe of `return_type` type
        Dataframe with 128 rows and the following columns:
        `['direction', 'strength', 'frequency']`.
    """
    return _get_dataset("wind", return_type=return_type)


def election(return_type="pandas"):
    """
    Each row represents voting results for an electoral district in the 2013 Montreal
    mayoral election.

    Parameters
    ----------
    return_type: {'pandas', 'polars', 'pyarrow', 'modin', 'cudf'}
        Type of the resulting dataframe

    Returns
    -------
    Dataframe of `return_type` type
        Dataframe with 58 rows and the following columns:
        `['district', 'Coderre', 'Bergeron', 'Joly', 'total', 'winner', 'result', 'district_id']`.
    """
    return _get_dataset("election", return_type=return_type)


def election_geojson():
    """
    Each feature represents an electoral district in the 2013 Montreal mayoral election.

    Returns
    -------
        A GeoJSON-formatted `dict` with 58 polygon or multi-polygon features whose `id`
        is an electoral district numerical ID and whose `district` property is the ID and
        district name.
    """
    import gzip
    import json
    import os

    path = os.path.join(
        os.path.dirname(os.path.dirname(__file__)),
        "package_data",
        "datasets",
        "election.geojson.gz",
    )
    with gzip.GzipFile(path, "r") as f:
        result = json.loads(f.read().decode("utf-8"))
    return result


def carshare(return_type="pandas"):
    """
    Each row represents the availability of car-sharing services near the centroid of a zone
    in Montreal over a month-long period.

    Parameters
    ----------
    return_type: {'pandas', 'polars', 'pyarrow', 'modin', 'cudf'}
        Type of the resulting dataframe

    Returns
    -------
    Dataframe of `return_type` type
        Dataframe` with 249 rows and the following columns:
        `['centroid_lat', 'centroid_lon', 'car_hours', 'peak_hour']`.
    """
    return _get_dataset("carshare", return_type=return_type)


def stocks(indexed=False, datetimes=False, return_type="pandas"):
    """
    Each row in this wide dataset represents closing prices from 6 tech stocks in 2018/2019.

    Parameters
    ----------
    indexed: bool
        Whether or not the 'date' column is used as the index and the column index
        is named 'company'. Applicable only if `return_type='pandas'`

    datetimes: bool
        Whether or not the 'date' column will be of datetime type

    return_type: {'pandas', 'polars', 'pyarrow', 'modin', 'cudf'}
        Type of the resulting dataframe

    Returns
    -------
    Dataframe of `return_type` type
        Dataframe with 100 rows and the following columns:
        `['date', 'GOOG', 'AAPL', 'AMZN', 'FB', 'NFLX', 'MSFT']`.
        If `indexed` is True, the 'date' column is used as the index and the column index
        is named 'company'
        If `datetimes` is True, the 'date' column will be a datetime column
    """
    if indexed and return_type not in BACKENDS_WITH_INDEX_SUPPORT:
        msg = f"Backend '{return_type}' does not support setting index"
        raise NotImplementedError(msg)

    df = nw.from_native(
        _get_dataset("stocks", return_type=return_type), eager_only=True
    ).with_columns(nw.col("date").cast(nw.String()))

    if datetimes:
        df = df.with_columns(nw.col("date").str.to_datetime())

    if indexed:  # then it must be pandas
        df = df.to_native().set_index("date")
        df.columns.name = "company"
        return df

    return df.to_native()


def experiment(indexed=False, return_type="pandas"):
    """
    Each row in this wide dataset represents the results of 100 simulated participants
    on three hypothetical experiments, along with their gender and control/treatment group.

    Parameters
    ----------
    indexed: bool
        If True, then the index is named "participant".
        Applicable only if `return_type='pandas'`

    return_type: {'pandas', 'polars', 'pyarrow', 'modin', 'cudf'}
        Type of the resulting dataframe

    Returns
    -------
    Dataframe of `return_type` type
        Dataframe with 100 rows and the following columns:
        `['experiment_1', 'experiment_2', 'experiment_3', 'gender', 'group']`.
        If `indexed` is True, the data frame index is named "participant"
    """

    if indexed and return_type not in BACKENDS_WITH_INDEX_SUPPORT:
        msg = f"Backend '{return_type}' does not support setting index"
        raise NotImplementedError(msg)

    df = nw.from_native(
        _get_dataset("experiment", return_type=return_type), eager_only=True
    )
    if indexed:  # then it must be pandas
        df = df.to_native()
        df.index.name = "participant"
        return df
    return df.to_native()


def medals_wide(indexed=False, return_type="pandas"):
    """
    This dataset represents the medal table for Olympic Short Track Speed Skating for the
    top three nations as of 2020.

    Parameters
    ----------
    indexed: bool
        Whether or not the 'nation' column is used as the index and the column index
        is named 'medal'. Applicable only if `return_type='pandas'`

    return_type: {'pandas', 'polars', 'pyarrow', 'modin', 'cudf'}
        Type of the resulting dataframe

    Returns
    -------
    Dataframe of `return_type` type
        Dataframe with 3 rows and the following columns:
        `['nation', 'gold', 'silver', 'bronze']`.
        If `indexed` is True, the 'nation' column is used as the index and the column index
        is named 'medal'
    """

    if indexed and return_type not in BACKENDS_WITH_INDEX_SUPPORT:
        msg = f"Backend '{return_type}' does not support setting index"
        raise NotImplementedError(msg)

    df = nw.from_native(
        _get_dataset("medals", return_type=return_type), eager_only=True
    )
    if indexed:  # then it must be pandas
        df = df.to_native().set_index("nation")
        df.columns.name = "medal"
        return df
    return df.to_native()


def medals_long(indexed=False, return_type="pandas"):
    """
    This dataset represents the medal table for Olympic Short Track Speed Skating for the
    top three nations as of 2020.

    Parameters
    ----------
    indexed: bool
        Whether or not the 'nation' column is used as the index.
        Applicable only if `return_type='pandas'`

    return_type: {'pandas', 'polars', 'pyarrow', 'modin', 'cudf'}
        Type of the resulting dataframe

    Returns
    -------
    Dataframe of `return_type` type
        Dataframe with 9 rows and the following columns: `['nation', 'medal', 'count']`.
        If `indexed` is True, the 'nation' column is used as the index.
    """

    if indexed and return_type not in BACKENDS_WITH_INDEX_SUPPORT:
        msg = f"Backend '{return_type}' does not support setting index"
        raise NotImplementedError(msg)

    df = nw.from_native(
        _get_dataset("medals", return_type=return_type), eager_only=True
    ).unpivot(
        index=["nation"],
        value_name="count",
        variable_name="medal",
    )
    if indexed:
        df = nw.maybe_set_index(df, "nation")
    return df.to_native()


def _get_dataset(d, return_type):
    """
    Loads the dataset using the specified backend.

    Notice that the available backends are 'pandas', 'polars', 'pyarrow' and they all have
    a `read_csv` function (pyarrow has it via pyarrow.csv). Therefore we can dynamically
    load the library using `importlib.import_module` and then call
    `backend.read_csv(filepath)`.

    Parameters
    ----------
    d: str
        Name of the dataset to load.

    return_type: {'pandas', 'polars', 'pyarrow', 'modin', 'cudf'}
        Type of the resulting dataframe

    Returns
    -------
    Dataframe of `return_type` type
    """
    filepath = os.path.join(
        os.path.dirname(os.path.dirname(__file__)),
        "package_data",
        "datasets",
        d + ".csv.gz",
    )

    if return_type not in AVAILABLE_BACKENDS:
        msg = (
            f"Unsupported return_type. Found {return_type}, expected one "
            f"of {AVAILABLE_BACKENDS}"
        )
        raise NotImplementedError(msg)

    try:
        if return_type == "pyarrow":
            module_to_load = "pyarrow.csv"
        elif return_type == "modin":
            module_to_load = "modin.pandas"
        else:
            module_to_load = return_type
        backend = import_module(module_to_load)
    except ModuleNotFoundError:
        msg = f"return_type={return_type}, but {return_type} is not installed"
        raise ModuleNotFoundError(msg)

    try:
        return backend.read_csv(filepath)
    except Exception as e:
        msg = f"Unable to read '{d}' dataset due to: {e}"
        raise Exception(msg).with_traceback(e.__traceback__)
