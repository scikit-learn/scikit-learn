"""
Built-in datasets for demonstration, educational and test purposes.
"""


def gapminder(datetimes=False, centroids=False, year=None, pretty_names=False):
    """
    Each row represents a country on a given year.

    https://www.gapminder.org/data/

    Returns:
        A `pandas.DataFrame` with 1704 rows and the following columns:
        `['country', 'continent', 'year', 'lifeExp', 'pop', 'gdpPercap',
        'iso_alpha', 'iso_num']`.
        If `datetimes` is True, the 'year' column will be a datetime column
        If `centroids` is True, two new columns are added: ['centroid_lat', 'centroid_lon']
        If `year` is an integer, the dataset will be filtered for that year
    """
    df = _get_dataset("gapminder")
    if year:
        df = df[df["year"] == year]
    if datetimes:
        df["year"] = (df["year"].astype(str) + "-01-01").astype("datetime64[ns]")
    if not centroids:
        df = df.drop(["centroid_lat", "centroid_lon"], axis=1)
    if pretty_names:
        df.rename(
            mapper=dict(
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
            ),
            axis="columns",
            inplace=True,
        )
    return df


def tips(pretty_names=False):
    """
    Each row represents a restaurant bill.

    https://vincentarelbundock.github.io/Rdatasets/doc/reshape2/tips.html

    Returns:
        A `pandas.DataFrame` with 244 rows and the following columns:
        `['total_bill', 'tip', 'sex', 'smoker', 'day', 'time', 'size']`."""

    df = _get_dataset("tips")
    if pretty_names:
        df.rename(
            mapper=dict(
                total_bill="Total Bill",
                tip="Tip",
                sex="Payer Gender",
                smoker="Smokers at Table",
                day="Day of Week",
                time="Meal",
                size="Party Size",
            ),
            axis="columns",
            inplace=True,
        )
    return df


def iris():
    """
    Each row represents a flower.

    https://en.wikipedia.org/wiki/Iris_flower_data_set

    Returns:
        A `pandas.DataFrame` with 150 rows and the following columns:
        `['sepal_length', 'sepal_width', 'petal_length', 'petal_width', 'species', 'species_id']`."""
    return _get_dataset("iris")


def wind():
    """
    Each row represents a level of wind intensity in a cardinal direction, and its frequency.

    Returns:
        A `pandas.DataFrame` with 128 rows and the following columns:
        `['direction', 'strength', 'frequency']`."""
    return _get_dataset("wind")


def election():
    """
    Each row represents voting results for an electoral district in the 2013 Montreal
    mayoral election.

    Returns:
        A `pandas.DataFrame` with 58 rows and the following columns:
        `['district', 'Coderre', 'Bergeron', 'Joly', 'total', 'winner', 'result', 'district_id']`."""
    return _get_dataset("election")


def election_geojson():
    """
    Each feature represents an electoral district in the 2013 Montreal mayoral election.

    Returns:
        A GeoJSON-formatted `dict` with 58 polygon or multi-polygon features whose `id`
        is an electoral district numerical ID and whose `district` property is the ID and
        district name."""
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


def carshare():
    """
    Each row represents the availability of car-sharing services near the centroid of a zone
    in Montreal over a month-long period.

    Returns:
        A `pandas.DataFrame` with 249 rows and the following columns:
        `['centroid_lat', 'centroid_lon', 'car_hours', 'peak_hour']`."""
    return _get_dataset("carshare")


def stocks(indexed=False, datetimes=False):
    """
    Each row in this wide dataset represents closing prices from 6 tech stocks in 2018/2019.

    Returns:
        A `pandas.DataFrame` with 100 rows and the following columns:
        `['date', 'GOOG', 'AAPL', 'AMZN', 'FB', 'NFLX', 'MSFT']`.
        If `indexed` is True, the 'date' column is used as the index and the column index
        If `datetimes` is True, the 'date' column will be a datetime column
        is named 'company'"""
    df = _get_dataset("stocks")
    if datetimes:
        df["date"] = df["date"].astype("datetime64[ns]")
    if indexed:
        df = df.set_index("date")
        df.columns.name = "company"
    return df


def experiment(indexed=False):
    """
    Each row in this wide dataset represents the results of 100 simulated participants
    on three hypothetical experiments, along with their gender and control/treatment group.


    Returns:
        A `pandas.DataFrame` with 100 rows and the following columns:
        `['experiment_1', 'experiment_2', 'experiment_3', 'gender', 'group']`.
        If `indexed` is True, the data frame index is named "participant" """
    df = _get_dataset("experiment")
    if indexed:
        df.index.name = "participant"
    return df


def medals_wide(indexed=False):
    """
    This dataset represents the medal table for Olympic Short Track Speed Skating for the
    top three nations as of 2020.

    Returns:
        A `pandas.DataFrame` with 3 rows and the following columns:
        `['nation', 'gold', 'silver', 'bronze']`.
        If `indexed` is True, the 'nation' column is used as the index and the column index
        is named 'medal'"""
    df = _get_dataset("medals")
    if indexed:
        df = df.set_index("nation")
        df.columns.name = "medal"
    return df


def medals_long(indexed=False):
    """
    This dataset represents the medal table for Olympic Short Track Speed Skating for the
    top three nations as of 2020.

    Returns:
        A `pandas.DataFrame` with 9 rows and the following columns:
        `['nation', 'medal', 'count']`.
        If `indexed` is True, the 'nation' column is used as the index."""
    df = _get_dataset("medals").melt(
        id_vars=["nation"], value_name="count", var_name="medal"
    )
    if indexed:
        df = df.set_index("nation")
    return df


def _get_dataset(d):
    import pandas
    import os

    return pandas.read_csv(
        os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            "package_data",
            "datasets",
            d + ".csv.gz",
        )
    )
