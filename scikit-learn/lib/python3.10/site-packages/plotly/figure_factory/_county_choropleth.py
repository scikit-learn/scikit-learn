import io
import numpy as np
import os
import pandas as pd
import warnings

from math import log, floor
from numbers import Number

from plotly import optional_imports
import plotly.colors as clrs
from plotly.figure_factory import utils
from plotly.exceptions import PlotlyError
import plotly.graph_objs as go

pd.options.mode.chained_assignment = None

shapely = optional_imports.get_module("shapely")
shapefile = optional_imports.get_module("shapefile")
gp = optional_imports.get_module("geopandas")
_plotly_geo = optional_imports.get_module("_plotly_geo")


def _create_us_counties_df(st_to_state_name_dict, state_to_st_dict):
    # URLS
    abs_dir_path = os.path.realpath(_plotly_geo.__file__)

    abs_plotly_geo_path = os.path.dirname(abs_dir_path)

    abs_package_data_dir_path = os.path.join(abs_plotly_geo_path, "package_data")

    shape_pre2010 = "gz_2010_us_050_00_500k.shp"
    shape_pre2010 = os.path.join(abs_package_data_dir_path, shape_pre2010)

    df_shape_pre2010 = gp.read_file(shape_pre2010)
    df_shape_pre2010["FIPS"] = df_shape_pre2010["STATE"] + df_shape_pre2010["COUNTY"]
    df_shape_pre2010["FIPS"] = pd.to_numeric(df_shape_pre2010["FIPS"])

    states_path = "cb_2016_us_state_500k.shp"
    states_path = os.path.join(abs_package_data_dir_path, states_path)

    df_state = gp.read_file(states_path)
    df_state = df_state[["STATEFP", "NAME", "geometry"]]
    df_state = df_state.rename(columns={"NAME": "STATE_NAME"})

    filenames = [
        "cb_2016_us_county_500k.dbf",
        "cb_2016_us_county_500k.shp",
        "cb_2016_us_county_500k.shx",
    ]

    for j in range(len(filenames)):
        filenames[j] = os.path.join(abs_package_data_dir_path, filenames[j])

    dbf = io.open(filenames[0], "rb")
    shp = io.open(filenames[1], "rb")
    shx = io.open(filenames[2], "rb")

    r = shapefile.Reader(shp=shp, shx=shx, dbf=dbf)

    attributes, geometry = [], []
    field_names = [field[0] for field in r.fields[1:]]
    for row in r.shapeRecords():
        geometry.append(shapely.geometry.shape(row.shape.__geo_interface__))
        attributes.append(dict(zip(field_names, row.record)))

    gdf = gp.GeoDataFrame(data=attributes, geometry=geometry)

    gdf["FIPS"] = gdf["STATEFP"] + gdf["COUNTYFP"]
    gdf["FIPS"] = pd.to_numeric(gdf["FIPS"])

    # add missing counties
    f = 46113
    singlerow = pd.DataFrame(
        [
            [
                st_to_state_name_dict["SD"],
                "SD",
                df_shape_pre2010[df_shape_pre2010["FIPS"] == f]["geometry"].iloc[0],
                df_shape_pre2010[df_shape_pre2010["FIPS"] == f]["FIPS"].iloc[0],
                "46",
                "Shannon",
            ]
        ],
        columns=["State", "ST", "geometry", "FIPS", "STATEFP", "NAME"],
        index=[max(gdf.index) + 1],
    )
    gdf = pd.concat([gdf, singlerow], sort=True)

    f = 51515
    singlerow = pd.DataFrame(
        [
            [
                st_to_state_name_dict["VA"],
                "VA",
                df_shape_pre2010[df_shape_pre2010["FIPS"] == f]["geometry"].iloc[0],
                df_shape_pre2010[df_shape_pre2010["FIPS"] == f]["FIPS"].iloc[0],
                "51",
                "Bedford City",
            ]
        ],
        columns=["State", "ST", "geometry", "FIPS", "STATEFP", "NAME"],
        index=[max(gdf.index) + 1],
    )
    gdf = pd.concat([gdf, singlerow], sort=True)

    f = 2270
    singlerow = pd.DataFrame(
        [
            [
                st_to_state_name_dict["AK"],
                "AK",
                df_shape_pre2010[df_shape_pre2010["FIPS"] == f]["geometry"].iloc[0],
                df_shape_pre2010[df_shape_pre2010["FIPS"] == f]["FIPS"].iloc[0],
                "02",
                "Wade Hampton",
            ]
        ],
        columns=["State", "ST", "geometry", "FIPS", "STATEFP", "NAME"],
        index=[max(gdf.index) + 1],
    )
    gdf = pd.concat([gdf, singlerow], sort=True)

    row_2198 = gdf[gdf["FIPS"] == 2198]
    row_2198.index = [max(gdf.index) + 1]
    row_2198.loc[row_2198.index[0], "FIPS"] = 2201
    row_2198.loc[row_2198.index[0], "STATEFP"] = "02"
    gdf = pd.concat([gdf, row_2198], sort=True)

    row_2105 = gdf[gdf["FIPS"] == 2105]
    row_2105.index = [max(gdf.index) + 1]
    row_2105.loc[row_2105.index[0], "FIPS"] = 2232
    row_2105.loc[row_2105.index[0], "STATEFP"] = "02"
    gdf = pd.concat([gdf, row_2105], sort=True)
    gdf = gdf.rename(columns={"NAME": "COUNTY_NAME"})

    gdf_reduced = gdf[["FIPS", "STATEFP", "COUNTY_NAME", "geometry"]]
    gdf_statefp = gdf_reduced.merge(df_state[["STATEFP", "STATE_NAME"]], on="STATEFP")

    ST = []
    for n in gdf_statefp["STATE_NAME"]:
        ST.append(state_to_st_dict[n])

    gdf_statefp["ST"] = ST
    return gdf_statefp, df_state


st_to_state_name_dict = {
    "AK": "Alaska",
    "AL": "Alabama",
    "AR": "Arkansas",
    "AZ": "Arizona",
    "CA": "California",
    "CO": "Colorado",
    "CT": "Connecticut",
    "DC": "District of Columbia",
    "DE": "Delaware",
    "FL": "Florida",
    "GA": "Georgia",
    "HI": "Hawaii",
    "IA": "Iowa",
    "ID": "Idaho",
    "IL": "Illinois",
    "IN": "Indiana",
    "KS": "Kansas",
    "KY": "Kentucky",
    "LA": "Louisiana",
    "MA": "Massachusetts",
    "MD": "Maryland",
    "ME": "Maine",
    "MI": "Michigan",
    "MN": "Minnesota",
    "MO": "Missouri",
    "MS": "Mississippi",
    "MT": "Montana",
    "NC": "North Carolina",
    "ND": "North Dakota",
    "NE": "Nebraska",
    "NH": "New Hampshire",
    "NJ": "New Jersey",
    "NM": "New Mexico",
    "NV": "Nevada",
    "NY": "New York",
    "OH": "Ohio",
    "OK": "Oklahoma",
    "OR": "Oregon",
    "PA": "Pennsylvania",
    "RI": "Rhode Island",
    "SC": "South Carolina",
    "SD": "South Dakota",
    "TN": "Tennessee",
    "TX": "Texas",
    "UT": "Utah",
    "VA": "Virginia",
    "VT": "Vermont",
    "WA": "Washington",
    "WI": "Wisconsin",
    "WV": "West Virginia",
    "WY": "Wyoming",
}

state_to_st_dict = {
    "Alabama": "AL",
    "Alaska": "AK",
    "American Samoa": "AS",
    "Arizona": "AZ",
    "Arkansas": "AR",
    "California": "CA",
    "Colorado": "CO",
    "Commonwealth of the Northern Mariana Islands": "MP",
    "Connecticut": "CT",
    "Delaware": "DE",
    "District of Columbia": "DC",
    "Florida": "FL",
    "Georgia": "GA",
    "Guam": "GU",
    "Hawaii": "HI",
    "Idaho": "ID",
    "Illinois": "IL",
    "Indiana": "IN",
    "Iowa": "IA",
    "Kansas": "KS",
    "Kentucky": "KY",
    "Louisiana": "LA",
    "Maine": "ME",
    "Maryland": "MD",
    "Massachusetts": "MA",
    "Michigan": "MI",
    "Minnesota": "MN",
    "Mississippi": "MS",
    "Missouri": "MO",
    "Montana": "MT",
    "Nebraska": "NE",
    "Nevada": "NV",
    "New Hampshire": "NH",
    "New Jersey": "NJ",
    "New Mexico": "NM",
    "New York": "NY",
    "North Carolina": "NC",
    "North Dakota": "ND",
    "Ohio": "OH",
    "Oklahoma": "OK",
    "Oregon": "OR",
    "Pennsylvania": "PA",
    "Puerto Rico": "",
    "Rhode Island": "RI",
    "South Carolina": "SC",
    "South Dakota": "SD",
    "Tennessee": "TN",
    "Texas": "TX",
    "United States Virgin Islands": "VI",
    "Utah": "UT",
    "Vermont": "VT",
    "Virginia": "VA",
    "Washington": "WA",
    "West Virginia": "WV",
    "Wisconsin": "WI",
    "Wyoming": "WY",
}

USA_XRANGE = [-125.0, -65.0]
USA_YRANGE = [25.0, 49.0]


def _human_format(number):
    units = ["", "K", "M", "G", "T", "P"]
    k = 1000.0
    magnitude = int(floor(log(number, k)))
    return "%.2f%s" % (number / k**magnitude, units[magnitude])


def _intervals_as_labels(array_of_intervals, round_legend_values, exponent_format):
    """
    Transform an number interval to a clean string for legend

    Example: [-inf, 30] to '< 30'
    """
    infs = [float("-inf"), float("inf")]
    string_intervals = []
    for interval in array_of_intervals:
        # round to 2nd decimal place
        if round_legend_values:
            rnd_interval = [
                (int(interval[i]) if interval[i] not in infs else interval[i])
                for i in range(2)
            ]
        else:
            rnd_interval = [round(interval[0], 2), round(interval[1], 2)]

        num0 = rnd_interval[0]
        num1 = rnd_interval[1]
        if exponent_format:
            if num0 not in infs:
                num0 = _human_format(num0)
            if num1 not in infs:
                num1 = _human_format(num1)
        else:
            if num0 not in infs:
                num0 = "{:,}".format(num0)
            if num1 not in infs:
                num1 = "{:,}".format(num1)

        if num0 == float("-inf"):
            as_str = "< {}".format(num1)
        elif num1 == float("inf"):
            as_str = "> {}".format(num0)
        else:
            as_str = "{} - {}".format(num0, num1)
        string_intervals.append(as_str)
    return string_intervals


def _calculations(
    df,
    fips,
    values,
    index,
    f,
    simplify_county,
    level,
    x_centroids,
    y_centroids,
    centroid_text,
    x_traces,
    y_traces,
    fips_polygon_map,
):
    # 0-pad FIPS code to ensure exactly 5 digits
    padded_f = str(f).zfill(5)
    if fips_polygon_map[f].type == "Polygon":
        x = fips_polygon_map[f].simplify(simplify_county).exterior.xy[0].tolist()
        y = fips_polygon_map[f].simplify(simplify_county).exterior.xy[1].tolist()

        x_c, y_c = fips_polygon_map[f].centroid.xy
        county_name_str = str(df[df["FIPS"] == f]["COUNTY_NAME"].iloc[0])
        state_name_str = str(df[df["FIPS"] == f]["STATE_NAME"].iloc[0])

        t_c = (
            "County: "
            + county_name_str
            + "<br>"
            + "State: "
            + state_name_str
            + "<br>"
            + "FIPS: "
            + padded_f
            + "<br>Value: "
            + str(values[index])
        )

        x_centroids.append(x_c[0])
        y_centroids.append(y_c[0])
        centroid_text.append(t_c)

        x_traces[level] = x_traces[level] + x + [np.nan]
        y_traces[level] = y_traces[level] + y + [np.nan]
    elif fips_polygon_map[f].type == "MultiPolygon":
        x = [
            poly.simplify(simplify_county).exterior.xy[0].tolist()
            for poly in fips_polygon_map[f].geoms
        ]
        y = [
            poly.simplify(simplify_county).exterior.xy[1].tolist()
            for poly in fips_polygon_map[f].geoms
        ]

        x_c = [poly.centroid.xy[0].tolist() for poly in fips_polygon_map[f].geoms]
        y_c = [poly.centroid.xy[1].tolist() for poly in fips_polygon_map[f].geoms]

        county_name_str = str(df[df["FIPS"] == f]["COUNTY_NAME"].iloc[0])
        state_name_str = str(df[df["FIPS"] == f]["STATE_NAME"].iloc[0])
        text = (
            "County: "
            + county_name_str
            + "<br>"
            + "State: "
            + state_name_str
            + "<br>"
            + "FIPS: "
            + padded_f
            + "<br>Value: "
            + str(values[index])
        )
        t_c = [text for poly in fips_polygon_map[f].geoms]
        x_centroids = x_c + x_centroids
        y_centroids = y_c + y_centroids
        centroid_text = t_c + centroid_text
        for x_y_idx in range(len(x)):
            x_traces[level] = x_traces[level] + x[x_y_idx] + [np.nan]
            y_traces[level] = y_traces[level] + y[x_y_idx] + [np.nan]

    return x_traces, y_traces, x_centroids, y_centroids, centroid_text


def create_choropleth(
    fips,
    values,
    scope=["usa"],
    binning_endpoints=None,
    colorscale=None,
    order=None,
    simplify_county=0.02,
    simplify_state=0.02,
    asp=None,
    show_hover=True,
    show_state_data=True,
    state_outline=None,
    county_outline=None,
    centroid_marker=None,
    round_legend_values=False,
    exponent_format=False,
    legend_title="",
    **layout_options,
):
    """
    **deprecated**, use instead
    :func:`plotly.express.choropleth` with custom GeoJSON.

    This function also requires `shapely`, `geopandas` and `plotly-geo` to be installed.

    Returns figure for county choropleth. Uses data from package_data.

    :param (list) fips: list of FIPS values which correspond to the con
        catination of state and county ids. An example is '01001'.
    :param (list) values: list of numbers/strings which correspond to the
        fips list. These are the values that will determine how the counties
        are colored.
    :param (list) scope: list of states and/or states abbreviations. Fits
        all states in the camera tightly. Selecting ['usa'] is the equivalent
        of appending all 50 states into your scope list. Selecting only 'usa'
        does not include 'Alaska', 'Puerto Rico', 'American Samoa',
        'Commonwealth of the Northern Mariana Islands', 'Guam',
        'United States Virgin Islands'. These must be added manually to the
        list.
        Default = ['usa']
    :param (list) binning_endpoints: ascending numbers which implicitly define
        real number intervals which are used as bins. The colorscale used must
        have the same number of colors as the number of bins and this will
        result in a categorical colormap.
    :param (list) colorscale: a list of colors with length equal to the
        number of categories of colors. The length must match either all
        unique numbers in the 'values' list or if endpoints is being used, the
        number of categories created by the endpoints.\n
        For example, if binning_endpoints = [4, 6, 8], then there are 4 bins:
        [-inf, 4), [4, 6), [6, 8), [8, inf)
    :param (list) order: a list of the unique categories (numbers/bins) in any
        desired order. This is helpful if you want to order string values to
        a chosen colorscale.
    :param (float) simplify_county: determines the simplification factor
        for the counties. The larger the number, the fewer vertices and edges
        each polygon has. See
        http://toblerity.org/shapely/manual.html#object.simplify for more
        information.
        Default = 0.02
    :param (float) simplify_state: simplifies the state outline polygon.
        See http://toblerity.org/shapely/manual.html#object.simplify for more
        information.
        Default = 0.02
    :param (float) asp: the width-to-height aspect ratio for the camera.
        Default = 2.5
    :param (bool) show_hover: show county hover and centroid info
    :param (bool) show_state_data: reveals state boundary lines
    :param (dict) state_outline: dict of attributes of the state outline
        including width and color. See
        https://plot.ly/python/reference/#scatter-marker-line for all valid
        params
    :param (dict) county_outline: dict of attributes of the county outline
        including width and color. See
        https://plot.ly/python/reference/#scatter-marker-line for all valid
        params
    :param (dict) centroid_marker: dict of attributes of the centroid marker.
        The centroid markers are invisible by default and appear visible on
        selection. See https://plot.ly/python/reference/#scatter-marker for
        all valid params
    :param (bool) round_legend_values: automatically round the numbers that
        appear in the legend to the nearest integer.
        Default = False
    :param (bool) exponent_format: if set to True, puts numbers in the K, M,
        B number format. For example 4000.0 becomes 4.0K
        Default = False
    :param (str) legend_title: title that appears above the legend
    :param **layout_options: a **kwargs argument for all layout parameters


    Example 1: Florida::

        import plotly.plotly as py
        import plotly.figure_factory as ff

        import numpy as np
        import pandas as pd

        df_sample = pd.read_csv(
            'https://raw.githubusercontent.com/plotly/datasets/master/minoritymajority.csv'
        )
        df_sample_r = df_sample[df_sample['STNAME'] == 'Florida']

        values = df_sample_r['TOT_POP'].tolist()
        fips = df_sample_r['FIPS'].tolist()

        binning_endpoints = list(np.mgrid[min(values):max(values):4j])
        colorscale = ["#030512","#1d1d3b","#323268","#3d4b94","#3e6ab0",
                    "#4989bc","#60a7c7","#85c5d3","#b7e0e4","#eafcfd"]
        fig = ff.create_choropleth(
            fips=fips, values=values, scope=['Florida'], show_state_data=True,
            colorscale=colorscale, binning_endpoints=binning_endpoints,
            round_legend_values=True, plot_bgcolor='rgb(229,229,229)',
            paper_bgcolor='rgb(229,229,229)', legend_title='Florida Population',
            county_outline={'color': 'rgb(255,255,255)', 'width': 0.5},
            exponent_format=True,
        )

    Example 2: New England::

        import plotly.figure_factory as ff

        import pandas as pd

        NE_states = ['Connecticut', 'Maine', 'Massachusetts',
                    'New Hampshire', 'Rhode Island']
        df_sample = pd.read_csv(
            'https://raw.githubusercontent.com/plotly/datasets/master/minoritymajority.csv'
        )
        df_sample_r = df_sample[df_sample['STNAME'].isin(NE_states)]
        colorscale = ['rgb(68.0, 1.0, 84.0)',
        'rgb(66.0, 64.0, 134.0)',
        'rgb(38.0, 130.0, 142.0)',
        'rgb(63.0, 188.0, 115.0)',
        'rgb(216.0, 226.0, 25.0)']

        values = df_sample_r['TOT_POP'].tolist()
        fips = df_sample_r['FIPS'].tolist()
        fig = ff.create_choropleth(
            fips=fips, values=values, scope=NE_states, show_state_data=True
        )
        fig.show()

    Example 3: California and Surrounding States::

        import plotly.figure_factory as ff

        import pandas as pd

        df_sample = pd.read_csv(
            'https://raw.githubusercontent.com/plotly/datasets/master/minoritymajority.csv'
        )
        df_sample_r = df_sample[df_sample['STNAME'] == 'California']

        values = df_sample_r['TOT_POP'].tolist()
        fips = df_sample_r['FIPS'].tolist()

        colorscale = [
            'rgb(193, 193, 193)',
            'rgb(239,239,239)',
            'rgb(195, 196, 222)',
            'rgb(144,148,194)',
            'rgb(101,104,168)',
            'rgb(65, 53, 132)'
        ]

        fig = ff.create_choropleth(
            fips=fips, values=values, colorscale=colorscale,
            scope=['CA', 'AZ', 'Nevada', 'Oregon', ' Idaho'],
            binning_endpoints=[14348, 63983, 134827, 426762, 2081313],
            county_outline={'color': 'rgb(255,255,255)', 'width': 0.5},
            legend_title='California Counties',
            title='California and Nearby States'
        )
        fig.show()

    Example 4: USA::

        import plotly.figure_factory as ff

        import numpy as np
        import pandas as pd

        df_sample = pd.read_csv(
            'https://raw.githubusercontent.com/plotly/datasets/master/laucnty16.csv'
        )
        df_sample['State FIPS Code'] = df_sample['State FIPS Code'].apply(
            lambda x: str(x).zfill(2)
        )
        df_sample['County FIPS Code'] = df_sample['County FIPS Code'].apply(
            lambda x: str(x).zfill(3)
        )
        df_sample['FIPS'] = (
            df_sample['State FIPS Code'] + df_sample['County FIPS Code']
        )

        binning_endpoints = list(np.linspace(1, 12, len(colorscale) - 1))
        colorscale = ["#f7fbff", "#ebf3fb", "#deebf7", "#d2e3f3", "#c6dbef",
                    "#b3d2e9", "#9ecae1", "#85bcdb", "#6baed6", "#57a0ce",
                    "#4292c6", "#3082be", "#2171b5", "#1361a9", "#08519c",
                    "#0b4083","#08306b"]
        fips = df_sample['FIPS']
        values = df_sample['Unemployment Rate (%)']
        fig = ff.create_choropleth(
            fips=fips, values=values, scope=['usa'],
            binning_endpoints=binning_endpoints, colorscale=colorscale,
            show_hover=True, centroid_marker={'opacity': 0},
            asp=2.9, title='USA by Unemployment %',
            legend_title='Unemployment %'
        )
        fig.show()
    """
    # ensure optional modules imported
    if not _plotly_geo:
        raise ValueError(
            """
The create_choropleth figure factory requires the plotly-geo package.
Install using pip with:

$ pip install plotly-geo

Or, install using conda with

$ conda install -c plotly plotly-geo
"""
        )

    if not gp or not shapefile or not shapely:
        raise ImportError(
            "geopandas, pyshp and shapely must be installed for this figure "
            "factory.\n\nRun the following commands to install the correct "
            "versions of the following modules:\n\n"
            "```\n"
            "$ pip install geopandas==0.3.0\n"
            "$ pip install pyshp==1.2.10\n"
            "$ pip install shapely==1.6.3\n"
            "```\n"
            "If you are using Windows, follow this post to properly "
            "install geopandas and dependencies:"
            "http://geoffboeing.com/2014/09/using-geopandas-windows/\n\n"
            "If you are using Anaconda, do not use PIP to install the "
            "packages above. Instead use conda to install them:\n\n"
            "```\n"
            "$ conda install plotly\n"
            "$ conda install geopandas\n"
            "```"
        )

    df, df_state = _create_us_counties_df(st_to_state_name_dict, state_to_st_dict)

    fips_polygon_map = dict(zip(df["FIPS"].tolist(), df["geometry"].tolist()))

    if not state_outline:
        state_outline = {"color": "rgb(240, 240, 240)", "width": 1}
    if not county_outline:
        county_outline = {"color": "rgb(0, 0, 0)", "width": 0}
    if not centroid_marker:
        centroid_marker = {"size": 3, "color": "white", "opacity": 1}

    # ensure centroid markers appear on selection
    if "opacity" not in centroid_marker:
        centroid_marker.update({"opacity": 1})

    if len(fips) != len(values):
        raise PlotlyError("fips and values must be the same length")

    # make fips, values into lists
    if isinstance(fips, pd.core.series.Series):
        fips = fips.tolist()
    if isinstance(values, pd.core.series.Series):
        values = values.tolist()

    # make fips numeric
    fips = map(lambda x: int(x), fips)

    if binning_endpoints:
        intervals = utils.endpts_to_intervals(binning_endpoints)
        LEVELS = _intervals_as_labels(intervals, round_legend_values, exponent_format)
    else:
        if not order:
            LEVELS = sorted(list(set(values)))
        else:
            # check if order is permutation
            # of unique color col values
            same_sets = sorted(list(set(values))) == set(order)
            no_duplicates = not any(order.count(x) > 1 for x in order)
            if same_sets and no_duplicates:
                LEVELS = order
            else:
                raise PlotlyError(
                    "if you are using a custom order of unique values from "
                    "your color column, you must: have all the unique values "
                    "in your order and have no duplicate items"
                )

    if not colorscale:
        colorscale = []
        viridis_colors = clrs.colorscale_to_colors(clrs.PLOTLY_SCALES["Viridis"])
        viridis_colors = clrs.color_parser(viridis_colors, clrs.hex_to_rgb)
        viridis_colors = clrs.color_parser(viridis_colors, clrs.label_rgb)
        viri_len = len(viridis_colors) + 1
        viri_intervals = utils.endpts_to_intervals(list(np.linspace(0, 1, viri_len)))[
            1:-1
        ]

        for L in np.linspace(0, 1, len(LEVELS)):
            for idx, inter in enumerate(viri_intervals):
                if L == 0:
                    break
                elif inter[0] < L <= inter[1]:
                    break

            intermed = (L - viri_intervals[idx][0]) / (
                viri_intervals[idx][1] - viri_intervals[idx][0]
            )

            float_color = clrs.find_intermediate_color(
                viridis_colors[idx], viridis_colors[idx], intermed, colortype="rgb"
            )

            # make R,G,B into int values
            float_color = clrs.unlabel_rgb(float_color)
            float_color = clrs.unconvert_from_RGB_255(float_color)
            int_rgb = clrs.convert_to_RGB_255(float_color)
            int_rgb = clrs.label_rgb(int_rgb)

            colorscale.append(int_rgb)

    if len(colorscale) < len(LEVELS):
        raise PlotlyError(
            "You have {} LEVELS. Your number of colors in 'colorscale' must "
            "be at least the number of LEVELS: {}. If you are "
            "using 'binning_endpoints' then 'colorscale' must have at "
            "least len(binning_endpoints) + 2 colors".format(
                len(LEVELS), min(LEVELS, LEVELS[:20])
            )
        )

    color_lookup = dict(zip(LEVELS, colorscale))
    x_traces = dict(zip(LEVELS, [[] for i in range(len(LEVELS))]))
    y_traces = dict(zip(LEVELS, [[] for i in range(len(LEVELS))]))

    # scope
    if isinstance(scope, str):
        raise PlotlyError("'scope' must be a list/tuple/sequence")

    scope_names = []
    extra_states = [
        "Alaska",
        "Commonwealth of the Northern Mariana Islands",
        "Puerto Rico",
        "Guam",
        "United States Virgin Islands",
        "American Samoa",
    ]
    for state in scope:
        if state.lower() == "usa":
            scope_names = df["STATE_NAME"].unique()
            scope_names = list(scope_names)
            for ex_st in extra_states:
                try:
                    scope_names.remove(ex_st)
                except ValueError:
                    pass
        else:
            if state in st_to_state_name_dict.keys():
                state = st_to_state_name_dict[state]
            scope_names.append(state)
    df_state = df_state[df_state["STATE_NAME"].isin(scope_names)]

    plot_data = []
    x_centroids = []
    y_centroids = []
    centroid_text = []
    fips_not_in_shapefile = []
    if not binning_endpoints:
        for index, f in enumerate(fips):
            level = values[index]
            try:
                fips_polygon_map[f].type

                (
                    x_traces,
                    y_traces,
                    x_centroids,
                    y_centroids,
                    centroid_text,
                ) = _calculations(
                    df,
                    fips,
                    values,
                    index,
                    f,
                    simplify_county,
                    level,
                    x_centroids,
                    y_centroids,
                    centroid_text,
                    x_traces,
                    y_traces,
                    fips_polygon_map,
                )
            except KeyError:
                fips_not_in_shapefile.append(f)

    else:
        for index, f in enumerate(fips):
            for j, inter in enumerate(intervals):
                if inter[0] < values[index] <= inter[1]:
                    break
            level = LEVELS[j]

            try:
                fips_polygon_map[f].type

                (
                    x_traces,
                    y_traces,
                    x_centroids,
                    y_centroids,
                    centroid_text,
                ) = _calculations(
                    df,
                    fips,
                    values,
                    index,
                    f,
                    simplify_county,
                    level,
                    x_centroids,
                    y_centroids,
                    centroid_text,
                    x_traces,
                    y_traces,
                    fips_polygon_map,
                )
            except KeyError:
                fips_not_in_shapefile.append(f)

    if len(fips_not_in_shapefile) > 0:
        msg = (
            "Unrecognized FIPS Values\n\nWhoops! It looks like you are "
            "trying to pass at least one FIPS value that is not in "
            "our shapefile of FIPS and data for the counties. Your "
            "choropleth will still show up but these counties cannot "
            "be shown.\nUnrecognized FIPS are: {}".format(fips_not_in_shapefile)
        )
        warnings.warn(msg)

    x_states = []
    y_states = []
    for index, row in df_state.iterrows():
        if df_state["geometry"][index].type == "Polygon":
            x = row.geometry.simplify(simplify_state).exterior.xy[0].tolist()
            y = row.geometry.simplify(simplify_state).exterior.xy[1].tolist()
            x_states = x_states + x
            y_states = y_states + y
        elif df_state["geometry"][index].type == "MultiPolygon":
            x = [
                poly.simplify(simplify_state).exterior.xy[0].tolist()
                for poly in df_state["geometry"][index].geoms
            ]
            y = [
                poly.simplify(simplify_state).exterior.xy[1].tolist()
                for poly in df_state["geometry"][index].geoms
            ]
            for segment in range(len(x)):
                x_states = x_states + x[segment]
                y_states = y_states + y[segment]
                x_states.append(np.nan)
                y_states.append(np.nan)
        x_states.append(np.nan)
        y_states.append(np.nan)

    for lev in LEVELS:
        county_data = dict(
            type="scatter",
            mode="lines",
            x=x_traces[lev],
            y=y_traces[lev],
            line=county_outline,
            fill="toself",
            fillcolor=color_lookup[lev],
            name=lev,
            hoverinfo="none",
        )
        plot_data.append(county_data)

    if show_hover:
        hover_points = dict(
            type="scatter",
            showlegend=False,
            legendgroup="centroids",
            x=x_centroids,
            y=y_centroids,
            text=centroid_text,
            name="US Counties",
            mode="markers",
            marker={"color": "white", "opacity": 0},
            hoverinfo="text",
        )
        centroids_on_select = dict(
            selected=dict(marker=centroid_marker),
            unselected=dict(marker=dict(opacity=0)),
        )
        hover_points.update(centroids_on_select)
        plot_data.append(hover_points)

    if show_state_data:
        state_data = dict(
            type="scatter",
            legendgroup="States",
            line=state_outline,
            x=x_states,
            y=y_states,
            hoverinfo="text",
            showlegend=False,
            mode="lines",
        )
        plot_data.append(state_data)

    DEFAULT_LAYOUT = dict(
        hovermode="closest",
        xaxis=dict(
            autorange=False,
            range=USA_XRANGE,
            showgrid=False,
            zeroline=False,
            fixedrange=True,
            showticklabels=False,
        ),
        yaxis=dict(
            autorange=False,
            range=USA_YRANGE,
            showgrid=False,
            zeroline=False,
            fixedrange=True,
            showticklabels=False,
        ),
        margin=dict(t=40, b=20, r=20, l=20),
        width=900,
        height=450,
        dragmode="select",
        legend=dict(traceorder="reversed", xanchor="right", yanchor="top", x=1, y=1),
        annotations=[],
    )
    fig = dict(data=plot_data, layout=DEFAULT_LAYOUT)
    fig["layout"].update(layout_options)
    fig["layout"]["annotations"].append(
        dict(
            x=1,
            y=1.05,
            xref="paper",
            yref="paper",
            xanchor="right",
            showarrow=False,
            text="<b>" + legend_title + "</b>",
        )
    )

    if len(scope) == 1 and scope[0].lower() == "usa":
        xaxis_range_low = -125.0
        xaxis_range_high = -55.0
        yaxis_range_low = 25.0
        yaxis_range_high = 49.0
    else:
        xaxis_range_low = float("inf")
        xaxis_range_high = float("-inf")
        yaxis_range_low = float("inf")
        yaxis_range_high = float("-inf")
        for trace in fig["data"]:
            if all(isinstance(n, Number) for n in trace["x"]):
                calc_x_min = min(trace["x"] or [float("inf")])
                calc_x_max = max(trace["x"] or [float("-inf")])
                if calc_x_min < xaxis_range_low:
                    xaxis_range_low = calc_x_min
                if calc_x_max > xaxis_range_high:
                    xaxis_range_high = calc_x_max
            if all(isinstance(n, Number) for n in trace["y"]):
                calc_y_min = min(trace["y"] or [float("inf")])
                calc_y_max = max(trace["y"] or [float("-inf")])
                if calc_y_min < yaxis_range_low:
                    yaxis_range_low = calc_y_min
                if calc_y_max > yaxis_range_high:
                    yaxis_range_high = calc_y_max

    # camera zoom
    fig["layout"]["xaxis"]["range"] = [xaxis_range_low, xaxis_range_high]
    fig["layout"]["yaxis"]["range"] = [yaxis_range_low, yaxis_range_high]

    # aspect ratio
    if asp is None:
        usa_x_range = USA_XRANGE[1] - USA_XRANGE[0]
        usa_y_range = USA_YRANGE[1] - USA_YRANGE[0]
        asp = usa_x_range / usa_y_range

    # based on your figure
    width = float(
        fig["layout"]["xaxis"]["range"][1] - fig["layout"]["xaxis"]["range"][0]
    )
    height = float(
        fig["layout"]["yaxis"]["range"][1] - fig["layout"]["yaxis"]["range"][0]
    )

    center = (
        sum(fig["layout"]["xaxis"]["range"]) / 2.0,
        sum(fig["layout"]["yaxis"]["range"]) / 2.0,
    )

    if height / width > (1 / asp):
        new_width = asp * height
        fig["layout"]["xaxis"]["range"][0] = center[0] - new_width * 0.5
        fig["layout"]["xaxis"]["range"][1] = center[0] + new_width * 0.5
    else:
        new_height = (1 / asp) * width
        fig["layout"]["yaxis"]["range"][0] = center[1] - new_height * 0.5
        fig["layout"]["yaxis"]["range"][1] = center[1] + new_height * 0.5

    return go.Figure(fig)
