"""
colors
=====

Functions that manipulate colors and arrays of colors.

-----
There are three basic types of color types: rgb, hex and tuple:

rgb - An rgb color is a string of the form 'rgb(a,b,c)' where a, b and c are
integers between 0 and 255 inclusive.

hex - A hex color is a string of the form '#xxxxxx' where each x is a
character that belongs to the set [0,1,2,3,4,5,6,7,8,9,a,b,c,d,e,f]. This is
just the set of characters used in the hexadecimal numeric system.

tuple - A tuple color is a 3-tuple of the form (a,b,c) where a, b and c are
floats between 0 and 1 inclusive.

-----
Colormaps and Colorscales:
A colormap or a colorscale is a correspondence between values - Pythonic
objects such as strings and floats - to colors.

There are typically two main types of colormaps that exist: numerical and
categorical colormaps.

Numerical:
----------
Numerical colormaps are used when the coloring column being used takes a
spectrum of values or numbers.

A classic example from the Plotly library:
```
rainbow_colorscale =  [
    [0, 'rgb(150,0,90)'], [0.125, 'rgb(0,0,200)'],
    [0.25, 'rgb(0,25,255)'], [0.375, 'rgb(0,152,255)'],
    [0.5, 'rgb(44,255,150)'], [0.625, 'rgb(151,255,0)'],
    [0.75, 'rgb(255,234,0)'], [0.875, 'rgb(255,111,0)'],
    [1, 'rgb(255,0,0)']
]
```

Notice that this colorscale is a list of lists with each inner list containing
a number and a color. These left hand numbers in the nested lists go from 0 to
1, and they are like pointers tell you when a number is mapped to a specific
color.

If you have a column of numbers `col_num` that you want to plot, and you know

```
min(col_num) = 0
max(col_num) = 100
```

then if you pull out the number `12.5` in the list and want to figure out what
color the corresponding chart element (bar, scatter plot, etc) is going to be,
you'll figure out that proportionally 12.5 to 100 is the same as 0.125 to 1.
So, the point will be mapped to 'rgb(0,0,200)'.

All other colors between the pinned values in a colorscale are linearly
interpolated.

Categorical:
------------
Alternatively, a categorical colormap is used to assign a specific value in a
color column to a specific color everytime it appears in the dataset.

A column of strings in a panadas.dataframe that is chosen to serve as the
color index would naturally use a categorical colormap. However, you can
choose to use a categorical colormap with a column of numbers.

Be careful! If you have a lot of unique numbers in your color column you will
end up with a colormap that is massive and may slow down graphing performance.
"""

import decimal
from numbers import Number

from _plotly_utils import exceptions


# Built-in qualitative color sequences and sequential,
# diverging and cyclical color scales.
#
# Initially ported over from plotly_express
from . import (  # noqa: F401
    qualitative,
    sequential,
    diverging,
    cyclical,
    cmocean,
    colorbrewer,
    carto,
    plotlyjs,
)

DEFAULT_PLOTLY_COLORS = [
    "rgb(31, 119, 180)",
    "rgb(255, 127, 14)",
    "rgb(44, 160, 44)",
    "rgb(214, 39, 40)",
    "rgb(148, 103, 189)",
    "rgb(140, 86, 75)",
    "rgb(227, 119, 194)",
    "rgb(127, 127, 127)",
    "rgb(188, 189, 34)",
    "rgb(23, 190, 207)",
]

PLOTLY_SCALES = {
    "Greys": [[0, "rgb(0,0,0)"], [1, "rgb(255,255,255)"]],
    "YlGnBu": [
        [0, "rgb(8,29,88)"],
        [0.125, "rgb(37,52,148)"],
        [0.25, "rgb(34,94,168)"],
        [0.375, "rgb(29,145,192)"],
        [0.5, "rgb(65,182,196)"],
        [0.625, "rgb(127,205,187)"],
        [0.75, "rgb(199,233,180)"],
        [0.875, "rgb(237,248,217)"],
        [1, "rgb(255,255,217)"],
    ],
    "Greens": [
        [0, "rgb(0,68,27)"],
        [0.125, "rgb(0,109,44)"],
        [0.25, "rgb(35,139,69)"],
        [0.375, "rgb(65,171,93)"],
        [0.5, "rgb(116,196,118)"],
        [0.625, "rgb(161,217,155)"],
        [0.75, "rgb(199,233,192)"],
        [0.875, "rgb(229,245,224)"],
        [1, "rgb(247,252,245)"],
    ],
    "YlOrRd": [
        [0, "rgb(128,0,38)"],
        [0.125, "rgb(189,0,38)"],
        [0.25, "rgb(227,26,28)"],
        [0.375, "rgb(252,78,42)"],
        [0.5, "rgb(253,141,60)"],
        [0.625, "rgb(254,178,76)"],
        [0.75, "rgb(254,217,118)"],
        [0.875, "rgb(255,237,160)"],
        [1, "rgb(255,255,204)"],
    ],
    "Bluered": [[0, "rgb(0,0,255)"], [1, "rgb(255,0,0)"]],
    # modified RdBu based on
    # www.sandia.gov/~kmorel/documents/ColorMaps/ColorMapsExpanded.pdf
    "RdBu": [
        [0, "rgb(5,10,172)"],
        [0.35, "rgb(106,137,247)"],
        [0.5, "rgb(190,190,190)"],
        [0.6, "rgb(220,170,132)"],
        [0.7, "rgb(230,145,90)"],
        [1, "rgb(178,10,28)"],
    ],
    # Scale for non-negative numeric values
    "Reds": [
        [0, "rgb(220,220,220)"],
        [0.2, "rgb(245,195,157)"],
        [0.4, "rgb(245,160,105)"],
        [1, "rgb(178,10,28)"],
    ],
    # Scale for non-positive numeric values
    "Blues": [
        [0, "rgb(5,10,172)"],
        [0.35, "rgb(40,60,190)"],
        [0.5, "rgb(70,100,245)"],
        [0.6, "rgb(90,120,245)"],
        [0.7, "rgb(106,137,247)"],
        [1, "rgb(220,220,220)"],
    ],
    "Picnic": [
        [0, "rgb(0,0,255)"],
        [0.1, "rgb(51,153,255)"],
        [0.2, "rgb(102,204,255)"],
        [0.3, "rgb(153,204,255)"],
        [0.4, "rgb(204,204,255)"],
        [0.5, "rgb(255,255,255)"],
        [0.6, "rgb(255,204,255)"],
        [0.7, "rgb(255,153,255)"],
        [0.8, "rgb(255,102,204)"],
        [0.9, "rgb(255,102,102)"],
        [1, "rgb(255,0,0)"],
    ],
    "Rainbow": [
        [0, "rgb(150,0,90)"],
        [0.125, "rgb(0,0,200)"],
        [0.25, "rgb(0,25,255)"],
        [0.375, "rgb(0,152,255)"],
        [0.5, "rgb(44,255,150)"],
        [0.625, "rgb(151,255,0)"],
        [0.75, "rgb(255,234,0)"],
        [0.875, "rgb(255,111,0)"],
        [1, "rgb(255,0,0)"],
    ],
    "Portland": [
        [0, "rgb(12,51,131)"],
        [0.25, "rgb(10,136,186)"],
        [0.5, "rgb(242,211,56)"],
        [0.75, "rgb(242,143,56)"],
        [1, "rgb(217,30,30)"],
    ],
    "Jet": [
        [0, "rgb(0,0,131)"],
        [0.125, "rgb(0,60,170)"],
        [0.375, "rgb(5,255,255)"],
        [0.625, "rgb(255,255,0)"],
        [0.875, "rgb(250,0,0)"],
        [1, "rgb(128,0,0)"],
    ],
    "Hot": [
        [0, "rgb(0,0,0)"],
        [0.3, "rgb(230,0,0)"],
        [0.6, "rgb(255,210,0)"],
        [1, "rgb(255,255,255)"],
    ],
    "Blackbody": [
        [0, "rgb(0,0,0)"],
        [0.2, "rgb(230,0,0)"],
        [0.4, "rgb(230,210,0)"],
        [0.7, "rgb(255,255,255)"],
        [1, "rgb(160,200,255)"],
    ],
    "Earth": [
        [0, "rgb(0,0,130)"],
        [0.1, "rgb(0,180,180)"],
        [0.2, "rgb(40,210,40)"],
        [0.4, "rgb(230,230,50)"],
        [0.6, "rgb(120,70,20)"],
        [1, "rgb(255,255,255)"],
    ],
    "Electric": [
        [0, "rgb(0,0,0)"],
        [0.15, "rgb(30,0,100)"],
        [0.4, "rgb(120,0,100)"],
        [0.6, "rgb(160,90,0)"],
        [0.8, "rgb(230,200,0)"],
        [1, "rgb(255,250,220)"],
    ],
    "Viridis": [
        [0, "#440154"],
        [0.06274509803921569, "#48186a"],
        [0.12549019607843137, "#472d7b"],
        [0.18823529411764706, "#424086"],
        [0.25098039215686274, "#3b528b"],
        [0.3137254901960784, "#33638d"],
        [0.3764705882352941, "#2c728e"],
        [0.4392156862745098, "#26828e"],
        [0.5019607843137255, "#21918c"],
        [0.5647058823529412, "#1fa088"],
        [0.6274509803921569, "#28ae80"],
        [0.6901960784313725, "#3fbc73"],
        [0.7529411764705882, "#5ec962"],
        [0.8156862745098039, "#84d44b"],
        [0.8784313725490196, "#addc30"],
        [0.9411764705882353, "#d8e219"],
        [1, "#fde725"],
    ],
    "Cividis": [
        [0.000000, "rgb(0,32,76)"],
        [0.058824, "rgb(0,42,102)"],
        [0.117647, "rgb(0,52,110)"],
        [0.176471, "rgb(39,63,108)"],
        [0.235294, "rgb(60,74,107)"],
        [0.294118, "rgb(76,85,107)"],
        [0.352941, "rgb(91,95,109)"],
        [0.411765, "rgb(104,106,112)"],
        [0.470588, "rgb(117,117,117)"],
        [0.529412, "rgb(131,129,120)"],
        [0.588235, "rgb(146,140,120)"],
        [0.647059, "rgb(161,152,118)"],
        [0.705882, "rgb(176,165,114)"],
        [0.764706, "rgb(192,177,109)"],
        [0.823529, "rgb(209,191,102)"],
        [0.882353, "rgb(225,204,92)"],
        [0.941176, "rgb(243,219,79)"],
        [1.000000, "rgb(255,233,69)"],
    ],
}


def color_parser(colors, function):
    """
    Takes color(s) and a function and applies the function on the color(s)

    In particular, this function identifies whether the given color object
    is an iterable or not and applies the given color-parsing function to
    the color or iterable of colors. If given an iterable, it will only be
    able to work with it if all items in the iterable are of the same type
    - rgb string, hex string or tuple
    """
    if isinstance(colors, str):
        return function(colors)

    if isinstance(colors, tuple) and isinstance(colors[0], Number):
        return function(colors)

    if hasattr(colors, "__iter__"):
        if isinstance(colors, tuple):
            new_color_tuple = tuple(function(item) for item in colors)
            return new_color_tuple

        else:
            new_color_list = [function(item) for item in colors]
            return new_color_list


def validate_colors(colors, colortype="tuple"):
    """
    Validates color(s) and returns a list of color(s) of a specified type
    """
    from numbers import Number

    if colors is None:
        colors = DEFAULT_PLOTLY_COLORS

    if isinstance(colors, str):
        if colors in PLOTLY_SCALES:
            colors_list = colorscale_to_colors(PLOTLY_SCALES[colors])
            # TODO: fix _gantt.py/_scatter.py so that they can accept the
            # actual colorscale and not just a list of the first and last
            # color in the plotly colorscale. In resolving this issue we
            # will be removing the immediate line below
            colors = [colors_list[0]] + [colors_list[-1]]
        elif "rgb" in colors or "#" in colors:
            colors = [colors]
        else:
            raise exceptions.PlotlyError(
                "If your colors variable is a string, it must be a "
                "Plotly scale, an rgb color or a hex color."
            )

    elif isinstance(colors, tuple):
        if isinstance(colors[0], Number):
            colors = [colors]
        else:
            colors = list(colors)

    # convert color elements in list to tuple color
    for j, each_color in enumerate(colors):
        if "rgb" in each_color:
            each_color = color_parser(each_color, unlabel_rgb)
            for value in each_color:
                if value > 255.0:
                    raise exceptions.PlotlyError(
                        "Whoops! The elements in your rgb colors "
                        "tuples cannot exceed 255.0."
                    )
            each_color = color_parser(each_color, unconvert_from_RGB_255)
            colors[j] = each_color

        if "#" in each_color:
            each_color = color_parser(each_color, hex_to_rgb)
            each_color = color_parser(each_color, unconvert_from_RGB_255)

            colors[j] = each_color

        if isinstance(each_color, tuple):
            for value in each_color:
                if value > 1.0:
                    raise exceptions.PlotlyError(
                        "Whoops! The elements in your colors tuples cannot exceed 1.0."
                    )
            colors[j] = each_color

    if colortype == "rgb" and not isinstance(colors, str):
        for j, each_color in enumerate(colors):
            rgb_color = color_parser(each_color, convert_to_RGB_255)
            colors[j] = color_parser(rgb_color, label_rgb)

    return colors


def validate_colors_dict(colors, colortype="tuple"):
    """
    Validates dictionary of color(s)
    """
    # validate each color element in the dictionary
    for key in colors:
        if "rgb" in colors[key]:
            colors[key] = color_parser(colors[key], unlabel_rgb)
            for value in colors[key]:
                if value > 255.0:
                    raise exceptions.PlotlyError(
                        "Whoops! The elements in your rgb colors "
                        "tuples cannot exceed 255.0."
                    )
            colors[key] = color_parser(colors[key], unconvert_from_RGB_255)

        if "#" in colors[key]:
            colors[key] = color_parser(colors[key], hex_to_rgb)
            colors[key] = color_parser(colors[key], unconvert_from_RGB_255)

        if isinstance(colors[key], tuple):
            for value in colors[key]:
                if value > 1.0:
                    raise exceptions.PlotlyError(
                        "Whoops! The elements in your colors tuples cannot exceed 1.0."
                    )

    if colortype == "rgb":
        for key in colors:
            colors[key] = color_parser(colors[key], convert_to_RGB_255)
            colors[key] = color_parser(colors[key], label_rgb)

    return colors


def convert_colors_to_same_type(
    colors,
    colortype="rgb",
    scale=None,
    return_default_colors=False,
    num_of_defualt_colors=2,
):
    """
    Converts color(s) to the specified color type

    Takes a single color or an iterable of colors, as well as a list of scale
    values, and outputs a 2-pair of the list of color(s) converted all to an
    rgb or tuple color type, aswell as the scale as the second element. If
    colors is a Plotly Scale name, then 'scale' will be forced to the scale
    from the respective colorscale and the colors in that colorscale will also
    be coverted to the selected colortype. If colors is None, then there is an
    option to return portion of the DEFAULT_PLOTLY_COLORS

    :param (str|tuple|list) colors: either a plotly scale name, an rgb or hex
        color, a color tuple or a list/tuple of colors
    :param (list) scale: see docs for validate_scale_values()

    :rtype (tuple) (colors_list, scale) if scale is None in the function call,
        then scale will remain None in the returned tuple
    """
    colors_list = []

    if colors is None and return_default_colors is True:
        colors_list = DEFAULT_PLOTLY_COLORS[0:num_of_defualt_colors]

    if isinstance(colors, str):
        if colors in PLOTLY_SCALES:
            colors_list = colorscale_to_colors(PLOTLY_SCALES[colors])
            if scale is None:
                scale = colorscale_to_scale(PLOTLY_SCALES[colors])

        elif "rgb" in colors or "#" in colors:
            colors_list = [colors]

    elif isinstance(colors, tuple):
        if isinstance(colors[0], Number):
            colors_list = [colors]
        else:
            colors_list = list(colors)

    elif isinstance(colors, list):
        colors_list = colors

    # validate scale
    if scale is not None:
        validate_scale_values(scale)

        if len(colors_list) != len(scale):
            raise exceptions.PlotlyError(
                "Make sure that the length of your scale matches the length "
                "of your list of colors which is {}.".format(len(colors_list))
            )

    # convert all colors to rgb
    for j, each_color in enumerate(colors_list):
        if "#" in each_color:
            each_color = color_parser(each_color, hex_to_rgb)
            each_color = color_parser(each_color, label_rgb)
            colors_list[j] = each_color

        elif isinstance(each_color, tuple):
            each_color = color_parser(each_color, convert_to_RGB_255)
            each_color = color_parser(each_color, label_rgb)
            colors_list[j] = each_color

    if colortype == "rgb":
        return (colors_list, scale)
    elif colortype == "tuple":
        for j, each_color in enumerate(colors_list):
            each_color = color_parser(each_color, unlabel_rgb)
            each_color = color_parser(each_color, unconvert_from_RGB_255)
            colors_list[j] = each_color
        return (colors_list, scale)
    else:
        raise exceptions.PlotlyError(
            "You must select either rgb or tuple for your colortype variable."
        )


def convert_dict_colors_to_same_type(colors_dict, colortype="rgb"):
    """
    Converts a colors in a dictionary of colors to the specified color type

    :param (dict) colors_dict: a dictionary whose values are single colors
    """
    for key in colors_dict:
        if "#" in colors_dict[key]:
            colors_dict[key] = color_parser(colors_dict[key], hex_to_rgb)
            colors_dict[key] = color_parser(colors_dict[key], label_rgb)

        elif isinstance(colors_dict[key], tuple):
            colors_dict[key] = color_parser(colors_dict[key], convert_to_RGB_255)
            colors_dict[key] = color_parser(colors_dict[key], label_rgb)

    if colortype == "rgb":
        return colors_dict
    elif colortype == "tuple":
        for key in colors_dict:
            colors_dict[key] = color_parser(colors_dict[key], unlabel_rgb)
            colors_dict[key] = color_parser(colors_dict[key], unconvert_from_RGB_255)
        return colors_dict
    else:
        raise exceptions.PlotlyError(
            "You must select either rgb or tuple for your colortype variable."
        )


def validate_scale_values(scale):
    """
    Validates scale values from a colorscale

    :param (list) scale: a strictly increasing list of floats that begins
        with 0 and ends with 1. Its usage derives from a colorscale which is
        a list of two-lists (a list with two elements) of the form
        [value, color] which are used to determine how interpolation weighting
        works between the colors in the colorscale. Therefore scale is just
        the extraction of these values from the two-lists in order
    """
    if len(scale) < 2:
        raise exceptions.PlotlyError(
            "You must input a list of scale values that has at least two values."
        )

    if (scale[0] != 0) or (scale[-1] != 1):
        raise exceptions.PlotlyError(
            "The first and last number in your scale must be 0.0 and 1.0 respectively."
        )

    if not all(x < y for x, y in zip(scale, scale[1:])):
        raise exceptions.PlotlyError(
            "'scale' must be a list that contains a strictly increasing "
            "sequence of numbers."
        )


def validate_colorscale(colorscale):
    """Validate the structure, scale values and colors of colorscale."""
    if not isinstance(colorscale, list):
        # TODO Write tests for these exceptions
        raise exceptions.PlotlyError("A valid colorscale must be a list.")
    if not all(isinstance(innerlist, list) for innerlist in colorscale):
        raise exceptions.PlotlyError("A valid colorscale must be a list of lists.")
    colorscale_colors = colorscale_to_colors(colorscale)
    scale_values = colorscale_to_scale(colorscale)

    validate_scale_values(scale_values)
    validate_colors(colorscale_colors)


def make_colorscale(colors, scale=None):
    """
    Makes a colorscale from a list of colors and a scale

    Takes a list of colors and scales and constructs a colorscale based
    on the colors in sequential order. If 'scale' is left empty, a linear-
    interpolated colorscale will be generated. If 'scale' is a specificed
    list, it must be the same legnth as colors and must contain all floats
    For documentation regarding to the form of the output, see
    https://plot.ly/python/reference/#mesh3d-colorscale

    :param (list) colors: a list of single colors
    """
    colorscale = []

    # validate minimum colors length of 2
    if len(colors) < 2:
        raise exceptions.PlotlyError(
            "You must input a list of colors that has at least two colors."
        )

    if scale is None:
        scale_incr = 1.0 / (len(colors) - 1)
        return [[i * scale_incr, color] for i, color in enumerate(colors)]

    else:
        if len(colors) != len(scale):
            raise exceptions.PlotlyError(
                "The length of colors and scale must be the same."
            )

        validate_scale_values(scale)

        colorscale = [list(tup) for tup in zip(scale, colors)]
        return colorscale


def find_intermediate_color(lowcolor, highcolor, intermed, colortype="tuple"):
    """
    Returns the color at a given distance between two colors

    This function takes two color tuples, where each element is between 0
    and 1, along with a value 0 < intermed < 1 and returns a color that is
    intermed-percent from lowcolor to highcolor. If colortype is set to 'rgb',
    the function will automatically convert the rgb type to a tuple, find the
    intermediate color and return it as an rgb color.
    """
    if colortype == "rgb":
        # convert to tuple color, eg. (1, 0.45, 0.7)
        lowcolor = unlabel_rgb(lowcolor)
        highcolor = unlabel_rgb(highcolor)

    diff_0 = float(highcolor[0] - lowcolor[0])
    diff_1 = float(highcolor[1] - lowcolor[1])
    diff_2 = float(highcolor[2] - lowcolor[2])

    inter_med_tuple = (
        lowcolor[0] + intermed * diff_0,
        lowcolor[1] + intermed * diff_1,
        lowcolor[2] + intermed * diff_2,
    )

    if colortype == "rgb":
        # back to an rgb string, e.g. rgb(30, 20, 10)
        inter_med_rgb = label_rgb(inter_med_tuple)
        return inter_med_rgb

    return inter_med_tuple


def unconvert_from_RGB_255(colors):
    """
    Return a tuple where each element gets divided by 255

    Takes a (list of) color tuple(s) where each element is between 0 and
    255. Returns the same tuples where each tuple element is normalized to
    a value between 0 and 1
    """
    return (colors[0] / (255.0), colors[1] / (255.0), colors[2] / (255.0))


def convert_to_RGB_255(colors):
    """
    Multiplies each element of a triplet by 255

    Each coordinate of the color tuple is rounded to the nearest float and
    then is turned into an integer. If a number is of the form x.5, then
    if x is odd, the number rounds up to (x+1). Otherwise, it rounds down
    to just x. This is the way rounding works in Python 3 and in current
    statistical analysis to avoid rounding bias

    :param (list) rgb_components: grabs the three R, G and B values to be
        returned as computed in the function
    """
    rgb_components = []

    for component in colors:
        rounded_num = decimal.Decimal(str(component * 255.0)).quantize(
            decimal.Decimal("1"), rounding=decimal.ROUND_HALF_EVEN
        )
        # convert rounded number to an integer from 'Decimal' form
        rounded_num = int(rounded_num)
        rgb_components.append(rounded_num)

    return (rgb_components[0], rgb_components[1], rgb_components[2])


def n_colors(lowcolor, highcolor, n_colors, colortype="tuple"):
    """
    Splits a low and high color into a list of n_colors colors in it

    Accepts two color tuples and returns a list of n_colors colors
    which form the intermediate colors between lowcolor and highcolor
    from linearly interpolating through RGB space. If colortype is 'rgb'
    the function will return a list of colors in the same form.
    """
    if colortype == "rgb":
        # convert to tuple
        lowcolor = unlabel_rgb(lowcolor)
        highcolor = unlabel_rgb(highcolor)

    diff_0 = float(highcolor[0] - lowcolor[0])
    incr_0 = diff_0 / (n_colors - 1)
    diff_1 = float(highcolor[1] - lowcolor[1])
    incr_1 = diff_1 / (n_colors - 1)
    diff_2 = float(highcolor[2] - lowcolor[2])
    incr_2 = diff_2 / (n_colors - 1)
    list_of_colors = []

    def _constrain_color(c):
        if c > 255.0:
            return 255.0
        elif c < 0.0:
            return 0.0
        else:
            return c

    for index in range(n_colors):
        new_tuple = (
            _constrain_color(lowcolor[0] + (index * incr_0)),
            _constrain_color(lowcolor[1] + (index * incr_1)),
            _constrain_color(lowcolor[2] + (index * incr_2)),
        )
        list_of_colors.append(new_tuple)

    if colortype == "rgb":
        # back to an rgb string
        list_of_colors = color_parser(list_of_colors, label_rgb)

    return list_of_colors


def label_rgb(colors):
    """
    Takes tuple (a, b, c) and returns an rgb color 'rgb(a, b, c)'
    """
    return "rgb(%s, %s, %s)" % (colors[0], colors[1], colors[2])


def unlabel_rgb(colors):
    """
    Takes rgb color(s) 'rgb(a, b, c)' and returns tuple(s) (a, b, c)

    This function takes either an 'rgb(a, b, c)' color or a list of
    such colors and returns the color tuples in tuple(s) (a, b, c)
    """
    str_vals = ""
    for index in range(len(colors)):
        try:
            float(colors[index])
            str_vals = str_vals + colors[index]
        except ValueError:
            if colors[index] == "," or colors[index] == ".":
                str_vals = str_vals + colors[index]

    str_vals = str_vals + ","
    numbers = []
    str_num = ""
    for char in str_vals:
        if char != ",":
            str_num = str_num + char
        else:
            numbers.append(float(str_num))
            str_num = ""
    return (numbers[0], numbers[1], numbers[2])


def hex_to_rgb(value):
    """
    Calculates rgb values from a hex color code.

    :param (string) value: Hex color string

    :rtype (tuple) (r_value, g_value, b_value): tuple of rgb values
    """
    value = value.lstrip("#")
    hex_total_length = len(value)
    rgb_section_length = hex_total_length // 3
    return tuple(
        int(value[i : i + rgb_section_length], 16)
        for i in range(0, hex_total_length, rgb_section_length)
    )


def colorscale_to_colors(colorscale):
    """
    Extracts the colors from colorscale as a list
    """
    color_list = []
    for item in colorscale:
        color_list.append(item[1])
    return color_list


def colorscale_to_scale(colorscale):
    """
    Extracts the interpolation scale values from colorscale as a list
    """
    scale_list = []
    for item in colorscale:
        scale_list.append(item[0])
    return scale_list


def convert_colorscale_to_rgb(colorscale):
    """
    Converts the colors in a colorscale to rgb colors

    A colorscale is an array of arrays, each with a numeric value as the
    first item and a color as the second. This function specifically is
    converting a colorscale with tuple colors (each coordinate between 0
    and 1) into a colorscale with the colors transformed into rgb colors
    """
    for color in colorscale:
        color[1] = convert_to_RGB_255(color[1])

    for color in colorscale:
        color[1] = label_rgb(color[1])
    return colorscale


def named_colorscales():
    """
    Returns lowercased names of built-in continuous colorscales.
    """
    from _plotly_utils.basevalidators import ColorscaleValidator

    return [c for c in ColorscaleValidator("", "").named_colorscales]


def get_colorscale(name):
    """
    Returns the colorscale for a given name. See `named_colorscales` for the
    built-in colorscales.
    """
    from _plotly_utils.basevalidators import ColorscaleValidator

    if not isinstance(name, str):
        raise exceptions.PlotlyError("Name argument have to be a string.")

    name = name.lower()
    if name[-2:] == "_r":
        should_reverse = True
        name = name[:-2]
    else:
        should_reverse = False

    if name in ColorscaleValidator("", "").named_colorscales:
        colorscale = ColorscaleValidator("", "").named_colorscales[name]
    else:
        raise exceptions.PlotlyError(f"Colorscale {name} is not a built-in scale.")

    if should_reverse:
        colorscale = colorscale[::-1]
    return make_colorscale(colorscale)


def sample_colorscale(colorscale, samplepoints, low=0.0, high=1.0, colortype="rgb"):
    """
    Samples a colorscale at specific points.

    Interpolates between colors in a colorscale to find the specific colors
    corresponding to the specified sample values. The colorscale can be specified
    as a list of `[scale, color]` pairs, as a list of colors, or as a named
    plotly colorscale. The samplepoints can be specefied as an iterable of specific
    points in the range [0.0, 1.0], or as an integer number of points which will
    be spaced equally between the low value (default 0.0) and the high value
    (default 1.0). The output is a list of colors, formatted according to the
    specified colortype.
    """
    from bisect import bisect_left

    try:
        validate_colorscale(colorscale)
    except exceptions.PlotlyError:
        if isinstance(colorscale, str):
            colorscale = get_colorscale(colorscale)
        else:
            colorscale = make_colorscale(colorscale)

    scale = colorscale_to_scale(colorscale)
    validate_scale_values(scale)
    colors = colorscale_to_colors(colorscale)
    colors = validate_colors(colors, colortype="tuple")

    if isinstance(samplepoints, int):
        samplepoints = [
            low + idx / (samplepoints - 1) * (high - low) for idx in range(samplepoints)
        ]
    elif isinstance(samplepoints, float):
        samplepoints = [samplepoints]

    sampled_colors = []
    for point in samplepoints:
        high = bisect_left(scale, point)
        low = high - 1
        interpolant = (point - scale[low]) / (scale[high] - scale[low])
        sampled_color = find_intermediate_color(colors[low], colors[high], interpolant)
        sampled_colors.append(sampled_color)
    return validate_colors(sampled_colors, colortype=colortype)
