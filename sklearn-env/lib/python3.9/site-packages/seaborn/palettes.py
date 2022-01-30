import colorsys
from itertools import cycle

import numpy as np
import matplotlib as mpl

from .external import husl

from .utils import desaturate, get_color_cycle
from .colors import xkcd_rgb, crayons


__all__ = ["color_palette", "hls_palette", "husl_palette", "mpl_palette",
           "dark_palette", "light_palette", "diverging_palette",
           "blend_palette", "xkcd_palette", "crayon_palette",
           "cubehelix_palette", "set_color_codes"]


SEABORN_PALETTES = dict(
    deep=["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3",
          "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD"],
    deep6=["#4C72B0", "#55A868", "#C44E52",
           "#8172B3", "#CCB974", "#64B5CD"],
    muted=["#4878D0", "#EE854A", "#6ACC64", "#D65F5F", "#956CB4",
           "#8C613C", "#DC7EC0", "#797979", "#D5BB67", "#82C6E2"],
    muted6=["#4878D0", "#6ACC64", "#D65F5F",
            "#956CB4", "#D5BB67", "#82C6E2"],
    pastel=["#A1C9F4", "#FFB482", "#8DE5A1", "#FF9F9B", "#D0BBFF",
            "#DEBB9B", "#FAB0E4", "#CFCFCF", "#FFFEA3", "#B9F2F0"],
    pastel6=["#A1C9F4", "#8DE5A1", "#FF9F9B",
             "#D0BBFF", "#FFFEA3", "#B9F2F0"],
    bright=["#023EFF", "#FF7C00", "#1AC938", "#E8000B", "#8B2BE2",
            "#9F4800", "#F14CC1", "#A3A3A3", "#FFC400", "#00D7FF"],
    bright6=["#023EFF", "#1AC938", "#E8000B",
             "#8B2BE2", "#FFC400", "#00D7FF"],
    dark=["#001C7F", "#B1400D", "#12711C", "#8C0800", "#591E71",
          "#592F0D", "#A23582", "#3C3C3C", "#B8850A", "#006374"],
    dark6=["#001C7F", "#12711C", "#8C0800",
           "#591E71", "#B8850A", "#006374"],
    colorblind=["#0173B2", "#DE8F05", "#029E73", "#D55E00", "#CC78BC",
                "#CA9161", "#FBAFE4", "#949494", "#ECE133", "#56B4E9"],
    colorblind6=["#0173B2", "#029E73", "#D55E00",
                 "#CC78BC", "#ECE133", "#56B4E9"]
)


MPL_QUAL_PALS = {
    "tab10": 10, "tab20": 20, "tab20b": 20, "tab20c": 20,
    "Set1": 9, "Set2": 8, "Set3": 12,
    "Accent": 8, "Paired": 12,
    "Pastel1": 9, "Pastel2": 8, "Dark2": 8,
}


QUAL_PALETTE_SIZES = MPL_QUAL_PALS.copy()
QUAL_PALETTE_SIZES.update({k: len(v) for k, v in SEABORN_PALETTES.items()})
QUAL_PALETTES = list(QUAL_PALETTE_SIZES.keys())


class _ColorPalette(list):
    """Set the color palette in a with statement, otherwise be a list."""
    def __enter__(self):
        """Open the context."""
        from .rcmod import set_palette
        self._orig_palette = color_palette()
        set_palette(self)
        return self

    def __exit__(self, *args):
        """Close the context."""
        from .rcmod import set_palette
        set_palette(self._orig_palette)

    def as_hex(self):
        """Return a color palette with hex codes instead of RGB values."""
        hex = [mpl.colors.rgb2hex(rgb) for rgb in self]
        return _ColorPalette(hex)

    def _repr_html_(self):
        """Rich display of the color palette in an HTML frontend."""
        s = 55
        n = len(self)
        html = f'<svg  width="{n * s}" height="{s}">'
        for i, c in enumerate(self.as_hex()):
            html += (
                f'<rect x="{i * s}" y="0" width="{s}" height="{s}" style="fill:{c};'
                'stroke-width:2;stroke:rgb(255,255,255)"/>'
            )
        html += '</svg>'
        return html


def color_palette(palette=None, n_colors=None, desat=None, as_cmap=False):
    """Return a list of colors or continuous colormap defining a palette.

    Possible ``palette`` values include:
        - Name of a seaborn palette (deep, muted, bright, pastel, dark, colorblind)
        - Name of matplotlib colormap
        - 'husl' or 'hls'
        - 'ch:<cubehelix arguments>'
        - 'light:<color>', 'dark:<color>', 'blend:<color>,<color>',
        - A sequence of colors in any format matplotlib accepts

    Calling this function with ``palette=None`` will return the current
    matplotlib color cycle.

    This function can also be used in a ``with`` statement to temporarily
    set the color cycle for a plot or set of plots.

    See the :ref:`tutorial <palette_tutorial>` for more information.

    Parameters
    ----------
    palette : None, string, or sequence, optional
        Name of palette or None to return current palette. If a sequence, input
        colors are used but possibly cycled and desaturated.
    n_colors : int, optional
        Number of colors in the palette. If ``None``, the default will depend
        on how ``palette`` is specified. Named palettes default to 6 colors,
        but grabbing the current palette or passing in a list of colors will
        not change the number of colors unless this is specified. Asking for
        more colors than exist in the palette will cause it to cycle. Ignored
        when ``as_cmap`` is True.
    desat : float, optional
        Proportion to desaturate each color by.
    as_cmap : bool
        If True, return a :class:`matplotlib.colors.Colormap`.

    Returns
    -------
    list of RGB tuples or :class:`matplotlib.colors.Colormap`

    See Also
    --------
    set_palette : Set the default color cycle for all plots.
    set_color_codes : Reassign color codes like ``"b"``, ``"g"``, etc. to
                      colors from one of the seaborn palettes.

    Examples
    --------

    .. include:: ../docstrings/color_palette.rst

    """
    if palette is None:
        palette = get_color_cycle()
        if n_colors is None:
            n_colors = len(palette)

    elif not isinstance(palette, str):
        palette = palette
        if n_colors is None:
            n_colors = len(palette)
    else:

        if n_colors is None:
            # Use all colors in a qualitative palette or 6 of another kind
            n_colors = QUAL_PALETTE_SIZES.get(palette, 6)

        if palette in SEABORN_PALETTES:
            # Named "seaborn variant" of matplotlib default color cycle
            palette = SEABORN_PALETTES[palette]

        elif palette == "hls":
            # Evenly spaced colors in cylindrical RGB space
            palette = hls_palette(n_colors, as_cmap=as_cmap)

        elif palette == "husl":
            # Evenly spaced colors in cylindrical Lab space
            palette = husl_palette(n_colors, as_cmap=as_cmap)

        elif palette.lower() == "jet":
            # Paternalism
            raise ValueError("No.")

        elif palette.startswith("ch:"):
            # Cubehelix palette with params specified in string
            args, kwargs = _parse_cubehelix_args(palette)
            palette = cubehelix_palette(n_colors, *args, **kwargs, as_cmap=as_cmap)

        elif palette.startswith("light:"):
            # light palette to color specified in string
            _, color = palette.split(":")
            reverse = color.endswith("_r")
            if reverse:
                color = color[:-2]
            palette = light_palette(color, n_colors, reverse=reverse, as_cmap=as_cmap)

        elif palette.startswith("dark:"):
            # light palette to color specified in string
            _, color = palette.split(":")
            reverse = color.endswith("_r")
            if reverse:
                color = color[:-2]
            palette = dark_palette(color, n_colors, reverse=reverse, as_cmap=as_cmap)

        elif palette.startswith("blend:"):
            # blend palette between colors specified in string
            _, colors = palette.split(":")
            colors = colors.split(",")
            palette = blend_palette(colors, n_colors, as_cmap=as_cmap)

        else:
            try:
                # Perhaps a named matplotlib colormap?
                palette = mpl_palette(palette, n_colors, as_cmap=as_cmap)
            except ValueError:
                raise ValueError("%s is not a valid palette name" % palette)

    if desat is not None:
        palette = [desaturate(c, desat) for c in palette]

    if not as_cmap:

        # Always return as many colors as we asked for
        pal_cycle = cycle(palette)
        palette = [next(pal_cycle) for _ in range(n_colors)]

        # Always return in r, g, b tuple format
        try:
            palette = map(mpl.colors.colorConverter.to_rgb, palette)
            palette = _ColorPalette(palette)
        except ValueError:
            raise ValueError(f"Could not generate a palette for {palette}")

    return palette


def hls_palette(n_colors=6, h=.01, l=.6, s=.65, as_cmap=False):  # noqa
    """Get a set of evenly spaced colors in HLS hue space.

    h, l, and s should be between 0 and 1

    Parameters
    ----------

    n_colors : int
        number of colors in the palette
    h : float
        first hue
    l : float
        lightness
    s : float
        saturation

    Returns
    -------
    list of RGB tuples or :class:`matplotlib.colors.Colormap`

    See Also
    --------
    husl_palette : Make a palette using evenly spaced hues in the HUSL system.

    Examples
    --------

    Create a palette of 10 colors with the default parameters:

    .. plot::
        :context: close-figs

        >>> import seaborn as sns; sns.set_theme()
        >>> sns.palplot(sns.hls_palette(10))

    Create a palette of 10 colors that begins at a different hue value:

    .. plot::
        :context: close-figs

        >>> sns.palplot(sns.hls_palette(10, h=.5))

    Create a palette of 10 colors that are darker than the default:

    .. plot::
        :context: close-figs

        >>> sns.palplot(sns.hls_palette(10, l=.4))

    Create a palette of 10 colors that are less saturated than the default:

    .. plot::
        :context: close-figs

        >>> sns.palplot(sns.hls_palette(10, s=.4))

    """
    if as_cmap:
        n_colors = 256
    hues = np.linspace(0, 1, int(n_colors) + 1)[:-1]
    hues += h
    hues %= 1
    hues -= hues.astype(int)
    palette = [colorsys.hls_to_rgb(h_i, l, s) for h_i in hues]
    if as_cmap:
        return mpl.colors.ListedColormap(palette, "hls")
    else:
        return _ColorPalette(palette)


def husl_palette(n_colors=6, h=.01, s=.9, l=.65, as_cmap=False):  # noqa
    """Get a set of evenly spaced colors in HUSL hue space.

    h, s, and l should be between 0 and 1

    Parameters
    ----------

    n_colors : int
        number of colors in the palette
    h : float
        first hue
    s : float
        saturation
    l : float
        lightness

    Returns
    -------
    list of RGB tuples or :class:`matplotlib.colors.Colormap`

    See Also
    --------
    hls_palette : Make a palette using evently spaced circular hues in the
                  HSL system.

    Examples
    --------

    Create a palette of 10 colors with the default parameters:

    .. plot::
        :context: close-figs

        >>> import seaborn as sns; sns.set_theme()
        >>> sns.palplot(sns.husl_palette(10))

    Create a palette of 10 colors that begins at a different hue value:

    .. plot::
        :context: close-figs

        >>> sns.palplot(sns.husl_palette(10, h=.5))

    Create a palette of 10 colors that are darker than the default:

    .. plot::
        :context: close-figs

        >>> sns.palplot(sns.husl_palette(10, l=.4))

    Create a palette of 10 colors that are less saturated than the default:

    .. plot::
        :context: close-figs

        >>> sns.palplot(sns.husl_palette(10, s=.4))

    """
    if as_cmap:
        n_colors = 256
    hues = np.linspace(0, 1, int(n_colors) + 1)[:-1]
    hues += h
    hues %= 1
    hues *= 359
    s *= 99
    l *= 99  # noqa
    palette = [_color_to_rgb((h_i, s, l), input="husl") for h_i in hues]
    if as_cmap:
        return mpl.colors.ListedColormap(palette, "hsl")
    else:
        return _ColorPalette(palette)


def mpl_palette(name, n_colors=6, as_cmap=False):
    """Return discrete colors from a matplotlib palette.

    Note that this handles the qualitative colorbrewer palettes
    properly, although if you ask for more colors than a particular
    qualitative palette can provide you will get fewer than you are
    expecting. In contrast, asking for qualitative color brewer palettes
    using :func:`color_palette` will return the expected number of colors,
    but they will cycle.

    If you are using the IPython notebook, you can also use the function
    :func:`choose_colorbrewer_palette` to interactively select palettes.

    Parameters
    ----------
    name : string
        Name of the palette. This should be a named matplotlib colormap.
    n_colors : int
        Number of discrete colors in the palette.

    Returns
    -------
    list of RGB tuples or :class:`matplotlib.colors.Colormap`

    Examples
    --------

    Create a qualitative colorbrewer palette with 8 colors:

    .. plot::
        :context: close-figs

        >>> import seaborn as sns; sns.set_theme()
        >>> sns.palplot(sns.mpl_palette("Set2", 8))

    Create a sequential colorbrewer palette:

    .. plot::
        :context: close-figs

        >>> sns.palplot(sns.mpl_palette("Blues"))

    Create a diverging palette:

    .. plot::
        :context: close-figs

        >>> sns.palplot(sns.mpl_palette("seismic", 8))

    Create a "dark" sequential palette:

    .. plot::
        :context: close-figs

        >>> sns.palplot(sns.mpl_palette("GnBu_d"))

    """
    if name.endswith("_d"):
        sub_name = name[:-2]
        if sub_name.endswith("_r"):
            reverse = True
            sub_name = sub_name[:-2]
        else:
            reverse = False
        pal = color_palette(sub_name, 2) + ["#333333"]
        if reverse:
            pal = pal[::-1]
        cmap = blend_palette(pal, n_colors, as_cmap=True)
    else:
        cmap = mpl.cm.get_cmap(name)

    if name in MPL_QUAL_PALS:
        bins = np.linspace(0, 1, MPL_QUAL_PALS[name])[:n_colors]
    else:
        bins = np.linspace(0, 1, int(n_colors) + 2)[1:-1]
    palette = list(map(tuple, cmap(bins)[:, :3]))

    if as_cmap:
        return cmap
    else:
        return _ColorPalette(palette)


def _color_to_rgb(color, input):
    """Add some more flexibility to color choices."""
    if input == "hls":
        color = colorsys.hls_to_rgb(*color)
    elif input == "husl":
        color = husl.husl_to_rgb(*color)
        color = tuple(np.clip(color, 0, 1))
    elif input == "xkcd":
        color = xkcd_rgb[color]

    return mpl.colors.to_rgb(color)


def dark_palette(color, n_colors=6, reverse=False, as_cmap=False, input="rgb"):
    """Make a sequential palette that blends from dark to ``color``.

    This kind of palette is good for data that range between relatively
    uninteresting low values and interesting high values.

    The ``color`` parameter can be specified in a number of ways, including
    all options for defining a color in matplotlib and several additional
    color spaces that are handled by seaborn. You can also use the database
    of named colors from the XKCD color survey.

    If you are using the IPython notebook, you can also choose this palette
    interactively with the :func:`choose_dark_palette` function.

    Parameters
    ----------
    color : base color for high values
        hex, rgb-tuple, or html color name
    n_colors : int, optional
        number of colors in the palette
    reverse : bool, optional
        if True, reverse the direction of the blend
    as_cmap : bool, optional
        If True, return a :class:`matplotlib.colors.Colormap`.
    input : {'rgb', 'hls', 'husl', xkcd'}
        Color space to interpret the input color. The first three options
        apply to tuple inputs and the latter applies to string inputs.

    Returns
    -------
    list of RGB tuples or :class:`matplotlib.colors.Colormap`

    See Also
    --------
    light_palette : Create a sequential palette with bright low values.
    diverging_palette : Create a diverging palette with two colors.

    Examples
    --------

    Generate a palette from an HTML color:

    .. plot::
        :context: close-figs

        >>> import seaborn as sns; sns.set_theme()
        >>> sns.palplot(sns.dark_palette("purple"))

    Generate a palette that decreases in lightness:

    .. plot::
        :context: close-figs

        >>> sns.palplot(sns.dark_palette("seagreen", reverse=True))

    Generate a palette from an HUSL-space seed:

    .. plot::
        :context: close-figs

        >>> sns.palplot(sns.dark_palette((260, 75, 60), input="husl"))

    Generate a colormap object:

    .. plot::
        :context: close-figs

        >>> from numpy import arange
        >>> x = arange(25).reshape(5, 5)
        >>> cmap = sns.dark_palette("#2ecc71", as_cmap=True)
        >>> ax = sns.heatmap(x, cmap=cmap)

    """
    rgb = _color_to_rgb(color, input)
    h, s, l = husl.rgb_to_husl(*rgb)
    gray_s, gray_l = .15 * s, 15
    gray = _color_to_rgb((h, gray_s, gray_l), input="husl")
    colors = [rgb, gray] if reverse else [gray, rgb]
    return blend_palette(colors, n_colors, as_cmap)


def light_palette(color, n_colors=6, reverse=False, as_cmap=False, input="rgb"):
    """Make a sequential palette that blends from light to ``color``.

    This kind of palette is good for data that range between relatively
    uninteresting low values and interesting high values.

    The ``color`` parameter can be specified in a number of ways, including
    all options for defining a color in matplotlib and several additional
    color spaces that are handled by seaborn. You can also use the database
    of named colors from the XKCD color survey.

    If you are using the IPython notebook, you can also choose this palette
    interactively with the :func:`choose_light_palette` function.

    Parameters
    ----------
    color : base color for high values
        hex code, html color name, or tuple in ``input`` space.
    n_colors : int, optional
        number of colors in the palette
    reverse : bool, optional
        if True, reverse the direction of the blend
    as_cmap : bool, optional
        If True, return a :class:`matplotlib.colors.Colormap`.
    input : {'rgb', 'hls', 'husl', xkcd'}
        Color space to interpret the input color. The first three options
        apply to tuple inputs and the latter applies to string inputs.

    Returns
    -------
    list of RGB tuples or :class:`matplotlib.colors.Colormap`

    See Also
    --------
    dark_palette : Create a sequential palette with dark low values.
    diverging_palette : Create a diverging palette with two colors.

    Examples
    --------

    Generate a palette from an HTML color:

    .. plot::
        :context: close-figs

        >>> import seaborn as sns; sns.set_theme()
        >>> sns.palplot(sns.light_palette("purple"))

    Generate a palette that increases in lightness:

    .. plot::
        :context: close-figs

        >>> sns.palplot(sns.light_palette("seagreen", reverse=True))

    Generate a palette from an HUSL-space seed:

    .. plot::
        :context: close-figs

        >>> sns.palplot(sns.light_palette((260, 75, 60), input="husl"))

    Generate a colormap object:

    .. plot::
        :context: close-figs

        >>> from numpy import arange
        >>> x = arange(25).reshape(5, 5)
        >>> cmap = sns.light_palette("#2ecc71", as_cmap=True)
        >>> ax = sns.heatmap(x, cmap=cmap)

    """
    rgb = _color_to_rgb(color, input)
    h, s, l = husl.rgb_to_husl(*rgb)
    gray_s, gray_l = .15 * s, 95
    gray = _color_to_rgb((h, gray_s, gray_l), input="husl")
    colors = [rgb, gray] if reverse else [gray, rgb]
    return blend_palette(colors, n_colors, as_cmap)


def diverging_palette(h_neg, h_pos, s=75, l=50, sep=1, n=6,  # noqa
                      center="light", as_cmap=False):
    """Make a diverging palette between two HUSL colors.

    If you are using the IPython notebook, you can also choose this palette
    interactively with the :func:`choose_diverging_palette` function.

    Parameters
    ----------
    h_neg, h_pos : float in [0, 359]
        Anchor hues for negative and positive extents of the map.
    s : float in [0, 100], optional
        Anchor saturation for both extents of the map.
    l : float in [0, 100], optional
        Anchor lightness for both extents of the map.
    sep : int, optional
        Size of the intermediate region.
    n : int, optional
        Number of colors in the palette (if not returning a cmap)
    center : {"light", "dark"}, optional
        Whether the center of the palette is light or dark
    as_cmap : bool, optional
        If True, return a :class:`matplotlib.colors.Colormap`.

    Returns
    -------
    list of RGB tuples or :class:`matplotlib.colors.Colormap`

    See Also
    --------
    dark_palette : Create a sequential palette with dark values.
    light_palette : Create a sequential palette with light values.

    Examples
    --------

    Generate a blue-white-red palette:

    .. plot::
        :context: close-figs

        >>> import seaborn as sns; sns.set_theme()
        >>> sns.palplot(sns.diverging_palette(240, 10, n=9))

    Generate a brighter green-white-purple palette:

    .. plot::
        :context: close-figs

        >>> sns.palplot(sns.diverging_palette(150, 275, s=80, l=55, n=9))

    Generate a blue-black-red palette:

    .. plot::
        :context: close-figs

        >>> sns.palplot(sns.diverging_palette(250, 15, s=75, l=40,
        ...                                   n=9, center="dark"))

    Generate a colormap object:

    .. plot::
        :context: close-figs

        >>> from numpy import arange
        >>> x = arange(25).reshape(5, 5)
        >>> cmap = sns.diverging_palette(220, 20, as_cmap=True)
        >>> ax = sns.heatmap(x, cmap=cmap)

    """
    palfunc = dict(dark=dark_palette, light=light_palette)[center]
    n_half = int(128 - (sep // 2))
    neg = palfunc((h_neg, s, l), n_half, reverse=True, input="husl")
    pos = palfunc((h_pos, s, l), n_half, input="husl")
    midpoint = dict(light=[(.95, .95, .95)], dark=[(.133, .133, .133)])[center]
    mid = midpoint * sep
    pal = blend_palette(np.concatenate([neg, mid, pos]), n, as_cmap=as_cmap)
    return pal


def blend_palette(colors, n_colors=6, as_cmap=False, input="rgb"):
    """Make a palette that blends between a list of colors.

    Parameters
    ----------
    colors : sequence of colors in various formats interpreted by ``input``
        hex code, html color name, or tuple in ``input`` space.
    n_colors : int, optional
        Number of colors in the palette.
    as_cmap : bool, optional
        If True, return a :class:`matplotlib.colors.Colormap`.

    Returns
    -------
    list of RGB tuples or :class:`matplotlib.colors.Colormap`

    """
    colors = [_color_to_rgb(color, input) for color in colors]
    name = "blend"
    pal = mpl.colors.LinearSegmentedColormap.from_list(name, colors)
    if not as_cmap:
        rgb_array = pal(np.linspace(0, 1, int(n_colors)))[:, :3]  # no alpha
        pal = _ColorPalette(map(tuple, rgb_array))
    return pal


def xkcd_palette(colors):
    """Make a palette with color names from the xkcd color survey.

    See xkcd for the full list of colors: https://xkcd.com/color/rgb/

    This is just a simple wrapper around the ``seaborn.xkcd_rgb`` dictionary.

    Parameters
    ----------
    colors : list of strings
        List of keys in the ``seaborn.xkcd_rgb`` dictionary.

    Returns
    -------
    palette : seaborn color palette
        Returns the list of colors as RGB tuples in an object that behaves like
        other seaborn color palettes.

    See Also
    --------
    crayon_palette : Make a palette with Crayola crayon colors.

    """
    palette = [xkcd_rgb[name] for name in colors]
    return color_palette(palette, len(palette))


def crayon_palette(colors):
    """Make a palette with color names from Crayola crayons.

    Colors are taken from here:
    https://en.wikipedia.org/wiki/List_of_Crayola_crayon_colors

    This is just a simple wrapper around the ``seaborn.crayons`` dictionary.

    Parameters
    ----------
    colors : list of strings
        List of keys in the ``seaborn.crayons`` dictionary.

    Returns
    -------
    palette : seaborn color palette
        Returns the list of colors as rgb tuples in an object that behaves like
        other seaborn color palettes.

    See Also
    --------
    xkcd_palette : Make a palette with named colors from the XKCD color survey.

    """
    palette = [crayons[name] for name in colors]
    return color_palette(palette, len(palette))


def cubehelix_palette(n_colors=6, start=0, rot=.4, gamma=1.0, hue=0.8,
                      light=.85, dark=.15, reverse=False, as_cmap=False):
    """Make a sequential palette from the cubehelix system.

    This produces a colormap with linearly-decreasing (or increasing)
    brightness. That means that information will be preserved if printed to
    black and white or viewed by someone who is colorblind.  "cubehelix" is
    also available as a matplotlib-based palette, but this function gives the
    user more control over the look of the palette and has a different set of
    defaults.

    In addition to using this function, it is also possible to generate a
    cubehelix palette generally in seaborn using a string-shorthand; see the
    example below.

    Parameters
    ----------
    n_colors : int
        Number of colors in the palette.
    start : float, 0 <= start <= 3
        The hue at the start of the helix.
    rot : float
        Rotations around the hue wheel over the range of the palette.
    gamma : float 0 <= gamma
        Gamma factor to emphasize darker (gamma < 1) or lighter (gamma > 1)
        colors.
    hue : float, 0 <= hue <= 1
        Saturation of the colors.
    dark : float 0 <= dark <= 1
        Intensity of the darkest color in the palette.
    light : float 0 <= light <= 1
        Intensity of the lightest color in the palette.
    reverse : bool
        If True, the palette will go from dark to light.
    as_cmap : bool
        If True, return a :class:`matplotlib.colors.Colormap`.

    Returns
    -------
    list of RGB tuples or :class:`matplotlib.colors.Colormap`

    See Also
    --------
    choose_cubehelix_palette : Launch an interactive widget to select cubehelix
                               palette parameters.
    dark_palette : Create a sequential palette with dark low values.
    light_palette : Create a sequential palette with bright low values.

    References
    ----------
    Green, D. A. (2011). "A colour scheme for the display of astronomical
    intensity images". Bulletin of the Astromical Society of India, Vol. 39,
    p. 289-295.

    Examples
    --------

    Generate the default palette:

    .. plot::
        :context: close-figs

        >>> import seaborn as sns; sns.set_theme()
        >>> sns.palplot(sns.cubehelix_palette())

    Rotate backwards from the same starting location:

    .. plot::
        :context: close-figs

        >>> sns.palplot(sns.cubehelix_palette(rot=-.4))

    Use a different starting point and shorter rotation:

    .. plot::
        :context: close-figs

        >>> sns.palplot(sns.cubehelix_palette(start=2.8, rot=.1))

    Reverse the direction of the lightness ramp:

    .. plot::
        :context: close-figs

        >>> sns.palplot(sns.cubehelix_palette(reverse=True))

    Generate a colormap object:

    .. plot::
        :context: close-figs

        >>> from numpy import arange
        >>> x = arange(25).reshape(5, 5)
        >>> cmap = sns.cubehelix_palette(as_cmap=True)
        >>> ax = sns.heatmap(x, cmap=cmap)

    Use the full lightness range:

    .. plot::
        :context: close-figs

        >>> cmap = sns.cubehelix_palette(dark=0, light=1, as_cmap=True)
        >>> ax = sns.heatmap(x, cmap=cmap)

    Use through the :func:`color_palette` interface:

    .. plot::
        :context: close-figs

        >>> sns.palplot(sns.color_palette("ch:2,r=.2,l=.6"))

    """
    def get_color_function(p0, p1):
        # Copied from matplotlib because it lives in private module
        def color(x):
            # Apply gamma factor to emphasise low or high intensity values
            xg = x ** gamma

            # Calculate amplitude and angle of deviation from the black
            # to white diagonal in the plane of constant
            # perceived intensity.
            a = hue * xg * (1 - xg) / 2

            phi = 2 * np.pi * (start / 3 + rot * x)

            return xg + a * (p0 * np.cos(phi) + p1 * np.sin(phi))
        return color

    cdict = {
        "red": get_color_function(-0.14861, 1.78277),
        "green": get_color_function(-0.29227, -0.90649),
        "blue": get_color_function(1.97294, 0.0),
    }

    cmap = mpl.colors.LinearSegmentedColormap("cubehelix", cdict)

    x = np.linspace(light, dark, int(n_colors))
    pal = cmap(x)[:, :3].tolist()
    if reverse:
        pal = pal[::-1]

    if as_cmap:
        x_256 = np.linspace(light, dark, 256)
        if reverse:
            x_256 = x_256[::-1]
        pal_256 = cmap(x_256)
        cmap = mpl.colors.ListedColormap(pal_256, "seaborn_cubehelix")
        return cmap
    else:
        return _ColorPalette(pal)


def _parse_cubehelix_args(argstr):
    """Turn stringified cubehelix params into args/kwargs."""

    if argstr.startswith("ch:"):
        argstr = argstr[3:]

    if argstr.endswith("_r"):
        reverse = True
        argstr = argstr[:-2]
    else:
        reverse = False

    if not argstr:
        return [], {"reverse": reverse}

    all_args = argstr.split(",")

    args = [float(a.strip(" ")) for a in all_args if "=" not in a]

    kwargs = [a.split("=") for a in all_args if "=" in a]
    kwargs = {k.strip(" "): float(v.strip(" ")) for k, v in kwargs}

    kwarg_map = dict(
        s="start", r="rot", g="gamma",
        h="hue", l="light", d="dark",  # noqa: E741
    )

    kwargs = {kwarg_map.get(k, k): v for k, v in kwargs.items()}

    if reverse:
        kwargs["reverse"] = True

    return args, kwargs


def set_color_codes(palette="deep"):
    """Change how matplotlib color shorthands are interpreted.

    Calling this will change how shorthand codes like "b" or "g"
    are interpreted by matplotlib in subsequent plots.

    Parameters
    ----------
    palette : {deep, muted, pastel, dark, bright, colorblind}
        Named seaborn palette to use as the source of colors.

    See Also
    --------
    set : Color codes can be set through the high-level seaborn style
          manager.
    set_palette : Color codes can also be set through the function that
                  sets the matplotlib color cycle.

    Examples
    --------

    Map matplotlib color codes to the default seaborn palette.

    .. plot::
        :context: close-figs

        >>> import matplotlib.pyplot as plt
        >>> import seaborn as sns; sns.set_theme()
        >>> sns.set_color_codes()
        >>> _ = plt.plot([0, 1], color="r")

    Use a different seaborn palette.

    .. plot::
        :context: close-figs

        >>> sns.set_color_codes("dark")
        >>> _ = plt.plot([0, 1], color="g")
        >>> _ = plt.plot([0, 2], color="m")

    """
    if palette == "reset":
        colors = [(0., 0., 1.), (0., .5, 0.), (1., 0., 0.), (.75, 0., .75),
                  (.75, .75, 0.), (0., .75, .75), (0., 0., 0.)]
    elif not isinstance(palette, str):
        err = "set_color_codes requires a named seaborn palette"
        raise TypeError(err)
    elif palette in SEABORN_PALETTES:
        if not palette.endswith("6"):
            palette = palette + "6"
        colors = SEABORN_PALETTES[palette] + [(.1, .1, .1)]
    else:
        err = "Cannot set colors with palette '{}'".format(palette)
        raise ValueError(err)

    for code, color in zip("bgrmyck", colors):
        rgb = mpl.colors.colorConverter.to_rgb(color)
        mpl.colors.colorConverter.colors[code] = rgb
        mpl.colors.colorConverter.cache[code] = rgb
