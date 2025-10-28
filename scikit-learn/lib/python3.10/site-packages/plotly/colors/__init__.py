# ruff: noqa: F405

"""For a list of colors available in `plotly.colors`, please see

* the `tutorial on discrete color sequences <https://plotly.com/python/discrete-color/#color-sequences-in-plotly-express>`_
* the `list of built-in continuous color scales <https://plotly.com/python/builtin-colorscales/>`_
* the `tutorial on continuous colors <https://plotly.com/python/colorscales/>`_

Color scales and sequences are available within the following namespaces

* cyclical
* diverging
* qualitative
* sequential
"""

from _plotly_utils.colors import *  # noqa: F403

__all__ = [
    "named_colorscales",
    "cyclical",
    "diverging",
    "sequential",
    "qualitative",
    "colorbrewer",
    "carto",
    "cmocean",
    "color_parser",
    "colorscale_to_colors",
    "colorscale_to_scale",
    "convert_colors_to_same_type",
    "convert_colorscale_to_rgb",
    "convert_dict_colors_to_same_type",
    "convert_to_RGB_255",
    "find_intermediate_color",
    "hex_to_rgb",
    "label_rgb",
    "make_colorscale",
    "n_colors",
    "sample_colorscale",
    "unconvert_from_RGB_255",
    "unlabel_rgb",
    "validate_colors",
    "validate_colors_dict",
    "validate_colorscale",
    "validate_scale_values",
    "plotlyjs",
    "DEFAULT_PLOTLY_COLORS",
    "PLOTLY_SCALES",
    "get_colorscale",
]
