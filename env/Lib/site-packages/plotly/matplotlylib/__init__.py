"""
matplotlylib
============

This module converts matplotlib figure objects into JSON structures which can
be understood and visualized by Plotly.

Most of the functionality should be accessed through the parent directory's
'tools' module or 'plotly' package.

"""
from plotly.matplotlylib.renderer import PlotlyRenderer
from plotly.matplotlylib.mplexporter import Exporter
