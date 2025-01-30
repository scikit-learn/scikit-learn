"""
Matplotlib Renderers
====================
This submodule contains renderer objects which define renderer behavior used
within the Exporter class.  The base renderer class is :class:`Renderer`, an
abstract base class
"""

from .base import Renderer
from .vega_renderer import VegaRenderer, fig_to_vega
from .vincent_renderer import VincentRenderer, fig_to_vincent
from .fake_renderer import FakeRenderer, FullFakeRenderer
