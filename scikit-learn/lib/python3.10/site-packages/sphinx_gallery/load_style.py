"""Only load CSS and modify html_static_path.

This should not be used at the same time as sphinx_gallery.gen_gallery.
"""

from . import __version__, glr_path_static
from .directives import ImageSg, imagesg_addnode


def config_inited(app, config):
    """Append path to packaged static files to `html_static_path`."""
    path = glr_path_static()
    if path not in config.html_static_path:
        config.html_static_path.append(path)
    app.add_css_file("sg_gallery.css")


def setup(app):
    """Sphinx setup."""
    app.require_sphinx("1.8")
    app.connect("config-inited", config_inited)
    app.add_directive("image-sg", ImageSg)
    imagesg_addnode(app)
    return {
        "parallel_read_safe": True,
        "version": __version__,
    }
