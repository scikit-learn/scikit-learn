# Default settings for image generation


class _Defaults(object):
    """
    Class to store default settings for image generation.
    """

    def __init__(self):
        self.default_format = "png"
        self.default_width = 700
        self.default_height = 500
        self.default_scale = 1
        self.mathjax = None
        self.topojson = None
        self.plotlyjs = None


defaults = _Defaults()
