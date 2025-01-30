from .base import Renderer


class FakeRenderer(Renderer):
    """
    Fake Renderer

    This is a fake renderer which simply outputs a text tree representing the
    elements found in the plot(s).  This is used in the unit tests for the
    package.

    Below are the methods your renderer must implement. You are free to do
    anything you wish within the renderer (i.e. build an XML or JSON
    representation, call an external API, etc.)  Here the renderer just
    builds a simple string representation for testing purposes.
    """

    def __init__(self):
        self.output = ""

    def open_figure(self, fig, props):
        self.output += "opening figure\n"

    def close_figure(self, fig):
        self.output += "closing figure\n"

    def open_axes(self, ax, props):
        self.output += "  opening axes\n"

    def close_axes(self, ax):
        self.output += "  closing axes\n"

    def open_legend(self, legend, props):
        self.output += "    opening legend\n"

    def close_legend(self, legend):
        self.output += "    closing legend\n"

    def draw_text(
        self, text, position, coordinates, style, text_type=None, mplobj=None
    ):
        self.output += "    draw text '{0}' {1}\n".format(text, text_type)

    def draw_path(
        self,
        data,
        coordinates,
        pathcodes,
        style,
        offset=None,
        offset_coordinates="data",
        mplobj=None,
    ):
        self.output += "    draw path with {0} vertices\n".format(data.shape[0])

    def draw_image(self, imdata, extent, coordinates, style, mplobj=None):
        self.output += "    draw image of size {0}\n".format(len(imdata))


class FullFakeRenderer(FakeRenderer):
    """
    Renderer with the full complement of methods.

    When the following are left undefined, they will be implemented via
    other methods in the class.  They can be defined explicitly for
    more efficient or specialized use within the renderer implementation.
    """

    def draw_line(self, data, coordinates, style, label, mplobj=None):
        self.output += "    draw line with {0} points\n".format(data.shape[0])

    def draw_markers(self, data, coordinates, style, label, mplobj=None):
        self.output += "    draw {0} markers\n".format(data.shape[0])

    def draw_path_collection(
        self,
        paths,
        path_coordinates,
        path_transforms,
        offsets,
        offset_coordinates,
        offset_order,
        styles,
        mplobj=None,
    ):
        self.output += "    draw path collection " "with {0} offsets\n".format(
            offsets.shape[0]
        )
