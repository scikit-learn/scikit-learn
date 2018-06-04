
import six

from .backend_cairo import cairo, FigureCanvasCairo, RendererCairo
from .backend_qt5 import QtCore, QtGui, _BackendQT5, FigureCanvasQT
from .qt_compat import QT_API


class FigureCanvasQTCairo(FigureCanvasQT, FigureCanvasCairo):
    def __init__(self, figure):
        super(FigureCanvasQTCairo, self).__init__(figure=figure)
        self._renderer = RendererCairo(self.figure.dpi)
        self._renderer.set_width_height(-1, -1)  # Invalid values.

    def draw(self):
        if hasattr(self._renderer.gc, "ctx"):
            self.figure.draw(self._renderer)
        super(FigureCanvasQTCairo, self).draw()

    def paintEvent(self, event):
        self._update_dpi()
        dpi_ratio = self._dpi_ratio
        width = dpi_ratio * self.width()
        height = dpi_ratio * self.height()
        if (width, height) != self._renderer.get_canvas_width_height():
            surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
            self._renderer.set_ctx_from_surface(surface)
            self._renderer.set_width_height(width, height)
            self.figure.draw(self._renderer)
        buf = self._renderer.gc.ctx.get_target().get_data()
        qimage = QtGui.QImage(buf, width, height,
                              QtGui.QImage.Format_ARGB32_Premultiplied)
        # Adjust the buf reference count to work around a memory leak bug in
        # QImage under PySide on Python 3.
        if QT_API == 'PySide' and six.PY3:
            import ctypes
            ctypes.c_long.from_address(id(buf)).value = 1
        if hasattr(qimage, 'setDevicePixelRatio'):
            # Not available on Qt4 or some older Qt5.
            qimage.setDevicePixelRatio(dpi_ratio)
        painter = QtGui.QPainter(self)
        painter.drawImage(0, 0, qimage)
        self._draw_rect_callback(painter)
        painter.end()


@_BackendQT5.export
class _BackendQT5Cairo(_BackendQT5):
    FigureCanvas = FigureCanvasQTCairo
