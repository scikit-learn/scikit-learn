"""
Displays Agg images in the browser, with interactivity
"""
# The WebAgg backend is divided into two modules:
#
# - `backend_webagg_core.py` contains code necessary to embed a WebAgg
#   plot inside of a web application, and communicate in an abstract
#   way over a web socket.
#
# - `backend_webagg.py` contains a concrete implementation of a basic
#   application, implemented with tornado.

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import datetime
import io
import json
import os
import warnings

import numpy as np
import tornado

from matplotlib.backends import backend_agg
from matplotlib.backend_bases import _Backend
from matplotlib import backend_bases
from matplotlib import _png


# http://www.cambiaresearch.com/articles/15/javascript-char-codes-key-codes
_SHIFT_LUT = {59: ':',
              61: '+',
              173: '_',
              186: ':',
              187: '+',
              188: '<',
              189: '_',
              190: '>',
              191: '?',
              192: '~',
              219: '{',
              220: '|',
              221: '}',
              222: '"'}

_LUT = {8: 'backspace',
        9: 'tab',
        13: 'enter',
        16: 'shift',
        17: 'control',
        18: 'alt',
        19: 'pause',
        20: 'caps',
        27: 'escape',
        32: ' ',
        33: 'pageup',
        34: 'pagedown',
        35: 'end',
        36: 'home',
        37: 'left',
        38: 'up',
        39: 'right',
        40: 'down',
        45: 'insert',
        46: 'delete',
        91: 'super',
        92: 'super',
        93: 'select',
        106: '*',
        107: '+',
        109: '-',
        110: '.',
        111: '/',
        144: 'num_lock',
        145: 'scroll_lock',
        186: ':',
        187: '=',
        188: ',',
        189: '-',
        190: '.',
        191: '/',
        192: '`',
        219: '[',
        220: '\\',
        221: ']',
        222: "'"}


def _handle_key(key):
    """Handle key codes"""
    code = int(key[key.index('k') + 1:])
    value = chr(code)
    # letter keys
    if code >= 65 and code <= 90:
        if 'shift+' in key:
            key = key.replace('shift+', '')
        else:
            value = value.lower()
    # number keys
    elif code >= 48 and code <= 57:
        if 'shift+' in key:
            value = ')!@#$%^&*('[int(value)]
            key = key.replace('shift+', '')
    # function keys
    elif code >= 112 and code <= 123:
        value = 'f%s' % (code - 111)
    # number pad keys
    elif code >= 96 and code <= 105:
        value = '%s' % (code - 96)
    # keys with shift alternatives
    elif code in _SHIFT_LUT and 'shift+' in key:
        key = key.replace('shift+', '')
        value = _SHIFT_LUT[code]
    elif code in _LUT:
        value = _LUT[code]
    key = key[:key.index('k')] + value
    return key


class FigureCanvasWebAggCore(backend_agg.FigureCanvasAgg):
    supports_blit = False

    def __init__(self, *args, **kwargs):
        backend_agg.FigureCanvasAgg.__init__(self, *args, **kwargs)

        # Set to True when the renderer contains data that is newer
        # than the PNG buffer.
        self._png_is_old = True

        # Set to True by the `refresh` message so that the next frame
        # sent to the clients will be a full frame.
        self._force_full = True

        # Store the current image mode so that at any point, clients can
        # request the information. This should be changed by calling
        # self.set_image_mode(mode) so that the notification can be given
        # to the connected clients.
        self._current_image_mode = 'full'

        # Store the DPI ratio of the browser.  This is the scaling that
        # occurs automatically for all images on a HiDPI display.
        self._dpi_ratio = 1

    def show(self):
        # show the figure window
        from matplotlib.pyplot import show
        show()

    def draw(self):
        renderer = self.get_renderer(cleared=True)

        self._png_is_old = True

        backend_agg.RendererAgg.lock.acquire()
        try:
            self.figure.draw(renderer)
        finally:
            backend_agg.RendererAgg.lock.release()
            # Swap the frames
            self.manager.refresh_all()

    def draw_idle(self):
        self.send_event("draw")

    def set_image_mode(self, mode):
        """
        Set the image mode for any subsequent images which will be sent
        to the clients. The modes may currently be either 'full' or 'diff'.

        Note: diff images may not contain transparency, therefore upon
        draw this mode may be changed if the resulting image has any
        transparent component.

        """
        if mode not in ['full', 'diff']:
            raise ValueError('image mode must be either full or diff.')
        if self._current_image_mode != mode:
            self._current_image_mode = mode
            self.handle_send_image_mode(None)

    def get_diff_image(self):
        if self._png_is_old:
            renderer = self.get_renderer()

            # The buffer is created as type uint32 so that entire
            # pixels can be compared in one numpy call, rather than
            # needing to compare each plane separately.
            buff = (np.frombuffer(renderer.buffer_rgba(), dtype=np.uint32)
                    .reshape((renderer.height, renderer.width)))

            # If any pixels have transparency, we need to force a full
            # draw as we cannot overlay new on top of old.
            pixels = buff.view(dtype=np.uint8).reshape(buff.shape + (4,))

            if self._force_full or np.any(pixels[:, :, 3] != 255):
                self.set_image_mode('full')
                output = buff
            else:
                self.set_image_mode('diff')
                last_buffer = (np.frombuffer(self._last_renderer.buffer_rgba(),
                                             dtype=np.uint32)
                               .reshape((renderer.height, renderer.width)))
                diff = buff != last_buffer
                output = np.where(diff, buff, 0)

            # TODO: We should write a new version of write_png that
            # handles the differencing inline
            buff = _png.write_png(
                output.view(dtype=np.uint8).reshape(output.shape + (4,)),
                None, compression=6, filter=_png.PNG_FILTER_NONE)

            # Swap the renderer frames
            self._renderer, self._last_renderer = (
                self._last_renderer, renderer)
            self._force_full = False
            self._png_is_old = False
            return buff

    def get_renderer(self, cleared=None):
        # Mirrors super.get_renderer, but caches the old one
        # so that we can do things such as produce a diff image
        # in get_diff_image
        _, _, w, h = self.figure.bbox.bounds
        w, h = int(w), int(h)
        key = w, h, self.figure.dpi
        try:
            self._lastKey, self._renderer
        except AttributeError:
            need_new_renderer = True
        else:
            need_new_renderer = (self._lastKey != key)

        if need_new_renderer:
            self._renderer = backend_agg.RendererAgg(
                w, h, self.figure.dpi)
            self._last_renderer = backend_agg.RendererAgg(
                w, h, self.figure.dpi)
            self._lastKey = key

        elif cleared:
            self._renderer.clear()

        return self._renderer

    def handle_event(self, event):
        e_type = event['type']
        handler = getattr(self, 'handle_{0}'.format(e_type),
                          self.handle_unknown_event)
        return handler(event)

    def handle_unknown_event(self, event):
        warnings.warn('Unhandled message type {0}. {1}'.format(
            event['type'], event))

    def handle_ack(self, event):
        # Network latency tends to decrease if traffic is flowing
        # in both directions.  Therefore, the browser sends back
        # an "ack" message after each image frame is received.
        # This could also be used as a simple sanity check in the
        # future, but for now the performance increase is enough
        # to justify it, even if the server does nothing with it.
        pass

    def handle_draw(self, event):
        self.draw()

    def _handle_mouse(self, event):
        x = event['x']
        y = event['y']
        y = self.get_renderer().height - y

        # Javascript button numbers and matplotlib button numbers are
        # off by 1
        button = event['button'] + 1

        # The right mouse button pops up a context menu, which
        # doesn't work very well, so use the middle mouse button
        # instead.  It doesn't seem that it's possible to disable
        # the context menu in recent versions of Chrome.  If this
        # is resolved, please also adjust the docstring in MouseEvent.
        if button == 2:
            button = 3

        e_type = event['type']
        guiEvent = event.get('guiEvent', None)
        if e_type == 'button_press':
            self.button_press_event(x, y, button, guiEvent=guiEvent)
        elif e_type == 'button_release':
            self.button_release_event(x, y, button, guiEvent=guiEvent)
        elif e_type == 'motion_notify':
            self.motion_notify_event(x, y, guiEvent=guiEvent)
        elif e_type == 'figure_enter':
            self.enter_notify_event(xy=(x, y), guiEvent=guiEvent)
        elif e_type == 'figure_leave':
            self.leave_notify_event()
        elif e_type == 'scroll':
            self.scroll_event(x, y, event['step'], guiEvent=guiEvent)
    handle_button_press = handle_button_release = handle_motion_notify = \
        handle_figure_enter = handle_figure_leave = handle_scroll = \
        _handle_mouse

    def _handle_key(self, event):
        key = _handle_key(event['key'])
        e_type = event['type']
        guiEvent = event.get('guiEvent', None)
        if e_type == 'key_press':
            self.key_press_event(key, guiEvent=guiEvent)
        elif e_type == 'key_release':
            self.key_release_event(key, guiEvent=guiEvent)
    handle_key_press = handle_key_release = _handle_key

    def handle_toolbar_button(self, event):
        # TODO: Be more suspicious of the input
        getattr(self.toolbar, event['name'])()

    def handle_refresh(self, event):
        figure_label = self.figure.get_label()
        if not figure_label:
            figure_label = "Figure {0}".format(self.manager.num)
        self.send_event('figure_label', label=figure_label)
        self._force_full = True
        self.draw_idle()

    def handle_resize(self, event):
        x, y = event.get('width', 800), event.get('height', 800)
        x, y = int(x) * self._dpi_ratio, int(y) * self._dpi_ratio
        fig = self.figure
        # An attempt at approximating the figure size in pixels.
        fig.set_size_inches(x / fig.dpi, y / fig.dpi, forward=False)

        _, _, w, h = self.figure.bbox.bounds
        # Acknowledge the resize, and force the viewer to update the
        # canvas size to the figure's new size (which is hopefully
        # identical or within a pixel or so).
        self._png_is_old = True
        self.manager.resize(w, h)
        self.resize_event()

    def handle_send_image_mode(self, event):
        # The client requests notification of what the current image mode is.
        self.send_event('image_mode', mode=self._current_image_mode)

    def handle_set_dpi_ratio(self, event):
        dpi_ratio = event.get('dpi_ratio', 1)
        if dpi_ratio != self._dpi_ratio:
            # We don't want to scale up the figure dpi more than once.
            if not hasattr(self.figure, '_original_dpi'):
                self.figure._original_dpi = self.figure.dpi
            self.figure.dpi = dpi_ratio * self.figure._original_dpi
            self._dpi_ratio = dpi_ratio
            self._force_full = True
            self.draw_idle()

    def send_event(self, event_type, **kwargs):
        self.manager._send_event(event_type, **kwargs)


_JQUERY_ICON_CLASSES = {
    'home': 'ui-icon ui-icon-home',
    'back': 'ui-icon ui-icon-circle-arrow-w',
    'forward': 'ui-icon ui-icon-circle-arrow-e',
    'zoom_to_rect': 'ui-icon ui-icon-search',
    'move': 'ui-icon ui-icon-arrow-4',
    'download': 'ui-icon ui-icon-disk',
    None: None,
}


class NavigationToolbar2WebAgg(backend_bases.NavigationToolbar2):

    # Use the standard toolbar items + download button
    toolitems = [(text, tooltip_text, _JQUERY_ICON_CLASSES[image_file],
                  name_of_method)
                 for text, tooltip_text, image_file, name_of_method
                 in (backend_bases.NavigationToolbar2.toolitems +
                     (('Download', 'Download plot', 'download', 'download'),))
                 if image_file in _JQUERY_ICON_CLASSES]

    def _init_toolbar(self):
        self.message = ''
        self.cursor = 0

    def set_message(self, message):
        if message != self.message:
            self.canvas.send_event("message", message=message)
        self.message = message

    def set_cursor(self, cursor):
        if cursor != self.cursor:
            self.canvas.send_event("cursor", cursor=cursor)
        self.cursor = cursor

    def draw_rubberband(self, event, x0, y0, x1, y1):
        self.canvas.send_event(
            "rubberband", x0=x0, y0=y0, x1=x1, y1=y1)

    def release_zoom(self, event):
        backend_bases.NavigationToolbar2.release_zoom(self, event)
        self.canvas.send_event(
            "rubberband", x0=-1, y0=-1, x1=-1, y1=-1)

    def save_figure(self, *args):
        """Save the current figure"""
        self.canvas.send_event('save')


class FigureManagerWebAgg(backend_bases.FigureManagerBase):
    ToolbarCls = NavigationToolbar2WebAgg

    def __init__(self, canvas, num):
        backend_bases.FigureManagerBase.__init__(self, canvas, num)

        self.web_sockets = set()

        self.toolbar = self._get_toolbar(canvas)

    def show(self):
        pass

    def _get_toolbar(self, canvas):
        toolbar = self.ToolbarCls(canvas)
        return toolbar

    def resize(self, w, h):
        self._send_event(
            'resize',
            size=(w / self.canvas._dpi_ratio, h / self.canvas._dpi_ratio))

    def set_window_title(self, title):
        self._send_event('figure_label', label=title)

    # The following methods are specific to FigureManagerWebAgg

    def add_web_socket(self, web_socket):
        assert hasattr(web_socket, 'send_binary')
        assert hasattr(web_socket, 'send_json')

        self.web_sockets.add(web_socket)

        _, _, w, h = self.canvas.figure.bbox.bounds
        self.resize(w, h)
        self._send_event('refresh')

    def remove_web_socket(self, web_socket):
        self.web_sockets.remove(web_socket)

    def handle_json(self, content):
        self.canvas.handle_event(content)

    def refresh_all(self):
        if self.web_sockets:
            diff = self.canvas.get_diff_image()
            if diff is not None:
                for s in self.web_sockets:
                    s.send_binary(diff)

    @classmethod
    def get_javascript(cls, stream=None):
        if stream is None:
            output = io.StringIO()
        else:
            output = stream

        with io.open(os.path.join(
                os.path.dirname(__file__),
                "web_backend", "js",
                "mpl.js"), encoding='utf8') as fd:
            output.write(fd.read())

        toolitems = []
        for name, tooltip, image, method in cls.ToolbarCls.toolitems:
            if name is None:
                toolitems.append(['', '', '', ''])
            else:
                toolitems.append([name, tooltip, image, method])
        output.write("mpl.toolbar_items = {0};\n\n".format(
            json.dumps(toolitems)))

        extensions = []
        for filetype, ext in sorted(FigureCanvasWebAggCore.
                                    get_supported_filetypes_grouped().
                                    items()):
            if not ext[0] == 'pgf':  # pgf does not support BytesIO
                extensions.append(ext[0])
        output.write("mpl.extensions = {0};\n\n".format(
            json.dumps(extensions)))

        output.write("mpl.default_extension = {0};".format(
            json.dumps(FigureCanvasWebAggCore.get_default_filetype())))

        if stream is None:
            return output.getvalue()

    @classmethod
    def get_static_file_path(cls):
        return os.path.join(os.path.dirname(__file__), 'web_backend')

    def _send_event(self, event_type, **kwargs):
        payload = {'type': event_type}
        payload.update(kwargs)
        for s in self.web_sockets:
            s.send_json(payload)


class TimerTornado(backend_bases.TimerBase):
    def _timer_start(self):
        self._timer_stop()
        if self._single:
            ioloop = tornado.ioloop.IOLoop.instance()
            self._timer = ioloop.add_timeout(
                datetime.timedelta(milliseconds=self.interval),
                self._on_timer)
        else:
            self._timer = tornado.ioloop.PeriodicCallback(
                self._on_timer,
                self.interval)
            self._timer.start()

    def _timer_stop(self):
        if self._timer is None:
            return
        elif self._single:
            ioloop = tornado.ioloop.IOLoop.instance()
            ioloop.remove_timeout(self._timer)
        else:
            self._timer.stop()

        self._timer = None

    def _timer_set_interval(self):
        # Only stop and restart it if the timer has already been started
        if self._timer is not None:
            self._timer_stop()
            self._timer_start()


@_Backend.export
class _BackendWebAggCoreAgg(_Backend):
    FigureCanvas = FigureCanvasWebAggCore
    FigureManager = FigureManagerWebAgg
