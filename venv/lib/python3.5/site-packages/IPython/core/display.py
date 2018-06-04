# -*- coding: utf-8 -*-
"""Top-level display functions for displaying object in different formats."""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.


from binascii import b2a_hex, b2a_base64, hexlify
import json
import mimetypes
import os
import struct
import sys
import warnings
from copy import deepcopy

from IPython.utils.py3compat import cast_unicode
from IPython.testing.skipdoctest import skip_doctest

__all__ = ['display', 'display_pretty', 'display_html', 'display_markdown',
'display_svg', 'display_png', 'display_jpeg', 'display_latex', 'display_json',
'display_javascript', 'display_pdf', 'DisplayObject', 'TextDisplayObject',
'Pretty', 'HTML', 'Markdown', 'Math', 'Latex', 'SVG', 'ProgressBar', 'JSON',
'GeoJSON', 'Javascript', 'Image', 'clear_output', 'set_matplotlib_formats',
'set_matplotlib_close', 'publish_display_data', 'update_display', 'DisplayHandle',
'Video']

#-----------------------------------------------------------------------------
# utility functions
#-----------------------------------------------------------------------------

def _safe_exists(path):
    """Check path, but don't let exceptions raise"""
    try:
        return os.path.exists(path)
    except Exception:
        return False

def _merge(d1, d2):
    """Like update, but merges sub-dicts instead of clobbering at the top level.

    Updates d1 in-place
    """

    if not isinstance(d2, dict) or not isinstance(d1, dict):
        return d2
    for key, value in d2.items():
        d1[key] = _merge(d1.get(key), value)
    return d1

def _display_mimetype(mimetype, objs, raw=False, metadata=None):
    """internal implementation of all display_foo methods

    Parameters
    ----------
    mimetype : str
        The mimetype to be published (e.g. 'image/png')
    objs : tuple of objects
        The Python objects to display, or if raw=True raw text data to
        display.
    raw : bool
        Are the data objects raw data or Python objects that need to be
        formatted before display? [default: False]
    metadata : dict (optional)
        Metadata to be associated with the specific mimetype output.
    """
    if metadata:
        metadata = {mimetype: metadata}
    if raw:
        # turn list of pngdata into list of { 'image/png': pngdata }
        objs = [ {mimetype: obj} for obj in objs ]
    display(*objs, raw=raw, metadata=metadata, include=[mimetype])

#-----------------------------------------------------------------------------
# Main functions
#-----------------------------------------------------------------------------

# use * to indicate transient is keyword-only
def publish_display_data(data, metadata=None, source=None, *, transient=None, **kwargs):
    """Publish data and metadata to all frontends.

    See the ``display_data`` message in the messaging documentation for
    more details about this message type.

    Keys of data and metadata can be any mime-type.

    Parameters
    ----------
    data : dict
        A dictionary having keys that are valid MIME types (like
        'text/plain' or 'image/svg+xml') and values that are the data for
        that MIME type. The data itself must be a JSON'able data
        structure. Minimally all data should have the 'text/plain' data,
        which can be displayed by all frontends. If more than the plain
        text is given, it is up to the frontend to decide which
        representation to use.
    metadata : dict
        A dictionary for metadata related to the data. This can contain
        arbitrary key, value pairs that frontends can use to interpret
        the data. mime-type keys matching those in data can be used
        to specify metadata about particular representations.
    source : str, deprecated
        Unused.
    transient : dict, keyword-only
        A dictionary of transient data, such as display_id.
        """
    from IPython.core.interactiveshell import InteractiveShell

    display_pub = InteractiveShell.instance().display_pub

    # only pass transient if supplied,
    # to avoid errors with older ipykernel.
    # TODO: We could check for ipykernel version and provide a detailed upgrade message.
    if transient:
        kwargs['transient'] = transient

    display_pub.publish(
        data=data,
        metadata=metadata,
        **kwargs
    )


def _new_id():
    """Generate a new random text id with urandom"""
    return b2a_hex(os.urandom(16)).decode('ascii')


def display(*objs, include=None, exclude=None, metadata=None, transient=None, display_id=None, **kwargs):
    """Display a Python object in all frontends.

    By default all representations will be computed and sent to the frontends.
    Frontends can decide which representation is used and how.

    In terminal IPython this will be similar to using :func:`print`, for use in richer
    frontends see Jupyter notebook examples with rich display logic.

    Parameters
    ----------
    objs : tuple of objects
        The Python objects to display.
    raw : bool, optional
        Are the objects to be displayed already mimetype-keyed dicts of raw display data,
        or Python objects that need to be formatted before display? [default: False]
    include : list, tuple or set, optional
        A list of format type strings (MIME types) to include in the
        format data dict. If this is set *only* the format types included
        in this list will be computed.
    exclude : list, tuple or set, optional
        A list of format type strings (MIME types) to exclude in the format
        data dict. If this is set all format types will be computed,
        except for those included in this argument.
    metadata : dict, optional
        A dictionary of metadata to associate with the output.
        mime-type keys in this dictionary will be associated with the individual
        representation formats, if they exist.
    transient : dict, optional
        A dictionary of transient data to associate with the output.
        Data in this dict should not be persisted to files (e.g. notebooks).
    display_id : str, bool optional
        Set an id for the display.
        This id can be used for updating this display area later via update_display.
        If given as `True`, generate a new `display_id`
    kwargs: additional keyword-args, optional
        Additional keyword-arguments are passed through to the display publisher.

    Returns
    -------

    handle: DisplayHandle
        Returns a handle on updatable displays for use with :func:`update_display`,
        if `display_id` is given. Returns :any:`None` if no `display_id` is given
        (default).

    Examples
    --------

    >>> class Json(object):
    ...     def __init__(self, json):
    ...         self.json = json
    ...     def _repr_pretty_(self, pp, cycle):
    ...         import json
    ...         pp.text(json.dumps(self.json, indent=2))
    ...     def __repr__(self):
    ...         return str(self.json)
    ...

    >>> d = Json({1:2, 3: {4:5}})

    >>> print(d)
    {1: 2, 3: {4: 5}}

    >>> display(d)
    {
      "1": 2,
      "3": {
        "4": 5
      }
    }

    >>> def int_formatter(integer, pp, cycle):
    ...     pp.text('I'*integer)

    >>> plain = get_ipython().display_formatter.formatters['text/plain']
    >>> plain.for_type(int, int_formatter)
    <function _repr_pprint at 0x...>
    >>> display(7-5)
    II

    >>> del plain.type_printers[int]
    >>> display(7-5)
    2

    See Also
    --------

    :func:`update_display`

    Notes
    -----

    In Python, objects can declare their textual representation using the
    `__repr__` method. IPython expands on this idea and allows objects to declare
    other, rich representations including:

      - HTML
      - JSON
      - PNG
      - JPEG
      - SVG
      - LaTeX

    A single object can declare some or all of these representations; all are
    handled by IPython's display system.

    The main idea of the first approach is that you have to implement special
    display methods when you define your class, one for each representation you
    want to use. Here is a list of the names of the special methods and the
    values they must return:

      - `_repr_html_`: return raw HTML as a string
      - `_repr_json_`: return a JSONable dict
      - `_repr_jpeg_`: return raw JPEG data
      - `_repr_png_`: return raw PNG data
      - `_repr_svg_`: return raw SVG data as a string
      - `_repr_latex_`: return LaTeX commands in a string surrounded by "$".
      - `_repr_mimebundle_`: return a full mimebundle containing the mapping
                             from all mimetypes to data.
                             Use this for any mime-type not listed above.

    When you are directly writing your own classes, you can adapt them for
    display in IPython by following the above approach. But in practice, you
    often need to work with existing classes that you can't easily modify.

    You can refer to the documentation on integrating with the display system in
    order to register custom formatters for already existing types
    (:ref:`integrating_rich_display`).

    .. versionadded:: 5.4 display available without import
    .. versionadded:: 6.1 display available without import

    Since IPython 5.4 and 6.1 :func:`display` is automatically made available to
    the user without import. If you are using display in a document that might
    be used in a pure python context or with older version of IPython, use the
    following import at the top of your file::

        from IPython.display import display

    """
    from IPython.core.interactiveshell import InteractiveShell
    
    if not InteractiveShell.initialized():
        # Directly print objects.
        print(*objs)
        return
    
    raw = kwargs.pop('raw', False)
    if transient is None:
        transient = {}
    if metadata is None:
        metadata={}
    if display_id:
        if display_id is True:
            display_id = _new_id()
        transient['display_id'] = display_id
    if kwargs.get('update') and 'display_id' not in transient:
        raise TypeError('display_id required for update_display')
    if transient:
        kwargs['transient'] = transient

    if not raw:
        format = InteractiveShell.instance().display_formatter.format

    for obj in objs:
        if raw:
            publish_display_data(data=obj, metadata=metadata, **kwargs)
        else:
            format_dict, md_dict = format(obj, include=include, exclude=exclude)
            if not format_dict:
                # nothing to display (e.g. _ipython_display_ took over)
                continue
            if metadata:
                # kwarg-specified metadata gets precedence
                _merge(md_dict, metadata)
            publish_display_data(data=format_dict, metadata=md_dict, **kwargs)
    if display_id:
        return DisplayHandle(display_id)


# use * for keyword-only display_id arg
def update_display(obj, *, display_id, **kwargs):
    """Update an existing display by id

    Parameters
    ----------

    obj:
        The object with which to update the display
    display_id: keyword-only
        The id of the display to update

    See Also
    --------

    :func:`display`
    """
    kwargs['update'] = True
    display(obj, display_id=display_id, **kwargs)


class DisplayHandle(object):
    """A handle on an updatable display

    Call `.update(obj)` to display a new object.

    Call `.display(obj`) to add a new instance of this display,
    and update existing instances.

    See Also
    --------

        :func:`display`, :func:`update_display`

    """

    def __init__(self, display_id=None):
        if display_id is None:
            display_id = _new_id()
        self.display_id = display_id

    def __repr__(self):
        return "<%s display_id=%s>" % (self.__class__.__name__, self.display_id)

    def display(self, obj, **kwargs):
        """Make a new display with my id, updating existing instances.
        
        Parameters
        ----------
        
        obj:
            object to display
        **kwargs:
            additional keyword arguments passed to display
        """
        display(obj, display_id=self.display_id, **kwargs)

    def update(self, obj, **kwargs):
        """Update existing displays with my id
        
        Parameters
        ----------
        
        obj:
            object to display
        **kwargs:
            additional keyword arguments passed to update_display
        """
        update_display(obj, display_id=self.display_id, **kwargs)


def display_pretty(*objs, **kwargs):
    """Display the pretty (default) representation of an object.

    Parameters
    ----------
    objs : tuple of objects
        The Python objects to display, or if raw=True raw text data to
        display.
    raw : bool
        Are the data objects raw data or Python objects that need to be
        formatted before display? [default: False]
    metadata : dict (optional)
        Metadata to be associated with the specific mimetype output.
    """
    _display_mimetype('text/plain', objs, **kwargs)


def display_html(*objs, **kwargs):
    """Display the HTML representation of an object.

    Note: If raw=False and the object does not have a HTML
    representation, no HTML will be shown.

    Parameters
    ----------
    objs : tuple of objects
        The Python objects to display, or if raw=True raw HTML data to
        display.
    raw : bool
        Are the data objects raw data or Python objects that need to be
        formatted before display? [default: False]
    metadata : dict (optional)
        Metadata to be associated with the specific mimetype output.
    """
    _display_mimetype('text/html', objs, **kwargs)


def display_markdown(*objs, **kwargs):
    """Displays the Markdown representation of an object.

    Parameters
    ----------
    objs : tuple of objects
        The Python objects to display, or if raw=True raw markdown data to
        display.
    raw : bool
        Are the data objects raw data or Python objects that need to be
        formatted before display? [default: False]
    metadata : dict (optional)
        Metadata to be associated with the specific mimetype output.
    """

    _display_mimetype('text/markdown', objs, **kwargs)


def display_svg(*objs, **kwargs):
    """Display the SVG representation of an object.

    Parameters
    ----------
    objs : tuple of objects
        The Python objects to display, or if raw=True raw svg data to
        display.
    raw : bool
        Are the data objects raw data or Python objects that need to be
        formatted before display? [default: False]
    metadata : dict (optional)
        Metadata to be associated with the specific mimetype output.
    """
    _display_mimetype('image/svg+xml', objs, **kwargs)


def display_png(*objs, **kwargs):
    """Display the PNG representation of an object.

    Parameters
    ----------
    objs : tuple of objects
        The Python objects to display, or if raw=True raw png data to
        display.
    raw : bool
        Are the data objects raw data or Python objects that need to be
        formatted before display? [default: False]
    metadata : dict (optional)
        Metadata to be associated with the specific mimetype output.
    """
    _display_mimetype('image/png', objs, **kwargs)


def display_jpeg(*objs, **kwargs):
    """Display the JPEG representation of an object.

    Parameters
    ----------
    objs : tuple of objects
        The Python objects to display, or if raw=True raw JPEG data to
        display.
    raw : bool
        Are the data objects raw data or Python objects that need to be
        formatted before display? [default: False]
    metadata : dict (optional)
        Metadata to be associated with the specific mimetype output.
    """
    _display_mimetype('image/jpeg', objs, **kwargs)


def display_latex(*objs, **kwargs):
    """Display the LaTeX representation of an object.

    Parameters
    ----------
    objs : tuple of objects
        The Python objects to display, or if raw=True raw latex data to
        display.
    raw : bool
        Are the data objects raw data or Python objects that need to be
        formatted before display? [default: False]
    metadata : dict (optional)
        Metadata to be associated with the specific mimetype output.
    """
    _display_mimetype('text/latex', objs, **kwargs)


def display_json(*objs, **kwargs):
    """Display the JSON representation of an object.

    Note that not many frontends support displaying JSON.

    Parameters
    ----------
    objs : tuple of objects
        The Python objects to display, or if raw=True raw json data to
        display.
    raw : bool
        Are the data objects raw data or Python objects that need to be
        formatted before display? [default: False]
    metadata : dict (optional)
        Metadata to be associated with the specific mimetype output.
    """
    _display_mimetype('application/json', objs, **kwargs)


def display_javascript(*objs, **kwargs):
    """Display the Javascript representation of an object.

    Parameters
    ----------
    objs : tuple of objects
        The Python objects to display, or if raw=True raw javascript data to
        display.
    raw : bool
        Are the data objects raw data or Python objects that need to be
        formatted before display? [default: False]
    metadata : dict (optional)
        Metadata to be associated with the specific mimetype output.
    """
    _display_mimetype('application/javascript', objs, **kwargs)


def display_pdf(*objs, **kwargs):
    """Display the PDF representation of an object.

    Parameters
    ----------
    objs : tuple of objects
        The Python objects to display, or if raw=True raw javascript data to
        display.
    raw : bool
        Are the data objects raw data or Python objects that need to be
        formatted before display? [default: False]
    metadata : dict (optional)
        Metadata to be associated with the specific mimetype output.
    """
    _display_mimetype('application/pdf', objs, **kwargs)


#-----------------------------------------------------------------------------
# Smart classes
#-----------------------------------------------------------------------------


class DisplayObject(object):
    """An object that wraps data to be displayed."""

    _read_flags = 'r'
    _show_mem_addr = False
    metadata = None

    def __init__(self, data=None, url=None, filename=None, metadata=None):
        """Create a display object given raw data.

        When this object is returned by an expression or passed to the
        display function, it will result in the data being displayed
        in the frontend. The MIME type of the data should match the
        subclasses used, so the Png subclass should be used for 'image/png'
        data. If the data is a URL, the data will first be downloaded
        and then displayed. If

        Parameters
        ----------
        data : unicode, str or bytes
            The raw data or a URL or file to load the data from
        url : unicode
            A URL to download the data from.
        filename : unicode
            Path to a local file to load the data from.
        metadata : dict
            Dict of metadata associated to be the object when displayed
        """
        if data is not None and isinstance(data, str):
            if data.startswith('http') and url is None:
                url = data
                filename = None
                data = None
            elif _safe_exists(data) and filename is None:
                url = None
                filename = data
                data = None

        self.data = data
        self.url = url
        self.filename = filename

        if metadata is not None:
            self.metadata = metadata
        elif self.metadata is None:
            self.metadata = {}

        self.reload()
        self._check_data()

    def __repr__(self):
        if not self._show_mem_addr:
            cls = self.__class__
            r = "<%s.%s object>" % (cls.__module__, cls.__name__)
        else:
            r = super(DisplayObject, self).__repr__()
        return r

    def _check_data(self):
        """Override in subclasses if there's something to check."""
        pass

    def _data_and_metadata(self):
        """shortcut for returning metadata with shape information, if defined"""
        if self.metadata:
            return self.data, deepcopy(self.metadata)
        else:
            return self.data

    def reload(self):
        """Reload the raw data from file or URL."""
        if self.filename is not None:
            with open(self.filename, self._read_flags) as f:
                self.data = f.read()
        elif self.url is not None:
            try:
                # Deferred import
                from urllib.request import urlopen
                response = urlopen(self.url)
                self.data = response.read()
                # extract encoding from header, if there is one:
                encoding = None
                for sub in response.headers['content-type'].split(';'):
                    sub = sub.strip()
                    if sub.startswith('charset'):
                        encoding = sub.split('=')[-1].strip()
                        break
                # decode data, if an encoding was specified
                if encoding:
                    self.data = self.data.decode(encoding, 'replace')
            except:
                self.data = None

class TextDisplayObject(DisplayObject):
    """Validate that display data is text"""
    def _check_data(self):
        if self.data is not None and not isinstance(self.data, str):
            raise TypeError("%s expects text, not %r" % (self.__class__.__name__, self.data))

class Pretty(TextDisplayObject):

    def _repr_pretty_(self, pp, cycle):
        return pp.text(self.data)


class HTML(TextDisplayObject):

    def _repr_html_(self):
        return self._data_and_metadata()

    def __html__(self):
        """
        This method exists to inform other HTML-using modules (e.g. Markupsafe,
        htmltag, etc) that this object is HTML and does not need things like
        special characters (<>&) escaped.
        """
        return self._repr_html_()


class Markdown(TextDisplayObject):

    def _repr_markdown_(self):
        return self._data_and_metadata()


class Math(TextDisplayObject):

    def _repr_latex_(self):
        s = "$$%s$$" % self.data.strip('$')
        if self.metadata:
            return s, deepcopy(self.metadata)
        else:
            return s


class Latex(TextDisplayObject):

    def _repr_latex_(self):
        return self._data_and_metadata()


class SVG(DisplayObject):

    _read_flags = 'rb'
    # wrap data in a property, which extracts the <svg> tag, discarding
    # document headers
    _data = None

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, svg):
        if svg is None:
            self._data = None
            return
        # parse into dom object
        from xml.dom import minidom
        x = minidom.parseString(svg)
        # get svg tag (should be 1)
        found_svg = x.getElementsByTagName('svg')
        if found_svg:
            svg = found_svg[0].toxml()
        else:
            # fallback on the input, trust the user
            # but this is probably an error.
            pass
        svg = cast_unicode(svg)
        self._data = svg
    
    def _repr_svg_(self):
        return self._data_and_metadata()

class ProgressBar(DisplayObject):
    """Progressbar supports displaying a progressbar like element 
    """
    def __init__(self, total):
        """Creates a new progressbar
        
        Parameters
        ----------
        total : int
            maximum size of the progressbar
        """
        self.total = total
        self._progress = 0
        self.html_width = '60ex'
        self.text_width = 60
        self._display_id = hexlify(os.urandom(8)).decode('ascii')

    def __repr__(self):
        fraction = self.progress / self.total
        filled = '=' * int(fraction * self.text_width)
        rest = ' ' * (self.text_width - len(filled))
        return '[{}{}] {}/{}'.format(
            filled, rest,
            self.progress, self.total,
        )

    def _repr_html_(self):
        return "<progress style='width:{}' max='{}' value='{}'></progress>".format(
            self.html_width, self.total, self.progress)

    def display(self):
        display(self, display_id=self._display_id)

    def update(self):
        display(self, display_id=self._display_id, update=True)

    @property
    def progress(self):
        return self._progress

    @progress.setter
    def progress(self, value):
        self._progress = value
        self.update()

    def __iter__(self):
        self.display()
        self._progress = -1 # First iteration is 0
        return self

    def __next__(self):
        """Returns current value and increments display by one."""
        self.progress += 1
        if self.progress < self.total:
            return self.progress
        else:
            raise StopIteration()

class JSON(DisplayObject):
    """JSON expects a JSON-able dict or list

    not an already-serialized JSON string.

    Scalar types (None, number, string) are not allowed, only dict or list containers.
    """
    # wrap data in a property, which warns about passing already-serialized JSON
    _data = None
    def __init__(self, data=None, url=None, filename=None, expanded=False, metadata=None, **kwargs):
        """Create a JSON display object given raw data.

        Parameters
        ----------
        data : dict or list
            JSON data to display. Not an already-serialized JSON string.
            Scalar types (None, number, string) are not allowed, only dict
            or list containers.
        url : unicode
            A URL to download the data from.
        filename : unicode
            Path to a local file to load the data from.
        expanded : boolean
            Metadata to control whether a JSON display component is expanded.
        metadata: dict
            Specify extra metadata to attach to the json display object.
        """
        self.metadata = {'expanded': expanded}
        if metadata:
            self.metadata.update(metadata)
        if kwargs:
            self.metadata.update(kwargs)
        super(JSON, self).__init__(data=data, url=url, filename=filename)

    def _check_data(self):
        if self.data is not None and not isinstance(self.data, (dict, list)):
            raise TypeError("%s expects JSONable dict or list, not %r" % (self.__class__.__name__, self.data))

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, data):
        if isinstance(data, str):
            if getattr(self, 'filename', None) is None:
                warnings.warn("JSON expects JSONable dict or list, not JSON strings")
            data = json.loads(data)
        self._data = data

    def _data_and_metadata(self):
        return self.data, self.metadata

    def _repr_json_(self):
        return self._data_and_metadata()

_css_t = """$("head").append($("<link/>").attr({
  rel:  "stylesheet",
  type: "text/css",
  href: "%s"
}));
"""

_lib_t1 = """$.getScript("%s", function () {
"""
_lib_t2 = """});
"""

class GeoJSON(JSON):
    """GeoJSON expects JSON-able dict

    not an already-serialized JSON string.

    Scalar types (None, number, string) are not allowed, only dict containers.
    """
    
    def __init__(self, *args, **kwargs):
        """Create a GeoJSON display object given raw data.

        Parameters
        ----------
        data : dict or list
            VegaLite data. Not an already-serialized JSON string.
            Scalar types (None, number, string) are not allowed, only dict
            or list containers.
        url_template : string
            Leaflet TileLayer URL template: http://leafletjs.com/reference.html#url-template
        layer_options : dict
            Leaflet TileLayer options: http://leafletjs.com/reference.html#tilelayer-options
        url : unicode
            A URL to download the data from.
        filename : unicode
            Path to a local file to load the data from.
        metadata: dict
            Specify extra metadata to attach to the json display object.

        Examples
        --------

        The following will display an interactive map of Mars with a point of
        interest on frontend that do support GeoJSON display.

            >>> from IPython.display import GeoJSON

            >>> GeoJSON(data={
            ...     "type": "Feature",
            ...     "geometry": {
            ...         "type": "Point",
            ...         "coordinates": [-81.327, 296.038]
            ...     }
            ... },
            ... url_template="http://s3-eu-west-1.amazonaws.com/whereonmars.cartodb.net/{basemap_id}/{z}/{x}/{y}.png",
            ... layer_options={
            ...     "basemap_id": "celestia_mars-shaded-16k_global",
            ...     "attribution" : "Celestia/praesepe",
            ...     "minZoom" : 0,
            ...     "maxZoom" : 18,
            ... })
            <IPython.core.display.GeoJSON object>

        In the terminal IPython, you will only see the text representation of
        the GeoJSON object.

        """
        
        super(GeoJSON, self).__init__(*args, **kwargs)


    def _ipython_display_(self):
        bundle = {
            'application/geo+json': self.data,
            'text/plain': '<IPython.display.GeoJSON object>'
        }
        metadata = {
            'application/geo+json': self.metadata
        }
        display(bundle, metadata=metadata, raw=True)

class Javascript(TextDisplayObject):

    def __init__(self, data=None, url=None, filename=None, lib=None, css=None):
        """Create a Javascript display object given raw data.

        When this object is returned by an expression or passed to the
        display function, it will result in the data being displayed
        in the frontend. If the data is a URL, the data will first be
        downloaded and then displayed.

        In the Notebook, the containing element will be available as `element`,
        and jQuery will be available.  Content appended to `element` will be
        visible in the output area.

        Parameters
        ----------
        data : unicode, str or bytes
            The Javascript source code or a URL to download it from.
        url : unicode
            A URL to download the data from.
        filename : unicode
            Path to a local file to load the data from.
        lib : list or str
            A sequence of Javascript library URLs to load asynchronously before
            running the source code. The full URLs of the libraries should
            be given. A single Javascript library URL can also be given as a
            string.
        css: : list or str
            A sequence of css files to load before running the source code.
            The full URLs of the css files should be given. A single css URL
            can also be given as a string.
        """
        if isinstance(lib, str):
            lib = [lib]
        elif lib is None:
            lib = []
        if isinstance(css, str):
            css = [css]
        elif css is None:
            css = []
        if not isinstance(lib, (list,tuple)):
            raise TypeError('expected sequence, got: %r' % lib)
        if not isinstance(css, (list,tuple)):
            raise TypeError('expected sequence, got: %r' % css)
        self.lib = lib
        self.css = css
        super(Javascript, self).__init__(data=data, url=url, filename=filename)

    def _repr_javascript_(self):
        r = ''
        for c in self.css:
            r += _css_t % c
        for l in self.lib:
            r += _lib_t1 % l
        r += self.data
        r += _lib_t2*len(self.lib)
        return r

# constants for identifying png/jpeg data
_PNG = b'\x89PNG\r\n\x1a\n'
_JPEG = b'\xff\xd8'

def _pngxy(data):
    """read the (width, height) from a PNG header"""
    ihdr = data.index(b'IHDR')
    # next 8 bytes are width/height
    return struct.unpack('>ii', data[ihdr+4:ihdr+12])

def _jpegxy(data):
    """read the (width, height) from a JPEG header"""
    # adapted from http://www.64lines.com/jpeg-width-height

    idx = 4
    while True:
        block_size = struct.unpack('>H', data[idx:idx+2])[0]
        idx = idx + block_size
        if data[idx:idx+2] == b'\xFF\xC0':
            # found Start of Frame
            iSOF = idx
            break
        else:
            # read another block
            idx += 2

    h, w = struct.unpack('>HH', data[iSOF+5:iSOF+9])
    return w, h

def _gifxy(data):
    """read the (width, height) from a GIF header"""
    return struct.unpack('<HH', data[6:10])


class Image(DisplayObject):

    _read_flags = 'rb'
    _FMT_JPEG = u'jpeg'
    _FMT_PNG = u'png'
    _FMT_GIF = u'gif'
    _ACCEPTABLE_EMBEDDINGS = [_FMT_JPEG, _FMT_PNG, _FMT_GIF]
    _MIMETYPES = {
        _FMT_PNG: 'image/png',
        _FMT_JPEG: 'image/jpeg',
        _FMT_GIF: 'image/gif',
    }

    def __init__(self, data=None, url=None, filename=None, format=None,
                 embed=None, width=None, height=None, retina=False,
                 unconfined=False, metadata=None):
        """Create a PNG/JPEG/GIF image object given raw data.

        When this object is returned by an input cell or passed to the
        display function, it will result in the image being displayed
        in the frontend.

        Parameters
        ----------
        data : unicode, str or bytes
            The raw image data or a URL or filename to load the data from.
            This always results in embedded image data.
        url : unicode
            A URL to download the data from. If you specify `url=`,
            the image data will not be embedded unless you also specify `embed=True`.
        filename : unicode
            Path to a local file to load the data from.
            Images from a file are always embedded.
        format : unicode
            The format of the image data (png/jpeg/jpg/gif). If a filename or URL is given
            for format will be inferred from the filename extension.
        embed : bool
            Should the image data be embedded using a data URI (True) or be
            loaded using an <img> tag. Set this to True if you want the image
            to be viewable later with no internet connection in the notebook.

            Default is `True`, unless the keyword argument `url` is set, then
            default value is `False`.

            Note that QtConsole is not able to display images if `embed` is set to `False`
        width : int
            Width in pixels to which to constrain the image in html
        height : int
            Height in pixels to which to constrain the image in html
        retina : bool
            Automatically set the width and height to half of the measured
            width and height.
            This only works for embedded images because it reads the width/height
            from image data.
            For non-embedded images, you can just set the desired display width
            and height directly.
        unconfined: bool
            Set unconfined=True to disable max-width confinement of the image.
        metadata: dict
            Specify extra metadata to attach to the image.

        Examples
        --------
        # embedded image data, works in qtconsole and notebook
        # when passed positionally, the first arg can be any of raw image data,
        # a URL, or a filename from which to load image data.
        # The result is always embedding image data for inline images.
        Image('http://www.google.fr/images/srpr/logo3w.png')
        Image('/path/to/image.jpg')
        Image(b'RAW_PNG_DATA...')

        # Specifying Image(url=...) does not embed the image data,
        # it only generates `<img>` tag with a link to the source.
        # This will not work in the qtconsole or offline.
        Image(url='http://www.google.fr/images/srpr/logo3w.png')

        """
        if filename is not None:
            ext = self._find_ext(filename)
        elif url is not None:
            ext = self._find_ext(url)
        elif data is None:
            raise ValueError("No image data found. Expecting filename, url, or data.")
        elif isinstance(data, str) and (
            data.startswith('http') or _safe_exists(data)
        ):
            ext = self._find_ext(data)
        else:
            ext = None

        if format is None:
            if ext is not None:
                if ext == u'jpg' or ext == u'jpeg':
                    format = self._FMT_JPEG
                elif ext == u'png':
                    format = self._FMT_PNG
                elif ext == u'gif':
                    format = self._FMT_GIF
                else:
                    format = ext.lower()
            elif isinstance(data, bytes):
                # infer image type from image data header,
                # only if format has not been specified.
                if data[:2] == _JPEG:
                    format = self._FMT_JPEG

        # failed to detect format, default png
        if format is None:
            format = self._FMT_PNG

        if format.lower() == 'jpg':
            # jpg->jpeg
            format = self._FMT_JPEG

        self.format = format.lower()
        self.embed = embed if embed is not None else (url is None)

        if self.embed and self.format not in self._ACCEPTABLE_EMBEDDINGS:
            raise ValueError("Cannot embed the '%s' image format" % (self.format))
        if self.embed:
            self._mimetype = self._MIMETYPES.get(self.format)

        self.width = width
        self.height = height
        self.retina = retina
        self.unconfined = unconfined
        super(Image, self).__init__(data=data, url=url, filename=filename, 
                metadata=metadata)

        if self.width is None and self.metadata.get('width', {}):
            self.width = metadata['width']

        if self.height is None and self.metadata.get('height', {}):
            self.height = metadata['height']

        if retina:
            self._retina_shape()


    def _retina_shape(self):
        """load pixel-doubled width and height from image data"""
        if not self.embed:
            return
        if self.format == self._FMT_PNG:
            w, h = _pngxy(self.data)
        elif self.format == self._FMT_JPEG:
            w, h = _jpegxy(self.data)
        elif self.format == self._FMT_GIF:
            w, h = _gifxy(self.data)
        else:
            # retina only supports png
            return
        self.width = w // 2
        self.height = h // 2

    def reload(self):
        """Reload the raw data from file or URL."""
        if self.embed:
            super(Image,self).reload()
            if self.retina:
                self._retina_shape()

    def _repr_html_(self):
        if not self.embed:
            width = height = klass = ''
            if self.width:
                width = ' width="%d"' % self.width
            if self.height:
                height = ' height="%d"' % self.height
            if self.unconfined:
                klass = ' class="unconfined"'
            return u'<img src="{url}"{width}{height}{klass}/>'.format(
                url=self.url,
                width=width,
                height=height,
                klass=klass,
            )

    def _repr_mimebundle_(self, include=None, exclude=None):
        """Return the image as a mimebundle

        Any new mimetype support should be implemented here.
        """
        if self.embed:
            mimetype = self._mimetype
            data, metadata = self._data_and_metadata(always_both=True)
            if metadata:
                metadata = {mimetype: metadata}
            return {mimetype: data}, metadata
        else:
            return {'text/html': self._repr_html_()}

    def _data_and_metadata(self, always_both=False):
        """shortcut for returning metadata with shape information, if defined"""
        b64_data = b2a_base64(self.data).decode('ascii')
        md = {}
        if self.metadata:
            md.update(self.metadata)
        if self.width:
            md['width'] = self.width
        if self.height:
            md['height'] = self.height
        if self.unconfined:
            md['unconfined'] = self.unconfined
        if md or always_both:
            return b64_data, md
        else:
            return b64_data

    def _repr_png_(self):
        if self.embed and self.format == self._FMT_PNG:
            return self._data_and_metadata()

    def _repr_jpeg_(self):
        if self.embed and self.format == self._FMT_JPEG:
            return self._data_and_metadata()

    def _find_ext(self, s):
        return s.split('.')[-1].lower()


class Video(DisplayObject):

    def __init__(self, data=None, url=None, filename=None, embed=False, mimetype=None):
        """Create a video object given raw data or an URL.

        When this object is returned by an input cell or passed to the
        display function, it will result in the video being displayed
        in the frontend.

        Parameters
        ----------
        data : unicode, str or bytes
            The raw video data or a URL or filename to load the data from.
            Raw data will require passing `embed=True`.
        url : unicode
            A URL for the video. If you specify `url=`,
            the image data will not be embedded.
        filename : unicode
            Path to a local file containing the video.
            Will be interpreted as a local URL unless `embed=True`.
        embed : bool
            Should the video be embedded using a data URI (True) or be
            loaded using a <video> tag (False).

            Since videos are large, embedding them should be avoided, if possible.
            You must confirm embedding as your intention by passing `embed=True`.

            Local files can be displayed with URLs without embedding the content, via::

                Video('./video.mp4')

        mimetype: unicode
            Specify the mimetype for embedded videos.
            Default will be guessed from file extension, if available.

        Examples
        --------

        Video('https://archive.org/download/Sita_Sings_the_Blues/Sita_Sings_the_Blues_small.mp4')
        Video('path/to/video.mp4')
        Video('path/to/video.mp4', embed=True)
        Video(b'raw-videodata', embed=True)
        """
        if url is None and isinstance(data, str) and data.startswith(('http:', 'https:')):
            url = data
            data = None
        elif os.path.exists(data):
            filename = data
            data = None

        if data and not embed:
            msg = ''.join([
                "To embed videos, you must pass embed=True ",
                "(this may make your notebook files huge)\n",
                "Consider passing Video(url='...')",
            ])
            raise ValueError(msg)

        self.mimetype = mimetype
        self.embed = embed
        super(Video, self).__init__(data=data, url=url, filename=filename)

    def _repr_html_(self):
        # External URLs and potentially local files are not embedded into the
        # notebook output.
        if not self.embed:
            url = self.url if self.url is not None else self.filename
            output = """<video src="{0}" controls>
      Your browser does not support the <code>video</code> element.
    </video>""".format(url)
            return output

        # Embedded videos are base64-encoded.
        mimetype = self.mimetype
        if self.filename is not None:
            if not mimetype:
                mimetype, _ = mimetypes.guess_type(self.filename)

            with open(self.filename, 'rb') as f:
                video = f.read()
        else:
            video = self.data
        if isinstance(video, str):
            # unicode input is already b64-encoded
            b64_video = video
        else:
            b64_video = b2a_base64(video).decode('ascii').rstrip()

        output = """<video controls>
 <source src="data:{0};base64,{1}" type="{0}">
 Your browser does not support the video tag.
 </video>""".format(mimetype, b64_video)
        return output

    def reload(self):
        # TODO
        pass


def clear_output(wait=False):
    """Clear the output of the current cell receiving output.

    Parameters
    ----------
    wait : bool [default: false]
        Wait to clear the output until new output is available to replace it."""
    from IPython.core.interactiveshell import InteractiveShell
    if InteractiveShell.initialized():
        InteractiveShell.instance().display_pub.clear_output(wait)
    else:
        print('\033[2K\r', end='')
        sys.stdout.flush()
        print('\033[2K\r', end='')
        sys.stderr.flush()


@skip_doctest
def set_matplotlib_formats(*formats, **kwargs):
    """Select figure formats for the inline backend. Optionally pass quality for JPEG.

    For example, this enables PNG and JPEG output with a JPEG quality of 90%::

        In [1]: set_matplotlib_formats('png', 'jpeg', quality=90)

    To set this in your config files use the following::

        c.InlineBackend.figure_formats = {'png', 'jpeg'}
        c.InlineBackend.print_figure_kwargs.update({'quality' : 90})

    Parameters
    ----------
    *formats : strs
        One or more figure formats to enable: 'png', 'retina', 'jpeg', 'svg', 'pdf'.
    **kwargs :
        Keyword args will be relayed to ``figure.canvas.print_figure``.
    """
    from IPython.core.interactiveshell import InteractiveShell
    from IPython.core.pylabtools import select_figure_formats
    # build kwargs, starting with InlineBackend config
    kw = {}
    from ipykernel.pylab.config import InlineBackend
    cfg = InlineBackend.instance()
    kw.update(cfg.print_figure_kwargs)
    kw.update(**kwargs)
    shell = InteractiveShell.instance()
    select_figure_formats(shell, formats, **kw)

@skip_doctest
def set_matplotlib_close(close=True):
    """Set whether the inline backend closes all figures automatically or not.

    By default, the inline backend used in the IPython Notebook will close all
    matplotlib figures automatically after each cell is run. This means that
    plots in different cells won't interfere. Sometimes, you may want to make
    a plot in one cell and then refine it in later cells. This can be accomplished
    by::

        In [1]: set_matplotlib_close(False)

    To set this in your config files use the following::

        c.InlineBackend.close_figures = False

    Parameters
    ----------
    close : bool
        Should all matplotlib figures be automatically closed after each cell is
        run?
    """
    from ipykernel.pylab.config import InlineBackend
    cfg = InlineBackend.instance()
    cfg.close_figures = close
