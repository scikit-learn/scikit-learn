import json
import decimal
import datetime
import warnings
from pathlib import Path

from plotly.io._utils import validate_coerce_fig_to_dict, validate_coerce_output_type
from _plotly_utils.optional_imports import get_module
from _plotly_utils.basevalidators import ImageUriValidator


# Orca configuration class
# ------------------------
class JsonConfig(object):
    _valid_engines = ("json", "orjson", "auto")

    def __init__(self):
        self._default_engine = "auto"

    @property
    def default_engine(self):
        return self._default_engine

    @default_engine.setter
    def default_engine(self, val):
        if val not in JsonConfig._valid_engines:
            raise ValueError(
                "Supported JSON engines include {valid}\n    Received {val}".format(
                    valid=JsonConfig._valid_engines, val=val
                )
            )

        if val == "orjson":
            self.validate_orjson()

        self._default_engine = val

    @classmethod
    def validate_orjson(cls):
        orjson = get_module("orjson")
        if orjson is None:
            raise ValueError("The orjson engine requires the orjson package")


config = JsonConfig()


def coerce_to_strict(const):
    """
    This is used to ultimately *encode* into strict JSON, see `encode`

    """
    # before python 2.7, 'true', 'false', 'null', were include here.
    if const in ("Infinity", "-Infinity", "NaN"):
        return None
    else:
        return const


_swap_json = (
    ("<", "\\u003c"),
    (">", "\\u003e"),
    ("/", "\\u002f"),
)
_swap_orjson = _swap_json + (
    ("\u2028", "\\u2028"),
    ("\u2029", "\\u2029"),
)


def _safe(json_str, _swap):
    out = json_str
    for unsafe_char, safe_char in _swap:
        if unsafe_char in out:
            out = out.replace(unsafe_char, safe_char)
    return out


def to_json_plotly(plotly_object, pretty=False, engine=None):
    """
    Convert a plotly/Dash object to a JSON string representation

    Parameters
    ----------
    plotly_object:
        A plotly/Dash object represented as a dict, graph_object, or Dash component

    pretty: bool (default False)
        True if JSON representation should be pretty-printed, False if
        representation should be as compact as possible.

    engine: str (default None)
        The JSON encoding engine to use. One of:
          - "json" for an engine based on the built-in Python json module
          - "orjson" for a faster engine that requires the orjson package
          - "auto" for the "orjson" engine if available, otherwise "json"
        If not specified, the default engine is set to the current value of
        plotly.io.json.config.default_engine.

    Returns
    -------
    str
        Representation of input object as a JSON string

    See Also
    --------
    to_json : Convert a plotly Figure to JSON with validation
    """
    orjson = get_module("orjson", should_load=True)

    # Determine json engine
    if engine is None:
        engine = config.default_engine

    if engine == "auto":
        if orjson is not None:
            engine = "orjson"
        else:
            engine = "json"
    elif engine not in ["orjson", "json"]:
        raise ValueError("Invalid json engine: %s" % engine)

    modules = {
        "sage_all": get_module("sage.all", should_load=False),
        "np": get_module("numpy", should_load=False),
        "pd": get_module("pandas", should_load=False),
        "image": get_module("PIL.Image", should_load=False),
    }

    # Dump to a JSON string and return
    # --------------------------------
    if engine == "json":
        opts = {}
        if pretty:
            opts["indent"] = 2
        else:
            # Remove all whitespace
            opts["separators"] = (",", ":")

        from _plotly_utils.utils import PlotlyJSONEncoder

        return _safe(
            json.dumps(plotly_object, cls=PlotlyJSONEncoder, **opts), _swap_json
        )
    elif engine == "orjson":
        JsonConfig.validate_orjson()
        opts = orjson.OPT_NON_STR_KEYS | orjson.OPT_SERIALIZE_NUMPY

        if pretty:
            opts |= orjson.OPT_INDENT_2

        # Plotly
        try:
            plotly_object = plotly_object.to_plotly_json()
        except AttributeError:
            pass

        # Try without cleaning
        try:
            return _safe(
                orjson.dumps(plotly_object, option=opts).decode("utf8"), _swap_orjson
            )
        except TypeError:
            pass

        cleaned = clean_to_json_compatible(
            plotly_object,
            numpy_allowed=True,
            datetime_allowed=True,
            modules=modules,
        )
        return _safe(orjson.dumps(cleaned, option=opts).decode("utf8"), _swap_orjson)


def to_json(fig, validate=True, pretty=False, remove_uids=True, engine=None):
    """
    Convert a figure to a JSON string representation

    Parameters
    ----------
    fig:
        Figure object or dict representing a figure

    validate: bool (default True)
        True if the figure should be validated before being converted to
        JSON, False otherwise.

    pretty: bool (default False)
        True if JSON representation should be pretty-printed, False if
        representation should be as compact as possible.

    remove_uids: bool (default True)
        True if trace UIDs should be omitted from the JSON representation

    engine: str (default None)
        The JSON encoding engine to use. One of:
          - "json" for an engine based on the built-in Python json module
          - "orjson" for a faster engine that requires the orjson package
          - "auto" for the "orjson" engine if available, otherwise "json"
        If not specified, the default engine is set to the current value of
        plotly.io.json.config.default_engine.

    Returns
    -------
    str
        Representation of figure as a JSON string

    See Also
    --------
    to_json_plotly : Convert an arbitrary plotly graph_object or Dash component to JSON
    """
    # Validate figure
    # ---------------
    fig_dict = validate_coerce_fig_to_dict(fig, validate)

    # Remove trace uid
    # ----------------
    if remove_uids:
        for trace in fig_dict.get("data", []):
            trace.pop("uid", None)

    return to_json_plotly(fig_dict, pretty=pretty, engine=engine)


def write_json(fig, file, validate=True, pretty=False, remove_uids=True, engine=None):
    """
    Convert a figure to JSON and write it to a file or writeable
    object.

    Note: A figure converted to JSON with one version of Plotly.py may not be compatible with another version.

    Parameters
    ----------
    fig:
        Figure object or dict representing a figure

    file: str or writeable
        A string representing a local file path or a writeable object
        (e.g. a pathlib.Path object or an open file descriptor)

    pretty: bool (default False)
        True if JSON representation should be pretty-printed, False if
        representation should be as compact as possible.

    remove_uids: bool (default True)
        True if trace UIDs should be omitted from the JSON representation

    engine: str (default None)
        The JSON encoding engine to use. One of:
          - "json" for an engine based on the built-in Python json module
          - "orjson" for a faster engine that requires the orjson package
          - "auto" for the "orjson" engine if available, otherwise "json"
        If not specified, the default engine is set to the current value of
        plotly.io.json.config.default_engine.
    Returns
    -------
    None
    """

    # Get JSON string
    # ---------------
    # Pass through validate argument and let to_json handle validation logic
    json_str = to_json(
        fig, validate=validate, pretty=pretty, remove_uids=remove_uids, engine=engine
    )

    # Try to cast `file` as a pathlib object `path`.
    # ----------------------------------------------
    if isinstance(file, str):
        # Use the standard Path constructor to make a pathlib object.
        path = Path(file)
    elif isinstance(file, Path):
        # `file` is already a Path object.
        path = file
    else:
        # We could not make a Path object out of file. Either `file` is an open file
        # descriptor with a `write()` method or it's an invalid object.
        path = None

    # Open file
    # ---------
    if path is None:
        # We previously failed to make sense of `file` as a pathlib object.
        # Attempt to write to `file` as an open file descriptor.
        try:
            file.write(json_str)
            return
        except AttributeError:
            pass
        raise ValueError(
            """
The 'file' argument '{file}' is not a string, pathlib.Path object, or file descriptor.
""".format(file=file)
        )
    else:
        # We previously succeeded in interpreting `file` as a pathlib object.
        # Now we can use `write_bytes()`.
        path.write_text(json_str)


def from_json_plotly(value, engine=None):
    """
    Parse JSON string using the specified JSON engine

    Parameters
    ----------
    value: str or bytes
        A JSON string or bytes object

    engine: str (default None)
        The JSON decoding engine to use. One of:
          - if "json", parse JSON using built in json module
          - if "orjson", parse using the faster orjson module, requires the orjson
            package
          - if "auto" use orjson module if available, otherwise use the json module

        If not specified, the default engine is set to the current value of
        plotly.io.json.config.default_engine.

    Returns
    -------
    dict

    See Also
    --------
    from_json_plotly : Parse JSON with plotly conventions into a dict
    """
    orjson = get_module("orjson", should_load=True)

    # Validate value
    # --------------
    if not isinstance(value, (str, bytes)):
        raise ValueError(
            """
from_json_plotly requires a string or bytes argument but received value of type {typ}
    Received value: {value}""".format(typ=type(value), value=value)
        )

    # Determine json engine
    if engine is None:
        engine = config.default_engine

    if engine == "auto":
        if orjson is not None:
            engine = "orjson"
        else:
            engine = "json"
    elif engine not in ["orjson", "json"]:
        raise ValueError("Invalid json engine: %s" % engine)

    if engine == "orjson":
        JsonConfig.validate_orjson()
        # orjson handles bytes input natively
        value_dict = orjson.loads(value)
    else:
        # decode bytes to str for built-in json module
        if isinstance(value, bytes):
            value = value.decode("utf-8")
        value_dict = json.loads(value)

    return value_dict


def from_json(value, output_type="Figure", skip_invalid=False, engine=None):
    """
    Construct a figure from a JSON string

    Parameters
    ----------
    value: str or bytes
        String or bytes object containing the JSON representation of a figure

    output_type: type or str (default 'Figure')
        The output figure type or type name.
        One of:  graph_objs.Figure, 'Figure', graph_objs.FigureWidget, 'FigureWidget'

    skip_invalid: bool (default False)
        False if invalid figure properties should result in an exception.
        True if invalid figure properties should be silently ignored.

    engine: str (default None)
        The JSON decoding engine to use. One of:
          - if "json", parse JSON using built in json module
          - if "orjson", parse using the faster orjson module, requires the orjson
            package
          - if "auto" use orjson module if available, otherwise use the json module

        If not specified, the default engine is set to the current value of
        plotly.io.json.config.default_engine.

    Raises
    ------
    ValueError
        if value is not a string, or if skip_invalid=False and value contains
        invalid figure properties

    Returns
    -------
    Figure or FigureWidget
    """

    # Decode JSON
    # -----------
    fig_dict = from_json_plotly(value, engine=engine)

    # Validate coerce output type
    # ---------------------------
    cls = validate_coerce_output_type(output_type)

    # Create and return figure
    # ------------------------
    fig = cls(fig_dict, skip_invalid=skip_invalid)
    return fig


def read_json(file, output_type="Figure", skip_invalid=False, engine=None):
    """
    Construct a figure from the JSON contents of a local file or readable
    Python object.

    Note: A figure converted to JSON with one version of Plotly.py may not be compatible with another version.

    Parameters
    ----------
    file: str or readable
       A string containing the path to a local file or a read-able Python
       object (e.g. a pathlib.Path object or an open file descriptor)

    output_type: type or str (default 'Figure')
        The output figure type or type name.
        One of:  graph_objs.Figure, 'Figure', graph_objs.FigureWidget, 'FigureWidget'

    skip_invalid: bool (default False)
        False if invalid figure properties should result in an exception.
        True if invalid figure properties should be silently ignored.

    engine: str (default None)
        The JSON decoding engine to use. One of:
          - if "json", parse JSON using built in json module
          - if "orjson", parse using the faster orjson module, requires the orjson
            package
          - if "auto" use orjson module if available, otherwise use the json module

        If not specified, the default engine is set to the current value of
        plotly.io.json.config.default_engine.

    Returns
    -------
    Figure or FigureWidget
    """

    # Try to cast `file` as a pathlib object `path`.
    if isinstance(file, str):
        # Use the standard Path constructor to make a pathlib object.
        path = Path(file)
    elif isinstance(file, Path):
        # `file` is already a Path object.
        path = file
    else:
        # We could not make a Path object out of file. Either `file` is an open file
        # descriptor with a `write()` method or it's an invalid object.
        path = None

    # Read file contents into JSON string
    # -----------------------------------
    if path is not None:
        json_str = path.read_text()
    else:
        json_str = file.read()

    # Construct and return figure
    # ---------------------------
    return from_json(
        json_str, skip_invalid=skip_invalid, output_type=output_type, engine=engine
    )


def clean_to_json_compatible(obj, **kwargs):
    # Try handling value as a scalar value that we have a conversion for.
    # Return immediately if we know we've hit a primitive value

    # Bail out fast for simple scalar types
    if isinstance(obj, (int, float, str)):
        return obj

    if isinstance(obj, dict):
        return {k: clean_to_json_compatible(v, **kwargs) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)):
        if obj:
            # Must process list recursively even though it may be slow
            return [clean_to_json_compatible(v, **kwargs) for v in obj]

    # unpack kwargs
    numpy_allowed = kwargs.get("numpy_allowed", False)
    datetime_allowed = kwargs.get("datetime_allowed", False)

    modules = kwargs.get("modules", {})
    sage_all = modules["sage_all"]
    np = modules["np"]
    pd = modules["pd"]
    image = modules["image"]

    # Sage
    if sage_all is not None:
        if obj in sage_all.RR:
            return float(obj)
        elif obj in sage_all.ZZ:
            return int(obj)

    # numpy
    if np is not None:
        if obj is np.ma.core.masked:
            return float("nan")
        elif isinstance(obj, np.ndarray):
            if numpy_allowed and obj.dtype.kind in ("b", "i", "u", "f"):
                return np.ascontiguousarray(obj)
            elif obj.dtype.kind == "M":
                # datetime64 array
                return np.datetime_as_string(obj).tolist()
            elif obj.dtype.kind == "U":
                return obj.tolist()
            elif obj.dtype.kind == "O":
                # Treat object array as a lists, continue processing
                obj = obj.tolist()
        elif isinstance(obj, np.datetime64):
            return str(obj)

    # pandas
    if pd is not None:
        if obj is pd.NaT or obj is pd.NA:
            return None
        elif isinstance(obj, (pd.Series, pd.DatetimeIndex)):
            if numpy_allowed and obj.dtype.kind in ("b", "i", "u", "f"):
                return np.ascontiguousarray(obj.values)
            elif obj.dtype.kind == "M":
                if isinstance(obj, pd.Series):
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", FutureWarning)
                        # Series.dt.to_pydatetime will return Index[object]
                        # https://github.com/pandas-dev/pandas/pull/52459
                        dt_values = np.array(obj.dt.to_pydatetime()).tolist()
                else:  # DatetimeIndex
                    dt_values = obj.to_pydatetime().tolist()

                if not datetime_allowed:
                    # Note: We don't need to handle dropping timezones here because
                    # numpy's datetime64 doesn't support them and pandas's tz_localize
                    # above drops them.
                    for i in range(len(dt_values)):
                        dt_values[i] = dt_values[i].isoformat()

                return dt_values

    # datetime and date
    try:
        # Need to drop timezone for scalar datetimes. Don't need to convert
        # to string since engine can do that
        obj = obj.to_pydatetime()
    except (TypeError, AttributeError):
        pass

    if not datetime_allowed:
        try:
            return obj.isoformat()
        except (TypeError, AttributeError):
            pass
    elif isinstance(obj, datetime.datetime):
        return obj

    # Try .tolist() convertible, do not recurse inside
    try:
        return obj.tolist()
    except AttributeError:
        pass

    # Do best we can with decimal
    if isinstance(obj, decimal.Decimal):
        return float(obj)

    # PIL
    if image is not None and isinstance(obj, image.Image):
        return ImageUriValidator.pil_image_to_uri(obj)

    # Plotly
    try:
        obj = obj.to_plotly_json()
    except AttributeError:
        pass

    # Recurse into lists and dictionaries
    if isinstance(obj, dict):
        return {k: clean_to_json_compatible(v, **kwargs) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)):
        if obj:
            # Must process list recursively even though it may be slow
            return [clean_to_json_compatible(v, **kwargs) for v in obj]

    return obj
