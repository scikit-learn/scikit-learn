import base64
import numbers
import textwrap
import uuid
from importlib import import_module
import copy
import io
import re
import sys
import narwhals.stable.v1 as nw

from _plotly_utils.optional_imports import get_module


# back-port of fullmatch from Py3.4+
def fullmatch(regex, string, flags=0):
    """Emulate python-3.4 re.fullmatch()."""
    if "pattern" in dir(regex):
        regex_string = regex.pattern
    else:
        regex_string = regex
    return re.match("(?:" + regex_string + r")\Z", string, flags=flags)


# Utility functions
# -----------------
def to_scalar_or_list(v):
    # Handle the case where 'v' is a non-native scalar-like type,
    # such as numpy.float32. Without this case, the object might be
    # considered numpy-convertable and therefore promoted to a
    # 0-dimensional array, but we instead want it converted to a
    # Python native scalar type ('float' in the example above).
    # We explicitly check if is has the 'item' method, which conventionally
    # converts these types to native scalars.
    np = get_module("numpy", should_load=False)
    pd = get_module("pandas", should_load=False)
    if np and np.isscalar(v) and hasattr(v, "item"):
        return v.item()
    if isinstance(v, (list, tuple)):
        return [to_scalar_or_list(e) for e in v]
    elif np and isinstance(v, np.ndarray):
        if v.ndim == 0:
            return v.item()
        return [to_scalar_or_list(e) for e in v]
    elif pd and isinstance(v, (pd.Series, pd.Index)):
        return [to_scalar_or_list(e) for e in v]
    elif is_numpy_convertable(v):
        return to_scalar_or_list(np.array(v))
    else:
        return v


def copy_to_readonly_numpy_array(v, kind=None, force_numeric=False):
    """
    Convert an array-like value into a read-only numpy array

    Parameters
    ----------
    v : array like
        Array like value (list, tuple, numpy array, pandas series, etc.)
    kind : str or tuple of str
        If specified, the numpy dtype kind (or kinds) that the array should
        have, or be converted to if possible.
        If not specified then let numpy infer the datatype
    force_numeric : bool
        If true, raise an exception if the resulting numpy array does not
        have a numeric dtype (i.e. dtype.kind not in ['u', 'i', 'f'])
    Returns
    -------
    np.ndarray
        Numpy array with the 'WRITEABLE' flag set to False
    """
    np = get_module("numpy")

    assert np is not None

    # ### Process kind ###
    if not kind:
        kind = ()
    elif isinstance(kind, str):
        kind = (kind,)

    first_kind = kind[0] if kind else None

    # u: unsigned int, i: signed int, f: float
    numeric_kinds = {"u", "i", "f"}
    kind_default_dtypes = {
        "u": "uint32",
        "i": "int32",
        "f": "float64",
        "O": "object",
    }

    # With `pass_through=True`, the original object will be returned if unable to convert
    # to a Narwhals DataFrame or Series.
    v = nw.from_native(v, allow_series=True, pass_through=True)

    if isinstance(v, nw.Series):
        if v.dtype == nw.Datetime and v.dtype.time_zone is not None:
            # Remove time zone so that local time is displayed
            v = v.dt.replace_time_zone(None).to_numpy()
        else:
            v = v.to_numpy()
    elif isinstance(v, nw.DataFrame):
        schema = v.schema
        overrides = {}
        for key, val in schema.items():
            if val == nw.Datetime and val.time_zone is not None:
                # Remove time zone so that local time is displayed
                overrides[key] = nw.col(key).dt.replace_time_zone(None)
        if overrides:
            v = v.with_columns(**overrides)
        v = v.to_numpy()

    if not isinstance(v, np.ndarray):
        # v has its own logic on how to convert itself into a numpy array
        if is_numpy_convertable(v):
            return copy_to_readonly_numpy_array(
                np.array(v), kind=kind, force_numeric=force_numeric
            )
        else:
            # v is not homogenous array
            v_list = [to_scalar_or_list(e) for e in v]

            # Lookup dtype for requested kind, if any
            dtype = kind_default_dtypes.get(first_kind, None)

            # construct new array from list
            new_v = np.array(v_list, order="C", dtype=dtype)
    elif v.dtype.kind in numeric_kinds:
        # v is a homogenous numeric array
        if kind and v.dtype.kind not in kind:
            # Kind(s) were specified and this array doesn't match
            # Convert to the default dtype for the first kind
            dtype = kind_default_dtypes.get(first_kind, None)
            new_v = np.ascontiguousarray(v.astype(dtype))
        else:
            # Either no kind was requested or requested kind is satisfied
            new_v = np.ascontiguousarray(v.copy())
    else:
        # v is a non-numeric homogenous array
        new_v = v.copy()

    # Handle force numeric param
    # --------------------------
    if force_numeric and new_v.dtype.kind not in numeric_kinds:
        raise ValueError(
            "Input value is not numeric and force_numeric parameter set to True"
        )

    if "U" not in kind:
        # Force non-numeric arrays to have object type
        # --------------------------------------------
        # Here we make sure that non-numeric arrays have the object
        # datatype. This works around cases like np.array([1, 2, '3']) where
        # numpy converts the integers to strings and returns array of dtype
        # '<U21'
        if new_v.dtype.kind not in ["u", "i", "f", "O", "M"]:
            new_v = np.array(v, dtype="object")

    # Set new array to be read-only
    # -----------------------------
    new_v.flags["WRITEABLE"] = False

    return new_v


def is_numpy_convertable(v):
    """
    Return whether a value is meaningfully convertable to a numpy array
    via 'numpy.array'
    """
    return hasattr(v, "__array__") or hasattr(v, "__array_interface__")


def is_homogeneous_array(v):
    """
    Return whether a value is considered to be a homogeneous array
    """
    np = get_module("numpy", should_load=False)
    pd = get_module("pandas", should_load=False)
    if (
        np
        and isinstance(v, np.ndarray)
        or (pd and isinstance(v, (pd.Series, pd.Index)))
        or (isinstance(v, nw.Series))
    ):
        return True
    if is_numpy_convertable(v):
        np = get_module("numpy", should_load=True)
        if np:
            v_numpy = np.array(v)
            # v is essentially a scalar and so shouldn't count as an array
            if v_numpy.shape == ():
                return False
            else:
                return True  # v_numpy.dtype.kind in ["u", "i", "f", "M", "U"]
    return False


def is_simple_array(v):
    """
    Return whether a value is considered to be an simple array
    """
    return isinstance(v, (list, tuple))


def is_array(v):
    """
    Return whether a value is considered to be an array
    """
    return is_simple_array(v) or is_homogeneous_array(v)


def type_str(v):
    """
    Return a type string of the form module.name for the input value v
    """
    if not isinstance(v, type):
        v = type(v)

    return "'{module}.{name}'".format(module=v.__module__, name=v.__name__)


def is_typed_array_spec(v):
    """
    Return whether a value is considered to be a typed array spec for plotly.js
    """
    return isinstance(v, dict) and "bdata" in v and "dtype" in v


def is_none_or_typed_array_spec(v):
    return v is None or is_typed_array_spec(v)


# Validators
# ----------
class BaseValidator(object):
    """
    Base class for all validator classes
    """

    def __init__(self, plotly_name, parent_name, role=None, **_):
        """
        Construct a validator instance

        Parameters
        ----------
        plotly_name : str
            Name of the property being validated
        parent_name : str
            Names of all of the ancestors of this property joined on '.'
            characters. e.g.
            plotly_name == 'range' and parent_name == 'layout.xaxis'
        role : str
            The role string for the property as specified in
            plot-schema.json
        """
        self.parent_name = parent_name
        self.plotly_name = plotly_name
        self.role = role
        self.array_ok = False

    def description(self):
        """
        Returns a string that describes the values that are acceptable
        to the validator

        Should start with:
            The '{plotly_name}' property is a...

        For consistancy, string should have leading 4-space indent
        """
        raise NotImplementedError()

    def raise_invalid_val(self, v, inds=None):
        """
        Helper method to raise an informative exception when an invalid
        value is passed to the validate_coerce method.

        Parameters
        ----------
        v :
            Value that was input to validate_coerce and could not be coerced
        inds: list of int or None (default)
            Indexes to display after property name. e.g. if self.plotly_name
            is 'prop' and inds=[2, 1] then the name in the validation error
            message will be 'prop[2][1]`
        Raises
        -------
        ValueError
        """
        name = self.plotly_name
        if inds:
            for i in inds:
                name += "[" + str(i) + "]"

        raise ValueError(
            """
    Invalid value of type {typ} received for the '{name}' property of {pname}
        Received value: {v}

{valid_clr_desc}""".format(
                name=name,
                pname=self.parent_name,
                typ=type_str(v),
                v=repr(v),
                valid_clr_desc=self.description(),
            )
        )

    def raise_invalid_elements(self, invalid_els):
        if invalid_els:
            raise ValueError(
                """
    Invalid element(s) received for the '{name}' property of {pname}
        Invalid elements include: {invalid}

{valid_clr_desc}""".format(
                    name=self.plotly_name,
                    pname=self.parent_name,
                    invalid=invalid_els[:10],
                    valid_clr_desc=self.description(),
                )
            )

    def validate_coerce(self, v):
        """
        Validate whether an input value is compatible with this property,
        and coerce the value to be compatible of possible.

        Parameters
        ----------
        v
            The input value to be validated

        Raises
        ------
        ValueError
            if `v` cannot be coerced into a compatible form

        Returns
        -------
        The input `v` in a form that's compatible with this property
        """
        raise NotImplementedError()

    def present(self, v):
        """
        Convert output value of a previous call to `validate_coerce` into a
        form suitable to be returned to the user on upon property
        access.

        Note: The value returned by present must be either immutable or an
        instance of BasePlotlyType, otherwise the value could be mutated by
        the user and we wouldn't get notified about the change.

        Parameters
        ----------
        v
            A value that was the ouput of a previous call the
            `validate_coerce` method on the same object

        Returns
        -------

        """
        if is_homogeneous_array(v):
            # Note: numpy array was already coerced into read-only form so
            # we don't need to copy it here.
            return v
        elif is_simple_array(v):
            return tuple(v)
        else:
            return v


class DataArrayValidator(BaseValidator):
    """
    "data_array": {
        "description": "An {array} of data. The value MUST be an
                        {array}, or we ignore it.",
        "requiredOpts": [],
        "otherOpts": [
            "dflt"
        ]
    },
    """

    def __init__(self, plotly_name, parent_name, **kwargs):
        super(DataArrayValidator, self).__init__(
            plotly_name=plotly_name, parent_name=parent_name, **kwargs
        )

        self.array_ok = True

    def description(self):
        return """\
    The '{plotly_name}' property is an array that may be specified as a tuple,
    list, numpy array, or pandas Series""".format(plotly_name=self.plotly_name)

    def validate_coerce(self, v):
        if is_none_or_typed_array_spec(v):
            pass
        elif is_homogeneous_array(v):
            v = copy_to_readonly_numpy_array(v)
        elif is_simple_array(v):
            v = to_scalar_or_list(v)
        else:
            self.raise_invalid_val(v)
        return v


class EnumeratedValidator(BaseValidator):
    """
    "enumerated": {
        "description": "Enumerated value type. The available values are
                        listed in `values`.",
        "requiredOpts": [
            "values"
        ],
        "otherOpts": [
            "dflt",
            "coerceNumber",
            "arrayOk"
        ]
    },
    """

    def __init__(
        self,
        plotly_name,
        parent_name,
        values,
        array_ok=False,
        coerce_number=False,
        **kwargs,
    ):
        super(EnumeratedValidator, self).__init__(
            plotly_name=plotly_name, parent_name=parent_name, **kwargs
        )

        # Save params
        # -----------
        self.values = values
        self.array_ok = array_ok
        # coerce_number is rarely used and not implemented
        self.coerce_number = coerce_number
        self.kwargs = kwargs

        # Handle regular expressions
        # --------------------------
        # Compiled regexs
        self.val_regexs = []

        # regex replacements that run before the matching regex
        # So far, this is only used to cast 'x1' -> 'x' for anchor-style
        # enumeration properties
        self.regex_replacements = []

        # Loop over enumeration values
        # ----------------------------
        # Look for regular expressions
        for v in self.values:
            if v and isinstance(v, str) and v[0] == "/" and v[-1] == "/" and len(v) > 1:
                # String is a regex with leading and trailing '/' character
                regex_str = v[1:-1]
                self.val_regexs.append(re.compile(regex_str))
                self.regex_replacements.append(
                    EnumeratedValidator.build_regex_replacement(regex_str)
                )
            else:
                self.val_regexs.append(None)
                self.regex_replacements.append(None)

    def __deepcopy__(self, memodict={}):
        """
        A custom deepcopy method is needed here because compiled regex
        objects don't support deepcopy
        """
        cls = self.__class__
        return cls(self.plotly_name, self.parent_name, values=self.values)

    @staticmethod
    def build_regex_replacement(regex_str):
        # Example: regex_str == r"^y([2-9]|[1-9][0-9]+)?$"
        #
        # When we see a regular expression like the one above, we want to
        # build regular expression replacement params that will remove a
        # suffix of 1 from the input string ('y1' -> 'y' in this example)
        #
        # Why?: Regular expressions like this one are used in enumeration
        # properties that refer to subplotids (e.g. layout.annotation.xref)
        # The regular expressions forbid suffixes of 1, like 'x1'. But we
        # want to accept 'x1' and coerce it into 'x'
        #
        # To be cautious, we only perform this conversion for enumerated
        # values that match the anchor-style regex
        match = re.match(
            r"\^(\w)\(\[2\-9\]\|\[1\-9\]\[0\-9\]\+\)\?\( domain\)\?\$", regex_str
        )

        if match:
            anchor_char = match.group(1)
            return "^" + anchor_char + "1$", anchor_char
        else:
            return None

    def perform_replacemenet(self, v):
        """
        Return v with any applicable regex replacements applied
        """
        if isinstance(v, str):
            for repl_args in self.regex_replacements:
                if repl_args:
                    v = re.sub(repl_args[0], repl_args[1], v)

        return v

    def description(self):
        # Separate regular values from regular expressions
        enum_vals = []
        enum_regexs = []
        for v, regex in zip(self.values, self.val_regexs):
            if regex is not None:
                enum_regexs.append(regex.pattern)
            else:
                enum_vals.append(v)
        desc = """\
    The '{name}' property is an enumeration that may be specified as:""".format(
            name=self.plotly_name
        )

        if enum_vals:
            enum_vals_str = "\n".join(
                textwrap.wrap(
                    repr(enum_vals),
                    initial_indent=" " * 12,
                    subsequent_indent=" " * 12,
                    break_on_hyphens=False,
                )
            )

            desc = (
                desc
                + """
      - One of the following enumeration values:
{enum_vals_str}""".format(enum_vals_str=enum_vals_str)
            )

        if enum_regexs:
            enum_regexs_str = "\n".join(
                textwrap.wrap(
                    repr(enum_regexs),
                    initial_indent=" " * 12,
                    subsequent_indent=" " * 12,
                    break_on_hyphens=False,
                )
            )

            desc = (
                desc
                + """
      - A string that matches one of the following regular expressions:
{enum_regexs_str}""".format(enum_regexs_str=enum_regexs_str)
            )

        if self.array_ok:
            desc = (
                desc
                + """
      - A tuple, list, or one-dimensional numpy array of the above"""
            )

        return desc

    def in_values(self, e):
        """
        Return whether a value matches one of the enumeration options
        """
        is_str = isinstance(e, str)
        for v, regex in zip(self.values, self.val_regexs):
            if is_str and regex:
                in_values = fullmatch(regex, e) is not None
                # in_values = regex.fullmatch(e) is not None
            else:
                in_values = e == v

            if in_values:
                return True

        return False

    def validate_coerce(self, v):
        if is_none_or_typed_array_spec(v):
            pass
        elif self.array_ok and is_array(v):
            v_replaced = [self.perform_replacemenet(v_el) for v_el in v]

            invalid_els = [e for e in v_replaced if (not self.in_values(e))]
            if invalid_els:
                self.raise_invalid_elements(invalid_els[:10])

            if is_homogeneous_array(v):
                v = copy_to_readonly_numpy_array(v)
            else:
                v = to_scalar_or_list(v)
        else:
            v = self.perform_replacemenet(v)
            if not self.in_values(v):
                self.raise_invalid_val(v)
        return v


class BooleanValidator(BaseValidator):
    """
    "boolean": {
        "description": "A boolean (true/false) value.",
        "requiredOpts": [],
        "otherOpts": [
            "dflt"
        ]
    },
    """

    def __init__(self, plotly_name, parent_name, **kwargs):
        super(BooleanValidator, self).__init__(
            plotly_name=plotly_name, parent_name=parent_name, **kwargs
        )

    def description(self):
        return """\
    The '{plotly_name}' property must be specified as a bool
    (either True, or False)""".format(plotly_name=self.plotly_name)

    def validate_coerce(self, v):
        if is_none_or_typed_array_spec(v):
            pass
        elif not isinstance(v, bool):
            self.raise_invalid_val(v)

        return v


class SrcValidator(BaseValidator):
    def __init__(self, plotly_name, parent_name, **kwargs):
        super(SrcValidator, self).__init__(
            plotly_name=plotly_name, parent_name=parent_name, **kwargs
        )

        self.chart_studio = get_module("chart_studio")

    def description(self):
        return """\
    The '{plotly_name}' property must be specified as a string or
    as a plotly.grid_objs.Column object""".format(plotly_name=self.plotly_name)

    def validate_coerce(self, v):
        if is_none_or_typed_array_spec(v):
            pass
        elif isinstance(v, str):
            pass
        elif self.chart_studio and isinstance(v, self.chart_studio.grid_objs.Column):
            # Convert to id string
            v = v.id
        else:
            self.raise_invalid_val(v)

        return v


class NumberValidator(BaseValidator):
    """
    "number": {
        "description": "A number or a numeric value (e.g. a number
                        inside a string). When applicable, values
                        greater (less) than `max` (`min`) are coerced to
                        the `dflt`.",
        "requiredOpts": [],
        "otherOpts": [
            "dflt",
            "min",
            "max",
            "arrayOk"
        ]
    },
    """

    def __init__(
        self, plotly_name, parent_name, min=None, max=None, array_ok=False, **kwargs
    ):
        super(NumberValidator, self).__init__(
            plotly_name=plotly_name, parent_name=parent_name, **kwargs
        )

        # Handle min
        if min is None and max is not None:
            # Max was specified, so make min -inf
            self.min_val = float("-inf")
        else:
            self.min_val = min

        # Handle max
        if max is None and min is not None:
            # Min was specified, so make min inf
            self.max_val = float("inf")
        else:
            self.max_val = max

        if min is not None or max is not None:
            self.has_min_max = True
        else:
            self.has_min_max = False

        self.array_ok = array_ok

    def description(self):
        desc = """\
    The '{plotly_name}' property is a number and may be specified as:""".format(
            plotly_name=self.plotly_name
        )

        if not self.has_min_max:
            desc = (
                desc
                + """
      - An int or float"""
            )

        else:
            desc = (
                desc
                + """
      - An int or float in the interval [{min_val}, {max_val}]""".format(
                    min_val=self.min_val, max_val=self.max_val
                )
            )

        if self.array_ok:
            desc = (
                desc
                + """
      - A tuple, list, or one-dimensional numpy array of the above"""
            )

        return desc

    def validate_coerce(self, v):
        if is_none_or_typed_array_spec(v):
            pass
        elif self.array_ok and is_homogeneous_array(v):
            np = get_module("numpy")
            try:
                v_array = copy_to_readonly_numpy_array(v, force_numeric=True)
            except (ValueError, TypeError, OverflowError):
                self.raise_invalid_val(v)

            # Check min/max
            if self.has_min_max:
                v_valid = np.logical_and(
                    self.min_val <= v_array, v_array <= self.max_val
                )

                if not np.all(v_valid):
                    # Grab up to the first 10 invalid values
                    v_invalid = np.logical_not(v_valid)
                    some_invalid_els = np.array(v, dtype="object")[v_invalid][
                        :10
                    ].tolist()

                    self.raise_invalid_elements(some_invalid_els)

            v = v_array  # Always numeric numpy array
        elif self.array_ok and is_simple_array(v):
            # Check numeric
            invalid_els = [e for e in v if not isinstance(e, numbers.Number)]

            if invalid_els:
                self.raise_invalid_elements(invalid_els[:10])

            # Check min/max
            if self.has_min_max:
                invalid_els = [e for e in v if not (self.min_val <= e <= self.max_val)]

                if invalid_els:
                    self.raise_invalid_elements(invalid_els[:10])

            v = to_scalar_or_list(v)
        else:
            # Check numeric
            if not isinstance(v, numbers.Number):
                self.raise_invalid_val(v)

            # Check min/max
            if self.has_min_max:
                if not (self.min_val <= v <= self.max_val):
                    self.raise_invalid_val(v)
        return v


class IntegerValidator(BaseValidator):
    """
    "integer": {
        "description": "An integer or an integer inside a string. When
                        applicable, values greater (less) than `max`
                        (`min`) are coerced to the `dflt`.",
        "requiredOpts": [],
        "otherOpts": [
            "dflt",
            "min",
            "max",
            "extras",
            "arrayOk"
        ]
    },
    """

    def __init__(
        self,
        plotly_name,
        parent_name,
        min=None,
        max=None,
        extras=None,
        array_ok=False,
        **kwargs,
    ):
        super(IntegerValidator, self).__init__(
            plotly_name=plotly_name, parent_name=parent_name, **kwargs
        )

        # Handle min
        if min is None and max is not None:
            # Max was specified, so make min -inf
            self.min_val = -sys.maxsize - 1
        else:
            self.min_val = min

        # Handle max
        if max is None and min is not None:
            # Min was specified, so make min inf
            self.max_val = sys.maxsize
        else:
            self.max_val = max

        if min is not None or max is not None:
            self.has_min_max = True
        else:
            self.has_min_max = False

        self.extras = extras if extras is not None else []
        self.array_ok = array_ok

    def description(self):
        desc = """\
    The '{plotly_name}' property is a integer and may be specified as:""".format(
            plotly_name=self.plotly_name
        )

        if not self.has_min_max:
            desc = (
                desc
                + """
      - An int (or float that will be cast to an int)"""
            )
        else:
            desc = desc + (
                """
      - An int (or float that will be cast to an int)
        in the interval [{min_val}, {max_val}]""".format(
                    min_val=self.min_val, max_val=self.max_val
                )
            )

        # Extras
        if self.extras:
            desc = desc + (
                """
        OR exactly one of {extras} (e.g. '{eg_extra}')"""
            ).format(extras=self.extras, eg_extra=self.extras[-1])

        if self.array_ok:
            desc = (
                desc
                + """
      - A tuple, list, or one-dimensional numpy array of the above"""
            )

        return desc

    def validate_coerce(self, v):
        if is_none_or_typed_array_spec(v):
            pass
        elif v in self.extras:
            return v
        elif self.array_ok and is_homogeneous_array(v):
            np = get_module("numpy")
            v_array = copy_to_readonly_numpy_array(
                v, kind=("i", "u"), force_numeric=True
            )

            if v_array.dtype.kind not in ["i", "u"]:
                self.raise_invalid_val(v)

            # Check min/max
            if self.has_min_max:
                v_valid = np.logical_and(
                    self.min_val <= v_array, v_array <= self.max_val
                )

                if not np.all(v_valid):
                    # Grab up to the first 10 invalid values
                    v_invalid = np.logical_not(v_valid)
                    some_invalid_els = np.array(v, dtype="object")[v_invalid][
                        :10
                    ].tolist()
                    self.raise_invalid_elements(some_invalid_els)

            v = v_array
        elif self.array_ok and is_simple_array(v):
            # Check integer type
            invalid_els = [
                e for e in v if not isinstance(e, int) and e not in self.extras
            ]

            if invalid_els:
                self.raise_invalid_elements(invalid_els[:10])

            # Check min/max
            if self.has_min_max:
                invalid_els = [
                    e
                    for e in v
                    if not (isinstance(e, int) and self.min_val <= e <= self.max_val)
                    and e not in self.extras
                ]

                if invalid_els:
                    self.raise_invalid_elements(invalid_els[:10])

            v = to_scalar_or_list(v)
        else:
            # Check int
            if not isinstance(v, int):
                # don't let int() cast strings to ints
                self.raise_invalid_val(v)

            # Check min/max
            if self.has_min_max:
                if not (self.min_val <= v <= self.max_val):
                    self.raise_invalid_val(v)

        return v


class StringValidator(BaseValidator):
    """
    "string": {
        "description": "A string value. Numbers are converted to strings
                        except for attributes with `strict` set to true.",
        "requiredOpts": [],
        "otherOpts": [
            "dflt",
            "noBlank",
            "strict",
            "arrayOk",
            "values"
        ]
    },
    """

    def __init__(
        self,
        plotly_name,
        parent_name,
        no_blank=False,
        strict=False,
        array_ok=False,
        values=None,
        **kwargs,
    ):
        super(StringValidator, self).__init__(
            plotly_name=plotly_name, parent_name=parent_name, **kwargs
        )
        self.no_blank = no_blank
        self.strict = strict
        self.array_ok = array_ok
        self.values = values

    @staticmethod
    def to_str_or_unicode_or_none(v):
        """
        Convert a value to a string if it's not None, a string,
        or a unicode (on Python 2).
        """
        if v is None or isinstance(v, str):
            return v
        else:
            return str(v)

    def description(self):
        desc = """\
    The '{plotly_name}' property is a string and must be specified as:""".format(
            plotly_name=self.plotly_name
        )

        if self.no_blank:
            desc = (
                desc
                + """
      - A non-empty string"""
            )
        elif self.values:
            valid_str = "\n".join(
                textwrap.wrap(
                    repr(self.values),
                    initial_indent=" " * 12,
                    subsequent_indent=" " * 12,
                    break_on_hyphens=False,
                )
            )

            desc = (
                desc
                + """
      - One of the following strings:
{valid_str}""".format(valid_str=valid_str)
            )
        else:
            desc = (
                desc
                + """
      - A string"""
            )

        if not self.strict:
            desc = (
                desc
                + """
      - A number that will be converted to a string"""
            )

        if self.array_ok:
            desc = (
                desc
                + """
      - A tuple, list, or one-dimensional numpy array of the above"""
            )

        return desc

    def validate_coerce(self, v):
        if is_none_or_typed_array_spec(v):
            pass
        elif self.array_ok and is_array(v):
            # If strict, make sure all elements are strings.
            if self.strict:
                invalid_els = [e for e in v if not isinstance(e, str)]
                if invalid_els:
                    self.raise_invalid_elements(invalid_els)

            if is_homogeneous_array(v):
                np = get_module("numpy")

                # If not strict, let numpy cast elements to strings
                v = copy_to_readonly_numpy_array(v, kind="U")

                # Check no_blank
                if self.no_blank:
                    invalid_els = v[v == ""][:10].tolist()
                    if invalid_els:
                        self.raise_invalid_elements(invalid_els)

                # Check values
                if self.values:
                    invalid_inds = np.logical_not(np.isin(v, self.values))
                    invalid_els = v[invalid_inds][:10].tolist()
                    if invalid_els:
                        self.raise_invalid_elements(invalid_els)
            elif is_simple_array(v):
                if not self.strict:
                    v = [StringValidator.to_str_or_unicode_or_none(e) for e in v]

                # Check no_blank
                if self.no_blank:
                    invalid_els = [e for e in v if e == ""]
                    if invalid_els:
                        self.raise_invalid_elements(invalid_els)

                # Check values
                if self.values:
                    invalid_els = [e for e in v if v not in self.values]
                    if invalid_els:
                        self.raise_invalid_elements(invalid_els)

                v = to_scalar_or_list(v)

        else:
            if self.strict:
                if not isinstance(v, str):
                    self.raise_invalid_val(v)
            else:
                if isinstance(v, str):
                    pass
                elif isinstance(v, (int, float)):
                    # Convert value to a string
                    v = str(v)
                else:
                    self.raise_invalid_val(v)

            if self.no_blank and len(v) == 0:
                self.raise_invalid_val(v)

            if self.values and v not in self.values:
                self.raise_invalid_val(v)

        return v


class ColorValidator(BaseValidator):
    """
    "color": {
        "description": "A string describing color. Supported formats:
                        - hex (e.g. '#d3d3d3')
                        - rgb (e.g. 'rgb(255, 0, 0)')
                        - rgba (e.g. 'rgb(255, 0, 0, 0.5)')
                        - hsl (e.g. 'hsl(0, 100%, 50%)')
                        - hsv (e.g. 'hsv(0, 100%, 100%)')
                        - named colors(full list:
                          http://www.w3.org/TR/css3-color/#svg-color)",
        "requiredOpts": [],
        "otherOpts": [
            "dflt",
            "arrayOk"
        ]
    },
    """

    re_hex = re.compile(r"#([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})")
    re_rgb_etc = re.compile(r"(rgb|hsl|hsv)a?\([\d.]+%?(,[\d.]+%?){2,3}\)")
    re_ddk = re.compile(r"var\(\-\-.*\)")

    named_colors = [
        "aliceblue",
        "antiquewhite",
        "aqua",
        "aquamarine",
        "azure",
        "beige",
        "bisque",
        "black",
        "blanchedalmond",
        "blue",
        "blueviolet",
        "brown",
        "burlywood",
        "cadetblue",
        "chartreuse",
        "chocolate",
        "coral",
        "cornflowerblue",
        "cornsilk",
        "crimson",
        "cyan",
        "darkblue",
        "darkcyan",
        "darkgoldenrod",
        "darkgray",
        "darkgrey",
        "darkgreen",
        "darkkhaki",
        "darkmagenta",
        "darkolivegreen",
        "darkorange",
        "darkorchid",
        "darkred",
        "darksalmon",
        "darkseagreen",
        "darkslateblue",
        "darkslategray",
        "darkslategrey",
        "darkturquoise",
        "darkviolet",
        "deeppink",
        "deepskyblue",
        "dimgray",
        "dimgrey",
        "dodgerblue",
        "firebrick",
        "floralwhite",
        "forestgreen",
        "fuchsia",
        "gainsboro",
        "ghostwhite",
        "gold",
        "goldenrod",
        "gray",
        "grey",
        "green",
        "greenyellow",
        "honeydew",
        "hotpink",
        "indianred",
        "indigo",
        "ivory",
        "khaki",
        "lavender",
        "lavenderblush",
        "lawngreen",
        "lemonchiffon",
        "lightblue",
        "lightcoral",
        "lightcyan",
        "lightgoldenrodyellow",
        "lightgray",
        "lightgrey",
        "lightgreen",
        "lightpink",
        "lightsalmon",
        "lightseagreen",
        "lightskyblue",
        "lightslategray",
        "lightslategrey",
        "lightsteelblue",
        "lightyellow",
        "lime",
        "limegreen",
        "linen",
        "magenta",
        "maroon",
        "mediumaquamarine",
        "mediumblue",
        "mediumorchid",
        "mediumpurple",
        "mediumseagreen",
        "mediumslateblue",
        "mediumspringgreen",
        "mediumturquoise",
        "mediumvioletred",
        "midnightblue",
        "mintcream",
        "mistyrose",
        "moccasin",
        "navajowhite",
        "navy",
        "oldlace",
        "olive",
        "olivedrab",
        "orange",
        "orangered",
        "orchid",
        "palegoldenrod",
        "palegreen",
        "paleturquoise",
        "palevioletred",
        "papayawhip",
        "peachpuff",
        "peru",
        "pink",
        "plum",
        "powderblue",
        "purple",
        "red",
        "rosybrown",
        "royalblue",
        "rebeccapurple",
        "saddlebrown",
        "salmon",
        "sandybrown",
        "seagreen",
        "seashell",
        "sienna",
        "silver",
        "skyblue",
        "slateblue",
        "slategray",
        "slategrey",
        "snow",
        "springgreen",
        "steelblue",
        "tan",
        "teal",
        "thistle",
        "tomato",
        "turquoise",
        "violet",
        "wheat",
        "white",
        "whitesmoke",
        "yellow",
        "yellowgreen",
    ]

    def __init__(
        self, plotly_name, parent_name, array_ok=False, colorscale_path=None, **kwargs
    ):
        super(ColorValidator, self).__init__(
            plotly_name=plotly_name, parent_name=parent_name, **kwargs
        )

        self.array_ok = array_ok

        # colorscale_path is the path to the colorscale associated with this
        # color property, or None if no such colorscale exists. Only colors
        # with an associated colorscale may take on numeric values
        self.colorscale_path = colorscale_path

    def numbers_allowed(self):
        return self.colorscale_path is not None

    def description(self):
        valid_color_description = """\
    The '{plotly_name}' property is a color and may be specified as:
      - A hex string (e.g. '#ff0000')
      - An rgb/rgba string (e.g. 'rgb(255,0,0)')
      - An hsl/hsla string (e.g. 'hsl(0,100%,50%)')
      - An hsv/hsva string (e.g. 'hsv(0,100%,100%)')
      - A named CSS color: see https://plotly.com/python/css-colors/ for a list""".format(
            plotly_name=self.plotly_name
        )

        if self.colorscale_path:
            valid_color_description = (
                valid_color_description
                + """
      - A number that will be interpreted as a color
        according to {colorscale_path}""".format(colorscale_path=self.colorscale_path)
            )

        if self.array_ok:
            valid_color_description = (
                valid_color_description
                + """
      - A list or array of any of the above"""
            )

        return valid_color_description

    def validate_coerce(self, v, should_raise=True):
        if is_none_or_typed_array_spec(v):
            pass
        elif self.array_ok and is_homogeneous_array(v):
            v = copy_to_readonly_numpy_array(v)
            if self.numbers_allowed() and v.dtype.kind in ["u", "i", "f"]:
                # Numbers are allowed and we have an array of numbers.
                # All good
                pass
            else:
                validated_v = [self.validate_coerce(e, should_raise=False) for e in v]

                invalid_els = self.find_invalid_els(v, validated_v)

                if invalid_els and should_raise:
                    self.raise_invalid_elements(invalid_els)

                # ### Check that elements have valid colors types ###
                elif self.numbers_allowed() or invalid_els:
                    v = copy_to_readonly_numpy_array(validated_v, kind="O")
                else:
                    v = copy_to_readonly_numpy_array(validated_v, kind="U")
        elif self.array_ok and is_simple_array(v):
            validated_v = [self.validate_coerce(e, should_raise=False) for e in v]

            invalid_els = self.find_invalid_els(v, validated_v)

            if invalid_els and should_raise:
                self.raise_invalid_elements(invalid_els)
            else:
                v = validated_v
        else:
            # Validate scalar color
            validated_v = self.vc_scalar(v)
            if validated_v is None and should_raise:
                self.raise_invalid_val(v)

            v = validated_v

        return v

    def find_invalid_els(self, orig, validated, invalid_els=None):
        """
        Helper method to find invalid elements in orig array.
        Elements are invalid if their corresponding element in
        the validated array is None.

        This method handles deeply nested list structures
        """
        if invalid_els is None:
            invalid_els = []

        for orig_el, validated_el in zip(orig, validated):
            if is_array(orig_el):
                self.find_invalid_els(orig_el, validated_el, invalid_els)
            else:
                if validated_el is None:
                    invalid_els.append(orig_el)

        return invalid_els

    def vc_scalar(self, v):
        """Helper to validate/coerce a scalar color"""
        return ColorValidator.perform_validate_coerce(
            v, allow_number=self.numbers_allowed()
        )

    @staticmethod
    def perform_validate_coerce(v, allow_number=None):
        """
        Validate, coerce, and return a single color value. If input cannot be
        coerced to a valid color then return None.

        Parameters
        ----------
        v : number or str
            Candidate color value

        allow_number : bool
            True if numbers are allowed as colors

        Returns
        -------
        number or str or None
        """

        if isinstance(v, numbers.Number) and allow_number:
            # If allow_numbers then any number is ok
            return v
        elif not isinstance(v, str):
            # If not allow_numbers then value must be a string
            return None
        else:
            # Remove spaces so regexes don't need to bother with them.
            v_normalized = v.replace(" ", "").lower()

            # if ColorValidator.re_hex.fullmatch(v_normalized):
            if fullmatch(ColorValidator.re_hex, v_normalized):
                # valid hex color (e.g. #f34ab3)
                return v
            elif fullmatch(ColorValidator.re_rgb_etc, v_normalized):
                # elif ColorValidator.re_rgb_etc.fullmatch(v_normalized):
                # Valid rgb(a), hsl(a), hsv(a) color
                # (e.g. rgba(10, 234, 200, 50%)
                return v
            elif fullmatch(ColorValidator.re_ddk, v_normalized):
                # Valid var(--*) DDK theme variable, inspired by CSS syntax
                # (e.g. var(--accent) )
                # DDK will crawl & eval var(-- colors for Graph theming
                return v
            elif v_normalized in ColorValidator.named_colors:
                # Valid named color (e.g. 'coral')
                return v
            else:
                # Not a valid color
                return None


class ColorlistValidator(BaseValidator):
    """
    "colorlist": {
      "description": "A list of colors. Must be an {array} containing
                      valid colors.",
      "requiredOpts": [],
      "otherOpts": [
        "dflt"
      ]
    }
    """

    def __init__(self, plotly_name, parent_name, **kwargs):
        super(ColorlistValidator, self).__init__(
            plotly_name=plotly_name, parent_name=parent_name, **kwargs
        )

    def description(self):
        return """\
    The '{plotly_name}' property is a colorlist that may be specified
    as a tuple, list, one-dimensional numpy array, or pandas Series of valid
    color strings""".format(plotly_name=self.plotly_name)

    def validate_coerce(self, v):
        if is_none_or_typed_array_spec(v):
            pass
        elif is_array(v):
            validated_v = [
                ColorValidator.perform_validate_coerce(e, allow_number=False) for e in v
            ]

            invalid_els = [
                el for el, validated_el in zip(v, validated_v) if validated_el is None
            ]
            if invalid_els:
                self.raise_invalid_elements(invalid_els)

            v = to_scalar_or_list(v)
        else:
            self.raise_invalid_val(v)
        return v


class ColorscaleValidator(BaseValidator):
    """
    "colorscale": {
        "description": "A Plotly colorscale either picked by a name:
                        (any of Greys, YlGnBu, Greens, YlOrRd, Bluered,
                        RdBu, Reds, Blues, Picnic, Rainbow, Portland,
                        Jet, Hot, Blackbody, Earth, Electric, Viridis)
                        customized as an {array} of 2-element {arrays}
                        where the first element is the normalized color
                        level value (starting at *0* and ending at *1*),
                        and the second item is a valid color string.",
        "requiredOpts": [],
        "otherOpts": [
            "dflt"
        ]
    },
    """

    def __init__(self, plotly_name, parent_name, **kwargs):
        super(ColorscaleValidator, self).__init__(
            plotly_name=plotly_name, parent_name=parent_name, **kwargs
        )

        # named colorscales initialized on first use
        self._named_colorscales = None

    @property
    def named_colorscales(self):
        if self._named_colorscales is None:
            import inspect
            import itertools
            from plotly import colors

            colorscale_members = itertools.chain(
                inspect.getmembers(colors.sequential),
                inspect.getmembers(colors.diverging),
                inspect.getmembers(colors.cyclical),
            )

            self._named_colorscales = {
                c[0].lower(): c[1]
                for c in colorscale_members
                if isinstance(c, tuple)
                and len(c) == 2
                and isinstance(c[0], str)
                and isinstance(c[1], list)
                and not c[0].endswith("_r")
                and not c[0].startswith("_")
            }

        return self._named_colorscales

    def description(self):
        colorscales_str = "\n".join(
            textwrap.wrap(
                repr(sorted(list(self.named_colorscales))),
                initial_indent=" " * 12,
                subsequent_indent=" " * 13,
                break_on_hyphens=False,
                width=80,
            )
        )

        desc = """\
    The '{plotly_name}' property is a colorscale and may be
    specified as:
      - A list of colors that will be spaced evenly to create the colorscale.
        Many predefined colorscale lists are included in the sequential, diverging,
        and cyclical modules in the plotly.colors package.
      - A list of 2-element lists where the first element is the
        normalized color level value (starting at 0 and ending at 1),
        and the second item is a valid color string.
        (e.g. [[0, 'green'], [0.5, 'red'], [1.0, 'rgb(0, 0, 255)']])
      - One of the following named colorscales:
{colorscales_str}.
        Appending '_r' to a named colorscale reverses it.
""".format(plotly_name=self.plotly_name, colorscales_str=colorscales_str)

        return desc

    def validate_coerce(self, v):
        v_valid = False

        if v is None:
            v_valid = True
        elif isinstance(v, str):
            v_lower = v.lower()
            if v_lower in self.named_colorscales:
                # Convert to color list
                v = self.named_colorscales[v_lower]
                v_valid = True
            elif v_lower.endswith("_r") and v_lower[:-2] in self.named_colorscales:
                v = self.named_colorscales[v_lower[:-2]][::-1]
                v_valid = True
            #
            if v_valid:
                # Convert to list of lists colorscale
                d = len(v) - 1
                v = [[(1.0 * i) / (1.0 * d), x] for i, x in enumerate(v)]

        elif is_array(v) and len(v) > 0:
            # If firset element is a string, treat as colorsequence
            if isinstance(v[0], str):
                invalid_els = [
                    e for e in v if ColorValidator.perform_validate_coerce(e) is None
                ]

                if len(invalid_els) == 0:
                    v_valid = True

                    # Convert to list of lists colorscale
                    d = len(v) - 1
                    v = [[(1.0 * i) / (1.0 * d), x] for i, x in enumerate(v)]
            else:
                invalid_els = [
                    e
                    for e in v
                    if (
                        not is_array(e)
                        or len(e) != 2
                        or not isinstance(e[0], numbers.Number)
                        or not (0 <= e[0] <= 1)
                        or not isinstance(e[1], str)
                        or ColorValidator.perform_validate_coerce(e[1]) is None
                    )
                ]

                if len(invalid_els) == 0:
                    v_valid = True

                    # Convert to list of lists
                    v = [
                        [e[0], ColorValidator.perform_validate_coerce(e[1])] for e in v
                    ]

        if not v_valid:
            self.raise_invalid_val(v)

        return v

    def present(self, v):
        # Return-type must be immutable
        if v is None:
            return None
        elif isinstance(v, str):
            return v
        else:
            return tuple([tuple(e) for e in v])


class AngleValidator(BaseValidator):
    """
    "angle": {
        "description": "A number (in degree) between -180 and 180.",
        "requiredOpts": [],
        "otherOpts": [
            "dflt",
            "arrayOk"
        ]
    },
    """

    def __init__(self, plotly_name, parent_name, array_ok=False, **kwargs):
        super(AngleValidator, self).__init__(
            plotly_name=plotly_name, parent_name=parent_name, **kwargs
        )
        self.array_ok = array_ok

    def description(self):
        desc = """\
    The '{plotly_name}' property is a angle (in degrees) that may be
    specified as a number between -180 and 180{array_ok}.
    Numeric values outside this range are converted to the equivalent value
    (e.g. 270 is converted to -90).
        """.format(
            plotly_name=self.plotly_name,
            array_ok=(
                ", or a list, numpy array or other iterable thereof"
                if self.array_ok
                else ""
            ),
        )

        return desc

    def validate_coerce(self, v):
        if is_none_or_typed_array_spec(v):
            pass
        elif self.array_ok and is_homogeneous_array(v):
            try:
                v_array = copy_to_readonly_numpy_array(v, force_numeric=True)
            except (ValueError, TypeError, OverflowError):
                self.raise_invalid_val(v)
            v = v_array  # Always numeric numpy array
            # Normalize v onto the interval [-180, 180)
            v = (v + 180) % 360 - 180
        elif self.array_ok and is_simple_array(v):
            # Check numeric
            invalid_els = [e for e in v if not isinstance(e, numbers.Number)]

            if invalid_els:
                self.raise_invalid_elements(invalid_els[:10])

            v = [(x + 180) % 360 - 180 for x in to_scalar_or_list(v)]
        elif not isinstance(v, numbers.Number):
            self.raise_invalid_val(v)
        else:
            # Normalize v onto the interval [-180, 180)
            v = (v + 180) % 360 - 180

        return v


class SubplotidValidator(BaseValidator):
    """
    "subplotid": {
        "description": "An id string of a subplot type (given by dflt),
                        optionally followed by an integer >1. e.g. if
                        dflt='geo', we can have 'geo', 'geo2', 'geo3',
                        ...",
        "requiredOpts": [
            "dflt"
        ],
        "otherOpts": [
            "regex"
        ]
    }
    """

    def __init__(self, plotly_name, parent_name, dflt=None, regex=None, **kwargs):
        if dflt is None and regex is None:
            raise ValueError("One or both of regex and deflt must be specified")

        super(SubplotidValidator, self).__init__(
            plotly_name=plotly_name, parent_name=parent_name, **kwargs
        )

        if dflt is not None:
            self.base = dflt
        else:
            # e.g. regex == '/^y([2-9]|[1-9][0-9]+)?$/'
            self.base = re.match(r"/\^(\w+)", regex).group(1)

        self.regex = self.base + r"(\d*)"

    def description(self):
        desc = """\
    The '{plotly_name}' property is an identifier of a particular
    subplot, of type '{base}', that may be specified as the string '{base}'
    optionally followed by an integer >= 1
    (e.g. '{base}', '{base}1', '{base}2', '{base}3', etc.)
        """.format(plotly_name=self.plotly_name, base=self.base)
        return desc

    def validate_coerce(self, v):
        if v is None:
            pass
        elif not isinstance(v, str):
            self.raise_invalid_val(v)
        else:
            # match = re.fullmatch(self.regex, v)
            match = fullmatch(self.regex, v)
            if not match:
                is_valid = False
            else:
                digit_str = match.group(1)
                if len(digit_str) > 0 and int(digit_str) == 0:
                    is_valid = False
                elif len(digit_str) > 0 and int(digit_str) == 1:
                    # Remove 1 suffix (e.g. x1 -> x)
                    v = self.base
                    is_valid = True
                else:
                    is_valid = True

            if not is_valid:
                self.raise_invalid_val(v)
        return v


class FlaglistValidator(BaseValidator):
    """
    "flaglist": {
        "description": "A string representing a combination of flags
                        (order does not matter here). Combine any of the
                        available `flags` with *+*.
                        (e.g. ('lines+markers')). Values in `extras`
                        cannot be combined.",
        "requiredOpts": [
            "flags"
        ],
        "otherOpts": [
            "dflt",
            "extras",
            "arrayOk"
        ]
    },
    """

    def __init__(
        self, plotly_name, parent_name, flags, extras=None, array_ok=False, **kwargs
    ):
        super(FlaglistValidator, self).__init__(
            plotly_name=plotly_name, parent_name=parent_name, **kwargs
        )
        self.flags = flags
        self.extras = extras if extras is not None else []
        self.array_ok = array_ok

    def description(self):
        desc = (
            """\
    The '{plotly_name}' property is a flaglist and may be specified
    as a string containing:"""
        ).format(plotly_name=self.plotly_name)

        # Flags
        desc = desc + (
            """
      - Any combination of {flags} joined with '+' characters
        (e.g. '{eg_flag}')"""
        ).format(flags=self.flags, eg_flag="+".join(self.flags[:2]))

        # Extras
        if self.extras:
            desc = desc + (
                """
        OR exactly one of {extras} (e.g. '{eg_extra}')"""
            ).format(extras=self.extras, eg_extra=self.extras[-1])

        if self.array_ok:
            desc = (
                desc
                + """
      - A list or array of the above"""
            )

        return desc

    def vc_scalar(self, v):
        if isinstance(v, str):
            v = v.strip()

        if v in self.extras:
            return v

        if not isinstance(v, str):
            return None

        # To be generous we accept flags separated on plus ('+'),
        # or comma (',') and we accept whitespace around the flags
        split_vals = [e.strip() for e in re.split("[,+]", v)]

        # Are all flags valid names?
        if all(f in self.flags for f in split_vals):
            return "+".join(split_vals)
        else:
            return None

    def validate_coerce(self, v):
        if is_none_or_typed_array_spec(v):
            pass
        elif self.array_ok and is_array(v):
            # Coerce individual strings
            validated_v = [self.vc_scalar(e) for e in v]

            invalid_els = [
                el for el, validated_el in zip(v, validated_v) if validated_el is None
            ]
            if invalid_els:
                self.raise_invalid_elements(invalid_els)

            if is_homogeneous_array(v):
                v = copy_to_readonly_numpy_array(validated_v, kind="U")
            else:
                v = to_scalar_or_list(v)
        else:
            validated_v = self.vc_scalar(v)
            if validated_v is None:
                self.raise_invalid_val(v)

            v = validated_v

        return v


class AnyValidator(BaseValidator):
    """
    "any": {
        "description": "Any type.",
        "requiredOpts": [],
        "otherOpts": [
            "dflt",
            "values",
            "arrayOk"
        ]
    },
    """

    def __init__(self, plotly_name, parent_name, values=None, array_ok=False, **kwargs):
        super(AnyValidator, self).__init__(
            plotly_name=plotly_name, parent_name=parent_name, **kwargs
        )
        self.values = values
        self.array_ok = array_ok

    def description(self):
        desc = """\
    The '{plotly_name}' property accepts values of any type
        """.format(plotly_name=self.plotly_name)
        return desc

    def validate_coerce(self, v):
        if is_none_or_typed_array_spec(v):
            pass
        elif self.array_ok and is_homogeneous_array(v):
            v = copy_to_readonly_numpy_array(v, kind="O")
        elif self.array_ok and is_simple_array(v):
            v = to_scalar_or_list(v)
        return v


class InfoArrayValidator(BaseValidator):
    """
    "info_array": {
        "description": "An {array} of plot information.",
        "requiredOpts": [
            "items"
        ],
        "otherOpts": [
            "dflt",
            "freeLength",
            "dimensions"
        ]
    }
    """

    def __init__(
        self,
        plotly_name,
        parent_name,
        items,
        free_length=None,
        dimensions=None,
        **kwargs,
    ):
        super(InfoArrayValidator, self).__init__(
            plotly_name=plotly_name, parent_name=parent_name, **kwargs
        )

        self.items = items
        self.dimensions = dimensions if dimensions else 1
        self.free_length = free_length

        # Instantiate validators for each info array element
        self.item_validators = []
        info_array_items = self.items if isinstance(self.items, list) else [self.items]

        for i, item in enumerate(info_array_items):
            element_name = "{name}[{i}]".format(name=plotly_name, i=i)
            item_validator = InfoArrayValidator.build_validator(
                item, element_name, parent_name
            )
            self.item_validators.append(item_validator)

    def description(self):
        # Cases
        #  1) self.items is array, self.dimensions is 1
        #       a) free_length=True
        #       b) free_length=False
        #  2) self.items is array, self.dimensions is 2
        #     (requires free_length=True)
        #  3) self.items is scalar (requires free_length=True)
        #       a) dimensions=1
        #       b) dimensions=2
        #
        # dimensions can be set to '1-2' to indicate the both are accepted
        #
        desc = """\
    The '{plotly_name}' property is an info array that may be specified as:\
""".format(plotly_name=self.plotly_name)

        if isinstance(self.items, list):
            # ### Case 1 ###
            if self.dimensions in (1, "1-2"):
                upto = " up to" if self.free_length and self.dimensions == 1 else ""
                desc += """

    * a list or tuple of{upto} {N} elements where:\
""".format(upto=upto, N=len(self.item_validators))

                for i, item_validator in enumerate(self.item_validators):
                    el_desc = item_validator.description().strip()
                    desc = (
                        desc
                        + """
({i}) {el_desc}""".format(i=i, el_desc=el_desc)
                    )

            # ### Case 2 ###
            if self.dimensions in ("1-2", 2):
                assert self.free_length

                desc += """

    * a 2D list where:"""
                for i, item_validator in enumerate(self.item_validators):
                    # Update name for 2d
                    orig_name = item_validator.plotly_name
                    item_validator.plotly_name = "{name}[i][{i}]".format(
                        name=self.plotly_name, i=i
                    )

                    el_desc = item_validator.description().strip()
                    desc = (
                        desc
                        + """
({i}) {el_desc}""".format(i=i, el_desc=el_desc)
                    )
                    item_validator.plotly_name = orig_name
        else:
            # ### Case 3 ###
            assert self.free_length
            item_validator = self.item_validators[0]
            orig_name = item_validator.plotly_name

            if self.dimensions in (1, "1-2"):
                item_validator.plotly_name = "{name}[i]".format(name=self.plotly_name)

                el_desc = item_validator.description().strip()

                desc += """
    * a list of elements where:
      {el_desc}
""".format(el_desc=el_desc)

            if self.dimensions in ("1-2", 2):
                item_validator.plotly_name = "{name}[i][j]".format(
                    name=self.plotly_name
                )

                el_desc = item_validator.description().strip()
                desc += """
    * a 2D list where:
      {el_desc}
""".format(el_desc=el_desc)

            item_validator.plotly_name = orig_name

        return desc

    @staticmethod
    def build_validator(validator_info, plotly_name, parent_name):
        datatype = validator_info["valType"]  # type: str
        validator_classname = datatype.title().replace("_", "") + "Validator"
        validator_class = eval(validator_classname)

        kwargs = {
            k: validator_info[k]
            for k in validator_info
            if k not in ["valType", "description", "role"]
        }

        return validator_class(
            plotly_name=plotly_name, parent_name=parent_name, **kwargs
        )

    def validate_element_with_indexed_name(self, val, validator, inds):
        """
        Helper to add indexes to a validator's name, call validate_coerce on
        a value, then restore the original validator name.

        This makes sure that if a validation error message is raised, the
        property name the user sees includes the index(es) of the offending
        element.

        Parameters
        ----------
        val:
            A value to be validated
        validator
            A validator
        inds
            List of one or more non-negative integers that represent the
            nested index of the value being validated
        Returns
        -------
        val
            validated value

        Raises
        ------
        ValueError
            if val fails validation
        """
        orig_name = validator.plotly_name
        new_name = self.plotly_name
        for i in inds:
            new_name += "[" + str(i) + "]"
        validator.plotly_name = new_name
        try:
            val = validator.validate_coerce(val)
        finally:
            validator.plotly_name = orig_name

        return val

    def validate_coerce(self, v):
        if is_none_or_typed_array_spec(v):
            return None
        elif not is_array(v):
            self.raise_invalid_val(v)

        # Save off original v value to use in error reporting
        orig_v = v

        # Convert everything into nested lists
        # This way we don't need to worry about nested numpy arrays
        v = to_scalar_or_list(v)

        is_v_2d = v and is_array(v[0])

        if is_v_2d and self.dimensions in ("1-2", 2):
            if is_array(self.items):
                # e.g. 2D list as parcoords.dimensions.constraintrange
                # check that all items are there for each nested element
                for i, row in enumerate(v):
                    # Check row length
                    if not is_array(row) or len(row) != len(self.items):
                        self.raise_invalid_val(orig_v[i], [i])

                    for j, validator in enumerate(self.item_validators):
                        row[j] = self.validate_element_with_indexed_name(
                            v[i][j], validator, [i, j]
                        )
            else:
                # e.g. 2D list as layout.grid.subplots
                # check that all elements match individual validator
                validator = self.item_validators[0]
                for i, row in enumerate(v):
                    if not is_array(row):
                        self.raise_invalid_val(orig_v[i], [i])

                    for j, el in enumerate(row):
                        row[j] = self.validate_element_with_indexed_name(
                            el, validator, [i, j]
                        )
        elif v and self.dimensions == 2:
            # e.g. 1D list passed as layout.grid.subplots
            self.raise_invalid_val(orig_v[0], [0])
        elif not is_array(self.items):
            # e.g. 1D list passed as layout.grid.xaxes
            validator = self.item_validators[0]
            for i, el in enumerate(v):
                v[i] = self.validate_element_with_indexed_name(el, validator, [i])

        elif not self.free_length and len(v) != len(self.item_validators):
            # e.g. 3 element list as layout.xaxis.range
            self.raise_invalid_val(orig_v)
        elif self.free_length and len(v) > len(self.item_validators):
            # e.g. 4 element list as layout.updatemenu.button.args
            self.raise_invalid_val(orig_v)
        else:
            # We have a 1D array of the correct length
            for i, (el, validator) in enumerate(zip(v, self.item_validators)):
                # Validate coerce elements
                v[i] = validator.validate_coerce(el)

        return v

    def present(self, v):
        if v is None:
            return None
        else:
            if (
                self.dimensions == 2
                or self.dimensions == "1-2"
                and v
                and is_array(v[0])
            ):
                # 2D case
                v = copy.deepcopy(v)
                for row in v:
                    for i, (el, validator) in enumerate(zip(row, self.item_validators)):
                        row[i] = validator.present(el)

                return tuple(tuple(row) for row in v)
            else:
                # 1D case
                v = copy.copy(v)
                # Call present on each of the item validators
                for i, (el, validator) in enumerate(zip(v, self.item_validators)):
                    # Validate coerce elements
                    v[i] = validator.present(el)

                # Return tuple form of
                return tuple(v)


class LiteralValidator(BaseValidator):
    """
    Validator for readonly literal values
    """

    def __init__(self, plotly_name, parent_name, val, **kwargs):
        super(LiteralValidator, self).__init__(
            plotly_name=plotly_name, parent_name=parent_name, **kwargs
        )
        self.val = val

    def validate_coerce(self, v):
        if v != self.val:
            raise ValueError(
                """\
    The '{plotly_name}' property of {parent_name} is read-only""".format(
                    plotly_name=self.plotly_name, parent_name=self.parent_name
                )
            )
        else:
            return v


class DashValidator(EnumeratedValidator):
    """
    Special case validator for handling dash properties that may be specified
    as lists of dash lengths.  These are not currently specified in the
    schema.

    "dash": {
        "valType": "string",
        "values": [
            "solid",
            "dot",
            "dash",
            "longdash",
            "dashdot",
            "longdashdot"
        ],
        "dflt": "solid",
        "role": "style",
        "editType": "style",
        "description": "Sets the dash style of lines. Set to a dash type
        string (*solid*, *dot*, *dash*, *longdash*, *dashdot*, or
        *longdashdot*) or a dash length list in px (eg *5px,10px,2px,2px*)."
    },
    """

    def __init__(self, plotly_name, parent_name, values, **kwargs):
        # Add regex to handle dash length lists
        dash_list_regex = r"/^\d+(\.\d+)?(px|%)?((,|\s)\s*\d+(\.\d+)?(px|%)?)*$/"

        values = values + [dash_list_regex]

        # Call EnumeratedValidator superclass
        super(DashValidator, self).__init__(
            plotly_name=plotly_name, parent_name=parent_name, values=values, **kwargs
        )

    def description(self):
        # Separate regular values from regular expressions
        enum_vals = []
        enum_regexs = []
        for v, regex in zip(self.values, self.val_regexs):
            if regex is not None:
                enum_regexs.append(regex.pattern)
            else:
                enum_vals.append(v)
        desc = """\
    The '{name}' property is an enumeration that may be specified as:""".format(
            name=self.plotly_name
        )

        if enum_vals:
            enum_vals_str = "\n".join(
                textwrap.wrap(
                    repr(enum_vals),
                    initial_indent=" " * 12,
                    subsequent_indent=" " * 12,
                    break_on_hyphens=False,
                    width=80,
                )
            )

            desc = (
                desc
                + """
      - One of the following dash styles:
{enum_vals_str}""".format(enum_vals_str=enum_vals_str)
            )

        desc = (
            desc
            + """
      - A string containing a dash length list in pixels or percentages
            (e.g. '5px 10px 2px 2px', '5, 10, 2, 2', '10% 20% 40%', etc.)
"""
        )
        return desc


class ImageUriValidator(BaseValidator):
    _PIL = None

    try:
        _PIL = import_module("PIL")
    except ImportError:
        pass

    def __init__(self, plotly_name, parent_name, **kwargs):
        super(ImageUriValidator, self).__init__(
            plotly_name=plotly_name, parent_name=parent_name, **kwargs
        )

    def description(self):
        desc = """\
    The '{plotly_name}' property is an image URI that may be specified as:
      - A remote image URI string
        (e.g. 'http://www.somewhere.com/image.png')
      - A data URI image string
        (e.g. 'data:image/png;base64,iVBORw0KGgoAAAANSU')
      - A PIL.Image.Image object which will be immediately converted
        to a data URI image string
        See http://pillow.readthedocs.io/en/latest/reference/Image.html
        """.format(plotly_name=self.plotly_name)
        return desc

    def validate_coerce(self, v):
        if v is None:
            pass
        elif isinstance(v, str):
            # Future possibilities:
            #   - Detect filesystem system paths and convert to URI
            #   - Validate either url or data uri
            pass
        elif self._PIL and isinstance(v, self._PIL.Image.Image):
            # Convert PIL image to png data uri string
            v = self.pil_image_to_uri(v)
        else:
            self.raise_invalid_val(v)

        return v

    @staticmethod
    def pil_image_to_uri(v):
        in_mem_file = io.BytesIO()
        v.save(in_mem_file, format="PNG")
        in_mem_file.seek(0)
        img_bytes = in_mem_file.read()
        base64_encoded_result_bytes = base64.b64encode(img_bytes)
        base64_encoded_result_str = base64_encoded_result_bytes.decode("ascii")
        v = "data:image/png;base64,{base64_encoded_result_str}".format(
            base64_encoded_result_str=base64_encoded_result_str
        )
        return v


class CompoundValidator(BaseValidator):
    def __init__(self, plotly_name, parent_name, data_class_str, data_docs, **kwargs):
        super(CompoundValidator, self).__init__(
            plotly_name=plotly_name, parent_name=parent_name, **kwargs
        )

        # Save element class string
        self.data_class_str = data_class_str
        self._data_class = None
        self.data_docs = data_docs
        self.module_str = CompoundValidator.compute_graph_obj_module_str(
            self.data_class_str, parent_name
        )

    @staticmethod
    def compute_graph_obj_module_str(data_class_str, parent_name):
        if parent_name == "frame" and data_class_str in ["Data", "Layout"]:
            # Special case. There are no graph_objs.frame.Data or
            # graph_objs.frame.Layout classes. These are remapped to
            # graph_objs.Data and graph_objs.Layout

            parent_parts = parent_name.split(".")
            module_str = ".".join(["plotly.graph_objs"] + parent_parts[1:])
        elif parent_name == "layout.template" and data_class_str == "Layout":
            # Remap template's layout to regular layout
            module_str = "plotly.graph_objs"
        elif "layout.template.data" in parent_name:
            # Remap template's traces to regular traces
            parent_name = parent_name.replace("layout.template.data.", "")
            if parent_name:
                module_str = "plotly.graph_objs." + parent_name
            else:
                module_str = "plotly.graph_objs"
        elif parent_name:
            module_str = "plotly.graph_objs." + parent_name
        else:
            module_str = "plotly.graph_objs"

        return module_str

    @property
    def data_class(self):
        if self._data_class is None:
            module = import_module(self.module_str)
            self._data_class = getattr(module, self.data_class_str)

        return self._data_class

    def description(self):
        desc = (
            """\
    The '{plotly_name}' property is an instance of {class_str}
    that may be specified as:
      - An instance of :class:`{module_str}.{class_str}`
      - A dict of string/value properties that will be passed
        to the {class_str} constructor"""
        ).format(
            plotly_name=self.plotly_name,
            class_str=self.data_class_str,
            module_str=self.module_str,
        )

        return desc

    def validate_coerce(self, v, skip_invalid=False, _validate=True):
        if v is None:
            v = self.data_class()

        elif isinstance(v, dict):
            v = self.data_class(v, skip_invalid=skip_invalid, _validate=_validate)

        elif isinstance(v, self.data_class):
            # Copy object
            v = self.data_class(v)
        else:
            if skip_invalid:
                v = self.data_class()
            else:
                self.raise_invalid_val(v)

        v._plotly_name = self.plotly_name
        return v

    def present(self, v):
        # Return compound object as-is
        return v


class TitleValidator(CompoundValidator):
    """
    This is a special validator to allow compound title properties
    (e.g. layout.title, layout.xaxis.title, etc.) to be set as strings
    or numbers.  These strings are mapped to the 'text' property of the
    compound validator.
    """

    def __init__(self, *args, **kwargs):
        super(TitleValidator, self).__init__(*args, **kwargs)

    def validate_coerce(self, v, skip_invalid=False):
        if isinstance(v, (str, int, float)):
            v = {"text": v}
        return super(TitleValidator, self).validate_coerce(v, skip_invalid=skip_invalid)


class CompoundArrayValidator(BaseValidator):
    def __init__(self, plotly_name, parent_name, data_class_str, data_docs, **kwargs):
        super(CompoundArrayValidator, self).__init__(
            plotly_name=plotly_name, parent_name=parent_name, **kwargs
        )

        # Save element class string
        self.data_class_str = data_class_str
        self._data_class = None

        self.data_docs = data_docs
        self.module_str = CompoundValidator.compute_graph_obj_module_str(
            self.data_class_str, parent_name
        )

    def description(self):
        desc = (
            """\
    The '{plotly_name}' property is a tuple of instances of
    {class_str} that may be specified as:
      - A list or tuple of instances of {module_str}.{class_str}
      - A list or tuple of dicts of string/value properties that
        will be passed to the {class_str} constructor"""
        ).format(
            plotly_name=self.plotly_name,
            class_str=self.data_class_str,
            module_str=self.module_str,
        )

        return desc

    @property
    def data_class(self):
        if self._data_class is None:
            module = import_module(self.module_str)
            self._data_class = getattr(module, self.data_class_str)

        return self._data_class

    def validate_coerce(self, v, skip_invalid=False):
        if v is None:
            v = []

        elif isinstance(v, (list, tuple)):
            res = []
            invalid_els = []
            for v_el in v:
                if isinstance(v_el, self.data_class):
                    res.append(self.data_class(v_el))
                elif isinstance(v_el, dict):
                    res.append(self.data_class(v_el, skip_invalid=skip_invalid))
                else:
                    if skip_invalid:
                        res.append(self.data_class())
                    else:
                        res.append(None)
                        invalid_els.append(v_el)

            if invalid_els:
                self.raise_invalid_elements(invalid_els)

            v = to_scalar_or_list(res)
        else:
            if skip_invalid:
                v = []
            else:
                self.raise_invalid_val(v)

        return v

    def present(self, v):
        # Return compound object as tuple
        return tuple(v)


class BaseDataValidator(BaseValidator):
    def __init__(
        self, class_strs_map, plotly_name, parent_name, set_uid=False, **kwargs
    ):
        super(BaseDataValidator, self).__init__(
            plotly_name=plotly_name, parent_name=parent_name, **kwargs
        )

        self.class_strs_map = class_strs_map
        self._class_map = {}
        self.set_uid = set_uid

    def description(self):
        trace_types = str(list(self.class_strs_map.keys()))

        trace_types_wrapped = "\n".join(
            textwrap.wrap(
                trace_types,
                initial_indent="            One of: ",
                subsequent_indent=" " * 21,
                width=79 - 12,
            )
        )

        desc = (
            """\
    The '{plotly_name}' property is a tuple of trace instances
    that may be specified as:
      - A list or tuple of trace instances
        (e.g. [Scatter(...), Bar(...)])
      - A single trace instance
        (e.g. Scatter(...), Bar(...), etc.)
      - A list or tuple of dicts of string/value properties where:
        - The 'type' property specifies the trace type
{trace_types}

        - All remaining properties are passed to the constructor of
          the specified trace type

        (e.g. [{{'type': 'scatter', ...}}, {{'type': 'bar, ...}}])"""
        ).format(plotly_name=self.plotly_name, trace_types=trace_types_wrapped)

        return desc

    def get_trace_class(self, trace_name):
        # Import trace classes
        if trace_name not in self._class_map:
            trace_module = import_module("plotly.graph_objs")
            trace_class_name = self.class_strs_map[trace_name]
            self._class_map[trace_name] = getattr(trace_module, trace_class_name)

        return self._class_map[trace_name]

    def validate_coerce(self, v, skip_invalid=False, _validate=True):
        from plotly.basedatatypes import BaseTraceType

        # Import Histogram2dcontour, this is the deprecated name of the
        # Histogram2dContour trace.
        from plotly.graph_objs import Histogram2dcontour

        if v is None:
            v = []
        else:
            if not isinstance(v, (list, tuple)):
                v = [v]

            res = []
            invalid_els = []
            for v_el in v:
                if isinstance(v_el, BaseTraceType):
                    if isinstance(v_el, Histogram2dcontour):
                        v_el = dict(type="histogram2dcontour", **v_el._props)
                    else:
                        v_el = v_el._props

                if isinstance(v_el, dict):
                    type_in_v_el = "type" in v_el
                    trace_type = v_el.pop("type", "scatter")

                    if trace_type not in self.class_strs_map:
                        if skip_invalid:
                            # Treat as scatter trace
                            trace = self.get_trace_class("scatter")(
                                skip_invalid=skip_invalid, _validate=_validate, **v_el
                            )
                            res.append(trace)
                        else:
                            res.append(None)
                            invalid_els.append(v_el)
                    else:
                        trace = self.get_trace_class(trace_type)(
                            skip_invalid=skip_invalid, _validate=_validate, **v_el
                        )
                        res.append(trace)

                    if type_in_v_el:
                        # Restore type in v_el
                        v_el["type"] = trace_type
                else:
                    if skip_invalid:
                        # Add empty scatter trace
                        trace = self.get_trace_class("scatter")()
                        res.append(trace)
                    else:
                        res.append(None)
                        invalid_els.append(v_el)

            if invalid_els:
                self.raise_invalid_elements(invalid_els)

            v = to_scalar_or_list(res)

            # Set new UIDs
            if self.set_uid:
                for trace in v:
                    trace.uid = str(uuid.uuid4())

        return v


class BaseTemplateValidator(CompoundValidator):
    def __init__(self, plotly_name, parent_name, data_class_str, data_docs, **kwargs):
        super(BaseTemplateValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            data_class_str=data_class_str,
            data_docs=data_docs,
            **kwargs,
        )

    def description(self):
        compound_description = super(BaseTemplateValidator, self).description()
        compound_description += """
      - The name of a registered template where current registered templates
        are stored in the plotly.io.templates configuration object. The names
        of all registered templates can be retrieved with:
            >>> import plotly.io as pio
            >>> list(pio.templates)  # doctest: +ELLIPSIS
            ['ggplot2', 'seaborn', 'simple_white', 'plotly', 'plotly_white', ...]

      - A string containing multiple registered template names, joined on '+'
        characters (e.g. 'template1+template2'). In this case the resulting
        template is computed by merging together the collection of registered
        templates"""

        return compound_description

    def validate_coerce(self, v, skip_invalid=False):
        import plotly.io as pio

        try:
            # Check if v is a template identifier
            # (could be any hashable object)
            if v in pio.templates:
                return copy.deepcopy(pio.templates[v])
            # Otherwise, if v is a string, check to see if it consists of
            # multiple template names joined on '+' characters
            elif isinstance(v, str):
                template_names = v.split("+")
                if all([name in pio.templates for name in template_names]):
                    return pio.templates.merge_templates(*template_names)

        except TypeError:
            # v is un-hashable
            pass

        # Check for empty template
        if v == {} or isinstance(v, self.data_class) and v.to_plotly_json() == {}:
            # Replace empty template with {'data': {'scatter': [{}]}} so that we can
            # tell the difference between an un-initialized template and a template
            # explicitly set to empty.
            return self.data_class(data_scatter=[{}])

        return super(BaseTemplateValidator, self).validate_coerce(
            v, skip_invalid=skip_invalid
        )
