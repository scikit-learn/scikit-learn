import warnings
from contextlib import nullcontext
from functools import cached_property, partial, wraps
from inspect import getsourcefile, getsourcelines

import pytest

import numpydoc.tests
from numpydoc import validate
from numpydoc.docscrape import get_doc_object
from numpydoc.validate import Validator

validate_one = validate.validate

ALL_CHECKS = set(validate.ERROR_MSGS.keys())


@pytest.mark.parametrize(
    ["checks", "expected"],
    [
        [{"all"}, ALL_CHECKS],
        [set(), set()],
        [{"EX01"}, {"EX01"}],
        [{"EX01", "SA01"}, {"EX01", "SA01"}],
        [{"all", "EX01", "SA01"}, ALL_CHECKS - {"EX01", "SA01"}],
        [{"all", "PR01"}, ALL_CHECKS - {"PR01"}],
    ],
)
def test_utils_get_validation_checks(checks, expected):
    """Ensure check selection is working."""
    assert validate.get_validation_checks(checks) == expected


@pytest.mark.parametrize(
    "checks",
    [
        {"every"},
        {None},
        {"SM10"},
        {"EX01", "SM10"},
    ],
)
def test_get_validation_checks_validity(checks):
    """Ensure that invalid checks are flagged."""
    with pytest.raises(ValueError, match="Unrecognized validation code"):
        _ = validate.get_validation_checks(checks)


class _DummyList(list):
    """Dummy list class to test validation."""


def test_no_file():
    """Test that validation can be done on functions made on the fly."""
    # Just a smoke test for now, <list> will have a None filename
    validate.validate("numpydoc.tests.test_validate._DummyList.clear")


@pytest.mark.parametrize(
    ["file_contents", "expected"],
    [
        ["class MyClass:\n    pass", {}],
        ["class MyClass:  # numpydoc ignore=EX01\n    pass", {1: ["EX01"]}],
        [
            "class MyClass:  # numpydoc ignore= EX01,SA01\n    pass",
            {1: ["EX01", "SA01"]},
        ],
        [
            "class MyClass:\n    def my_method():  # numpydoc ignore:EX01\n        pass",
            {2: ["EX01"]},
        ],
        [
            "class MyClass:\n    def my_method():  # numpydoc ignore: EX01,PR01\n        pass",
            {2: ["EX01", "PR01"]},
        ],
        [
            "class MyClass:  # numpydoc ignore=GL08\n    def my_method():  # numpydoc ignore:EX01,PR01\n        pass",
            {1: ["GL08"], 2: ["EX01", "PR01"]},
        ],
    ],
)
def test_extract_ignore_validation_comments(tmp_path, file_contents, expected):
    """Test that extraction of validation ignore comments is working."""
    filepath = tmp_path / "ignore_comments.py"
    with open(filepath, "w") as file:
        file.write(file_contents)
    assert validate.extract_ignore_validation_comments(filepath) == expected


@pytest.mark.parametrize(
    "assumed_encoding",
    (
        pytest.param("utf-8", id="utf8_codec"),
        pytest.param("cp1252", id="cp1252_codec"),
    ),
)
@pytest.mark.parametrize(
    ("classname", "actual_encoding"),
    (
        pytest.param("MÿClass", "cp1252", id="cp1252_file"),
        pytest.param("My\u0081Class", "utf-8", id="utf8_file"),
    ),
)
def test_encodings(tmp_path, classname, actual_encoding, assumed_encoding):
    """Test handling of different source file encodings."""
    # write file as bytes with `actual_encoding`
    filepath = tmp_path / "ignore_comments.py"
    file_contents = f"class {classname}:\n    pass"
    with open(filepath, "wb") as file:
        file.write(file_contents.encode(actual_encoding))
    # this should fail on the ÿ in MÿClass. It represents the (presumed rare) case where
    # a user's editor saved the source file in cp1252 (or anything other than utf-8).
    if actual_encoding == "cp1252" and assumed_encoding == "utf-8":
        context = partial(
            pytest.raises,
            UnicodeDecodeError,
            match="can't decode byte 0xff in position 7: invalid start byte",
        )
    # this is the more likely case: file was utf-8 encoded, but Python on Windows uses
    # the system codepage to read the file. This case is fixed by numpy/numpydoc#510
    elif actual_encoding == "utf-8" and assumed_encoding == "cp1252":
        context = partial(
            pytest.raises,
            UnicodeDecodeError,
            match="can't decode byte 0x81 in position 9: character maps to <undefined>",
        )
    else:
        context = nullcontext
    with context():
        result = validate.extract_ignore_validation_comments(filepath, assumed_encoding)
        assert result == {}


class GoodDocStrings:
    """
    Collection of good doc strings.

    This class contains a lot of docstrings that should pass the validation
    script without any errors.

    See Also
    --------
    AnotherClass : With its description.

    Examples
    --------
    >>> result = 1 + 1
    """

    def one_liner(self):
        """Allow one liner docstrings (including quotes)."""
        # This should fail, but not because of the position of the quotes

    def plot(self, kind, color="blue", **kwargs):
        """
        Generate a plot.

        Render the data in the Series as a matplotlib plot of the
        specified kind.

        Parameters
        ----------
        kind : str
            Kind of matplotlib plot, e.g.::

                "foo"

        color : str, default 'blue'
            Color name or rgb code.
        **kwargs
            These parameters will be passed to the matplotlib plotting
            function.

        See Also
        --------
        related : Something related.

        Examples
        --------
        >>> result = 1 + 1
        """

    def swap(self, arr, i, j, *args, **kwargs):
        """
        Swap two indices on an array.

        The extended summary can be multiple paragraphs, but just one
        is enough to pass the validation.

        Parameters
        ----------
        arr : list
            The list having indexes swapped.
        i, j : int
            The indexes being swapped.
        *args, **kwargs
            Extraneous parameters are being permitted.

        See Also
        --------
        related : Something related.

        Examples
        --------
        >>> result = 1 + 1
        """

    def sample(self):
        """
        Generate and return a random number.

        The value is sampled from a continuous uniform distribution between
        0 and 1.

        Returns
        -------
        float
            Random number generated.

            - Make sure you set a seed for reproducibility

        See Also
        --------
        related : Something related.

        Examples
        --------
        >>> result = 1 + 1
        """

    def random_letters(self):
        """
        Generate and return a sequence of random letters.

        The length of the returned string is also random, and is also
        returned.

        Returns
        -------
        length : int
            Length of the returned string.
        letters : str
            String of random letters.

            .. versionadded:: 0.1

        See Also
        --------
        related : Something related.

        Examples
        --------
        >>> result = 1 + 1
        """

    def sample_values(self):
        """
        Generate an infinite sequence of random numbers.

        The values are sampled from a continuous uniform distribution between
        0 and 1.

        Yields
        ------
        float
            Random number generated.

        See Also
        --------
        related : Something related.

        Examples
        --------
        >>> result = 1 + 1
        """

    def head(self):
        """
        Return the first 5 elements of the Series.

        This function is mainly useful to preview the values of the
        Series without displaying the whole of it.

        Returns
        -------
        int
            Subset of the original series with the 5 first values.

        See Also
        --------
        Series.tail : Return the last 5 elements of the Series.
        Series.iloc : Return a slice of the elements in the Series,
            which can also be used to return the first or last n.

        Examples
        --------
        >>> 1 + 1
        2
        """
        return 1

    def head1(self, n=5):
        """
        Return the first elements of the Series.

        This function is mainly useful to preview the values of the
        Series without displaying the whole of it.

        Parameters
        ----------
        n : int
            Number of values to return.

        Returns
        -------
        int
            Subset of the original series with the n first values.

        See Also
        --------
        tail : Return the last n elements of the Series.

        Examples
        --------
        >>> s = 10
        >>> s
        10

        With the `n` parameter, we can change the number of returned rows:

        >>> s + 1
        11
        """
        return 1

    def summary_starts_with_number(self, n=5):
        """
        2nd rule of summaries should allow this.

        3 Starting the summary with a number instead of a capital letter.
        Also in parameters, returns, see also...

        Parameters
        ----------
        n : int
            4 Number of values to return.

        Returns
        -------
        int
            5 Subset of the original series with the n first values.

        See Also
        --------
        tail : 6 Return the last n elements of the Series.

        Examples
        --------
        >>> s = 10
        >>> s
        10

        7 With the `n` parameter, we can change the number of returned rows:

        >>> s + 1
        11
        """
        return 1

    def contains(self, pat, case=True, na=float("NaN")):
        """
        Return whether each value contains `pat`.

        In this case, we are illustrating how to use sections, even
        if the example is simple enough and does not require them.

        Parameters
        ----------
        pat : str
            Pattern to check for within each element.
        case : bool, default True
            Whether check should be done with case sensitivity.
        na : object, default np.nan
            Fill value for missing data.

        See Also
        --------
        related : Something related.

        Examples
        --------
        >>> s = 25
        >>> s
        25

        **Case sensitivity**

        With `case_sensitive` set to `False` we can match `a` with both
        `a` and `A`:

        >>> s + 1
        26

        **Missing values**

        We can fill missing values in the output using the `na` parameter:

        >>> s * 2
        50
        """

    def mode(self, axis, numeric_only):
        """
        Ensure reST directives don't affect checks for leading periods.

        The extended summary can be multiple paragraphs, but just one
        is enough to pass the validation.

        Parameters
        ----------
        axis : str
            Sentence ending in period, followed by single directive.

            .. versionchanged:: 0.1.2

        numeric_only : bool
            Sentence ending in period, followed by multiple directives.

            .. versionadded:: 0.1.2
            .. deprecated:: 0.00.0
                A multiline description,
                which spans another line.

        See Also
        --------
        related : Something related.

        Examples
        --------
        >>> result = 1 + 1
        """

    def good_imports(self):
        """
        Ensure import other than numpy and pandas are fine.

        The extended summary can be multiple paragraphs, but just one
        is enough to pass the validation.

        See Also
        --------
        related : Something related.

        Examples
        --------
        This example does not import pandas or import numpy.
        >>> import datetime
        >>> datetime.MAXYEAR
        9999
        """

    def no_returns(self):
        """
        Say hello and have no returns.

        The extended summary can be multiple paragraphs, but just one
        is enough to pass the validation.

        See Also
        --------
        related : Something related.

        Examples
        --------
        >>> result = 1 + 1
        """

    def empty_returns(self):
        """
        Say hello and always return None.

        Since this function never returns a value, this
        docstring doesn't need a return section.

        See Also
        --------
        related : Something related.

        Examples
        --------
        >>> result = 1 + 1
        """

        def say_hello():
            return "Hello World!"

        say_hello()
        if True:
            return
        else:
            return

    def warnings(self):
        """
        Do one thing.

        Sometimes, this function does other things.

        Warnings
        --------
        This function may produce side effects when some condition
        is met.

        See Also
        --------
        related : Something related.

        Examples
        --------
        >>> result = 1 + 1
        """

    def multiple_variables_on_one_line(self, matrix, a, b, i, j):
        """
        Swap two values in a matrix.

        The extended summary can be multiple paragraphs, but just one
        is enough to pass the validation.

        Parameters
        ----------
        matrix : list of list
            A double list that represents a matrix.
        a, b : int
            The indices of the first value.
        i, j : int
            The indices of the second value.

        See Also
        --------
        related : Something related.

        Examples
        --------
        >>> result = 1 + 1
        """

    def other_parameters(self, param1, param2):
        """
        Ensure "Other Parameters" are recognized.

        The second parameter is used infrequently, so it is put in the
        "Other Parameters" section.

        Parameters
        ----------
        param1 : bool
            Description of commonly used parameter.

        Other Parameters
        ----------------
        param2 : str
            Description of infrequently used parameter.

        See Also
        --------
        related : Something related.

        Examples
        --------
        >>> result = 1 + 1
        """

    def valid_options_in_parameter_description_sets(self, bar):
        """
        Ensure a PR06 error is not raised when type is member of a set.

        Literal keywords like 'integer' are valid when specified in a set of
        valid options for a keyword parameter.

        Parameters
        ----------
        bar : {'integer', 'boolean'}
            The literal values of 'integer' and 'boolean' are part of an
            options set and thus should not be subject to PR06 warnings.

        See Also
        --------
        related : Something related.

        Examples
        --------
        >>> result = 1 + 1
        """

    def parameters_with_trailing_underscores(self, str_):
        r"""
        Ensure PR01 and PR02 errors are not raised with trailing underscores.

        Parameters with trailing underscores need to be escaped to render
        properly in the documentation since trailing underscores are used to
        create links. Doing so without also handling the change in the validation
        logic makes it impossible to both pass validation and render correctly.

        Parameters
        ----------
        str\_ : str
           Some text.

        See Also
        --------
        related : Something related.

        Examples
        --------
        >>> result = 1 + 1
        """

    def parameter_with_wrong_types_as_substrings(self, a, b, c, d, e, f):
        r"""
        Ensure PR06 doesn't fail when non-preferable types are substrings.

        While PR06 checks for parameter types which contain non-preferable type
        names like integer (int), string (str), and boolean (bool), PR06 should
        not fail if those types are used only as susbtrings in, for example,
        custom type names.

        Parameters
        ----------
        a : Myint
           Some text.
        b : intClass
           Some text.
        c : Mystring
           Some text.
        d : stringClass
           Some text.
        e : Mybool
           Some text.
        f : boolClass
           Some text.

        See Also
        --------
        related : Something related.

        Examples
        --------
        >>> result = 1 + 1
        """


class BadGenericDocStrings:
    """Everything here has a bad docstring"""

    def func(self):
        """Some function.

        With several mistakes in the docstring.

        It has a blank like after the signature `def func():`.

        The text 'Some function' should go in the line after the
        opening quotes of the docstring, not in the same line.

        There is a blank line between the docstring and the first line
        of code `foo = 1`.

        The closing quotes should be in the next line, not in this one."""

        foo = 1
        bar = 2
        return foo + bar

    def astype(self, dtype):
        """
        Casts Series type.

        Verb in third-person of the present simple, should be infinitive.
        """

    def astype1(self, dtype):
        """
        Method to cast Series type.

        Does not start with verb.
        """

    def astype2(self, dtype):
        """
        Cast Series type

        Missing dot at the end.
        """

    def astype3(self, dtype):
        """
        Cast Series type from its current type to the new type defined in
        the parameter dtype.

        Summary is too verbose and doesn't fit in a single line.
        """

    def two_linebreaks_between_sections(self, foo):
        """
        Test linebreaks message GL03.

        Note 2 blank lines before parameters section.


        Parameters
        ----------
        foo : str
            Description of foo parameter.
        """

    def linebreak_at_end_of_docstring(self, foo):
        """
        Test linebreaks message GL03.

        Note extra blank line at end of docstring.

        Parameters
        ----------
        foo : str
            Description of foo parameter.

        """

    def plot(self, kind, **kwargs):
        """
        Generate a plot.

        Render the data in the Series as a matplotlib plot of the
        specified kind.

        Note the blank line between the parameters title and the first
        parameter. Also, note that after the name of the parameter `kind`
        and before the colon, a space is missing.

        Also, note that the parameter descriptions do not start with a
        capital letter, and do not finish with a dot.

        Finally, the `**kwargs` parameter is missing.

        Parameters
        ----------

        kind: str
            kind of matplotlib plot
        """

    def unknown_section(self):
        """
        This section has an unknown section title.

        Unknown Section
        ---------------
        This should raise an error in the validation.
        """

    def sections_in_wrong_order(self):
        """
        This docstring has the sections in the wrong order.

        Parameters
        ----------
        name : str
            This section is in the right position.

        Examples
        --------
        >>> print("So far Examples is good, as it goes before Parameters")
        So far Examples is good, as it goes before Parameters

        See Also
        --------
        function : This should generate an error, as See Also needs to go
            before Examples.
        """

    def deprecation_in_wrong_order(self):
        """
        This docstring has the deprecation warning in the wrong order.

        This is the extended summary. The correct order should be
        summary, deprecation warning, extended summary.

        .. deprecated:: 1.0
            This should generate an error as it needs to go before
            extended summary.
        """

    def method_wo_docstrings(self):
        pass

    def directives_without_two_colons(self, first, second):
        """
        Ensure reST directives have trailing colons.

        Parameters
        ----------
        first : str
            Sentence ending in period, followed by single directive w/o colons.

            .. versionchanged 0.1.2

        second : bool
            Sentence ending in period, followed by multiple directives w/o
            colons.

            .. versionadded 0.1.2
            .. deprecated 0.00.0

        """


class WarnGenericFormat:
    """
    Those contains things that _may_ be incorrect formatting.
    """

    def too_short_header_underline(self, a, b):
        """
        The header line is too short.

        Parameters
        ------
        a, b : int
            Foo bar baz.
        """


class BadSummaries:
    def no_summary(self):
        """
        Returns
        -------
        int
            Always one.
        """

    def heading_whitespaces(self):
        """
           Summary with heading whitespaces.

        Returns
        -------
        int
            Always one.
        """

    def wrong_line(self):
        """Quotes are on the wrong line.

        Both opening and closing."""

    def no_punctuation(self):
        """
        Has the right line but forgets punctuation
        """

    def no_capitalization(self):
        """
        provides a lowercase summary.
        """

    def no_infinitive(self):
        """
        Started with a verb that is not infinitive.
        """

    def multi_line(self):
        """
        Extends beyond one line
        which is not correct.
        """

    def two_paragraph_multi_line(self):
        """
        Extends beyond one line
        which is not correct.

        Extends beyond one line, which in itself is correct but the
        previous short summary should still be an issue.
        """


class BadParameters:
    """
    Everything here has a problem with its Parameters section.
    """

    def no_type(self, value):
        """
        Lacks the type.

        Parameters
        ----------
        value
            A parameter without type.
        """

    def type_with_period(self, value):
        """
        Has period after type.

        Parameters
        ----------
        value : str.
            A parameter type should not finish with period.
        """

    def no_description(self, value):
        """
        Lacks the description.

        Parameters
        ----------
        value : str
        """

    def missing_params(self, kind, **kwargs):
        """
        Lacks kwargs in Parameters.

        Parameters
        ----------
        kind : str
            Foo bar baz.
        """

    def bad_colon_spacing(self, kind):
        """
        Has bad spacing in the type line.

        Parameters
        ----------
        kind: str
            Needs a space after kind.
        """

    def no_description_period(self, kind):
        """
        Forgets to add a period to the description.

        Parameters
        ----------
        kind : str
           Doesn't end with a dot
        """

    def no_description_period_with_directive(self, kind):
        """
        Forgets to add a period, and also includes a directive.

        Parameters
        ----------
        kind : str
           Doesn't end with a dot

           .. versionadded:: 0.00.0
        """

    def no_description_period_with_directives(self, kind):
        """
        Forgets to add a period, and also includes multiple directives.

        Parameters
        ----------
        kind : str
           Doesn't end with a dot

           .. versionchanged:: 0.00.0
           .. deprecated:: 0.00.0
        """

    def parameter_capitalization(self, kind):
        """
        Forgets to capitalize the description.

        Parameters
        ----------
        kind : str
           this is not capitalized.
        """

    def blank_lines(self, kind):
        """
        Adds a blank line after the section header.

        Parameters
        ----------

        kind : str
            Foo bar baz.
        """

    def integer_parameter(self, kind):
        """
        Uses integer instead of int.

        Parameters
        ----------
        kind : integer
            Foo bar baz.
        """

    def string_parameter(self, kind):
        """
        Uses string instead of str.

        Parameters
        ----------
        kind : string
            Foo bar baz.
        """

    def boolean_parameter(self, kind):
        """
        Uses boolean instead of bool.

        Parameters
        ----------
        kind : boolean
            Foo bar baz.
        """

    def list_incorrect_parameter_type(self, kind):
        """
        Uses list of boolean instead of list of bool.

        Parameters
        ----------
        kind : list of boolean, integer, float or string
            Foo bar baz.
        """

    def bad_parameter_spacing(self, a, b):
        """
        The parameters on the same line have an extra space between them.

        Parameters
        ----------
        a,  b : int
            Foo bar baz.
        """


class BadReturns:
    def return_not_documented(self):
        """
        Lacks section for Returns
        """
        return "Hello world!"

    def yield_not_documented(self):
        """
        Lacks section for Yields
        """
        yield "Hello world!"

    def no_type(self):
        """
        Returns documented but without type.

        Returns
        -------
        Some value.
        """
        return "Hello world!"

    def no_description(self):
        """
        Provides type but no description.

        Returns
        -------
        str
        """
        return "Hello world!"

    def no_punctuation(self):
        """
        Provides type and description but no period.

        Returns
        -------
        str
           A nice greeting
        """
        return "Hello world!"

    def named_single_return(self):
        """
        Provides name but returns only one value.

        Returns
        -------
        s : str
           A nice greeting.
        """
        return "Hello world!"

    def no_capitalization(self):
        """
        Forgets capitalization in return values description.

        Returns
        -------
        foo : str
           The first returned string.
        bar : str
           the second returned string.
        """
        return "Hello", "World!"

    def no_period_multi(self):
        """
        Forgets period in return values description.

        Returns
        -------
        foo : str
           The first returned string
        bar : str
           The second returned string.
        """
        return "Hello", "World!"


class BadSeeAlso:
    def no_desc(self):
        """
        Return the first 5 elements of the Series.

        See Also
        --------
        Series.tail
        """

    def desc_no_period(self):
        """
        Return the first 5 elements of the Series.

        See Also
        --------
        Series.tail : Return the last 5 elements of the Series.
        Series.iloc : Return a slice of the elements in the Series,
            which can also be used to return the first or last n
        """

    def desc_first_letter_lowercase(self):
        """
        Return the first 5 elements of the Series.

        See Also
        --------
        Series.tail : return the last 5 elements of the Series.
        Series.iloc : Return a slice of the elements in the Series,
            which can also be used to return the first or last n.
        """

    def prefix_pandas(self):
        """
        Have `pandas` prefix in See Also section.

        See Also
        --------
        pandas.Series.rename : Alter Series index labels or name.
        DataFrame.head : The first `n` rows of the caller object.
        """


class BadExamples:
    def missing_whitespace_around_arithmetic_operator(self):
        """
        Examples
        --------
        >>> 2 + 5
        7
        """

    def indentation_is_not_a_multiple_of_four(self):
        """
        Examples
        --------
        >>> if 2 + 5:
        ...     pass
        """

    def missing_whitespace_after_comma(self):
        """
        Examples
        --------
        >>> import datetime
        >>> value = datetime.date(2019, 1, 1)
        """


class TestValidator:
    def _import_path(self, klass=None, func=None):
        """
        Build the required import path for tests in this module.

        Parameters
        ----------
        klass : str
            Class name of object in module.
        func : str
            Function name of object in module.

        Returns
        -------
        str
            Import path of specified object in this module
        """
        base_path = "numpydoc.tests.test_validate"

        if klass:
            base_path = ".".join([base_path, klass])

        if func:
            base_path = ".".join([base_path, func])

        return base_path

    def test_one_liner(self, capsys):
        result = validate_one(
            self._import_path(klass="GoodDocStrings", func="one_liner")
        )
        errors = " ".join(err[1] for err in result["errors"])
        assert (
            "should start in the line immediately after the opening quotes"
            not in errors
        )
        assert "should be placed in the line after the last text" not in errors

    def test_good_class(self, capsys):
        errors = validate_one(self._import_path(klass="GoodDocStrings"))["errors"]
        assert isinstance(errors, list)
        assert not errors

    @pytest.mark.parametrize(
        "func",
        [
            "plot",
            "swap",
            "sample",
            "random_letters",
            "sample_values",
            "head",
            "head1",
            "summary_starts_with_number",
            "contains",
            "mode",
            "good_imports",
            "no_returns",
            "empty_returns",
            "multiple_variables_on_one_line",
            "other_parameters",
            "warnings",
            "valid_options_in_parameter_description_sets",
            "parameters_with_trailing_underscores",
            "parameter_with_wrong_types_as_substrings",
        ],
    )
    def test_good_functions(self, capsys, func):
        errors = validate_one(self._import_path(klass="GoodDocStrings", func=func))[
            "errors"
        ]
        assert isinstance(errors, list)
        assert not errors

    def test_bad_class(self, capsys):
        errors = validate_one(self._import_path(klass="BadGenericDocStrings"))["errors"]
        assert isinstance(errors, list)
        assert errors

    @pytest.mark.parametrize(
        "func",
        [
            "too_short_header_underline",
        ],
    )
    def test_bad_generic_functions(self, capsys, func):
        with pytest.warns(UserWarning):
            errors = validate_one(
                self._import_path(klass="WarnGenericFormat", func=func)
            )
        assert "is too short" in w.msg

    @pytest.mark.parametrize(
        "func",
        [
            "func",
            "astype",
            "astype1",
            "astype2",
            "astype3",
            "plot",
            "directives_without_two_colons",
        ],
    )
    def test_bad_generic_functions(self, capsys, func):
        errors = validate_one(
            self._import_path(klass="BadGenericDocStrings", func=func)
        )["errors"]
        assert isinstance(errors, list)
        assert errors

    @pytest.mark.parametrize(
        "klass,func,msgs",
        [
            # See Also tests
            (
                "BadGenericDocStrings",
                "unknown_section",
                ('Found unknown section "Unknown Section".',),
            ),
            (
                "BadGenericDocStrings",
                "sections_in_wrong_order",
                (
                    "Sections are in the wrong order. Correct order is: Parameters, "
                    "See Also, Examples",
                ),
            ),
            (
                "BadGenericDocStrings",
                "deprecation_in_wrong_order",
                ("Deprecation warning should precede extended summary",),
            ),
            (
                "BadGenericDocStrings",
                "directives_without_two_colons",
                (
                    "reST directives ['versionchanged', 'versionadded', "
                    "'deprecated'] must be followed by two colons",
                ),
            ),
            (
                "BadSeeAlso",
                "no_desc",
                ('Missing description for See Also "Series.tail" reference',),
            ),
            (
                "BadSeeAlso",
                "desc_no_period",
                ('Missing period at end of description for See Also "Series.iloc"',),
            ),
            (
                "BadSeeAlso",
                "desc_first_letter_lowercase",
                ('should be capitalized for See Also "Series.tail"',),
            ),
            # Summary tests
            (
                "BadSummaries",
                "no_summary",
                ("No summary found",),
            ),
            (
                "BadSummaries",
                "heading_whitespaces",
                ("Summary contains heading whitespaces",),
            ),
            (
                "BadSummaries",
                "wrong_line",
                (
                    "should start in the line immediately after the opening quotes",
                    "should be placed in the line after the last text",
                ),
            ),
            ("BadSummaries", "no_punctuation", ("Summary does not end with a period",)),
            (
                "BadSummaries",
                "no_capitalization",
                ("Summary does not start with a capital letter",),
            ),
            (
                "BadSummaries",
                "no_capitalization",
                ("Summary must start with infinitive verb",),
            ),
            ("BadSummaries", "multi_line", ("Summary should fit in a single line",)),
            (
                "BadSummaries",
                "two_paragraph_multi_line",
                ("Summary should fit in a single line",),
            ),
            # Parameters tests
            (
                "BadParameters",
                "no_type",
                ('Parameter "value" has no type',),
            ),
            (
                "BadParameters",
                "type_with_period",
                ('Parameter "value" type should not finish with "."',),
            ),
            (
                "BadParameters",
                "no_description",
                ('Parameter "value" has no description',),
            ),
            (
                "BadParameters",
                "missing_params",
                ("Parameters {'**kwargs'} not documented",),
            ),
            (
                "BadParameters",
                "bad_colon_spacing",
                (
                    'Parameter "kind" requires a space before the colon '
                    "separating the parameter name and type",
                ),
            ),
            (
                "BadParameters",
                "no_description_period",
                ('Parameter "kind" description should finish with "."',),
            ),
            (
                "BadParameters",
                "no_description_period_with_directive",
                ('Parameter "kind" description should finish with "."',),
            ),
            (
                "BadParameters",
                "parameter_capitalization",
                ('Parameter "kind" description should start with a capital letter',),
            ),
            (
                "BadParameters",
                "integer_parameter",
                ('Parameter "kind" type should use "int" instead of "integer"',),
            ),
            (
                "BadParameters",
                "string_parameter",
                ('Parameter "kind" type should use "str" instead of "string"',),
            ),
            (
                "BadParameters",
                "boolean_parameter",
                ('Parameter "kind" type should use "bool" instead of "boolean"',),
            ),
            (
                "BadParameters",
                "list_incorrect_parameter_type",
                ('Parameter "kind" type should use "bool" instead of "boolean"',),
            ),
            (
                "BadParameters",
                "list_incorrect_parameter_type",
                ('Parameter "kind" type should use "int" instead of "integer"',),
            ),
            (
                "BadParameters",
                "list_incorrect_parameter_type",
                ('Parameter "kind" type should use "str" instead of "string"',),
            ),
            (
                "BadParameters",
                "bad_parameter_spacing",
                ("Parameters {'b'} not documented", "Unknown parameters {' b'}"),
            ),
            pytest.param(
                "BadParameters",
                "blank_lines",
                ("No error yet?",),
                marks=pytest.mark.xfail,
            ),
            # Returns tests
            ("BadReturns", "return_not_documented", ("No Returns section found",)),
            ("BadReturns", "yield_not_documented", ("No Yields section found",)),
            pytest.param("BadReturns", "no_type", ("foo",), marks=pytest.mark.xfail),
            ("BadReturns", "no_description", ("Return value has no description",)),
            (
                "BadReturns",
                "no_punctuation",
                ('Return value description should finish with "."',),
            ),
            (
                "BadReturns",
                "named_single_return",
                (
                    "The first line of the Returns section should contain only the "
                    "type, unless multiple values are being returned",
                ),
            ),
            (
                "BadReturns",
                "no_capitalization",
                ("Return value description should start with a capital letter",),
            ),
            (
                "BadReturns",
                "no_period_multi",
                ('Return value description should finish with "."',),
            ),
            (
                "BadGenericDocStrings",
                "method_wo_docstrings",
                ("The object does not have a docstring",),
            ),
            (
                "BadGenericDocStrings",
                "two_linebreaks_between_sections",
                (
                    "Double line break found; please use only one blank line to "
                    "separate sections or paragraphs, and do not leave blank lines "
                    "at the end of docstrings",
                ),
            ),
            (
                "BadGenericDocStrings",
                "linebreak_at_end_of_docstring",
                (
                    "Double line break found; please use only one blank line to "
                    "separate sections or paragraphs, and do not leave blank lines "
                    "at the end of docstrings",
                ),
            ),
        ],
    )
    def test_bad_docstrings(self, capsys, klass, func, msgs):
        with warnings.catch_warnings(record=True) as w:
            result = validate_one(self._import_path(klass=klass, func=func))
        if len(w):
            assert all("Unknown section" in str(ww.message) for ww in w)
        for msg in msgs:
            assert msg in " ".join(err[1] for err in result["errors"])


def decorator(x):
    """Test decorator."""
    return x


@decorator
@decorator
class DecoratorClass:
    """
    Class and methods with decorators.

    * `DecoratorClass` has two decorators.
    * `DecoratorClass.test_no_decorator` has no decorator.
    * `DecoratorClass.test_property` has a `@property` decorator.
    * `DecoratorClass.test_cached_property` has a `@cached_property` decorator.
    * `DecoratorClass.test_three_decorators` has three decorators.

    `Validator.source_file_def_line` should return the `def` or `class` line number, not
    the line of the first decorator.
    """

    def test_no_decorator(self):
        """Test method without decorators."""

    @property
    def test_property(self):
        """Test property method."""

    @cached_property
    def test_cached_property(self):
        """Test property method."""

    @decorator
    @decorator
    @decorator
    def test_three_decorators(self):
        """Test method with three decorators."""


class TestValidatorClass:
    @pytest.mark.parametrize("invalid_name", ["unknown_mod", "unknown_mod.MyClass"])
    def test_raises_for_invalid_module_name(self, invalid_name):
        msg = f'No module can be imported from "{invalid_name}"'
        with pytest.raises(ImportError, match=msg):
            numpydoc.validate.Validator._load_obj(invalid_name)

    @pytest.mark.parametrize(
        "invalid_name", ["datetime.BadClassName", "datetime.bad_method_name"]
    )
    def test_raises_for_invalid_attribute_name(self, invalid_name):
        name_components = invalid_name.split(".")
        obj_name, invalid_attr_name = name_components[-2], name_components[-1]
        msg = f"'{obj_name}' has no attribute '{invalid_attr_name}'"
        with pytest.raises(AttributeError, match=msg):
            numpydoc.validate.Validator._load_obj(invalid_name)

    @pytest.mark.parametrize(
        ["decorated_obj", "def_line"],
        [
            [
                "numpydoc.tests.test_validate.DecoratorClass",
                getsourcelines(DecoratorClass)[-1] + 2,
            ],
            [
                "numpydoc.tests.test_validate.DecoratorClass.test_no_decorator",
                getsourcelines(DecoratorClass.test_no_decorator)[-1],
            ],
            [
                "numpydoc.tests.test_validate.DecoratorClass.test_property",
                getsourcelines(DecoratorClass.test_property.fget)[-1] + 1,
            ],
            [
                "numpydoc.tests.test_validate.DecoratorClass.test_cached_property",
                getsourcelines(DecoratorClass.test_cached_property.func)[-1] + 1,
            ],
            [
                "numpydoc.tests.test_validate.DecoratorClass.test_three_decorators",
                getsourcelines(DecoratorClass.test_three_decorators)[-1] + 3,
            ],
        ],
    )
    def test_source_file_def_line_with_decorators(self, decorated_obj, def_line):
        doc = numpydoc.validate.Validator(
            numpydoc.docscrape.get_doc_object(
                numpydoc.validate.Validator._load_obj(decorated_obj)
            )
        )
        assert doc.source_file_def_line == def_line

    @pytest.mark.parametrize(
        ["property", "file_name"],
        [
            [
                "numpydoc.tests.test_validate.DecoratorClass.test_property",
                getsourcefile(DecoratorClass.test_property.fget),
            ],
            [
                "numpydoc.tests.test_validate.DecoratorClass.test_cached_property",
                getsourcefile(DecoratorClass.test_cached_property.func),
            ],
        ],
    )
    def test_source_file_name_with_properties(self, property, file_name):
        doc = numpydoc.validate.Validator(
            numpydoc.docscrape.get_doc_object(
                numpydoc.validate.Validator._load_obj(property)
            )
        )
        assert doc.source_file_name == file_name


def test_is_generator_validation_with_decorator():
    """Ensure that the check for a Yields section when an object is a generator
    (YD01) works with decorated generators."""

    def tinsel(f):
        @wraps(f)
        def wrapper(*args, **kwargs):
            return f(*args, **kwargs)

        return wrapper

    def foo():
        """A simple generator"""
        yield from range(10)

    @tinsel
    def bar():
        """Generator wrapped once"""
        yield from range(10)

    @tinsel
    @tinsel
    @tinsel
    def baz():
        """Generator wrapped multiple times"""
        yield from range(10)

    # foo without wrapper is a generator
    v = Validator(get_doc_object(foo))
    assert v.is_generator_function

    # Wrapped once
    v = Validator(get_doc_object(bar))
    assert v.is_generator_function

    # Wrapped multiple times
    v = Validator(get_doc_object(baz))
    assert v.is_generator_function
