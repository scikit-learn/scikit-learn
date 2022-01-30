import pytest
import numpydoc.validate
import numpydoc.tests


validate_one = numpydoc.validate.validate


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
        pass

    def plot(self, kind, color="blue", **kwargs):
        """
        Generate a plot.

        Render the data in the Series as a matplotlib plot of the
        specified kind.

        Parameters
        ----------
        kind : str
            Kind of matplotlib plot, e.g.::

                'foo'

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
        pass

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
        pass

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
        pass

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
        pass

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
        pass

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

    def contains(self, pat, case=True, na=float('NaN')):
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
        pass

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
        pass

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
        pass

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
        pass

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
            return None

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
        pass

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
        pass

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
        pass

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


class BadGenericDocStrings:
    """Everything here has a bad docstring
    """

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
        pass

    def astype1(self, dtype):
        """
        Method to cast Series type.

        Does not start with verb.
        """
        pass

    def astype2(self, dtype):
        """
        Cast Series type

        Missing dot at the end.
        """
        pass

    def astype3(self, dtype):
        """
        Cast Series type from its current type to the new type defined in
        the parameter dtype.

        Summary is too verbose and doesn't fit in a single line.
        """
        pass

    def two_linebreaks_between_sections(self, foo):
        """
        Test linebreaks message GL03.

        Note 2 blank lines before parameters section.


        Parameters
        ----------
        foo : str
            Description of foo parameter.
        """
        pass

    def linebreak_at_end_of_docstring(self, foo):
        """
        Test linebreaks message GL03.

        Note extra blank line at end of docstring.

        Parameters
        ----------
        foo : str
            Description of foo parameter.

        """
        pass

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
        pass

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
        >>> print('So far Examples is good, as it goes before Parameters')
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
        pass

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
        pass


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
        pass

    def no_punctuation(self):
        """
        Has the right line but forgets punctuation
        """
        pass

    def no_capitalization(self):
        """
        provides a lowercase summary.
        """
        pass

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
        pass

    def integer_parameter(self, kind):
        """
        Uses integer instead of int.

        Parameters
        ----------
        kind : integer
            Foo bar baz.
        """
        pass

    def string_parameter(self, kind):
        """
        Uses string instead of str.

        Parameters
        ----------
        kind : string
            Foo bar baz.
        """
        pass

    def boolean_parameter(self, kind):
        """
        Uses boolean instead of bool.

        Parameters
        ----------
        kind : boolean
            Foo bar baz.
        """
        pass

    def list_incorrect_parameter_type(self, kind):
        """
        Uses list of boolean instead of list of bool.

        Parameters
        ----------
        kind : list of boolean, integer, float or string
            Foo bar baz.
        """
        pass

    def bad_parameter_spacing(self, a, b):
        """
        The parameters on the same line have an extra space between them.

        Parameters
        ----------
        a,  b : int
            Foo bar baz.
        """
        pass


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
        pass

    def desc_no_period(self):
        """
        Return the first 5 elements of the Series.

        See Also
        --------
        Series.tail : Return the last 5 elements of the Series.
        Series.iloc : Return a slice of the elements in the Series,
            which can also be used to return the first or last n
        """
        pass

    def desc_first_letter_lowercase(self):
        """
        Return the first 5 elements of the Series.

        See Also
        --------
        Series.tail : return the last 5 elements of the Series.
        Series.iloc : Return a slice of the elements in the Series,
            which can also be used to return the first or last n.
        """
        pass

    def prefix_pandas(self):
        """
        Have `pandas` prefix in See Also section.

        See Also
        --------
        pandas.Series.rename : Alter Series index labels or name.
        DataFrame.head : The first `n` rows of the caller object.
        """
        pass


class BadExamples:
    def missing_whitespace_around_arithmetic_operator(self):
        """
        Examples
        --------
        >>> 2+5
        7
        """
        pass

    def indentation_is_not_a_multiple_of_four(self):
        """
        Examples
        --------
        >>> if 2 + 5:
        ...   pass
        """
        pass

    def missing_whitespace_after_comma(self):
        """
        Examples
        --------
        >>> import datetime
        >>> value = datetime.date(2019,1,1)
        """
        pass


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
        result = validate_one(self._import_path(klass="GoodDocStrings", func='one_liner'))
        errors = " ".join(err[1] for err in result["errors"])
        assert 'should start in the line immediately after the opening quotes' not in errors
        assert 'should be placed in the line after the last text' not in errors

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
                self._import_path(klass="WarnGenericFormat", func=func)  # noqa:F821
                )
        assert 'is too short' in w.msg


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
            self._import_path(klass="BadGenericDocStrings", func=func)  # noqa:F821
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
                ("should start in the line immediately after the opening quotes",
                 "should be placed in the line after the last text"),
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
        with pytest.warns(None) as w:
            result = validate_one(self._import_path(klass=klass, func=func))
        if len(w):
            assert all('Unknown section' in str(ww.message) for ww in w)
        for msg in msgs:
            assert msg in " ".join(err[1] for err in result["errors"])


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
