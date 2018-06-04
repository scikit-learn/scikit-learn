"""Parsed source code checkers for docstring violations."""

import ast
import string
import sys
import tokenize as tk
from itertools import takewhile
from re import compile as re
from collections import namedtuple

from . import violations
from .config import IllegalConfiguration
from .parser import (Package, Module, Class, NestedClass, Definition, AllError,
                     Method, Function, NestedFunction, Parser, StringIO,
                     ParseError)
from .utils import log, is_blank, pairwise
from .wordlists import IMPERATIVE_VERBS, IMPERATIVE_BLACKLIST, stem


__all__ = ('check', )


# If possible (python >= 3.2) use tokenize.open to open files, so PEP 263
# encoding markers are interpreted.
try:
    tokenize_open = tk.open
except AttributeError:
    tokenize_open = open


def check_for(kind, terminal=False):
    def decorator(f):
        f._check_for = kind
        f._terminal = terminal
        return f
    return decorator


class ConventionChecker(object):
    """Checker for PEP 257 and numpy conventions.

    D10x: Missing docstrings
    D20x: Whitespace issues
    D30x: Docstring formatting
    D40x: Docstring content issues

    """

    SECTION_NAMES = ['Short Summary',
                     'Extended Summary',
                     'Parameters',
                     'Returns',
                     'Yields',
                     'Other Parameters',
                     'Raises',
                     'See Also',
                     'Notes',
                     'References',
                     'Examples',
                     'Attributes',
                     'Methods']

    def check_source(self, source, filename, ignore_decorators=None):
        module = parse(StringIO(source), filename)
        for definition in module:
            for this_check in self.checks:
                terminate = False
                if isinstance(definition, this_check._check_for):
                    skipping_all = (definition.skipped_error_codes == 'all')
                    decorator_skip = ignore_decorators is not None and any(
                        len(ignore_decorators.findall(dec.name)) > 0
                        for dec in definition.decorators)
                    if not skipping_all and not decorator_skip:
                        error = this_check(self, definition,
                                           definition.docstring)
                    else:
                        error = None
                    errors = error if hasattr(error, '__iter__') else [error]
                    for error in errors:
                        if error is not None and error.code not in \
                                definition.skipped_error_codes:
                            partition = this_check.__doc__.partition('.\n')
                            message, _, explanation = partition
                            error.set_context(explanation=explanation,
                                              definition=definition)
                            yield error
                            if this_check._terminal:
                                terminate = True
                                break
                if terminate:
                    break

    @property
    def checks(self):
        all = [this_check for this_check in vars(type(self)).values()
               if hasattr(this_check, '_check_for')]
        return sorted(all, key=lambda this_check: not this_check._terminal)

    @check_for(Definition, terminal=True)
    def check_docstring_missing(self, definition, docstring):
        """D10{0,1,2,3}: Public definitions should have docstrings.

        All modules should normally have docstrings.  [...] all functions and
        classes exported by a module should also have docstrings. Public
        methods (including the __init__ constructor) should also have
        docstrings.

        Note: Public (exported) definitions are either those with names listed
              in __all__ variable (if present), or those that do not start
              with a single underscore.

        """
        if (not docstring and definition.is_public or
                docstring and is_blank(ast.literal_eval(docstring))):
            codes = {Module: violations.D100,
                     Class: violations.D101,
                     NestedClass: violations.D106,
                     Method: (lambda: violations.D105() if definition.is_magic
                              else (violations.D107() if definition.is_init
                              else violations.D102())),
                     Function: violations.D103,
                     NestedFunction: violations.D103,
                     Package: violations.D104}
            return codes[type(definition)]()

    @check_for(Definition)
    def check_one_liners(self, definition, docstring):
        """D200: One-liner docstrings should fit on one line with quotes.

        The closing quotes are on the same line as the opening quotes.
        This looks better for one-liners.

        """
        if docstring:
            lines = ast.literal_eval(docstring).split('\n')
            if len(lines) > 1:
                non_empty_lines = sum(1 for l in lines if not is_blank(l))
                if non_empty_lines == 1:
                    return violations.D200(len(lines))

    @check_for(Function)
    def check_no_blank_before(self, function, docstring):  # def
        """D20{1,2}: No blank lines allowed around function/method docstring.

        There's no blank line either before or after the docstring.

        """
        if docstring:
            before, _, after = function.source.partition(docstring)
            blanks_before = list(map(is_blank, before.split('\n')[:-1]))
            blanks_after = list(map(is_blank, after.split('\n')[1:]))
            blanks_before_count = sum(takewhile(bool, reversed(blanks_before)))
            blanks_after_count = sum(takewhile(bool, blanks_after))
            if blanks_before_count != 0:
                yield violations.D201(blanks_before_count)
            if not all(blanks_after) and blanks_after_count != 0:
                yield violations.D202(blanks_after_count)

    @check_for(Class)
    def check_blank_before_after_class(self, class_, docstring):
        """D20{3,4}: Class docstring should have 1 blank line around them.

        Insert a blank line before and after all docstrings (one-line or
        multi-line) that document a class -- generally speaking, the class's
        methods are separated from each other by a single blank line, and the
        docstring needs to be offset from the first method by a blank line;
        for symmetry, put a blank line between the class header and the
        docstring.

        """
        # NOTE: this gives false-positive in this case
        # class Foo:
        #
        #     """Docstring."""
        #
        #
        # # comment here
        # def foo(): pass
        if docstring:
            before, _, after = class_.source.partition(docstring)
            blanks_before = list(map(is_blank, before.split('\n')[:-1]))
            blanks_after = list(map(is_blank, after.split('\n')[1:]))
            blanks_before_count = sum(takewhile(bool, reversed(blanks_before)))
            blanks_after_count = sum(takewhile(bool, blanks_after))
            if blanks_before_count != 0:
                yield violations.D211(blanks_before_count)
            if blanks_before_count != 1:
                yield violations.D203(blanks_before_count)
            if not all(blanks_after) and blanks_after_count != 1:
                yield violations.D204(blanks_after_count)

    @check_for(Definition)
    def check_blank_after_summary(self, definition, docstring):
        """D205: Put one blank line between summary line and description.

        Multi-line docstrings consist of a summary line just like a one-line
        docstring, followed by a blank line, followed by a more elaborate
        description. The summary line may be used by automatic indexing tools;
        it is important that it fits on one line and is separated from the
        rest of the docstring by a blank line.

        """
        if docstring:
            lines = ast.literal_eval(docstring).strip().split('\n')
            if len(lines) > 1:
                post_summary_blanks = list(map(is_blank, lines[1:]))
                blanks_count = sum(takewhile(bool, post_summary_blanks))
                if blanks_count != 1:
                    return violations.D205(blanks_count)

    @staticmethod
    def _get_docstring_indent(definition, docstring):
        """Return the indentation of the docstring's opening quotes."""
        before_docstring, _, _ = definition.source.partition(docstring)
        _, _, indent = before_docstring.rpartition('\n')
        return indent

    @check_for(Definition)
    def check_indent(self, definition, docstring):
        """D20{6,7,8}: The entire docstring should be indented same as code.

        The entire docstring is indented the same as the quotes at its
        first line.

        """
        if docstring:
            indent = self._get_docstring_indent(definition, docstring)
            lines = docstring.split('\n')
            if len(lines) > 1:
                lines = lines[1:]  # First line does not need indent.
                indents = [leading_space(l) for l in lines if not is_blank(l)]
                if set(' \t') == set(''.join(indents) + indent):
                    yield violations.D206()
                if (len(indents) > 1 and min(indents[:-1]) > indent or
                        indents[-1] > indent):
                    yield violations.D208()
                if min(indents) < indent:
                    yield violations.D207()

    @check_for(Definition)
    def check_newline_after_last_paragraph(self, definition, docstring):
        """D209: Put multi-line docstring closing quotes on separate line.

        Unless the entire docstring fits on a line, place the closing
        quotes on a line by themselves.

        """
        if docstring:
            lines = [l for l in ast.literal_eval(docstring).split('\n')
                     if not is_blank(l)]
            if len(lines) > 1:
                if docstring.split("\n")[-1].strip() not in ['"""', "'''"]:
                    return violations.D209()

    @check_for(Definition)
    def check_surrounding_whitespaces(self, definition, docstring):
        """D210: No whitespaces allowed surrounding docstring text."""
        if docstring:
            lines = ast.literal_eval(docstring).split('\n')
            if lines[0].startswith(' ') or \
                    len(lines) == 1 and lines[0].endswith(' '):
                return violations.D210()

    @check_for(Definition)
    def check_multi_line_summary_start(self, definition, docstring):
        """D21{2,3}: Multi-line docstring summary style check.

        A multi-line docstring summary should start either at the first,
        or separately at the second line of a docstring.

        """
        if docstring:
            start_triple = [
                '"""', "'''",
                'u"""', "u'''",
                'r"""', "r'''",
                'ur"""', "ur'''"
            ]

            lines = ast.literal_eval(docstring).split('\n')
            if len(lines) > 1:
                first = docstring.split("\n")[0].strip().lower()
                if first in start_triple:
                    return violations.D212()
                else:
                    return violations.D213()

    @check_for(Definition)
    def check_triple_double_quotes(self, definition, docstring):
        r'''D300: Use """triple double quotes""".

        For consistency, always use """triple double quotes""" around
        docstrings. Use r"""raw triple double quotes""" if you use any
        backslashes in your docstrings. For Unicode docstrings, use
        u"""Unicode triple-quoted strings""".

        Note: Exception to this is made if the docstring contains
              """ quotes in its body.

        '''
        if docstring:
            if '"""' in ast.literal_eval(docstring):
                # Allow ''' quotes if docstring contains """, because
                # otherwise """ quotes could not be expressed inside
                # docstring. Not in PEP 257.
                regex = re(r"[uU]?[rR]?'''[^'].*")
            else:
                regex = re(r'[uU]?[rR]?"""[^"].*')

            if not regex.match(docstring):
                illegal_matcher = re(r"""[uU]?[rR]?("+|'+).*""")
                illegal_quotes = illegal_matcher.match(docstring).group(1)
                return violations.D300(illegal_quotes)

    @check_for(Definition)
    def check_backslashes(self, definition, docstring):
        r'''D301: Use r""" if any backslashes in a docstring.

        Use r"""raw triple double quotes""" if you use any backslashes
        (\) in your docstrings.

        '''
        # Just check that docstring is raw, check_triple_double_quotes
        # ensures the correct quotes.
        if docstring and '\\' in docstring and not docstring.startswith(
                ('r', 'ur')):
            return violations.D301()

    @check_for(Definition)
    def check_unicode_docstring(self, definition, docstring):
        r'''D302: Use u""" for docstrings with Unicode.

        For Unicode docstrings, use u"""Unicode triple-quoted strings""".

        '''
        if 'unicode_literals' in definition.module.future_imports:
            return

        # Just check that docstring is unicode, check_triple_double_quotes
        # ensures the correct quotes.
        if docstring and sys.version_info[0] <= 2:
            if not is_ascii(docstring) and not docstring.startswith(
                    ('u', 'ur')):
                return violations.D302()

    @check_for(Definition)
    def check_ends_with_period(self, definition, docstring):
        """D400: First line should end with a period.

        The [first line of a] docstring is a phrase ending in a period.

        """
        if docstring:
            summary_line = ast.literal_eval(docstring).strip().split('\n')[0]
            if not summary_line.endswith('.'):
                return violations.D400(summary_line[-1])

    @check_for(Function)
    def check_imperative_mood(self, function, docstring):  # def context
        """D401: First line should be in imperative mood: 'Do', not 'Does'.

        [Docstring] prescribes the function or method's effect as a command:
        ("Do this", "Return that"), not as a description; e.g. don't write
        "Returns the pathname ...".

        """
        if docstring and not function.is_test:
            stripped = ast.literal_eval(docstring).strip()
            if stripped:
                first_word = stripped.split()[0]
                check_word = first_word.lower()

                if check_word in IMPERATIVE_BLACKLIST:
                    return violations.D401b(first_word)

                try:
                    correct_form = IMPERATIVE_VERBS.get(stem(check_word))
                except UnicodeDecodeError:
                    # This is raised when the docstring contains unicode
                    # characters in the first word, but is not a unicode
                    # string. In which case D302 will be reported. Ignoring.
                    return

                if correct_form and correct_form != check_word:
                    return violations.D401(
                        correct_form.capitalize(),
                        first_word
                    )

    @check_for(Function)
    def check_no_signature(self, function, docstring):  # def context
        """D402: First line should not be function's or method's "signature".

        The one-line docstring should NOT be a "signature" reiterating the
        function/method parameters (which can be obtained by introspection).

        """
        if docstring:
            first_line = ast.literal_eval(docstring).strip().split('\n')[0]
            if function.name + '(' in first_line.replace(' ', ''):
                return violations.D402()

    @check_for(Function)
    def check_capitalized(self, function, docstring):
        """D403: First word of the first line should be properly capitalized.

        The [first line of a] docstring is a phrase ending in a period.

        """
        if docstring:
            first_word = ast.literal_eval(docstring).split()[0]
            if first_word == first_word.upper():
                return
            for char in first_word:
                if char not in string.ascii_letters and char != "'":
                    return
            if first_word != first_word.capitalize():
                return violations.D403(first_word.capitalize(), first_word)

    @check_for(Definition)
    def check_starts_with_this(self, function, docstring):
        """D404: First word of the docstring should not be `This`.

        Docstrings should use short, simple language. They should not begin
        with "This class is [..]" or "This module contains [..]".

        """
        if docstring:
            first_word = ast.literal_eval(docstring).split()[0]
            if first_word.lower() == 'this':
                return violations.D404()

    @staticmethod
    def _get_leading_words(line):
        """Return any leading set of words from `line`.

        For example, if `line` is "  Hello world!!!", returns "Hello world".
        """
        result = re("[A-Za-z ]+").match(line.strip())
        if result is not None:
            return result.group()

    @staticmethod
    def _is_a_docstring_section(context):
        """Check if the suspected context is really a section header.

        Lets have a look at the following example docstring:
            '''Title.

            Some part of the docstring that specifies what the function
            returns. <----- Not a real section name. It has a suffix and the
                            previous line is not empty and does not end with
                            a punctuation sign.

            This is another line in the docstring. It describes stuff,
            but we forgot to add a blank line between it and the section name.
            Returns  <----- A real section name. The previous line ends with
            -------         a period, therefore it is in a new
                            grammatical context.
            Bla.

            '''

        To make sure this is really a section we check these conditions:
            * There's no suffix to the section name.
            * The previous line ends with punctuation.
            * The previous line is empty.

        If one of the conditions is true, we will consider the line as
        a section name.
        """
        section_name_suffix = context.line.lstrip(context.section_name).strip()

        punctuation = [',', ';', '.', '-', '\\', '/', ']', '}', ')']
        prev_line_ends_with_punctuation = \
            any(context.previous_line.strip().endswith(x) for x in punctuation)

        return (is_blank(section_name_suffix) or
                prev_line_ends_with_punctuation or
                is_blank(context.previous_line))

    @classmethod
    def _check_section_underline(cls, section_name, context, indentation):
        """D4{07,08,09,12}, D215: Section underline checks.

        Check for correct formatting for docstring sections. Checks that:
            * The line that follows the section name contains
              dashes (D40{7,8}).
            * The amount of dashes is equal to the length of the section
              name (D409).
            * The section's content does not begin in the line that follows
              the section header (D412).
            * The indentation of the dashed line is equal to the docstring's
              indentation (D215).
        """
        blank_lines_after_header = 0

        for line in context.following_lines:
            if not is_blank(line):
                break
            blank_lines_after_header += 1
        else:
            # There are only blank lines after the header.
            yield violations.D407(section_name)
            return

        non_empty_line = context.following_lines[blank_lines_after_header]
        dash_line_found = ''.join(set(non_empty_line.strip())) == '-'

        if not dash_line_found:
            yield violations.D407(section_name)
            if blank_lines_after_header > 0:
                yield violations.D412(section_name)
        else:
            if blank_lines_after_header > 0:
                yield violations.D408(section_name)

            if non_empty_line.strip() != "-" * len(section_name):
                yield violations.D409(len(section_name),
                                      section_name,
                                      len(non_empty_line.strip()))

            if leading_space(non_empty_line) > indentation:
                yield violations.D215(section_name)

            line_after_dashes_index = blank_lines_after_header + 1
            # If the line index after the dashes is in range (perhaps we have
            # a header + underline followed by another section header).
            if line_after_dashes_index < len(context.following_lines):
                line_after_dashes = \
                    context.following_lines[line_after_dashes_index]
                if is_blank(line_after_dashes):
                    rest_of_lines = \
                        context.following_lines[line_after_dashes_index:]
                    if not is_blank(''.join(rest_of_lines)):
                        yield violations.D412(section_name)
                    else:
                        yield violations.D414(section_name)
            else:
                yield violations.D414(section_name)

    @classmethod
    def _check_section(cls, docstring, definition, context):
        """D4{05,06,10,11,13}, D214: Section name checks.

        Check for valid section names. Checks that:
            * The section name is properly capitalized (D405).
            * The section is not over-indented (D214).
            * The section name has no superfluous suffix to it (D406).
            * There's a blank line after the section (D410, D413).
            * There's a blank line before the section (D411).

        Also yields all the errors from `_check_section_underline`.
        """
        capitalized_section = context.section_name.title()
        indentation = cls._get_docstring_indent(definition, docstring)

        if (context.section_name not in cls.SECTION_NAMES and
                capitalized_section in cls.SECTION_NAMES):
            yield violations.D405(capitalized_section, context.section_name)

        if leading_space(context.line) > indentation:
            yield violations.D214(capitalized_section)

        suffix = context.line.strip().lstrip(context.section_name)
        if suffix:
            yield violations.D406(capitalized_section, context.line.strip())

        if (not context.following_lines or
                not is_blank(context.following_lines[-1])):
            if context.is_last_section:
                yield violations.D413(capitalized_section)
            else:
                yield violations.D410(capitalized_section)

        if not is_blank(context.previous_line):
            yield violations.D411(capitalized_section)

        for err in cls._check_section_underline(capitalized_section,
                                                context,
                                                indentation):
            yield err

    @check_for(Definition)
    def check_docstring_sections(self, definition, docstring):
        """D21{4,5}, D4{05,06,07,08,09,10}: Docstring sections checks.

        Check the general format of a sectioned docstring:
            '''This is my one-liner.

            Short Summary
            -------------
            This is my summary.

            Returns
            -------
            None.

            '''

        Section names appear in `SECTION_NAMES`.
        """
        if not docstring:
            return

        lines = docstring.split("\n")
        if len(lines) < 2:
            return

        lower_section_names = [s.lower() for s in self.SECTION_NAMES]

        def _suspected_as_section(_line):
            result = self._get_leading_words(_line.lower())
            return result in lower_section_names

        # Finding our suspects.
        suspected_section_indices = [i for i, line in enumerate(lines) if
                                     _suspected_as_section(line)]

        SectionContext = namedtuple('SectionContext', ('section_name',
                                                       'previous_line',
                                                       'line',
                                                       'following_lines',
                                                       'original_index',
                                                       'is_last_section'))

        # First - create a list of possible contexts. Note that the
        # `following_linex` member is until the end of the docstring.
        contexts = (SectionContext(self._get_leading_words(lines[i].strip()),
                                   lines[i - 1],
                                   lines[i],
                                   lines[i + 1:],
                                   i,
                                   False)
                    for i in suspected_section_indices)

        # Now that we have manageable objects - rule out false positives.
        contexts = (c for c in contexts if self._is_a_docstring_section(c))

        # Now we shall trim the `following lines` field to only reach the
        # next section name.
        for a, b in pairwise(contexts, None):
            end = -1 if b is None else b.original_index
            new_ctx = SectionContext(a.section_name,
                                     a.previous_line,
                                     a.line,
                                     lines[a.original_index + 1:end],
                                     a.original_index,
                                     b is None)
            for err in self._check_section(docstring, definition, new_ctx):
                yield err


parse = Parser()


def check(filenames, select=None, ignore=None, ignore_decorators=None):
    """Generate docstring errors that exist in `filenames` iterable.

    By default, the PEP-257 convention is checked. To specifically define the
    set of error codes to check for, supply either `select` or `ignore` (but
    not both). In either case, the parameter should be a collection of error
    code strings, e.g., {'D100', 'D404'}.

    When supplying `select`, only specified error codes will be reported.
    When supplying `ignore`, all error codes which were not specified will be
    reported.

    Note that ignored error code refer to the entire set of possible
    error codes, which is larger than just the PEP-257 convention. To your
    convenience, you may use `pydocstyle.violations.conventions.pep257` as
    a base set to add or remove errors from.

    Examples
    ---------
    >>> check(['pydocstyle.py'])
    <generator object check at 0x...>

    >>> check(['pydocstyle.py'], select=['D100'])
    <generator object check at 0x...>

    >>> check(['pydocstyle.py'], ignore=conventions.pep257 - {'D100'})
    <generator object check at 0x...>

    """
    if select is not None and ignore is not None:
        raise IllegalConfiguration('Cannot pass both select and ignore. '
                                   'They are mutually exclusive.')
    elif select is not None:
        checked_codes = select
    elif ignore is not None:
        checked_codes = list(set(violations.ErrorRegistry.get_error_codes()) -
                             set(ignore))
    else:
        checked_codes = violations.conventions.pep257

    for filename in filenames:
        log.info('Checking file %s.', filename)
        try:
            with tokenize_open(filename) as file:
                source = file.read()
            for error in ConventionChecker().check_source(source, filename,
                                                          ignore_decorators):
                code = getattr(error, 'code', None)
                if code in checked_codes:
                    yield error
        except (EnvironmentError, AllError, ParseError) as error:
            log.warning('Error in file %s: %s', filename, error)
            yield error
        except tk.TokenError:
            yield SyntaxError('invalid syntax in file %s' % filename)


def is_ascii(string):
    return all(ord(char) < 128 for char in string)


def leading_space(string):
    return re('\s*').match(string).group()
