#!/usr/bin/env python
"""
Analyze docstrings to detect errors.

Call ``validate(object_name_to_validate)`` to get a dictionary
with all the detected errors.
"""

import ast
import collections
import functools
import importlib
import inspect
import os
import pydoc
import re
import textwrap
import tokenize
from copy import deepcopy
from typing import Any, Dict, List, Optional, Set

from .docscrape import get_doc_object

DIRECTIVES = ["versionadded", "versionchanged", "deprecated"]
DIRECTIVE_PATTERN = re.compile(
    r"^\s*\.\. ({})(?!::)".format("|".join(DIRECTIVES)), re.IGNORECASE | re.MULTILINE
)
ALLOWED_SECTIONS = [
    "Parameters",
    "Attributes",
    "Methods",
    "Returns",
    "Yields",
    "Other Parameters",
    "Raises",
    "Warns",
    "Warnings",
    "See Also",
    "Notes",
    "References",
    "Examples",
]
# NOTE: The following comment is a sentinel for embedding in the docs - do not
# modify/remove
# start-err-msg
ERROR_MSGS = {
    "GL01": "Docstring text (summary) should start in the line immediately "
    "after the opening quotes (not in the same line, or leaving a "
    "blank line in between)",
    "GL02": "Closing quotes should be placed in the line after the last text "
    "in the docstring (do not close the quotes in the same line as "
    "the text, or leave a blank line between the last text and the "
    "quotes)",
    "GL03": "Double line break found; please use only one blank line to "
    "separate sections or paragraphs, and do not leave blank lines "
    "at the end of docstrings",
    "GL05": 'Tabs found at the start of line "{line_with_tabs}", please use '
    "whitespace only",
    "GL06": 'Found unknown section "{section}". Allowed sections are: '
    "{allowed_sections}",
    "GL07": "Sections are in the wrong order. Correct order is: {correct_sections}",
    "GL08": "The object does not have a docstring",
    "GL09": "Deprecation warning should precede extended summary",
    "GL10": "reST directives {directives} must be followed by two colons",
    "SS01": "No summary found (a short summary in a single line should be "
    "present at the beginning of the docstring)",
    "SS02": "Summary does not start with a capital letter",
    "SS03": "Summary does not end with a period",
    "SS04": "Summary contains heading whitespaces",
    "SS05": "Summary must start with infinitive verb, not third person "
    '(e.g. use "Generate" instead of "Generates")',
    "SS06": "Summary should fit in a single line",
    "ES01": "No extended summary found",
    "PR01": "Parameters {missing_params} not documented",
    "PR02": "Unknown parameters {unknown_params}",
    "PR03": "Wrong parameters order. Actual: {actual_params}. "
    "Documented: {documented_params}",
    "PR04": 'Parameter "{param_name}" has no type',
    "PR05": 'Parameter "{param_name}" type should not finish with "."',
    "PR06": 'Parameter "{param_name}" type should use "{right_type}" instead '
    'of "{wrong_type}"',
    "PR07": 'Parameter "{param_name}" has no description',
    "PR08": 'Parameter "{param_name}" description should start with a '
    "capital letter",
    "PR09": 'Parameter "{param_name}" description should finish with "."',
    "PR10": 'Parameter "{param_name}" requires a space before the colon '
    "separating the parameter name and type",
    "RT01": "No Returns section found",
    "RT02": "The first line of the Returns section should contain only the "
    "type, unless multiple values are being returned",
    "RT03": "Return value has no description",
    "RT04": "Return value description should start with a capital letter",
    "RT05": 'Return value description should finish with "."',
    "YD01": "No Yields section found",
    "SA01": "See Also section not found",
    "SA02": "Missing period at end of description for See Also "
    '"{reference_name}" reference',
    "SA03": "Description should be capitalized for See Also "
    '"{reference_name}" reference',
    "SA04": 'Missing description for See Also "{reference_name}" reference',
    "EX01": "No examples section found",
}
# end-err-msg
# NOTE: The above comment is a sentinel for embedding in the docs - do not
# modify/remove

# Ignore these when evaluating end-of-line-"." checks
IGNORE_STARTS = (" ", "* ", "- ")

# inline comments that can suppress individual checks per line
IGNORE_COMMENT_PATTERN = re.compile("(?:.* numpydoc ignore[=|:] ?)(.+)")


def _unwrap(obj):
    """Iteratively traverse obj.__wrapped__ until first non-wrapped obj found."""
    while hasattr(obj, "__wrapped__"):
        obj = obj.__wrapped__
    return obj


# This function gets called once per function/method to be validated.
# We have to balance memory usage with performance here. It shouldn't be too
# bad to store these `dict`s (they should be rare), but to be safe let's keep
# the limit low-ish. This was set by looking at scipy, numpy, matplotlib,
# and pandas, and they had between ~500 and ~1300 .py files as of 2023-08-16.
@functools.lru_cache(maxsize=2000)
def extract_ignore_validation_comments(
    filepath: Optional[os.PathLike],
    encoding: str = "utf-8",
) -> Dict[int, List[str]]:
    """
    Extract inline comments indicating certain validation checks should be ignored.

    Parameters
    ----------
    filepath : os.PathLike
        Path to the file being inspected.

    Returns
    -------
    dict[int, list[str]]
        Mapping of line number to a list of checks to ignore.
    """
    numpydoc_ignore_comments = {}
    try:
        file = open(filepath, encoding=encoding)
    except (OSError, TypeError):  # can be None, nonexistent, or unreadable
        return numpydoc_ignore_comments
    with file:
        last_declaration = 1
        declarations = ["def", "class"]
        for token in tokenize.generate_tokens(file.readline):
            if token.type == tokenize.NAME and token.string in declarations:
                last_declaration = token.start[0]
            if token.type == tokenize.COMMENT:
                match = re.match(IGNORE_COMMENT_PATTERN, token.string)
                if match:
                    rules = match.group(1).split(",")
                    numpydoc_ignore_comments[last_declaration] = rules
    return numpydoc_ignore_comments


def get_validation_checks(validation_checks: Set[str]) -> Set[str]:
    """
    Get the set of validation checks to report on.

    Parameters
    ----------
    validation_checks : set[str]
        A set of validation checks to report on. If the set is ``{"all"}``,
        all checks will be reported. If the set contains just specific checks,
        only those will be reported on. If the set contains both ``"all"`` and
        specific checks, all checks except those included in the set will be
        reported on.

    Returns
    -------
    set[str]
        The set of validation checks to report on.
    """
    valid_error_codes = set(ERROR_MSGS.keys())
    if "all" in validation_checks:
        block = deepcopy(validation_checks)
        validation_checks = valid_error_codes - block

    # Ensure that the validation check set contains only valid error codes
    invalid_error_codes = validation_checks - valid_error_codes
    if invalid_error_codes:
        raise ValueError(
            f"Unrecognized validation code(s) in numpydoc_validation_checks "
            f"config value: {invalid_error_codes}"
        )

    return validation_checks


def error(code, **kwargs):
    """
    Return a tuple with the error code and the message with variables replaced.

    This is syntactic sugar so instead of:
    - `('PR02', ERROR_MSGS['PR02'].format(doctest_log=log))`

    We can simply use:
    - `error('PR02', doctest_log=log)`

    Parameters
    ----------
    code : str
        Error code.
    **kwargs
        Values for the variables in the error messages

    Returns
    -------
    code : str
        Error code.
    message : str
        Error message with variables replaced.
    """
    return code, ERROR_MSGS[code].format(**kwargs)


class Validator:
    # TODO Can all this class be merged into NumpyDocString?
    def __init__(self, doc_object):
        self.doc = doc_object
        self.obj = self.doc._obj
        self.code_obj = inspect.unwrap(self.obj)
        self.raw_doc = self.obj.__doc__ or ""
        self.clean_doc = pydoc.getdoc(self.obj)

    @property
    def name(self):
        return ".".join([self.obj.__module__, self.obj.__name__])

    @staticmethod
    def _load_obj(name):
        """
        Import Python object from its name as string.

        Parameters
        ----------
        name : str
            Object name to import (e.g. pandas.Series.str.upper)

        Returns
        -------
        object
            Python object that can be a class, method, function...

        Examples
        --------
        >>> Validator._load_obj("datetime.datetime")
        <class 'datetime.datetime'>
        """
        for maxsplit in range(name.count(".") + 1):
            module, *func_parts = name.rsplit(".", maxsplit)
            try:
                obj = importlib.import_module(module)
            except ImportError:
                pass
            else:
                break
        else:
            raise ImportError(f'No module can be imported from "{name}"')

        for part in func_parts:
            obj = getattr(obj, part)
        return obj

    @property
    def type(self):
        return type(self.obj).__name__

    @property
    def is_function_or_method(self):
        return inspect.isfunction(self.obj)

    @property
    def is_generator_function(self):
        return inspect.isgeneratorfunction(_unwrap(self.obj))

    @property
    def source_file_name(self):
        """
        File name where the object is implemented (e.g. pandas/core/frame.py).
        """
        try:
            if isinstance(self.code_obj, property):
                fname = inspect.getsourcefile(self.code_obj.fget)
            elif isinstance(self.code_obj, functools.cached_property):
                fname = inspect.getsourcefile(self.code_obj.func)
            else:
                fname = inspect.getsourcefile(self.code_obj)

        except TypeError:
            # In some cases the object is something complex like a cython
            # object that can't be easily introspected. And it's better to
            # return the source code file of the object as None, than crash
            pass
        else:
            return fname

    # When calling validate, files are parsed twice
    @staticmethod
    @functools.lru_cache(maxsize=4000)
    def _getsourcelines(obj: Any):
        return inspect.getsourcelines(obj)

    @property
    def source_file_def_line(self):
        """
        Number of line where the object is defined in its file.
        """
        try:
            if isinstance(self.code_obj, property):
                sourcelines = self._getsourcelines(self.code_obj.fget)
            elif isinstance(self.code_obj, functools.cached_property):
                sourcelines = self._getsourcelines(self.code_obj.func)
            else:
                sourcelines = self._getsourcelines(self.code_obj)
            # getsourcelines will return the line of the first decorator found for the
            # current function. We have to find the def declaration after that.
            def_line = next(
                i
                for i, x in enumerate(
                    re.match("^ *(def|class) ", s) for s in sourcelines[0]
                )
                if x is not None
            )
            return sourcelines[-1] + def_line
        except (OSError, TypeError):
            # In some cases the object is something complex like a cython
            # object that can't be easily introspected. And it's better to
            # return the line number as None, than crash
            pass

    @property
    def start_blank_lines(self):
        i = None
        if self.raw_doc:
            for i, row in enumerate(self.raw_doc.split("\n")):
                if row.strip():
                    break
        return i

    @property
    def end_blank_lines(self):
        i = None
        if self.raw_doc:
            for i, row in enumerate(reversed(self.raw_doc.split("\n"))):
                if row.strip():
                    break
        return i

    @property
    def double_blank_lines(self):
        prev = True
        for row in self.raw_doc.split("\n"):
            if not prev and not row.strip():
                return True
            prev = row.strip()
        return False

    @property
    def section_titles(self):
        sections = []
        self.doc._doc.reset()
        while not self.doc._doc.eof():
            content = self.doc._read_to_next_section()
            if (
                len(content) > 1
                and len(content[0]) == len(content[1])
                and set(content[1]) == {"-"}
            ):
                sections.append(content[0])
        return sections

    @property
    def summary(self):
        return " ".join(self.doc["Summary"])

    @property
    def num_summary_lines(self):
        return len(self.doc["Summary"])

    @property
    def extended_summary(self):
        if not self.doc["Extended Summary"] and len(self.doc["Summary"]) > 1:
            return " ".join(self.doc["Summary"])
        return " ".join(self.doc["Extended Summary"])

    def _doc_parameters(self, sections):
        parameters = collections.OrderedDict()
        for section in sections:
            for names, type_, desc in self.doc[section]:
                for name in names.split(", "):
                    parameters[name] = (type_, desc)
        return parameters

    @property
    def doc_parameters(self):
        return self._doc_parameters(["Parameters"])

    @property
    def doc_other_parameters(self):
        return self._doc_parameters(["Other Parameters"])

    @property
    def doc_all_parameters(self):
        return self._doc_parameters(["Parameters", "Other Parameters"])

    @property
    def signature_parameters(self):
        def add_stars(param_name, info):
            """
            Add stars to *args and **kwargs parameters
            """
            if info.kind == inspect.Parameter.VAR_POSITIONAL:
                return f"*{param_name}"
            elif info.kind == inspect.Parameter.VAR_KEYWORD:
                return f"**{param_name}"
            else:
                return param_name

        if inspect.isclass(self.obj):
            if hasattr(self.obj, "_accessors") and (
                self.name.split(".")[-1] in self.obj._accessors
            ):
                # accessor classes have a signature but don't want to show this
                return tuple()
        try:
            sig = inspect.signature(self.obj)
        except (TypeError, ValueError):
            # Some objects, mainly in C extensions do not support introspection
            # of the signature
            return tuple()

        params = tuple(
            add_stars(parameter, sig.parameters[parameter])
            for parameter in sig.parameters
        )
        if params and params[0] in ("self", "cls"):
            return params[1:]
        return params

    @property
    def parameter_mismatches(self):
        errs = []
        signature_params = self.signature_parameters
        all_params = tuple(param.replace("\\", "") for param in self.doc_all_parameters)
        missing = set(signature_params) - set(all_params)
        if missing:
            errs.append(error("PR01", missing_params=str(missing)))
        extra = set(all_params) - set(signature_params)
        if extra:
            errs.append(error("PR02", unknown_params=str(extra)))
        if (
            not missing
            and not extra
            and signature_params != all_params
            and not (not signature_params and not all_params)
        ):
            errs.append(
                error(
                    "PR03", actual_params=signature_params, documented_params=all_params
                )
            )

        return errs

    @property
    def directives_without_two_colons(self):
        return DIRECTIVE_PATTERN.findall(self.raw_doc)

    def parameter_type(self, param):
        return self.doc_all_parameters[param][0]

    @property
    def see_also(self):
        result = collections.OrderedDict()
        for funcs, desc in self.doc["See Also"]:
            for func, _ in funcs:
                result[func] = "".join(desc)

        return result

    @property
    def examples(self):
        return self.doc["Examples"]

    @property
    def returns(self):
        return self.doc["Returns"]

    @property
    def yields(self):
        return self.doc["Yields"]

    @property
    def method_source(self):
        try:
            source = inspect.getsource(self.obj)
        except TypeError:
            return ""
        return textwrap.dedent(source)

    @property
    def method_returns_something(self):
        """
        Check if the docstrings method can return something.

        Bare returns, returns valued None and returns from nested functions are
        disconsidered.

        Returns
        -------
        bool
            Whether the docstrings method can return something.
        """

        def get_returns_not_on_nested_functions(node):
            returns = [node] if isinstance(node, ast.Return) else []
            for child in ast.iter_child_nodes(node):
                # Ignore nested functions and its subtrees.
                if not isinstance(child, ast.FunctionDef):
                    child_returns = get_returns_not_on_nested_functions(child)
                    returns.extend(child_returns)
            return returns

        tree = ast.parse(self.method_source).body
        if tree:
            returns = get_returns_not_on_nested_functions(tree[0])
            return_values = [r.value for r in returns]
            # Replace Constant nodes valued None for None.
            for i, v in enumerate(return_values):
                if isinstance(v, ast.Constant) and v.value is None:
                    return_values[i] = None
            return any(return_values)
        else:
            return False

    @property
    def deprecated(self):
        return ".. deprecated:: " in (self.summary + self.extended_summary)


def _check_desc(desc, code_no_desc, code_no_upper, code_no_period, **kwargs):
    # Find and strip out any sphinx directives
    desc = "\n".join(desc)
    for directive in DIRECTIVES:
        full_directive = f".. {directive}"
        if full_directive in desc:
            # Only retain any description before the directive
            desc = desc[: desc.index(full_directive)].rstrip("\n")
    desc = desc.split("\n")

    errs = list()
    if not "".join(desc):
        errs.append(error(code_no_desc, **kwargs))
    else:
        if desc[0][0].isalpha() and not desc[0][0].isupper():
            errs.append(error(code_no_upper, **kwargs))
        # Not ending in "." is only an error if the last bit is not
        # indented (e.g., quote or code block)
        if not desc[-1].endswith(".") and not desc[-1].startswith(IGNORE_STARTS):
            errs.append(error(code_no_period, **kwargs))
    return errs


def validate(obj_name, validator_cls=None, **validator_kwargs):
    """
    Validate the docstring.

    Parameters
    ----------
    obj_name : str
        The name of the object whose docstring will be evaluated, e.g.
        'pandas.read_csv'. The string must include the full, unabbreviated
        package/module names, i.e. 'pandas.read_csv', not 'pd.read_csv' or
        'read_csv'.
    validator_cls : Validator, optional
        The Validator class to use. By default, :class:`Validator` will be used.
    **validator_kwargs
        Keyword arguments to pass to ``validator_cls`` upon initialization.
        Note that ``obj_name`` will be passed as a named argument as well.

    Returns
    -------
    dict
        A dictionary containing all the information obtained from validating
        the docstring.

    Notes
    -----
    The errors codes are defined as:
    - First two characters: Section where the error happens:
       * GL: Global (no section, like section ordering errors)
       * SS: Short summary
       * ES: Extended summary
       * PR: Parameters
       * RT: Returns
       * YD: Yields
       * RS: Raises
       * WN: Warns
       * SA: See Also
       * NT: Notes
       * RF: References
       * EX: Examples
    - Last two characters: Numeric error code inside the section

    For example, PR02 is the second codified error in the Parameters section
    (which in this case is assigned to the error when unknown parameters are documented).

    The error codes, their corresponding error messages, and the details on how
    they are validated, are not documented more than in the source code of this
    function.
    """
    if not validator_cls:
        if isinstance(obj_name, str):
            doc = Validator(get_doc_object(Validator._load_obj(obj_name)))
        else:
            doc = Validator(obj_name)
    else:
        doc = validator_cls(obj_name=obj_name, **validator_kwargs)

    # lineno is only 0 if we have a module docstring in the file, and we are
    # validating that, so we change to 1 for readability of the output
    ignore_validation_comments = extract_ignore_validation_comments(
        doc.source_file_name
    ).get(doc.source_file_def_line or 1, [])

    errs = []
    if not doc.raw_doc:
        if "GL08" not in ignore_validation_comments:
            errs.append(error("GL08"))
        return {
            "type": doc.type,
            "docstring": doc.clean_doc,
            "deprecated": doc.deprecated,
            "file": doc.source_file_name,
            "file_line": doc.source_file_def_line,
            "errors": errs,
            "examples_errors": "",
        }

    if doc.start_blank_lines != 1 and "\n" in doc.raw_doc:
        errs.append(error("GL01"))
    if doc.end_blank_lines != 1 and "\n" in doc.raw_doc:
        errs.append(error("GL02"))
    if doc.double_blank_lines:
        errs.append(error("GL03"))
    for line in doc.raw_doc.splitlines():
        if re.match("^ *\t", line):
            errs.append(error("GL05", line_with_tabs=line.lstrip()))

    unexpected_sections = [
        section for section in doc.section_titles if section not in ALLOWED_SECTIONS
    ]
    for section in unexpected_sections:
        errs.append(
            error("GL06", section=section, allowed_sections=", ".join(ALLOWED_SECTIONS))
        )

    correct_order = [
        section for section in ALLOWED_SECTIONS if section in doc.section_titles
    ]
    if correct_order != doc.section_titles:
        errs.append(error("GL07", correct_sections=", ".join(correct_order)))

    if doc.deprecated and not doc.extended_summary.startswith(".. deprecated:: "):
        errs.append(error("GL09"))

    directives_without_two_colons = doc.directives_without_two_colons
    if directives_without_two_colons:
        errs.append(error("GL10", directives=directives_without_two_colons))

    if not doc.summary:
        errs.append(error("SS01"))
    else:
        if doc.summary[0].isalpha() and not doc.summary[0].isupper():
            errs.append(error("SS02"))
        if doc.summary[-1] != ".":
            errs.append(error("SS03"))
        if doc.summary != doc.summary.lstrip():
            errs.append(error("SS04"))
        elif doc.is_function_or_method and doc.summary.split(" ")[0][-1] == "s":
            errs.append(error("SS05"))
        if doc.num_summary_lines > 1:
            errs.append(error("SS06"))

    if not doc.extended_summary:
        errs.append(("ES01", "No extended summary found"))

    # PR01: Parameters not documented
    # PR02: Unknown parameters
    # PR03: Wrong parameters order
    errs += doc.parameter_mismatches

    for param, kind_desc in doc.doc_all_parameters.items():
        if not param.startswith("*"):  # Check can ignore var / kwargs
            if not doc.parameter_type(param):
                if ":" in param:
                    errs.append(error("PR10", param_name=param.split(":")[0]))
                else:
                    errs.append(error("PR04", param_name=param))
            else:
                if doc.parameter_type(param)[-1] == ".":
                    errs.append(error("PR05", param_name=param))
                # skip common_type_error checks when the param type is a set of
                # options
                if "{" in doc.parameter_type(param):
                    continue
                common_type_errors = [
                    ("integer", "int"),
                    ("boolean", "bool"),
                    ("string", "str"),
                ]
                for wrong_type, right_type in common_type_errors:
                    if wrong_type in set(re.split(r"\W", doc.parameter_type(param))):
                        errs.append(
                            error(
                                "PR06",
                                param_name=param,
                                right_type=right_type,
                                wrong_type=wrong_type,
                            )
                        )
        errs.extend(_check_desc(kind_desc[1], "PR07", "PR08", "PR09", param_name=param))

    if doc.is_function_or_method:
        if not doc.returns:
            if doc.method_returns_something:
                errs.append(error("RT01"))
        else:
            if len(doc.returns) == 1 and doc.returns[0].name:
                errs.append(error("RT02"))
            for name_or_type, type_, desc in doc.returns:
                errs.extend(_check_desc(desc, "RT03", "RT04", "RT05"))

        if not doc.yields and doc.is_generator_function:
            errs.append(error("YD01"))

    if not doc.see_also:
        errs.append(error("SA01"))
    else:
        for rel_name, rel_desc in doc.see_also.items():
            if rel_desc:
                if not rel_desc.endswith("."):
                    errs.append(error("SA02", reference_name=rel_name))
                if rel_desc[0].isalpha() and not rel_desc[0].isupper():
                    errs.append(error("SA03", reference_name=rel_name))
            else:
                errs.append(error("SA04", reference_name=rel_name))

    if not doc.examples:
        errs.append(error("EX01"))

    errs = [err for err in errs if err[0] not in ignore_validation_comments]

    return {
        "type": doc.type,
        "docstring": doc.clean_doc,
        "deprecated": doc.deprecated,
        "file": doc.source_file_name,
        "file_line": doc.source_file_def_line,
        "errors": errs,
    }
