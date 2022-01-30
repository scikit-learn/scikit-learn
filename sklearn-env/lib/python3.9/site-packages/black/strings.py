"""
Simple formatting on strings. Further string formatting code is in trans.py.
"""

import regex as re
import sys
from typing import List, Pattern


STRING_PREFIX_CHARS = "furbFURB"  # All possible string prefix characters.


def sub_twice(regex: Pattern[str], replacement: str, original: str) -> str:
    """Replace `regex` with `replacement` twice on `original`.

    This is used by string normalization to perform replaces on
    overlapping matches.
    """
    return regex.sub(replacement, regex.sub(replacement, original))


def has_triple_quotes(string: str) -> bool:
    """
    Returns:
        True iff @string starts with three quotation characters.
    """
    raw_string = string.lstrip(STRING_PREFIX_CHARS)
    return raw_string[:3] in {'"""', "'''"}


def lines_with_leading_tabs_expanded(s: str) -> List[str]:
    """
    Splits string into lines and expands only leading tabs (following the normal
    Python rules)
    """
    lines = []
    for line in s.splitlines():
        # Find the index of the first non-whitespace character after a string of
        # whitespace that includes at least one tab
        match = re.match(r"\s*\t+\s*(\S)", line)
        if match:
            first_non_whitespace_idx = match.start(1)

            lines.append(
                line[:first_non_whitespace_idx].expandtabs()
                + line[first_non_whitespace_idx:]
            )
        else:
            lines.append(line)
    return lines


def fix_docstring(docstring: str, prefix: str) -> str:
    # https://www.python.org/dev/peps/pep-0257/#handling-docstring-indentation
    if not docstring:
        return ""
    lines = lines_with_leading_tabs_expanded(docstring)
    # Determine minimum indentation (first line doesn't count):
    indent = sys.maxsize
    for line in lines[1:]:
        stripped = line.lstrip()
        if stripped:
            indent = min(indent, len(line) - len(stripped))
    # Remove indentation (first line is special):
    trimmed = [lines[0].strip()]
    if indent < sys.maxsize:
        last_line_idx = len(lines) - 2
        for i, line in enumerate(lines[1:]):
            stripped_line = line[indent:].rstrip()
            if stripped_line or i == last_line_idx:
                trimmed.append(prefix + stripped_line)
            else:
                trimmed.append("")
    return "\n".join(trimmed)


def get_string_prefix(string: str) -> str:
    """
    Pre-conditions:
        * assert_is_leaf_string(@string)

    Returns:
        @string's prefix (e.g. '', 'r', 'f', or 'rf').
    """
    assert_is_leaf_string(string)

    prefix = ""
    prefix_idx = 0
    while string[prefix_idx] in STRING_PREFIX_CHARS:
        prefix += string[prefix_idx]
        prefix_idx += 1

    return prefix


def assert_is_leaf_string(string: str) -> None:
    """
    Checks the pre-condition that @string has the format that you would expect
    of `leaf.value` where `leaf` is some Leaf such that `leaf.type ==
    token.STRING`. A more precise description of the pre-conditions that are
    checked are listed below.

    Pre-conditions:
        * @string starts with either ', ", <prefix>', or <prefix>" where
        `set(<prefix>)` is some subset of `set(STRING_PREFIX_CHARS)`.
        * @string ends with a quote character (' or ").

    Raises:
        AssertionError(...) if the pre-conditions listed above are not
        satisfied.
    """
    dquote_idx = string.find('"')
    squote_idx = string.find("'")
    if -1 in [dquote_idx, squote_idx]:
        quote_idx = max(dquote_idx, squote_idx)
    else:
        quote_idx = min(squote_idx, dquote_idx)

    assert (
        0 <= quote_idx < len(string) - 1
    ), f"{string!r} is missing a starting quote character (' or \")."
    assert string[-1] in (
        "'",
        '"',
    ), f"{string!r} is missing an ending quote character (' or \")."
    assert set(string[:quote_idx]).issubset(
        set(STRING_PREFIX_CHARS)
    ), f"{set(string[:quote_idx])} is NOT a subset of {set(STRING_PREFIX_CHARS)}."


def normalize_string_prefix(s: str, remove_u_prefix: bool = False) -> str:
    """Make all string prefixes lowercase.

    If remove_u_prefix is given, also removes any u prefix from the string.
    """
    match = re.match(r"^([" + STRING_PREFIX_CHARS + r"]*)(.*)$", s, re.DOTALL)
    assert match is not None, f"failed to match string {s!r}"
    orig_prefix = match.group(1)
    new_prefix = orig_prefix.replace("F", "f").replace("B", "b").replace("U", "u")
    if remove_u_prefix:
        new_prefix = new_prefix.replace("u", "")
    return f"{new_prefix}{match.group(2)}"


def normalize_string_quotes(s: str) -> str:
    """Prefer double quotes but only if it doesn't cause more escaping.

    Adds or removes backslashes as appropriate. Doesn't parse and fix
    strings nested in f-strings.
    """
    value = s.lstrip(STRING_PREFIX_CHARS)
    if value[:3] == '"""':
        return s

    elif value[:3] == "'''":
        orig_quote = "'''"
        new_quote = '"""'
    elif value[0] == '"':
        orig_quote = '"'
        new_quote = "'"
    else:
        orig_quote = "'"
        new_quote = '"'
    first_quote_pos = s.find(orig_quote)
    if first_quote_pos == -1:
        return s  # There's an internal error

    prefix = s[:first_quote_pos]
    unescaped_new_quote = re.compile(rf"(([^\\]|^)(\\\\)*){new_quote}")
    escaped_new_quote = re.compile(rf"([^\\]|^)\\((?:\\\\)*){new_quote}")
    escaped_orig_quote = re.compile(rf"([^\\]|^)\\((?:\\\\)*){orig_quote}")
    body = s[first_quote_pos + len(orig_quote) : -len(orig_quote)]
    if "r" in prefix.casefold():
        if unescaped_new_quote.search(body):
            # There's at least one unescaped new_quote in this raw string
            # so converting is impossible
            return s

        # Do not introduce or remove backslashes in raw strings
        new_body = body
    else:
        # remove unnecessary escapes
        new_body = sub_twice(escaped_new_quote, rf"\1\2{new_quote}", body)
        if body != new_body:
            # Consider the string without unnecessary escapes as the original
            body = new_body
            s = f"{prefix}{orig_quote}{body}{orig_quote}"
        new_body = sub_twice(escaped_orig_quote, rf"\1\2{orig_quote}", new_body)
        new_body = sub_twice(unescaped_new_quote, rf"\1\\{new_quote}", new_body)
    if "f" in prefix.casefold():
        matches = re.findall(
            r"""
            (?:[^{]|^)\{  # start of the string or a non-{ followed by a single {
                ([^{].*?)  # contents of the brackets except if begins with {{
            \}(?:[^}]|$)  # A } followed by end of the string or a non-}
            """,
            new_body,
            re.VERBOSE,
        )
        for m in matches:
            if "\\" in str(m):
                # Do not introduce backslashes in interpolated expressions
                return s

    if new_quote == '"""' and new_body[-1:] == '"':
        # edge case:
        new_body = new_body[:-1] + '\\"'
    orig_escape_count = body.count("\\")
    new_escape_count = new_body.count("\\")
    if new_escape_count > orig_escape_count:
        return s  # Do not introduce more escaping

    if new_escape_count == orig_escape_count and orig_quote == '"':
        return s  # Prefer double quotes

    return f"{prefix}{new_quote}{new_body}{new_quote}"
