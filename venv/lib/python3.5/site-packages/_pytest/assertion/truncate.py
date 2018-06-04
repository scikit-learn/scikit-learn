"""
Utilities for truncating assertion output.

Current default behaviour is to truncate assertion explanations at
~8 terminal lines, unless running in "-vv" mode or running on CI.
"""
from __future__ import absolute_import, division, print_function
import os

import six


DEFAULT_MAX_LINES = 8
DEFAULT_MAX_CHARS = 8 * 80
USAGE_MSG = "use '-vv' to show"


def truncate_if_required(explanation, item, max_length=None):
    """
    Truncate this assertion explanation if the given test item is eligible.
    """
    if _should_truncate_item(item):
        return _truncate_explanation(explanation)
    return explanation


def _should_truncate_item(item):
    """
    Whether or not this test item is eligible for truncation.
    """
    verbose = item.config.option.verbose
    return verbose < 2 and not _running_on_ci()


def _running_on_ci():
    """Check if we're currently running on a CI system."""
    env_vars = ['CI', 'BUILD_NUMBER']
    return any(var in os.environ for var in env_vars)


def _truncate_explanation(input_lines, max_lines=None, max_chars=None):
    """
    Truncate given list of strings that makes up the assertion explanation.

    Truncates to either 8 lines, or 640 characters - whichever the input reaches
    first. The remaining lines will be replaced by a usage message.
    """

    if max_lines is None:
        max_lines = DEFAULT_MAX_LINES
    if max_chars is None:
        max_chars = DEFAULT_MAX_CHARS

    # Check if truncation required
    input_char_count = len("".join(input_lines))
    if len(input_lines) <= max_lines and input_char_count <= max_chars:
        return input_lines

    # Truncate first to max_lines, and then truncate to max_chars if max_chars
    # is exceeded.
    truncated_explanation = input_lines[:max_lines]
    truncated_explanation = _truncate_by_char_count(truncated_explanation, max_chars)

    # Add ellipsis to final line
    truncated_explanation[-1] = truncated_explanation[-1] + "..."

    # Append useful message to explanation
    truncated_line_count = len(input_lines) - len(truncated_explanation)
    truncated_line_count += 1  # Account for the part-truncated final line
    msg = '...Full output truncated'
    if truncated_line_count == 1:
        msg += ' ({0} line hidden)'.format(truncated_line_count)
    else:
        msg += ' ({0} lines hidden)'.format(truncated_line_count)
    msg += ", {0}" .format(USAGE_MSG)
    truncated_explanation.extend([
        six.text_type(""),
        six.text_type(msg),
    ])
    return truncated_explanation


def _truncate_by_char_count(input_lines, max_chars):
    # Check if truncation required
    if len("".join(input_lines)) <= max_chars:
        return input_lines

    # Find point at which input length exceeds total allowed length
    iterated_char_count = 0
    for iterated_index, input_line in enumerate(input_lines):
        if iterated_char_count + len(input_line) > max_chars:
            break
        iterated_char_count += len(input_line)

    # Create truncated explanation with modified final line
    truncated_result = input_lines[:iterated_index]
    final_line = input_lines[iterated_index]
    if final_line:
        final_line_truncate_point = max_chars - iterated_char_count
        final_line = final_line[:final_line_truncate_point]
    truncated_result.append(final_line)
    return truncated_result
