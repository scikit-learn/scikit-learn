"""Utilities for assertion debugging"""
from __future__ import absolute_import, division, print_function
import pprint

import _pytest._code
import py
import six
from ..compat import Sequence

u = six.text_type

# The _reprcompare attribute on the util module is used by the new assertion
# interpretation code and assertion rewriter to detect this plugin was
# loaded and in turn call the hooks defined here as part of the
# DebugInterpreter.
_reprcompare = None


# the re-encoding is needed for python2 repr
# with non-ascii characters (see issue 877 and 1379)
def ecu(s):
    try:
        return u(s, 'utf-8', 'replace')
    except TypeError:
        return s


def format_explanation(explanation):
    """This formats an explanation

    Normally all embedded newlines are escaped, however there are
    three exceptions: \n{, \n} and \n~.  The first two are intended
    cover nested explanations, see function and attribute explanations
    for examples (.visit_Call(), visit_Attribute()).  The last one is
    for when one explanation needs to span multiple lines, e.g. when
    displaying diffs.
    """
    explanation = ecu(explanation)
    lines = _split_explanation(explanation)
    result = _format_lines(lines)
    return u('\n').join(result)


def _split_explanation(explanation):
    """Return a list of individual lines in the explanation

    This will return a list of lines split on '\n{', '\n}' and '\n~'.
    Any other newlines will be escaped and appear in the line as the
    literal '\n' characters.
    """
    raw_lines = (explanation or u('')).split('\n')
    lines = [raw_lines[0]]
    for values in raw_lines[1:]:
        if values and values[0] in ['{', '}', '~', '>']:
            lines.append(values)
        else:
            lines[-1] += '\\n' + values
    return lines


def _format_lines(lines):
    """Format the individual lines

    This will replace the '{', '}' and '~' characters of our mini
    formatting language with the proper 'where ...', 'and ...' and ' +
    ...' text, taking care of indentation along the way.

    Return a list of formatted lines.
    """
    result = lines[:1]
    stack = [0]
    stackcnt = [0]
    for line in lines[1:]:
        if line.startswith('{'):
            if stackcnt[-1]:
                s = u('and   ')
            else:
                s = u('where ')
            stack.append(len(result))
            stackcnt[-1] += 1
            stackcnt.append(0)
            result.append(u(' +') + u('  ') * (len(stack) - 1) + s + line[1:])
        elif line.startswith('}'):
            stack.pop()
            stackcnt.pop()
            result[stack[-1]] += line[1:]
        else:
            assert line[0] in ['~', '>']
            stack[-1] += 1
            indent = len(stack) if line.startswith('~') else len(stack) - 1
            result.append(u('  ') * indent + line[1:])
    assert len(stack) == 1
    return result


# Provide basestring in python3
try:
    basestring = basestring
except NameError:
    basestring = str


def assertrepr_compare(config, op, left, right):
    """Return specialised explanations for some operators/operands"""
    width = 80 - 15 - len(op) - 2  # 15 chars indentation, 1 space around op
    left_repr = py.io.saferepr(left, maxsize=int(width // 2))
    right_repr = py.io.saferepr(right, maxsize=width - len(left_repr))

    summary = u('%s %s %s') % (ecu(left_repr), op, ecu(right_repr))

    def issequence(x):
        return isinstance(x, Sequence) and not isinstance(x, basestring)

    def istext(x):
        return isinstance(x, basestring)

    def isdict(x):
        return isinstance(x, dict)

    def isset(x):
        return isinstance(x, (set, frozenset))

    def isiterable(obj):
        try:
            iter(obj)
            return not istext(obj)
        except TypeError:
            return False

    verbose = config.getoption('verbose')
    explanation = None
    try:
        if op == '==':
            if istext(left) and istext(right):
                explanation = _diff_text(left, right, verbose)
            else:
                if issequence(left) and issequence(right):
                    explanation = _compare_eq_sequence(left, right, verbose)
                elif isset(left) and isset(right):
                    explanation = _compare_eq_set(left, right, verbose)
                elif isdict(left) and isdict(right):
                    explanation = _compare_eq_dict(left, right, verbose)
                if isiterable(left) and isiterable(right):
                    expl = _compare_eq_iterable(left, right, verbose)
                    if explanation is not None:
                        explanation.extend(expl)
                    else:
                        explanation = expl
        elif op == 'not in':
            if istext(left) and istext(right):
                explanation = _notin_text(left, right, verbose)
    except Exception:
        explanation = [
            u('(pytest_assertion plugin: representation of details failed.  '
              'Probably an object has a faulty __repr__.)'),
            u(_pytest._code.ExceptionInfo())]

    if not explanation:
        return None

    return [summary] + explanation


def _diff_text(left, right, verbose=False):
    """Return the explanation for the diff between text or bytes

    Unless --verbose is used this will skip leading and trailing
    characters which are identical to keep the diff minimal.

    If the input are bytes they will be safely converted to text.
    """
    from difflib import ndiff
    explanation = []
    if isinstance(left, six.binary_type):
        left = u(repr(left)[1:-1]).replace(r'\n', '\n')
    if isinstance(right, six.binary_type):
        right = u(repr(right)[1:-1]).replace(r'\n', '\n')
    if not verbose:
        i = 0  # just in case left or right has zero length
        for i in range(min(len(left), len(right))):
            if left[i] != right[i]:
                break
        if i > 42:
            i -= 10                 # Provide some context
            explanation = [u('Skipping %s identical leading '
                             'characters in diff, use -v to show') % i]
            left = left[i:]
            right = right[i:]
        if len(left) == len(right):
            for i in range(len(left)):
                if left[-i] != right[-i]:
                    break
            if i > 42:
                i -= 10     # Provide some context
                explanation += [u('Skipping %s identical trailing '
                                  'characters in diff, use -v to show') % i]
                left = left[:-i]
                right = right[:-i]
    keepends = True
    explanation += [line.strip('\n')
                    for line in ndiff(left.splitlines(keepends),
                                      right.splitlines(keepends))]
    return explanation


def _compare_eq_iterable(left, right, verbose=False):
    if not verbose:
        return [u('Use -v to get the full diff')]
    # dynamic import to speedup pytest
    import difflib

    try:
        left_formatting = pprint.pformat(left).splitlines()
        right_formatting = pprint.pformat(right).splitlines()
        explanation = [u('Full diff:')]
    except Exception:
        # hack: PrettyPrinter.pformat() in python 2 fails when formatting items that can't be sorted(), ie, calling
        # sorted() on a list would raise. See issue #718.
        # As a workaround, the full diff is generated by using the repr() string of each item of each container.
        left_formatting = sorted(repr(x) for x in left)
        right_formatting = sorted(repr(x) for x in right)
        explanation = [u('Full diff (fallback to calling repr on each item):')]
    explanation.extend(line.strip() for line in difflib.ndiff(left_formatting, right_formatting))
    return explanation


def _compare_eq_sequence(left, right, verbose=False):
    explanation = []
    for i in range(min(len(left), len(right))):
        if left[i] != right[i]:
            explanation += [u('At index %s diff: %r != %r')
                            % (i, left[i], right[i])]
            break
    if len(left) > len(right):
        explanation += [u('Left contains more items, first extra item: %s')
                        % py.io.saferepr(left[len(right)],)]
    elif len(left) < len(right):
        explanation += [
            u('Right contains more items, first extra item: %s') %
            py.io.saferepr(right[len(left)],)]
    return explanation


def _compare_eq_set(left, right, verbose=False):
    explanation = []
    diff_left = left - right
    diff_right = right - left
    if diff_left:
        explanation.append(u('Extra items in the left set:'))
        for item in diff_left:
            explanation.append(py.io.saferepr(item))
    if diff_right:
        explanation.append(u('Extra items in the right set:'))
        for item in diff_right:
            explanation.append(py.io.saferepr(item))
    return explanation


def _compare_eq_dict(left, right, verbose=False):
    explanation = []
    common = set(left).intersection(set(right))
    same = dict((k, left[k]) for k in common if left[k] == right[k])
    if same and verbose < 2:
        explanation += [u('Omitting %s identical items, use -vv to show') %
                        len(same)]
    elif same:
        explanation += [u('Common items:')]
        explanation += pprint.pformat(same).splitlines()
    diff = set(k for k in common if left[k] != right[k])
    if diff:
        explanation += [u('Differing items:')]
        for k in diff:
            explanation += [py.io.saferepr({k: left[k]}) + ' != ' +
                            py.io.saferepr({k: right[k]})]
    extra_left = set(left) - set(right)
    if extra_left:
        explanation.append(u('Left contains more items:'))
        explanation.extend(pprint.pformat(
            dict((k, left[k]) for k in extra_left)).splitlines())
    extra_right = set(right) - set(left)
    if extra_right:
        explanation.append(u('Right contains more items:'))
        explanation.extend(pprint.pformat(
            dict((k, right[k]) for k in extra_right)).splitlines())
    return explanation


def _notin_text(term, text, verbose=False):
    index = text.find(term)
    head = text[:index]
    tail = text[index + len(term):]
    correct_text = head + tail
    diff = _diff_text(correct_text, text, verbose)
    newdiff = [u('%s is contained here:') % py.io.saferepr(term, maxsize=42)]
    for line in diff:
        if line.startswith(u('Skipping')):
            continue
        if line.startswith(u('- ')):
            continue
        if line.startswith(u('+ ')):
            newdiff.append(u('  ') + line[2:])
        else:
            newdiff.append(line)
    return newdiff
