r"""
Compiler for a regular grammar.

Example usage::

    # Create and compile grammar.
    p = compile('add \s+ (?P<var1>[^\s]+)  \s+  (?P<var2>[^\s]+)')

    # Match input string.
    m = p.match('add 23 432')

    # Get variables.
    m.variables().get('var1')  # Returns "23"
    m.variables().get('var2')  # Returns "432"


Partial matches are possible::

    # Create and compile grammar.
    p = compile('''
        # Operators with two arguments.
        ((?P<operator1>[^\s]+)  \s+ (?P<var1>[^\s]+)  \s+  (?P<var2>[^\s]+)) |

        # Operators with only one arguments.
        ((?P<operator2>[^\s]+)  \s+ (?P<var1>[^\s]+))
    ''')

    # Match partial input string.
    m = p.match_prefix('add 23')

    # Get variables. (Notice that both operator1 and operator2 contain the
    # value "add".) This is because our input is incomplete, and we don't know
    # yet in which rule of the regex we we'll end up. It could also be that
    # `operator1` and `operator2` have a different autocompleter and we want to
    # call all possible autocompleters that would result in valid input.)
    m.variables().get('var1')  # Returns "23"
    m.variables().get('operator1')  # Returns "add"
    m.variables().get('operator2')  # Returns "add"

"""
from __future__ import unicode_literals
import re

from six.moves import range
from .regex_parser import Any, Sequence, Regex, Variable, Repeat, Lookahead
from .regex_parser import parse_regex, tokenize_regex

__all__ = (
    'compile',
)


# Name of the named group in the regex, matching trailing input.
# (Trailing input is when the input contains characters after the end of the
# expression has been matched.)
_INVALID_TRAILING_INPUT = 'invalid_trailing'


class _CompiledGrammar(object):
    """
    Compiles a grammar. This will take the parse tree of a regular expression
    and compile the grammar.

    :param root_node: :class~`.regex_parser.Node` instance.
    :param escape_funcs: `dict` mapping variable names to escape callables.
    :param unescape_funcs: `dict` mapping variable names to unescape callables.
    """
    def __init__(self, root_node, escape_funcs=None, unescape_funcs=None):
        self.root_node = root_node
        self.escape_funcs = escape_funcs or {}
        self.unescape_funcs = unescape_funcs or {}

        #: Dictionary that will map the redex names to Node instances.
        self._group_names_to_nodes = {}  # Maps regex group names to varnames.
        counter = [0]

        def create_group_func(node):
            name = 'n%s' % counter[0]
            self._group_names_to_nodes[name] = node.varname
            counter[0] += 1
            return name

        # Compile regex strings.
        self._re_pattern = '^%s$' % self._transform(root_node, create_group_func)
        self._re_prefix_patterns = list(self._transform_prefix(root_node, create_group_func))

        # Compile the regex itself.
        flags = re.DOTALL  # Note that we don't need re.MULTILINE! (^ and $
                           # still represent the start and end of input text.)
        self._re = re.compile(self._re_pattern, flags)
        self._re_prefix = [re.compile(t, flags) for t in self._re_prefix_patterns]

        # We compile one more set of regexes, similar to `_re_prefix`, but accept any trailing
        # input. This will ensure that we can still highlight the input correctly, even when the
        # input contains some additional characters at the end that don't match the grammar.)
        self._re_prefix_with_trailing_input = [
            re.compile(r'(?:%s)(?P<%s>.*?)$' % (t.rstrip('$'), _INVALID_TRAILING_INPUT), flags)
            for t in self._re_prefix_patterns]

    def escape(self, varname, value):
        """
        Escape `value` to fit in the place of this variable into the grammar.
        """
        f = self.escape_funcs.get(varname)
        return f(value) if f else value

    def unescape(self, varname, value):
        """
        Unescape `value`.
        """
        f = self.unescape_funcs.get(varname)
        return f(value) if f else value

    @classmethod
    def _transform(cls, root_node, create_group_func):
        """
        Turn a :class:`Node` object into a regular expression.

        :param root_node: The :class:`Node` instance for which we generate the grammar.
        :param create_group_func: A callable which takes a `Node` and returns the next
            free name for this node.
        """
        def transform(node):
            # Turn `Any` into an OR.
            if isinstance(node, Any):
                return '(?:%s)' % '|'.join(transform(c) for c in node.children)

            # Concatenate a `Sequence`
            elif isinstance(node, Sequence):
                return ''.join(transform(c) for c in node.children)

            # For Regex and Lookahead nodes, just insert them literally.
            elif isinstance(node, Regex):
                return node.regex

            elif isinstance(node, Lookahead):
                before = ('(?!' if node.negative else '(=')
                return before + transform(node.childnode) + ')'

            # A `Variable` wraps the children into a named group.
            elif isinstance(node, Variable):
                return '(?P<%s>%s)' % (create_group_func(node), transform(node.childnode))

            # `Repeat`.
            elif isinstance(node, Repeat):
                return '(?:%s){%i,%s}%s' % (
                    transform(node.childnode), node.min_repeat,
                    ('' if node.max_repeat is None else str(node.max_repeat)),
                    ('' if node.greedy else '?')
                )
            else:
                raise TypeError('Got %r' % (node, ))

        return transform(root_node)

    @classmethod
    def _transform_prefix(cls, root_node, create_group_func):
        """
        Yield all the regular expressions matching a prefix of the grammar
        defined by the `Node` instance.

        This can yield multiple expressions, because in the case of on OR
        operation in the grammar, we can have another outcome depending on
        which clause would appear first. E.g. "(A|B)C" is not the same as
        "(B|A)C" because the regex engine is lazy and takes the first match.
        However, because we the current input is actually a prefix of the
        grammar which meight not yet contain the data for "C", we need to know
        both intermediate states, in order to call the appropriate
        autocompletion for both cases.

        :param root_node: The :class:`Node` instance for which we generate the grammar.
        :param create_group_func: A callable which takes a `Node` and returns the next
            free name for this node.
        """
        def transform(node):
            # Generate regexes for all permutations of this OR. Each node
            # should be in front once.
            if isinstance(node, Any):
                for c in node.children:
                    for r in transform(c):
                        yield '(?:%s)?' % r

            # For a sequence. We can either have a match for the sequence
            # of all the children, or for an exact match of the first X
            # children, followed by a partial match of the next children.
            elif isinstance(node, Sequence):
                for i in range(len(node.children)):
                    a = [cls._transform(c, create_group_func) for c in node.children[:i]]
                    for c in transform(node.children[i]):
                        yield '(?:%s)' % (''.join(a) + c)

            elif isinstance(node, Regex):
                yield '(?:%s)?' % node.regex

            elif isinstance(node, Lookahead):
                if node.negative:
                    yield '(?!%s)' % cls._transform(node.childnode, create_group_func)
                else:
                    # Not sure what the correct semantics are in this case.
                    # (Probably it's not worth implementing this.)
                    raise Exception('Positive lookahead not yet supported.')

            elif isinstance(node, Variable):
                # (Note that we should not append a '?' here. the 'transform'
                # method will already recursively do that.)
                for c in transform(node.childnode):
                    yield '(?P<%s>%s)' % (create_group_func(node), c)

            elif isinstance(node, Repeat):
                # If we have a repetition of 8 times. That would mean that the
                # current input could have for instance 7 times a complete
                # match, followed by a partial match.
                prefix = cls._transform(node.childnode, create_group_func)

                for c in transform(node.childnode):
                    if node.max_repeat:
                        repeat_sign = '{,%i}' % (node.max_repeat - 1)
                    else:
                        repeat_sign = '*'
                    yield '(?:%s)%s%s(?:%s)?' % (
                        prefix,
                        repeat_sign,
                        ('' if node.greedy else '?'),
                        c)

            else:
                raise TypeError('Got %r' % node)

        for r in transform(root_node):
            yield '^%s$' % r

    def match(self, string):
        """
        Match the string with the grammar.
        Returns a :class:`Match` instance or `None` when the input doesn't match the grammar.

        :param string: The input string.
        """
        m = self._re.match(string)

        if m:
            return Match(string, [(self._re, m)], self._group_names_to_nodes, self.unescape_funcs)

    def match_prefix(self, string):
        """
        Do a partial match of the string with the grammar. The returned
        :class:`Match` instance can contain multiple representations of the
        match. This will never return `None`. If it doesn't match at all, the "trailing input"
        part will capture all of the input.

        :param string: The input string.
        """
        # First try to match using `_re_prefix`. If nothing is found, use the patterns that
        # also accept trailing characters.
        for patterns in [self._re_prefix, self._re_prefix_with_trailing_input]:
            matches = [(r, r.match(string)) for r in patterns]
            matches = [(r, m) for r, m in matches if m]

            if matches != []:
                return Match(string, matches, self._group_names_to_nodes, self.unescape_funcs)


class Match(object):
    """
    :param string: The input string.
    :param re_matches: List of (compiled_re_pattern, re_match) tuples.
    :param group_names_to_nodes: Dictionary mapping all the re group names to the matching Node instances.
    """
    def __init__(self, string, re_matches, group_names_to_nodes, unescape_funcs):
        self.string = string
        self._re_matches = re_matches
        self._group_names_to_nodes = group_names_to_nodes
        self._unescape_funcs = unescape_funcs

    def _nodes_to_regs(self):
        """
        Return a list of (varname, reg) tuples.
        """
        def get_tuples():
            for r, re_match in self._re_matches:
                for group_name, group_index in r.groupindex.items():
                    if group_name != _INVALID_TRAILING_INPUT:
                        reg = re_match.regs[group_index]
                        node = self._group_names_to_nodes[group_name]
                        yield (node, reg)

        return list(get_tuples())

    def _nodes_to_values(self):
        """
        Returns list of list of (Node, string_value) tuples.
        """
        def is_none(slice):
            return slice[0] == -1 and slice[1] == -1

        def get(slice):
            return self.string[slice[0]:slice[1]]

        return [(varname, get(slice), slice) for varname, slice in self._nodes_to_regs() if not is_none(slice)]

    def _unescape(self, varname, value):
        unwrapper = self._unescape_funcs.get(varname)
        return unwrapper(value) if unwrapper else value

    def variables(self):
        """
        Returns :class:`Variables` instance.
        """
        return Variables([(k, self._unescape(k, v), sl) for k, v, sl in self._nodes_to_values()])

    def trailing_input(self):
        """
        Get the `MatchVariable` instance, representing trailing input, if there is any.
        "Trailing input" is input at the end that does not match the grammar anymore, but
        when this is removed from the end of the input, the input would be a valid string.
        """
        slices = []

        # Find all regex group for the name _INVALID_TRAILING_INPUT.
        for r, re_match in self._re_matches:
            for group_name, group_index in r.groupindex.items():
                if group_name == _INVALID_TRAILING_INPUT:
                    slices.append(re_match.regs[group_index])

        # Take the smallest part. (Smaller trailing text means that a larger input has
        # been matched, so that is better.)
        if slices:
            slice = [max(i[0] for i in slices), max(i[1] for i in slices)]
            value = self.string[slice[0]:slice[1]]
            return MatchVariable('<trailing_input>', value, slice)

    def end_nodes(self):
        """
        Yields `MatchVariable` instances for all the nodes having their end
        position at the end of the input string.
        """
        for varname, reg in self._nodes_to_regs():
            # If this part goes until the end of the input string.
            if reg[1] == len(self.string):
                value = self._unescape(varname, self.string[reg[0]: reg[1]])
                yield MatchVariable(varname, value, (reg[0], reg[1]))


class Variables(object):
    def __init__(self, tuples):
        #: List of (varname, value, slice) tuples.
        self._tuples = tuples

    def __repr__(self):
        return '%s(%s)' % (
            self.__class__.__name__, ', '.join('%s=%r' % (k, v) for k, v, _ in self._tuples))

    def get(self, key, default=None):
        items = self.getall(key)
        return items[0] if items else default

    def getall(self, key):
        return [v for k, v, _ in self._tuples if k == key]

    def __getitem__(self, key):
        return self.get(key)

    def __iter__(self):
        """
        Yield `MatchVariable` instances.
        """
        for varname, value, slice in self._tuples:
            yield MatchVariable(varname, value, slice)


class MatchVariable(object):
    """
    Represents a match of a variable in the grammar.

    :param varname: (string) Name of the variable.
    :param value: (string) Value of this variable.
    :param slice: (start, stop) tuple, indicating the position of this variable
                  in the input string.
    """
    def __init__(self, varname, value, slice):
        self.varname = varname
        self.value = value
        self.slice = slice

        self.start = self.slice[0]
        self.stop = self.slice[1]

    def __repr__(self):
        return '%s(%r, %r)' % (self.__class__.__name__, self.varname, self.value)


def compile(expression, escape_funcs=None, unescape_funcs=None):
    """
    Compile grammar (given as regex string), returning a `CompiledGrammar`
    instance.
    """
    return _compile_from_parse_tree(
        parse_regex(tokenize_regex(expression)),
        escape_funcs=escape_funcs,
        unescape_funcs=unescape_funcs)


def _compile_from_parse_tree(root_node, *a, **kw):
    """
    Compile grammar (given as parse tree), returning a `CompiledGrammar`
    instance.
    """
    return _CompiledGrammar(root_node, *a, **kw)
