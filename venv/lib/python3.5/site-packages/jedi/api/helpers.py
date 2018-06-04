"""
Helpers for the API
"""
import re
from collections import namedtuple
from textwrap import dedent

from parso.python.parser import Parser
from parso.python import tree

from jedi._compatibility import u
from jedi.evaluate.syntax_tree import eval_atom
from jedi.evaluate.helpers import evaluate_call_of_leaf
from jedi.evaluate.compiled import get_string_context_set
from jedi.evaluate.base_context import ContextSet
from jedi.cache import call_signature_time_cache


CompletionParts = namedtuple('CompletionParts', ['path', 'has_dot', 'name'])


def sorted_definitions(defs):
    # Note: `or ''` below is required because `module_path` could be
    return sorted(defs, key=lambda x: (x.module_path or '', x.line or 0, x.column or 0))


def get_on_completion_name(module_node, lines, position):
    leaf = module_node.get_leaf_for_position(position)
    if leaf is None or leaf.type in ('string', 'error_leaf'):
        # Completions inside strings are a bit special, we need to parse the
        # string. The same is true for comments and error_leafs.
        line = lines[position[0] - 1]
        # The first step of completions is to get the name
        return re.search(r'(?!\d)\w+$|$', line[:position[1]]).group(0)
    elif leaf.type not in ('name', 'keyword'):
        return ''

    return leaf.value[:position[1] - leaf.start_pos[1]]


def _get_code(code_lines, start_pos, end_pos):
    # Get relevant lines.
    lines = code_lines[start_pos[0] - 1:end_pos[0]]
    # Remove the parts at the end of the line.
    lines[-1] = lines[-1][:end_pos[1]]
    # Remove first line indentation.
    lines[0] = lines[0][start_pos[1]:]
    return ''.join(lines)


class OnErrorLeaf(Exception):
    @property
    def error_leaf(self):
        return self.args[0]


def _get_code_for_stack(code_lines, module_node, position):
    leaf = module_node.get_leaf_for_position(position, include_prefixes=True)
    # It might happen that we're on whitespace or on a comment. This means
    # that we would not get the right leaf.
    if leaf.start_pos >= position:
        # If we're not on a comment simply get the previous leaf and proceed.
        leaf = leaf.get_previous_leaf()
        if leaf is None:
            return u('')  # At the beginning of the file.

    is_after_newline = leaf.type == 'newline'
    while leaf.type == 'newline':
        leaf = leaf.get_previous_leaf()
        if leaf is None:
            return u('')

    if leaf.type == 'error_leaf' or leaf.type == 'string':
        if leaf.start_pos[0] < position[0]:
            # On a different line, we just begin anew.
            return u('')

        # Error leafs cannot be parsed, completion in strings is also
        # impossible.
        raise OnErrorLeaf(leaf)
    else:
        user_stmt = leaf
        while True:
            if user_stmt.parent.type in ('file_input', 'suite', 'simple_stmt'):
                break
            user_stmt = user_stmt.parent

        if is_after_newline:
            if user_stmt.start_pos[1] > position[1]:
                # This means that it's actually a dedent and that means that we
                # start without context (part of a suite).
                return u('')

        # This is basically getting the relevant lines.
        return _get_code(code_lines, user_stmt.get_start_pos_of_prefix(), position)


def get_stack_at_position(grammar, code_lines, module_node, pos):
    """
    Returns the possible node names (e.g. import_from, xor_test or yield_stmt).
    """
    class EndMarkerReached(Exception):
        pass

    def tokenize_without_endmarker(code):
        # TODO This is for now not an official parso API that exists purely
        #   for Jedi.
        tokens = grammar._tokenize(code)
        for token_ in tokens:
            if token_.string == safeword:
                raise EndMarkerReached()
            elif token_.prefix.endswith(safeword):
                # This happens with comments.
                raise EndMarkerReached()
            else:
                yield token_

    # The code might be indedented, just remove it.
    code = dedent(_get_code_for_stack(code_lines, module_node, pos))
    # We use a word to tell Jedi when we have reached the start of the
    # completion.
    # Use Z as a prefix because it's not part of a number suffix.
    safeword = 'ZZZ_USER_WANTS_TO_COMPLETE_HERE_WITH_JEDI'
    code = code + ' ' + safeword

    p = Parser(grammar._pgen_grammar, error_recovery=True)
    try:
        p.parse(tokens=tokenize_without_endmarker(code))
    except EndMarkerReached:
        return Stack(p.pgen_parser.stack)
    raise SystemError("This really shouldn't happen. There's a bug in Jedi.")


class Stack(list):
    def get_node_names(self, grammar):
        for dfa, state, (node_number, nodes) in self:
            yield grammar.number2symbol[node_number]

    def get_nodes(self):
        for dfa, state, (node_number, nodes) in self:
            for node in nodes:
                yield node


def get_possible_completion_types(pgen_grammar, stack):
    def add_results(label_index):
        try:
            grammar_labels.append(inversed_tokens[label_index])
        except KeyError:
            try:
                keywords.append(inversed_keywords[label_index])
            except KeyError:
                t, v = pgen_grammar.labels[label_index]
                assert t >= 256
                # See if it's a symbol and if we're in its first set
                inversed_keywords
                itsdfa = pgen_grammar.dfas[t]
                itsstates, itsfirst = itsdfa
                for first_label_index in itsfirst.keys():
                    add_results(first_label_index)

    inversed_keywords = dict((v, k) for k, v in pgen_grammar.keywords.items())
    inversed_tokens = dict((v, k) for k, v in pgen_grammar.tokens.items())

    keywords = []
    grammar_labels = []

    def scan_stack(index):
        dfa, state, node = stack[index]
        states, first = dfa
        arcs = states[state]

        for label_index, new_state in arcs:
            if label_index == 0:
                # An accepting state, check the stack below.
                scan_stack(index - 1)
            else:
                add_results(label_index)

    scan_stack(-1)

    return keywords, grammar_labels


def evaluate_goto_definition(evaluator, context, leaf):
    if leaf.type == 'name':
        # In case of a name we can just use goto_definition which does all the
        # magic itself.
        return evaluator.goto_definitions(context, leaf)

    parent = leaf.parent
    if parent.type == 'atom':
        return context.eval_node(leaf.parent)
    elif parent.type == 'trailer':
        return evaluate_call_of_leaf(context, leaf)
    elif isinstance(leaf, tree.Literal):
        return eval_atom(context, leaf)
    elif leaf.type in ('fstring_string', 'fstring_start', 'fstring_end'):
        return get_string_context_set(evaluator)
    return []


CallSignatureDetails = namedtuple(
    'CallSignatureDetails',
    ['bracket_leaf', 'call_index', 'keyword_name_str']
)


def _get_index_and_key(nodes, position):
    """
    Returns the amount of commas and the keyword argument string.
    """
    nodes_before = [c for c in nodes if c.start_pos < position]
    if nodes_before[-1].type == 'arglist':
        nodes_before = [c for c in nodes_before[-1].children if c.start_pos < position]

    key_str = None

    if nodes_before:
        last = nodes_before[-1]
        if last.type == 'argument' and last.children[1].end_pos <= position:
            # Checked if the argument
            key_str = last.children[0].value
        elif last == '=':
            key_str = nodes_before[-2].value

    return nodes_before.count(','), key_str


def _get_call_signature_details_from_error_node(node, position):
    for index, element in reversed(list(enumerate(node.children))):
        # `index > 0` means that it's a trailer and not an atom.
        if element == '(' and element.end_pos <= position and index > 0:
            # It's an error node, we don't want to match too much, just
            # until the parentheses is enough.
            children = node.children[index:]
            name = element.get_previous_leaf()
            if name is None:
                continue
            if name.type == 'name' or name.parent.type in ('trailer', 'atom'):
                return CallSignatureDetails(
                    element,
                    *_get_index_and_key(children, position)
                )


def get_call_signature_details(module, position):
    leaf = module.get_leaf_for_position(position, include_prefixes=True)
    if leaf.start_pos >= position:
        # Whitespace / comments after the leaf count towards the previous leaf.
        leaf = leaf.get_previous_leaf()
        if leaf is None:
            return None

    if leaf == ')':
        if leaf.end_pos == position:
            leaf = leaf.get_next_leaf()

    # Now that we know where we are in the syntax tree, we start to look at
    # parents for possible function definitions.
    node = leaf.parent
    while node is not None:
        if node.type in ('funcdef', 'classdef'):
            # Don't show call signatures if there's stuff before it that just
            # makes it feel strange to have a call signature.
            return None

        for n in node.children[::-1]:
            if n.start_pos < position and n.type == 'error_node':
                result = _get_call_signature_details_from_error_node(n, position)
                if result is not None:
                    return result

        if node.type == 'trailer' and node.children[0] == '(':
            leaf = node.get_previous_leaf()
            if leaf is None:
                return None
            return CallSignatureDetails(
                node.children[0], *_get_index_and_key(node.children, position))

        node = node.parent

    return None


@call_signature_time_cache("call_signatures_validity")
def cache_call_signatures(evaluator, context, bracket_leaf, code_lines, user_pos):
    """This function calculates the cache key."""
    line_index = user_pos[0] - 1

    before_cursor = code_lines[line_index][:user_pos[1]]
    other_lines = code_lines[bracket_leaf.start_pos[0]:line_index]
    whole = ''.join(other_lines + [before_cursor])
    before_bracket = re.match(r'.*\(', whole, re.DOTALL)

    module_path = context.get_root_context().py__file__()
    if module_path is None:
        yield None  # Don't cache!
    else:
        yield (module_path, before_bracket, bracket_leaf.start_pos)
    yield evaluate_goto_definition(
        evaluator,
        context,
        bracket_leaf.get_previous_leaf()
    )
