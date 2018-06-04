"""
Functions evaluating the syntax tree.
"""
import copy

from parso.python import tree

from jedi._compatibility import force_unicode, unicode
from jedi import debug
from jedi import parser_utils
from jedi.evaluate.base_context import ContextSet, NO_CONTEXTS, ContextualizedNode, \
    ContextualizedName, iterator_to_context_set, iterate_contexts
from jedi.evaluate import compiled
from jedi.evaluate import pep0484
from jedi.evaluate import recursion
from jedi.evaluate import helpers
from jedi.evaluate import analysis
from jedi.evaluate import imports
from jedi.evaluate import arguments
from jedi.evaluate.pep0484 import _evaluate_for_annotation
from jedi.evaluate.context import ClassContext, FunctionContext
from jedi.evaluate.context import iterable
from jedi.evaluate.context import TreeInstance, CompiledInstance
from jedi.evaluate.finder import NameFinder
from jedi.evaluate.helpers import is_string, is_literal, is_number, is_compiled
from jedi.evaluate.compiled.access import COMPARISON_OPERATORS


def _limit_context_infers(func):
    """
    This is for now the way how we limit type inference going wild. There are
    other ways to ensure recursion limits as well. This is mostly necessary
    because of instance (self) access that can be quite tricky to limit.

    I'm still not sure this is the way to go, but it looks okay for now and we
    can still go anther way in the future. Tests are there. ~ dave
    """
    def wrapper(context, *args, **kwargs):
        n = context.tree_node
        evaluator = context.evaluator
        try:
            evaluator.inferred_element_counts[n] += 1
            if evaluator.inferred_element_counts[n] > 300:
                debug.warning('In context %s there were too many inferences.', n)
                return NO_CONTEXTS
        except KeyError:
            evaluator.inferred_element_counts[n] = 1
        return func(context, *args, **kwargs)

    return wrapper


def _py__stop_iteration_returns(generators):
    results = ContextSet()
    for generator in generators:
        try:
            method = generator.py__stop_iteration_returns
        except AttributeError:
            debug.warning('%s is not actually a generator', generator)
        else:
            results |= method()
    return results


@debug.increase_indent
@_limit_context_infers
def eval_node(context, element):
    debug.dbg('eval_node %s@%s', element, element.start_pos)
    evaluator = context.evaluator
    typ = element.type
    if typ in ('name', 'number', 'string', 'atom', 'strings'):
        return eval_atom(context, element)
    elif typ == 'keyword':
        # For False/True/None
        if element.value in ('False', 'True', 'None'):
            return ContextSet(compiled.builtin_from_name(evaluator, element.value))
        # else: print e.g. could be evaluated like this in Python 2.7
        return NO_CONTEXTS
    elif typ == 'lambdef':
        return ContextSet(FunctionContext(evaluator, context, element))
    elif typ == 'expr_stmt':
        return eval_expr_stmt(context, element)
    elif typ in ('power', 'atom_expr'):
        first_child = element.children[0]
        children = element.children[1:]
        had_await = False
        if first_child.type == 'keyword' and first_child.value == 'await':
            had_await = True
            first_child = children.pop(0)

        context_set = eval_atom(context, first_child)
        for trailer in children:
            if trailer == '**':  # has a power operation.
                right = context.eval_node(children[1])
                context_set = _eval_comparison(
                    evaluator,
                    context,
                    context_set,
                    trailer,
                    right
                )
                break
            context_set = eval_trailer(context, context_set, trailer)

        if had_await:
            await_context_set = context_set.py__getattribute__(u"__await__")
            if not await_context_set:
                debug.warning('Tried to run py__await__ on context %s', context)
            context_set = ContextSet()
            return _py__stop_iteration_returns(await_context_set.execute_evaluated())
        return context_set
    elif typ in ('testlist_star_expr', 'testlist',):
        # The implicit tuple in statements.
        return ContextSet(iterable.SequenceLiteralContext(evaluator, context, element))
    elif typ in ('not_test', 'factor'):
        context_set = context.eval_node(element.children[-1])
        for operator in element.children[:-1]:
            context_set = eval_factor(context_set, operator)
        return context_set
    elif typ == 'test':
        # `x if foo else y` case.
        return (context.eval_node(element.children[0]) |
                context.eval_node(element.children[-1]))
    elif typ == 'operator':
        # Must be an ellipsis, other operators are not evaluated.
        # In Python 2 ellipsis is coded as three single dot tokens, not
        # as one token 3 dot token.
        if element.value not in ('.', '...'):
            origin = element.parent
            raise AssertionError("unhandled operator %s in %s " % (repr(element.value), origin))
        return ContextSet(compiled.builtin_from_name(evaluator, u'Ellipsis'))
    elif typ == 'dotted_name':
        context_set = eval_atom(context, element.children[0])
        for next_name in element.children[2::2]:
            # TODO add search_global=True?
            context_set = context_set.py__getattribute__(next_name, name_context=context)
        return context_set
    elif typ == 'eval_input':
        return eval_node(context, element.children[0])
    elif typ == 'annassign':
        return pep0484._evaluate_for_annotation(context, element.children[1])
    elif typ == 'yield_expr':
        if len(element.children) and element.children[1].type == 'yield_arg':
            # Implies that it's a yield from.
            element = element.children[1].children[1]
            generators = context.eval_node(element)
            return _py__stop_iteration_returns(generators)

        # Generator.send() is not implemented.
        return NO_CONTEXTS
    else:
        return eval_or_test(context, element)


def eval_trailer(context, base_contexts, trailer):
    trailer_op, node = trailer.children[:2]
    if node == ')':  # `arglist` is optional.
        node = None

    if trailer_op == '[':
        trailer_op, node, _ = trailer.children

        # TODO It's kind of stupid to cast this from a context set to a set.
        foo = set(base_contexts)
        # special case: PEP0484 typing module, see
        # https://github.com/davidhalter/jedi/issues/663
        result = ContextSet()
        for typ in list(foo):
            if isinstance(typ, (ClassContext, TreeInstance)):
                typing_module_types = pep0484.py__getitem__(context, typ, node)
                if typing_module_types is not None:
                    foo.remove(typ)
                    result |= typing_module_types

        return result | base_contexts.get_item(
            eval_subscript_list(context.evaluator, context, node),
            ContextualizedNode(context, trailer)
        )
    else:
        debug.dbg('eval_trailer: %s in %s', trailer, base_contexts)
        if trailer_op == '.':
            return base_contexts.py__getattribute__(
                name_context=context,
                name_or_str=node
            )
        else:
            assert trailer_op == '(', 'trailer_op is actually %s' % trailer_op
            args = arguments.TreeArguments(context.evaluator, context, node, trailer)
            return base_contexts.execute(args)


def eval_atom(context, atom):
    """
    Basically to process ``atom`` nodes. The parser sometimes doesn't
    generate the node (because it has just one child). In that case an atom
    might be a name or a literal as well.
    """
    if atom.type == 'name':
        # This is the first global lookup.
        stmt = tree.search_ancestor(
            atom, 'expr_stmt', 'lambdef'
        ) or atom
        if stmt.type == 'lambdef':
            stmt = atom
        return context.py__getattribute__(
            name_or_str=atom,
            position=stmt.start_pos,
            search_global=True
        )

    elif isinstance(atom, tree.Literal):
        string = context.evaluator.compiled_subprocess.safe_literal_eval(atom.value)
        return ContextSet(compiled.create_simple_object(context.evaluator, string))
    elif atom.type == 'strings':
        # Will be multiple string.
        context_set = eval_atom(context, atom.children[0])
        for string in atom.children[1:]:
            right = eval_atom(context, string)
            context_set = _eval_comparison(context.evaluator, context, context_set, u'+', right)
        return context_set
    else:
        c = atom.children
        # Parentheses without commas are not tuples.
        if c[0] == '(' and not len(c) == 2 \
                and not(c[1].type == 'testlist_comp' and
                        len(c[1].children) > 1):
            return context.eval_node(c[1])

        try:
            comp_for = c[1].children[1]
        except (IndexError, AttributeError):
            pass
        else:
            if comp_for == ':':
                # Dict comprehensions have a colon at the 3rd index.
                try:
                    comp_for = c[1].children[3]
                except IndexError:
                    pass

            if comp_for.type == 'comp_for':
                return ContextSet(iterable.comprehension_from_atom(
                    context.evaluator, context, atom
                ))

        # It's a dict/list/tuple literal.
        array_node = c[1]
        try:
            array_node_c = array_node.children
        except AttributeError:
            array_node_c = []
        if c[0] == '{' and (array_node == '}' or ':' in array_node_c):
            context = iterable.DictLiteralContext(context.evaluator, context, atom)
        else:
            context = iterable.SequenceLiteralContext(context.evaluator, context, atom)
        return ContextSet(context)


@_limit_context_infers
def eval_expr_stmt(context, stmt, seek_name=None):
    with recursion.execution_allowed(context.evaluator, stmt) as allowed:
        # Here we allow list/set to recurse under certain conditions. To make
        # it possible to resolve stuff like list(set(list(x))), this is
        # necessary.
        if not allowed and context.get_root_context() == context.evaluator.builtins_module:
            try:
                instance = context.instance
            except AttributeError:
                pass
            else:
                if instance.name.string_name in ('list', 'set'):
                    c = instance.get_first_non_keyword_argument_contexts()
                    if instance not in c:
                        allowed = True

        if allowed:
            return _eval_expr_stmt(context, stmt, seek_name)
    return NO_CONTEXTS


@debug.increase_indent
def _eval_expr_stmt(context, stmt, seek_name=None):
    """
    The starting point of the completion. A statement always owns a call
    list, which are the calls, that a statement does. In case multiple
    names are defined in the statement, `seek_name` returns the result for
    this name.

    :param stmt: A `tree.ExprStmt`.
    """
    debug.dbg('eval_expr_stmt %s (%s)', stmt, seek_name)
    rhs = stmt.get_rhs()
    context_set = context.eval_node(rhs)

    if seek_name:
        c_node = ContextualizedName(context, seek_name)
        context_set = check_tuple_assignments(context.evaluator, c_node, context_set)

    first_operator = next(stmt.yield_operators(), None)
    if first_operator not in ('=', None) and first_operator.type == 'operator':
        # `=` is always the last character in aug assignments -> -1
        operator = copy.copy(first_operator)
        operator.value = operator.value[:-1]
        name = stmt.get_defined_names()[0].value
        left = context.py__getattribute__(
            name, position=stmt.start_pos, search_global=True)

        for_stmt = tree.search_ancestor(stmt, 'for_stmt')
        if for_stmt is not None and for_stmt.type == 'for_stmt' and context_set \
                and parser_utils.for_stmt_defines_one_name(for_stmt):
            # Iterate through result and add the values, that's possible
            # only in for loops without clutter, because they are
            # predictable. Also only do it, if the variable is not a tuple.
            node = for_stmt.get_testlist()
            cn = ContextualizedNode(context, node)
            ordered = list(cn.infer().iterate(cn))

            for lazy_context in ordered:
                dct = {for_stmt.children[1].value: lazy_context.infer()}
                with helpers.predefine_names(context, for_stmt, dct):
                    t = context.eval_node(rhs)
                    left = _eval_comparison(context.evaluator, context, left, operator, t)
            context_set = left
        else:
            context_set = _eval_comparison(context.evaluator, context, left, operator, context_set)
    debug.dbg('eval_expr_stmt result %s', context_set)
    return context_set


def eval_or_test(context, or_test):
    iterator = iter(or_test.children)
    types = context.eval_node(next(iterator))
    for operator in iterator:
        right = next(iterator)
        if operator.type == 'comp_op':  # not in / is not
            operator = ' '.join(c.value for c in operator.children)

        # handle lazy evaluation of and/or here.
        if operator in ('and', 'or'):
            left_bools = set(left.py__bool__() for left in types)
            if left_bools == {True}:
                if operator == 'and':
                    types = context.eval_node(right)
            elif left_bools == {False}:
                if operator != 'and':
                    types = context.eval_node(right)
            # Otherwise continue, because of uncertainty.
        else:
            types = _eval_comparison(context.evaluator, context, types, operator,
                                     context.eval_node(right))
    debug.dbg('eval_or_test types %s', types)
    return types


@iterator_to_context_set
def eval_factor(context_set, operator):
    """
    Calculates `+`, `-`, `~` and `not` prefixes.
    """
    for context in context_set:
        if operator == '-':
            if is_number(context):
                yield context.negate()
        elif operator == 'not':
            value = context.py__bool__()
            if value is None:  # Uncertainty.
                return
            yield compiled.create_simple_object(context.evaluator, not value)
        else:
            yield context


def _literals_to_types(evaluator, result):
    # Changes literals ('a', 1, 1.0, etc) to its type instances (str(),
    # int(), float(), etc).
    new_result = NO_CONTEXTS
    for typ in result:
        if is_literal(typ):
            # Literals are only valid as long as the operations are
            # correct. Otherwise add a value-free instance.
            cls = compiled.builtin_from_name(evaluator, typ.name.string_name)
            new_result |= cls.execute_evaluated()
        else:
            new_result |= ContextSet(typ)
    return new_result


def _eval_comparison(evaluator, context, left_contexts, operator, right_contexts):
    if not left_contexts or not right_contexts:
        # illegal slices e.g. cause left/right_result to be None
        result = (left_contexts or NO_CONTEXTS) | (right_contexts or NO_CONTEXTS)
        return _literals_to_types(evaluator, result)
    else:
        # I don't think there's a reasonable chance that a string
        # operation is still correct, once we pass something like six
        # objects.
        if len(left_contexts) * len(right_contexts) > 6:
            return _literals_to_types(evaluator, left_contexts | right_contexts)
        else:
            return ContextSet.from_sets(
                _eval_comparison_part(evaluator, context, left, operator, right)
                for left in left_contexts
                for right in right_contexts
            )


def _is_tuple(context):
    return isinstance(context, iterable.Sequence) and context.array_type == 'tuple'


def _is_list(context):
    return isinstance(context, iterable.Sequence) and context.array_type == 'list'


def _bool_to_context(evaluator, bool_):
    return compiled.builtin_from_name(evaluator, force_unicode(str(bool_)))


def _eval_comparison_part(evaluator, context, left, operator, right):
    l_is_num = is_number(left)
    r_is_num = is_number(right)
    if isinstance(operator, unicode):
        str_operator = operator
    else:
        str_operator = force_unicode(str(operator.value))

    if str_operator == '*':
        # for iterables, ignore * operations
        if isinstance(left, iterable.Sequence) or is_string(left):
            return ContextSet(left)
        elif isinstance(right, iterable.Sequence) or is_string(right):
            return ContextSet(right)
    elif str_operator == '+':
        if l_is_num and r_is_num or is_string(left) and is_string(right):
            return ContextSet(left.execute_operation(right, str_operator))
        elif _is_tuple(left) and _is_tuple(right) or _is_list(left) and _is_list(right):
            return ContextSet(iterable.MergedArray(evaluator, (left, right)))
    elif str_operator == '-':
        if l_is_num and r_is_num:
            return ContextSet(left.execute_operation(right, str_operator))
    elif str_operator == '%':
        # With strings and numbers the left type typically remains. Except for
        # `int() % float()`.
        return ContextSet(left)
    elif str_operator in COMPARISON_OPERATORS:
        if is_compiled(left) and is_compiled(right):
            # Possible, because the return is not an option. Just compare.
            try:
                return ContextSet(left.execute_operation(right, str_operator))
            except TypeError:
                # Could be True or False.
                pass
        else:
            if str_operator in ('is', '!=', '==', 'is not'):
                operation = COMPARISON_OPERATORS[str_operator]
                bool_ = operation(left, right)
                return ContextSet(_bool_to_context(evaluator, bool_))

        return ContextSet(_bool_to_context(evaluator, True), _bool_to_context(evaluator, False))
    elif str_operator == 'in':
        return NO_CONTEXTS

    def check(obj):
        """Checks if a Jedi object is either a float or an int."""
        return isinstance(obj, CompiledInstance) and \
            obj.name.string_name in ('int', 'float')

    # Static analysis, one is a number, the other one is not.
    if str_operator in ('+', '-') and l_is_num != r_is_num \
            and not (check(left) or check(right)):
        message = "TypeError: unsupported operand type(s) for +: %s and %s"
        analysis.add(context, 'type-error-operation', operator,
                     message % (left, right))

    return ContextSet(left, right)


def _remove_statements(evaluator, context, stmt, name):
    """
    This is the part where statements are being stripped.

    Due to lazy evaluation, statements like a = func; b = a; b() have to be
    evaluated.
    """
    pep0484_contexts = \
        pep0484.find_type_from_comment_hint_assign(context, stmt, name)
    if pep0484_contexts:
        return pep0484_contexts

    return eval_expr_stmt(context, stmt, seek_name=name)


def tree_name_to_contexts(evaluator, context, tree_name):

    context_set = ContextSet()
    module_node = context.get_root_context().tree_node
    if module_node is not None:
        names = module_node.get_used_names().get(tree_name.value, [])
        for name in names:
            expr_stmt = name.parent

            correct_scope = parser_utils.get_parent_scope(name) == context.tree_node

            if expr_stmt.type == "expr_stmt" and expr_stmt.children[1].type == "annassign" and correct_scope:
                context_set |= _evaluate_for_annotation(context, expr_stmt.children[1].children[1])

    if context_set:
        return context_set

    types = []
    node = tree_name.get_definition(import_name_always=True)
    if node is None:
        node = tree_name.parent
        if node.type == 'global_stmt':
            context = evaluator.create_context(context, tree_name)
            finder = NameFinder(evaluator, context, context, tree_name.value)
            filters = finder.get_filters(search_global=True)
            # For global_stmt lookups, we only need the first possible scope,
            # which means the function itself.
            filters = [next(filters)]
            return finder.find(filters, attribute_lookup=False)
        elif node.type not in ('import_from', 'import_name'):
            raise ValueError("Should not happen. type: %s", node.type)

    typ = node.type
    if typ == 'for_stmt':
        types = pep0484.find_type_from_comment_hint_for(context, node, tree_name)
        if types:
            return types
    if typ == 'with_stmt':
        types = pep0484.find_type_from_comment_hint_with(context, node, tree_name)
        if types:
            return types

    if typ in ('for_stmt', 'comp_for'):
        try:
            types = context.predefined_names[node][tree_name.value]
        except KeyError:
            cn = ContextualizedNode(context, node.children[3])
            for_types = iterate_contexts(
                cn.infer(),
                contextualized_node=cn,
                is_async=node.parent.type == 'async_stmt',
            )
            c_node = ContextualizedName(context, tree_name)
            types = check_tuple_assignments(evaluator, c_node, for_types)
    elif typ == 'expr_stmt':
        types = _remove_statements(evaluator, context, node, tree_name)
    elif typ == 'with_stmt':
        context_managers = context.eval_node(node.get_test_node_from_name(tree_name))
        enter_methods = context_managers.py__getattribute__(u'__enter__')
        return enter_methods.execute_evaluated()
    elif typ in ('import_from', 'import_name'):
        types = imports.infer_import(context, tree_name)
    elif typ in ('funcdef', 'classdef'):
        types = _apply_decorators(context, node)
    elif typ == 'try_stmt':
        # TODO an exception can also be a tuple. Check for those.
        # TODO check for types that are not classes and add it to
        # the static analysis report.
        exceptions = context.eval_node(tree_name.get_previous_sibling().get_previous_sibling())
        types = exceptions.execute_evaluated()
    else:
        raise ValueError("Should not happen. type: %s" % typ)
    return types


def _apply_decorators(context, node):
    """
    Returns the function, that should to be executed in the end.
    This is also the places where the decorators are processed.
    """
    if node.type == 'classdef':
        decoratee_context = ClassContext(
            context.evaluator,
            parent_context=context,
            classdef=node
        )
    else:
        decoratee_context = FunctionContext(
            context.evaluator,
            parent_context=context,
            funcdef=node
        )
    initial = values = ContextSet(decoratee_context)
    for dec in reversed(node.get_decorators()):
        debug.dbg('decorator: %s %s', dec, values)
        dec_values = context.eval_node(dec.children[1])
        trailer_nodes = dec.children[2:-1]
        if trailer_nodes:
            # Create a trailer and evaluate it.
            trailer = tree.PythonNode('trailer', trailer_nodes)
            trailer.parent = dec
            dec_values = eval_trailer(context, dec_values, trailer)

        if not len(dec_values):
            debug.warning('decorator not found: %s on %s', dec, node)
            return initial

        values = dec_values.execute(arguments.ValuesArguments([values]))
        if not len(values):
            debug.warning('not possible to resolve wrappers found %s', node)
            return initial

        debug.dbg('decorator end %s', values)
    return values


def check_tuple_assignments(evaluator, contextualized_name, context_set):
    """
    Checks if tuples are assigned.
    """
    lazy_context = None
    for index, node in contextualized_name.assignment_indexes():
        cn = ContextualizedNode(contextualized_name.context, node)
        iterated = context_set.iterate(cn)
        for _ in range(index + 1):
            try:
                lazy_context = next(iterated)
            except StopIteration:
                # We could do this with the default param in next. But this
                # would allow this loop to run for a very long time if the
                # index number is high. Therefore break if the loop is
                # finished.
                return ContextSet()
        context_set = lazy_context.infer()
    return context_set


def eval_subscript_list(evaluator, context, index):
    """
    Handles slices in subscript nodes.
    """
    if index == ':':
        # Like array[:]
        return ContextSet(iterable.Slice(context, None, None, None))

    elif index.type == 'subscript' and not index.children[0] == '.':
        # subscript basically implies a slice operation, except for Python 2's
        # Ellipsis.
        # e.g. array[:3]
        result = []
        for el in index.children:
            if el == ':':
                if not result:
                    result.append(None)
            elif el.type == 'sliceop':
                if len(el.children) == 2:
                    result.append(el.children[1])
            else:
                result.append(el)
        result += [None] * (3 - len(result))

        return ContextSet(iterable.Slice(context, *result))
    elif index.type == 'subscriptlist':
        return NO_CONTEXTS

    # No slices
    return context.eval_node(index)
