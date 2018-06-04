"""
Contexts are the "values" that Python would return. However Contexts are at the
same time also the "contexts" that a user is currently sitting in.

A ContextSet is typically used to specify the return of a function or any other
static analysis operation. In jedi there are always multiple returns and not
just one.
"""
from parso.python.tree import ExprStmt, CompFor

from jedi import debug
from jedi._compatibility import Python3Method, zip_longest, unicode
from jedi.parser_utils import clean_scope_docstring, get_doc_with_call_signature
from jedi.common import BaseContextSet, BaseContext


class Context(BaseContext):
    """
    Should be defined, otherwise the API returns empty types.
    """

    predefined_names = {}
    tree_node = None
    """
    To be defined by subclasses.
    """

    @property
    def api_type(self):
        # By default just lower name of the class. Can and should be
        # overwritten.
        return self.__class__.__name__.lower()

    @debug.increase_indent
    def execute(self, arguments):
        """
        In contrast to py__call__ this function is always available.

        `hasattr(x, py__call__)` can also be checked to see if a context is
        executable.
        """
        if self.evaluator.is_analysis:
            arguments.eval_all()

        debug.dbg('execute: %s %s', self, arguments)
        from jedi.evaluate import stdlib
        try:
            # Some stdlib functions like super(), namedtuple(), etc. have been
            # hard-coded in Jedi to support them.
            return stdlib.execute(self.evaluator, self, arguments)
        except stdlib.NotInStdLib:
            pass

        try:
            func = self.py__call__
        except AttributeError:
            debug.warning("no execution possible %s", self)
            return NO_CONTEXTS
        else:
            context_set = func(arguments)
            debug.dbg('execute result: %s in %s', context_set, self)
            return context_set

        return self.evaluator.execute(self, arguments)

    def execute_evaluated(self, *value_list):
        """
        Execute a function with already executed arguments.
        """
        from jedi.evaluate.arguments import ValuesArguments
        arguments = ValuesArguments([ContextSet(value) for value in value_list])
        return self.execute(arguments)

    def iterate(self, contextualized_node=None, is_async=False):
        debug.dbg('iterate %s', self)
        try:
            if is_async:
                iter_method = self.py__aiter__
            else:
                iter_method = self.py__iter__
        except AttributeError:
            if contextualized_node is not None:
                from jedi.evaluate import analysis
                analysis.add(
                    contextualized_node.context,
                    'type-error-not-iterable',
                    contextualized_node.node,
                    message="TypeError: '%s' object is not iterable" % self)
            return iter([])
        else:
            return iter_method()

    def get_item(self, index_contexts, contextualized_node):
        from jedi.evaluate.compiled import CompiledObject
        from jedi.evaluate.context.iterable import Slice, Sequence
        result = ContextSet()

        for index in index_contexts:
            if isinstance(index, Slice):
                index = index.obj
            if isinstance(index, CompiledObject):
                try:
                    index = index.get_safe_value()
                except ValueError:
                    pass

            if type(index) not in (float, int, str, unicode, slice, bytes):
                # If the index is not clearly defined, we have to get all the
                # possiblities.
                if isinstance(self, Sequence) and self.array_type == 'dict':
                    result |= self.dict_values()
                else:
                    result |= iterate_contexts(ContextSet(self))
                continue

            # The actual getitem call.
            try:
                getitem = self.py__getitem__
            except AttributeError:
                from jedi.evaluate import analysis
                # TODO this context is probably not right.
                analysis.add(
                    contextualized_node.context,
                    'type-error-not-subscriptable',
                    contextualized_node.node,
                    message="TypeError: '%s' object is not subscriptable" % self
                )
            else:
                try:
                    result |= getitem(index)
                except IndexError:
                    result |= iterate_contexts(ContextSet(self))
                except KeyError:
                    # Must be a dict. Lists don't raise KeyErrors.
                    result |= self.dict_values()
        return result

    def eval_node(self, node):
        return self.evaluator.eval_element(self, node)

    @Python3Method
    def py__getattribute__(self, name_or_str, name_context=None, position=None,
                           search_global=False, is_goto=False,
                           analysis_errors=True):
        """
        :param position: Position of the last statement -> tuple of line, column
        """
        if name_context is None:
            name_context = self
        from jedi.evaluate import finder
        f = finder.NameFinder(self.evaluator, self, name_context, name_or_str,
                              position, analysis_errors=analysis_errors)
        filters = f.get_filters(search_global)
        if is_goto:
            return f.filter_name(filters)
        return f.find(filters, attribute_lookup=not search_global)

    def create_context(self, node, node_is_context=False, node_is_object=False):
        return self.evaluator.create_context(self, node, node_is_context, node_is_object)

    def is_class(self):
        return False

    def py__bool__(self):
        """
        Since Wrapper is a super class for classes, functions and modules,
        the return value will always be true.
        """
        return True

    def py__doc__(self, include_call_signature=False):
        try:
            self.tree_node.get_doc_node
        except AttributeError:
            return ''
        else:
            if include_call_signature:
                return get_doc_with_call_signature(self.tree_node)
            else:
                return clean_scope_docstring(self.tree_node)
        return None


def iterate_contexts(contexts, contextualized_node=None, is_async=False):
    """
    Calls `iterate`, on all contexts but ignores the ordering and just returns
    all contexts that the iterate functions yield.
    """
    return ContextSet.from_sets(
        lazy_context.infer()
        for lazy_context in contexts.iterate(contextualized_node, is_async=is_async)
    )


class TreeContext(Context):
    def __init__(self, evaluator, parent_context=None):
        super(TreeContext, self).__init__(evaluator, parent_context)
        self.predefined_names = {}

    def __repr__(self):
        return '<%s: %s>' % (self.__class__.__name__, self.tree_node)


class ContextualizedNode(object):
    def __init__(self, context, node):
        self.context = context
        self.node = node

    def get_root_context(self):
        return self.context.get_root_context()

    def infer(self):
        return self.context.eval_node(self.node)


class ContextualizedName(ContextualizedNode):
    # TODO merge with TreeNameDefinition?!
    @property
    def name(self):
        return self.node

    def assignment_indexes(self):
        """
        Returns an array of tuple(int, node) of the indexes that are used in
        tuple assignments.

        For example if the name is ``y`` in the following code::

            x, (y, z) = 2, ''

        would result in ``[(1, xyz_node), (0, yz_node)]``.
        """
        indexes = []
        node = self.node.parent
        compare = self.node
        while node is not None:
            if node.type in ('testlist', 'testlist_comp', 'testlist_star_expr', 'exprlist'):
                for i, child in enumerate(node.children):
                    if child == compare:
                        indexes.insert(0, (int(i / 2), node))
                        break
                else:
                    raise LookupError("Couldn't find the assignment.")
            elif isinstance(node, (ExprStmt, CompFor)):
                break

            compare = node
            node = node.parent
        return indexes


class ContextSet(BaseContextSet):
    def py__class__(self):
        return ContextSet.from_iterable(c.py__class__() for c in self._set)

    def iterate(self, contextualized_node=None, is_async=False):
        from jedi.evaluate.lazy_context import get_merged_lazy_context
        type_iters = [c.iterate(contextualized_node, is_async=is_async) for c in self._set]
        for lazy_contexts in zip_longest(*type_iters):
            yield get_merged_lazy_context(
                [l for l in lazy_contexts if l is not None]
            )


NO_CONTEXTS = ContextSet()


def iterator_to_context_set(func):
    def wrapper(*args, **kwargs):
        return ContextSet.from_iterable(func(*args, **kwargs))

    return wrapper
