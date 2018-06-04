from jedi.evaluate.base_context import ContextSet, NO_CONTEXTS

class AbstractLazyContext(object):
    def __init__(self, data):
        self.data = data

    def __repr__(self):
        return '<%s: %s>' % (self.__class__.__name__, self.data)

    def infer(self):
        raise NotImplementedError


class LazyKnownContext(AbstractLazyContext):
    """data is a context."""
    def infer(self):
        return ContextSet(self.data)


class LazyKnownContexts(AbstractLazyContext):
    """data is a ContextSet."""
    def infer(self):
        return self.data


class LazyUnknownContext(AbstractLazyContext):
    def __init__(self):
        super(LazyUnknownContext, self).__init__(None)

    def infer(self):
        return NO_CONTEXTS


class LazyTreeContext(AbstractLazyContext):
    def __init__(self, context, node):
        super(LazyTreeContext, self).__init__(node)
        self._context = context
        # We need to save the predefined names. It's an unfortunate side effect
        # that needs to be tracked otherwise results will be wrong.
        self._predefined_names = dict(context.predefined_names)

    def infer(self):
        old, self._context.predefined_names = \
            self._context.predefined_names, self._predefined_names
        try:
            return self._context.eval_node(self.data)
        finally:
            self._context.predefined_names = old


def get_merged_lazy_context(lazy_contexts):
    if len(lazy_contexts) > 1:
        return MergedLazyContexts(lazy_contexts)
    else:
        return lazy_contexts[0]


class MergedLazyContexts(AbstractLazyContext):
    """data is a list of lazy contexts."""
    def infer(self):
        return ContextSet.from_sets(l.infer() for l in self.data)
