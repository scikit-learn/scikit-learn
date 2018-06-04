from jedi.evaluate.filters import publish_method, BuiltinOverwrite
from jedi.evaluate.base_context import ContextSet


class AsyncBase(BuiltinOverwrite):
    def __init__(self, evaluator, func_execution_context):
        super(AsyncBase, self).__init__(evaluator)
        self.func_execution_context = func_execution_context

    @property
    def name(self):
        return self.get_object().name

    def __repr__(self):
        return "<%s of %s>" % (type(self).__name__, self.func_execution_context)


class Coroutine(AsyncBase):
    special_object_identifier = u'COROUTINE'

    @publish_method('__await__')
    def _await(self):
        return ContextSet(CoroutineWrapper(self.evaluator, self.func_execution_context))


class CoroutineWrapper(AsyncBase):
    special_object_identifier = u'COROUTINE_WRAPPER'

    def py__stop_iteration_returns(self):
        return self.func_execution_context.get_return_values()


class AsyncGenerator(AsyncBase):
    """Handling of `yield` functions."""
    special_object_identifier = u'ASYNC_GENERATOR'

    def py__aiter__(self):
        return self.func_execution_context.get_yield_lazy_contexts(is_async=True)
