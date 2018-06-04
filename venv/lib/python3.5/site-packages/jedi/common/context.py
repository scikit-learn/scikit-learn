class BaseContext(object):
    def __init__(self, evaluator, parent_context=None):
        self.evaluator = evaluator
        self.parent_context = parent_context

    def get_root_context(self):
        context = self
        while True:
            if context.parent_context is None:
                return context
            context = context.parent_context


class BaseContextSet(object):
    def __init__(self, *args):
        self._set = set(args)

    @classmethod
    def from_iterable(cls, iterable):
        return cls.from_set(set(iterable))

    @classmethod
    def from_set(cls, set_):
        self = cls()
        self._set = set_
        return self

    @classmethod
    def from_sets(cls, sets):
        """
        Used to work with an iterable of set.
        """
        aggregated = set()
        sets = list(sets)
        for set_ in sets:
            if isinstance(set_, BaseContextSet):
                aggregated |= set_._set
            else:
                aggregated |= set_
        return cls.from_set(aggregated)

    def __or__(self, other):
        return type(self).from_set(self._set | other._set)

    def __iter__(self):
        for element in self._set:
            yield element

    def __bool__(self):
        return bool(self._set)

    def __len__(self):
        return len(self._set)

    def __repr__(self):
        return '%s(%s)' % (self.__class__.__name__, ', '.join(str(s) for s in self._set))

    def filter(self, filter_func):
        return type(self).from_iterable(filter(filter_func, self._set))

    def __getattr__(self, name):
        def mapper(*args, **kwargs):
            return type(self).from_sets(
                getattr(context, name)(*args, **kwargs)
                for context in self._set
            )
        return mapper
