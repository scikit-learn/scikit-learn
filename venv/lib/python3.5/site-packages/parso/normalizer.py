from contextlib import contextmanager

from parso._compatibility import use_metaclass


class _NormalizerMeta(type):
    def __new__(cls, name, bases, dct):
        new_cls = type.__new__(cls, name, bases, dct)
        new_cls.rule_value_classes = {}
        new_cls.rule_type_classes = {}
        return new_cls


class Normalizer(use_metaclass(_NormalizerMeta)):
    def __init__(self, grammar, config):
        self.grammar = grammar
        self._config = config
        self.issues = []

        self._rule_type_instances = self._instantiate_rules('rule_type_classes')
        self._rule_value_instances = self._instantiate_rules('rule_value_classes')

    def _instantiate_rules(self, attr):
        dct = {}
        for base in type(self).mro():
            rules_map = getattr(base, attr, {})
            for type_, rule_classes in rules_map.items():
                new = [rule_cls(self) for rule_cls in rule_classes]
                dct.setdefault(type_, []).extend(new)
        return dct

    def walk(self, node):
        self.initialize(node)
        value = self.visit(node)
        self.finalize()
        return value

    def visit(self, node):
        try:
            children = node.children
        except AttributeError:
            return self.visit_leaf(node)
        else:
           with self.visit_node(node):
               return ''.join(self.visit(child) for child in children)

    @contextmanager
    def visit_node(self, node):
        self._check_type_rules(node)
        yield

    def _check_type_rules(self, node):
        for rule in self._rule_type_instances.get(node.type, []):
            rule.feed_node(node)

    def visit_leaf(self, leaf):
        self._check_type_rules(leaf)

        for rule in self._rule_value_instances.get(leaf.value, []):
            rule.feed_node(leaf)

        return leaf.prefix + leaf.value

    def initialize(self, node):
        pass

    def finalize(self):
        pass

    def add_issue(self, node, code, message):
        issue = Issue(node, code, message)
        if issue not in self.issues:
            self.issues.append(issue)
        return True

    @classmethod
    def register_rule(cls, **kwargs):
        """
        Use it as a class decorator::

            normalizer = Normalizer('grammar', 'config')
            @normalizer.register_rule(value='foo')
            class MyRule(Rule):
                error_code = 42
        """
        return cls._register_rule(**kwargs)

    @classmethod
    def _register_rule(cls, value=None, values=(), type=None, types=()):
        values = list(values)
        types = list(types)
        if value is not None:
            values.append(value)
        if type is not None:
            types.append(type)

        if not values and not types:
            raise ValueError("You must register at least something.")

        def decorator(rule_cls):
            for v in values:
                cls.rule_value_classes.setdefault(v, []).append(rule_cls)
            for t in types:
                cls.rule_type_classes.setdefault(t, []).append(rule_cls)
            return rule_cls

        return decorator


class NormalizerConfig(object):
    normalizer_class = Normalizer

    def create_normalizer(self, grammar):
        if self.normalizer_class is None:
            return None

        return self.normalizer_class(grammar, self)


class Issue(object):
    def __init__(self, node, code, message):
        self._node = node
        self.code = code
        """
        An integer code that stands for the type of error.
        """
        self.message = message
        """
        A message (string) for the issue.
        """
        self.start_pos = node.start_pos
        """
        The start position position of the error as a tuple (line, column). As
        always in |parso| the first line is 1 and the first column 0.
        """

    def __eq__(self, other):
        return self.start_pos == other.start_pos and self.code == other.code

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.code, self.start_pos))

    def __repr__(self):
        return '<%s: %s>' % (self.__class__.__name__, self.code)



class Rule(object):
    code = None
    message = None

    def __init__(self, normalizer):
        self._normalizer = normalizer

    def is_issue(self, node):
        raise NotImplementedError()

    def get_node(self, node):
        return node

    def _get_message(self, message):
        if message is None:
            message = self.message
            if message is None:
                raise ValueError("The message on the class is not set.")
        return message

    def add_issue(self, node, code=None, message=None):
        if code is None:
            code = self.code
            if code is None:
                raise ValueError("The error code on the class is not set.")

        message = self._get_message(message)

        self._normalizer.add_issue(node, code, message)

    def feed_node(self, node):
        if self.is_issue(node):
            issue_node = self.get_node(node)
            self.add_issue(issue_node)
