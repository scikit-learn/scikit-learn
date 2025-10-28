# Nodes for structural pattern matching.
#
# In a separate file because they're unlikely to be useful for much else.

from .Nodes import Node, StatNode, ErrorNode
from .Errors import error


class MatchNode(StatNode):
    """
    subject  ExprNode    The expression to be matched
    cases    [MatchCaseNode]  list of cases
    """

    child_attrs = ["subject", "cases"]

    def validate_irrefutable(self):
        found_irrefutable_case = None
        for case in self.cases:
            if isinstance(case, ErrorNode):
                # This validation happens before error nodes have been
                # transformed into actual errors, so we need to ignore them
                continue
            if found_irrefutable_case:
                error(
                    found_irrefutable_case.pos,
                    f"{found_irrefutable_case.pattern.irrefutable_message()} makes remaining patterns unreachable"
                )
                break
            if case.is_irrefutable():
                found_irrefutable_case = case
            case.validate_irrefutable()

    def analyse_expressions(self, env):
        error(self.pos, "Structural pattern match is not yet implemented")
        return self


class MatchCaseNode(Node):
    """
    pattern    PatternNode
    body       StatListNode
    guard      ExprNode or None
    """

    child_attrs = ["pattern", "body", "guard"]

    def is_irrefutable(self):
        if isinstance(self.pattern, ErrorNode):
            return True  # value doesn't really matter
        return self.pattern.is_irrefutable() and not self.guard

    def validate_targets(self):
        if isinstance(self.pattern, ErrorNode):
            return
        self.pattern.get_targets()

    def validate_irrefutable(self):
        if isinstance(self.pattern, ErrorNode):
            return
        self.pattern.validate_irrefutable()


class PatternNode(Node):
    """
    PatternNode is not an expression because
    it does several things (evaluating a boolean expression,
    assignment of targets), and they need to be done at different
    times.

    as_targets   [NameNode]    any target assign by "as"
    """

    child_attrs = ["as_targets"]

    def __init__(self, pos, **kwds):
        if "as_targets" not in kwds:
            kwds["as_targets"] = []
        super(PatternNode, self).__init__(pos, **kwds)

    def is_irrefutable(self):
        return False

    def get_targets(self):
        targets = self.get_main_pattern_targets()
        for target in self.as_targets:
            self.add_target_to_targets(targets, target.name)
        return targets

    def update_targets_with_targets(self, targets, other_targets):
        for name in targets.intersection(other_targets):
            error(self.pos, f"multiple assignments to name '{name}' in pattern")
        targets.update(other_targets)

    def add_target_to_targets(self, targets, target):
        if target in targets:
            error(self.pos, f"multiple assignments to name '{target}' in pattern")
        targets.add(target)

    def get_main_pattern_targets(self):
        # exclude "as" target
        raise NotImplementedError

    def validate_irrefutable(self):
        for attr in self.child_attrs:
            child = getattr(self, attr)
            if child is not None and isinstance(child, PatternNode):
                child.validate_irrefutable()


class MatchValuePatternNode(PatternNode):
    """
    value   ExprNode        # todo be more specific
    is_is_check   bool     Picks "is" or equality check
    """

    child_attrs = PatternNode.child_attrs + ["value"]
    is_is_check = False

    def get_main_pattern_targets(self):
        return set()


class MatchAndAssignPatternNode(PatternNode):
    """
    target   NameNode or None  the target to assign to (None = wildcard)
    is_star  bool
    """

    target = None
    is_star = False

    child_atts = PatternNode.child_attrs + ["target"]

    def is_irrefutable(self):
        return not self.is_star

    def irrefutable_message(self):
        if self.target:
            return "name capture '%s'" % self.target.name
        else:
            return "wildcard"

    def get_main_pattern_targets(self):
        if self.target:
            return {self.target.name}
        else:
            return set()


class OrPatternNode(PatternNode):
    """
    alternatives   list of PatternNodes
    """

    child_attrs = PatternNode.child_attrs + ["alternatives"]

    def get_first_irrefutable(self):
        for alternative in self.alternatives:
            if alternative.is_irrefutable():
                return alternative
        return None

    def is_irrefutable(self):
        return self.get_first_irrefutable() is not None

    def irrefutable_message(self):
        return self.get_first_irrefutable().irrefutable_message()

    def get_main_pattern_targets(self):
        child_targets = None
        for alternative in self.alternatives:
            alternative_targets = alternative.get_targets()
            if child_targets is not None and child_targets != alternative_targets:
                error(self.pos, "alternative patterns bind different names")
            child_targets = alternative_targets
        return child_targets

    def validate_irrefutable(self):
        super(OrPatternNode, self).validate_irrefutable()
        found_irrefutable_case = None
        for alternative in self.alternatives:
            if found_irrefutable_case:
                error(
                    found_irrefutable_case.pos,
                    f"{found_irrefutable_case.irrefutable_message()} makes remaining patterns unreachable"
                )
                break
            if alternative.is_irrefutable():
                found_irrefutable_case = alternative
            alternative.validate_irrefutable()


class MatchSequencePatternNode(PatternNode):
    """
    patterns   list of PatternNodes
    """

    child_attrs = PatternNode.child_attrs + ["patterns"]

    def get_main_pattern_targets(self):
        targets = set()
        for pattern in self.patterns:
            self.update_targets_with_targets(targets, pattern.get_targets())
        return targets


class MatchMappingPatternNode(PatternNode):
    """
    keys   list of NameNodes
    value_patterns  list of PatternNodes of equal length to keys
    double_star_capture_target  NameNode or None
    """

    keys = []
    value_patterns = []
    double_star_capture_target = None

    child_attrs = PatternNode.child_attrs + [
        "keys",
        "value_patterns",
        "double_star_capture_target",
    ]

    def get_main_pattern_targets(self):
        targets = set()
        for pattern in self.value_patterns:
            self.update_targets_with_targets(targets, pattern.get_targets())
        if self.double_star_capture_target:
            self.add_target_to_targets(targets, self.double_star_capture_target.name)
        return targets


class ClassPatternNode(PatternNode):
    """
    class_  NameNode or AttributeNode
    positional_patterns  list of PatternNodes
    keyword_pattern_names    list of NameNodes
    keyword_pattern_patterns    list of PatternNodes
                                (same length as keyword_pattern_names)
    """

    class_ = None
    positional_patterns = []
    keyword_pattern_names = []
    keyword_pattern_patterns = []

    child_attrs = PatternNode.child_attrs + [
        "class_",
        "positional_patterns",
        "keyword_pattern_names",
        "keyword_pattern_patterns",
    ]

    def get_main_pattern_targets(self):
        targets = set()
        for pattern in self.positional_patterns + self.keyword_pattern_patterns:
            self.update_targets_with_targets(targets, pattern.get_targets())
        return targets
