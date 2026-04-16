""" StaticExpressions gathers constant expression that involve types.  """

from pythran.passmanager import NodeAnalysis


class HasStaticExpression(NodeAnalysis):
    def __init__(self):
        self.result = False

    def visit_Attribute(self, node):
        self.generic_visit(node)
        self.result |= node.attr == 'is_none'


class StaticExpressions(NodeAnalysis):

    """Identify constant expressions."""

    ResultType = set
    def __init__(self):
        super().__init__()
        self.constant_expressions = set()

    def add(self, node):
        self.result.add(node)
        return True

    def not_add(self, _):
        return False

    def match_all(self, *args):
        assert len(args) > 1, "at least two arguments"
        static = False
        const = True
        for value in args:
            if self.visit(value):
                static = True
            else:
                const &= value in self.constant_expressions
        return static and const

    def visit_BoolOp(self, node):
        return self.match_all(*node.values) and self.add(node)

    def visit_BinOp(self, node):
        return self.match_all(node.left, node.right) and self.add(node)

    def visit_UnaryOp(self, node):
        return self.visit(node.operand) and self.add(node)

    def visit_IfExp(self, node):
        return (self.match_all(node.test, node.body, node.orelse)
                and self.add(node))

    def visit_Compare(self, node):
        return self.match_all(node.left, *node.comparators) and self.add(node)

    def visit_Call(self, node):
        return self.visit(node.func)and self.add(node)  # very limited

    def visit_Attribute(self, node):
        return node.attr in ('is_none', 'isinstance')

    def visit_Constant(self, node):
        self.constant_expressions.add(node)

    visit_Subscript = not_add
    visit_Name = not_add
    visit_Dict = not_add
    visit_List = not_add
    visit_Tuple = not_add
    visit_Set = not_add
    visit_Slice = not_add
    visit_Index = not_add
