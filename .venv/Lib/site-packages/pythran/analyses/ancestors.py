"""
Ancestors computes the ancestors of each node
"""

from pythran.passmanager import ModuleAnalysis

class Ancestors(ModuleAnalysis):
    '''
    Associate each node with the list of its ancestors

    Based on the tree view of the AST: each node has the Module as parent.
    The result of this analysis is a dictionary with nodes as key,
    and list of nodes as values.
    '''

    ResultType = dict

    def __init__(self):
        super().__init__()
        self.current = tuple()

    def generic_visit(self, node):
        self.result[node] = current = self.current
        self.current += node,
        super(Ancestors, self).generic_visit(node)
        self.current = current

    visit = generic_visit


class AncestorsWithBody(Ancestors):

    # Overload the visit method set from Ancestors
    visit = ModuleAnalysis.visit

    def visit_metadata(self, node):
        if hasattr(node, 'metadata'):
            self.generic_visit(node.metadata)

    def visit_body(self, body):
        body_as_tuple = tuple(body)
        self.result[body_as_tuple] = current = self.current
        self.current += body_as_tuple,
        for stmt in body:
            self.generic_visit(stmt)
        self.current = current

    def visit_If(self, node):
        self.result[node] = current = self.current
        self.current += node,
        self.generic_visit(node.test)
        self.visit_metadata(node)
        self.visit_body(node.body)
        self.visit_body(node.orelse)
        self.current = current

    def visit_While(self, node):
        self.result[node] = current = self.current
        self.current += node,
        self.generic_visit(node.test)
        self.visit_metadata(node)
        self.visit_body(node.body)
        self.visit_body(node.orelse)
        self.current = current

    def visit_For(self, node):
        self.result[node] = current = self.current
        self.current += node,
        self.generic_visit(node.target)
        self.generic_visit(node.iter)
        self.visit_metadata(node)
        self.visit_body(node.body)
        self.visit_body(node.orelse)
        self.current = current

    def visit_Try(self, node):
        self.result[node] = current = self.current
        self.current += node,
        self.visit_metadata(node)
        self.visit_body(node.body)
        for handler in node.handlers:
            self.generic_visit(handler)
        self.visit_body(node.orelse)
        self.visit_body(node.finalbody)
        self.current = current
