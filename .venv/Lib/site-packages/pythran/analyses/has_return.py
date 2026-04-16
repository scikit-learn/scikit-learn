"""
HasReturn detects if there's a return or yield statement
HasBreak detects if there's a break statement
HasContinue detects if there's a continue statement
"""

from pythran.passmanager import NodeAnalysis


class HasReturn(NodeAnalysis):

    ResultType = bool

    def visit_Return(self, _):
        self.result = True

    def visit_Yield(self, _):
        self.result = True


class HasBreak(NodeAnalysis):

    ResultType = bool

    def visit_For(self, _):
        return

    visit_While = visit_For

    def visit_Break(self, _):
        self.result = True


class HasContinue(NodeAnalysis):

    ResultType = bool

    def visit_For(self, _):
        return

    visit_While = visit_For

    def visit_Continue(self, _):
        self.result = True
