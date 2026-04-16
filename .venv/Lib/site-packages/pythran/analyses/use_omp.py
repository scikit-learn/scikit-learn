"""
UseOMP detects if a function use OpenMP
"""

from pythran.passmanager import FunctionAnalysis


class UseOMP(FunctionAnalysis):
    """Detects if a function use openMP"""
    ResultType = bool

    def visit_OMPDirective(self, _):
        self.result = True
