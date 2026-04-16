""" Inlinable list function that may be inlined. """

from pythran.passmanager import ModuleAnalysis
from pythran.analyses import Identifiers
from pythran.analyses.pure_expressions import PureExpressions
import pythran.metadata as metadata

import gast as ast
import copy


class Inlinable(ModuleAnalysis[PureExpressions]):

    """ Determine set of inlinable function.

    A function can be inlined if it has only one statement and doesn't
    recurse on itself.
    """

    ResultType = dict

    def visit_FunctionDef(self, node):
        """ Determine this function definition can be inlined. """
        non_pass = [i for i, s in enumerate(node.body)
                    if not isinstance(s, ast.Pass)]

        if len(non_pass) != 1:
            return

        sindex = non_pass[0]
        sbody = node.body[sindex]
        if not isinstance(sbody, (ast.Call, ast.Return)):
            return

        # only consider static return if they are pure
        if metadata.get(sbody, metadata.StaticReturn):
            if sbody not in self.pure_expressions:
                return

        ids = self.gather(Identifiers, sbody)
        # FIXME : It marks "not inlinable" def foo(foo): return foo
        if node.name not in ids:
            self.result[node.name] = copy.deepcopy(node)
            self.result[node.name].body = [self.result[node.name].body[sindex]]
