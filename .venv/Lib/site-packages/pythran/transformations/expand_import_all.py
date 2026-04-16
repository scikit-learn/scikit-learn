""" ExpandImportAll replaces import * by all their modules. """

from pythran.passmanager import Transformation
from pythran.tables import MODULES

import gast as ast


class ExpandImportAll(Transformation):
    '''
    Expands all import when '*' detected

    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> node = ast.parse("from math import *")
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(ExpandImportAll, node)
    >>> print(pm.dump(backend.Python, node))
    from math import acos, acosh, asin, asinh, atan, atan2, atanh, ceil, \
copysign, cos, cosh, degrees, e, erf, erfc, exp, expm1, fabs, factorial, \
floor, fmod, frexp, gamma, hypot, isinf, isnan, ldexp, lgamma, log, log10, \
log1p, modf, pi, pow, radians, sin, sinh, sqrt, tan, tanh, trunc
    '''

    def visit_ImportFrom(self, node):
        for alias in node.names:
            if alias.name == '*':
                self.update = True
                node.names.pop()
                node.names.extend(ast.alias(fname, None)
                                  for fname in sorted(MODULES[node.module]))
        return node
