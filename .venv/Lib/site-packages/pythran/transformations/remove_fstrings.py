"""Turns f-strings to format syntax with modulus"""

import gast as ast

from pythran.passmanager import Transformation
from pythran.syntax import PythranSyntaxError


class RemoveFStrings(Transformation, ast.NodeTransformer):
    """Turns f-strings to format syntax with modulus

    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> node = ast.parse("f'a = {1+1:4d}; b = {b:s};'")
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(RemoveFStrings, node)
    >>> print(pm.dump(backend.Python, node))
    ('a = %4d; b = %s;' % ((1 + 1), b))
    """

    def visit_JoinedStr(self, node):
        if len(node.values) == 1 and not isinstance(
            node.values[0], ast.FormattedValue
        ):
            # f-strings with no reference to variable (like `f"bar"`, see #1767)
            return node.values[0]

        if not any(
            isinstance(value, ast.FormattedValue) for value in node.values
        ):
            # nothing to do (not a f-string)
            return node

        base_str = ""
        elements = []
        for value in node.values:
            if isinstance(value, ast.Constant):
                base_str += value.value.replace("%", "%%")
            elif isinstance(value, ast.FormattedValue):
                base_str += "%"
                if value.format_spec is None:
                    raise PythranSyntaxError(
                        "f-strings without format specifier not supported", value
                    )
                base_str += value.format_spec.values[0].value
                elements.append(value.value)
            else:
                raise NotImplementedError

        return ast.BinOp(
            left=ast.Constant(value=base_str, kind=None),
            op=ast.Mod(),
            right=ast.Tuple(elts=elements, ctx=ast.Load()),
        )
