""" LoopFullUnrolling fully unrolls loops with static bounds. """

from pythran import metadata
from pythran.analyses import HasBreak, HasContinue, NodeCount
from pythran.openmp import OMPDirective
from pythran.conversion import to_ast
from pythran.passmanager import Transformation

from copy import deepcopy
import gast as ast


class LoopFullUnrolling(Transformation):
    '''
    Fully unroll loops with static bounds

    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> node = ast.parse('for j in [1,2,3]: i += j')
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(LoopFullUnrolling, node)
    >>> print(pm.dump(backend.Python, node))
    j = 1
    i += j
    j = 2
    i += j
    j = 3
    i += j

    >>> node = ast.parse('for j in (a,b): i += j')
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(LoopFullUnrolling, node)
    >>> print(pm.dump(backend.Python, node))
    j = a
    i += j
    j = b
    i += j

    >>> node = ast.parse('for j in {1}: i += j')
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(LoopFullUnrolling, node)
    >>> print(pm.dump(backend.Python, node))
    j = 1
    i += j

    >>> node = ast.parse('for j in builtins.enumerate("1"): j')
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(LoopFullUnrolling, node)
    >>> print(pm.dump(backend.Python, node))
    j = (0, '1')
    j
    '''

    MAX_NODE_COUNT = 4096

    def visit_For(self, node):
        # if the user added some OpenMP directive, trust him and no unroll
        if metadata.get(node, OMPDirective):
            return node  # don't visit children because of collapse

        # first unroll children if needed or possible
        self.generic_visit(node)

        # a break or continue in the loop prevents unrolling too
        has_break = any(self.gather(HasBreak, n)
                        for n in node.body)
        has_cont = any(self.gather(HasContinue, n)
                       for n in node.body)

        if has_break or has_cont:
            return node

        # do not unroll too much to prevent code growth
        node_count = self.gather(NodeCount, node)

        def unroll(elt, body):
            return [ast.Assign([deepcopy(node.target)], elt, None)] + body

        def dc(body, i, n):
            if i == n - 1:
                return body
            else:
                return deepcopy(body)

        def getrange(n):
            return getattr(getattr(n, 'func', None), 'attr', None)

        if isinstance(node.iter, (ast.Tuple, ast.List)):
            elts_count = len(node.iter.elts)
            total_count = node_count * elts_count
            issmall = total_count < LoopFullUnrolling.MAX_NODE_COUNT
            if issmall:
                self.update = True
                return sum([unroll(elt, dc(node.body, i, elts_count))
                            for i, elt in enumerate(node.iter.elts)], [])
        ast.fix_missing_locations(node.iter)
        code = compile(ast.gast_to_ast(ast.Expression(node.iter)),
                       '<loop unrolling>', 'eval')
        try:
            values = list(eval(code, {'builtins': __import__('builtins')}))
        except Exception:
            return node

        values_count = len(values)
        total_count = node_count * values_count
        issmall = total_count < LoopFullUnrolling.MAX_NODE_COUNT
        if issmall:
            try:
                new_node = sum([unroll(to_ast(elt),
                                       dc(node.body, i, values_count))
                                for i, elt in enumerate(values)], [])
                self.update = True
                return new_node
            except Exception:
                return node
        return node
