import ast
import gast


def _generate_translators(to):

    class Translator(ast.NodeTransformer):

        def _visit(self, node):
            if isinstance(node, ast.AST):
                return self.visit(node)
            elif isinstance(node, list):
                return [self._visit(n) for n in node]
            else:
                return node

        def generic_visit(self, node):
            cls = type(node).__name__
            try:
                new_node = getattr(to, cls)()
            except AttributeError:
                # handle nodes that are not part of the AST
                return

            for field in node._fields:
                setattr(new_node, field, self._visit(getattr(node, field)))

            for attr in node._attributes:
                try:
                    setattr(new_node, attr, getattr(node, attr))
                except AttributeError:
                    pass
            return new_node

    return Translator


AstToGAst = _generate_translators(gast)

GAstToAst = _generate_translators(ast)
