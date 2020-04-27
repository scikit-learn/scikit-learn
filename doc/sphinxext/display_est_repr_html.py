"""
Primarily used to display the html output of `_repr_html_` of estimators
"""
from sphinx.util.docutils import SphinxDirective
from contextlib import redirect_stderr, redirect_stdout
from docutils import nodes
from io import StringIO


class DisplayEstimatorRepr(SphinxDirective):
    """Execute Python and runs `_repr_html_` on the last element on the code
    block. The last element in the code block should be an estimator with
    support for `_repr_html_`.
    """

    has_content = True
    required_arguments = 0
    optional_arguments = 0

    def execute(self, code, format_str):
        code_parts = code.split('\n')
        final_output = code_parts[-1]
        final_est = final_output.lstrip(' ')
        n_whitespace = len(final_output) - len(final_est)
        code_parts[-1] = " " * n_whitespace + format_str.format(final_est)
        code = '\n'.join(code_parts)

        output, err = StringIO(), StringIO()
        with redirect_stdout(output), redirect_stderr(err):
            exec(code)

        return f"{output.getvalue()}{err.getvalue()}"

    def run(self):
        output = []
        code = "\n".join(self.content)
        repr_html = self.execute(
            code, format_str='print({}._repr_mimebundle_()["text/html"])')

        input_code = nodes.literal_block(rawsource=code, text=code)
        input_code['language'] = 'python'
        output.append(input_code)

        repr_html = f"<p>{repr_html}</p>"
        html_node = nodes.raw('', repr_html, format='html')
        output.append(html_node)

        if self.env.app.builder.name == 'latex':
            code_results_latex = r"""
            \begin{sphinxadmonition}{note}{Note:}
            The HTML output of this code snippet can only been seen on the HTML
            version of the documentation. The following is a string
            representation.
            \end{sphinxadmonition}
            """
            latex_node = nodes.raw('', code_results_latex, format='latex')
            output.append(latex_node)

            str_repr = self.execute(code, format_str='print(repr({}))')
            str_repr_node = nodes.literal_block(rawsource=str_repr,
                                                text=str_repr)
            str_repr_node['language'] = 'python'
            output.append(str_repr_node)

        return output


def setup(app):
    app.add_directive('display_estimator_repr_html', DisplayEstimatorRepr)
    return {'parallel_read_safe': True, 'parallel_write_safe': True}
