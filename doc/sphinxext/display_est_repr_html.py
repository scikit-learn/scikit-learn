"""
Primary used to display the html output of `_repr_html_` of estimators
"""
import sys
from docutils.parsers.rst import Directive
from docutils import nodes
from io import StringIO


class DisplayEstimatorRepr(Directive):
    "Execute Python code and includes stdout as HTML"

    has_content = True
    required_arguments = 0
    optional_arguments = 0

    def execute(self, code):
        code_parts = code.split('\n')
        final_output = code_parts[-1]
        code_parts[-1] = f'print({final_output}._repr_html_())'
        code = '\n'.join(code_parts)
        orig_stdout, orig_stderr = sys.stdout, sys.stderr

        output, err = StringIO(), StringIO()

        sys.stdout, sys.stderr = output, err
        exec(code)
        sys.stdout, sys.stderr = orig_stdout, orig_stderr

        return f"{output.getvalue()}{err.getvalue()}"

    def run(self):
        output = []
        code = "\n".join(self.content)
        code_results = self.execute(code)

        input_code = nodes.literal_block(code, code)
        input_code['language'] = 'python'
        output.append(input_code)

        html_node = nodes.raw('', code_results, format='html')
        output.append(html_node)

        code_results_latex = r"""
        \begin{sphinxadmonition}{note}{Note:}
        The HTML output of this code snippet can only been seen on the HTML
        version of the documentation
        \end{sphinxadmonition}
        """
        latex_node = nodes.raw('', code_results_latex, format='latex')
        output.append(latex_node)
        return output


def setup(app):
    app.add_directive('display_estimator_repr_html', DisplayEstimatorRepr)

    return {
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
