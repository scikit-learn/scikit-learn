"""
Primary used to display the html output of `_repr_html_` of estimators
"""
import sys
from docutils.parsers.rst import Directive
from docutils import nodes
from io import StringIO
from sphinx import addnodes


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

        return "".join(['<div style="font-size: 1.1em;display: inline-block">',
                        output.getvalue(), err.getvalue(), "</div>"])

    def run(self):
        output = []
        code = "\n".join(self.content)
        code_results = self.execute(code)

        input_code = nodes.literal_block(code, code)
        input_code['language'] = 'python'
        output.append(input_code)

        onlynode_html = addnodes.only(expr='html')
        onlynode_html += nodes.raw('', code_results, format='html')
        output.append(onlynode_html)

        onlynode_latex = addnodes.only(expr='latex')
        onlynode_latex += nodes.raw('', code_results, format='html')
        onlynode_latex += nodes.note('The HTML output of this code snippet '
                                     'can only been seen on the HTML version '
                                     'of the docs.')
        output.append(onlynode_latex)

        return output


def setup(app):
    app.add_directive('display_estimator_repr_html', DisplayEstimatorRepr)
