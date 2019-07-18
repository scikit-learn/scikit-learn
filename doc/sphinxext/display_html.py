"""
Primary used to display the html output  `sklearn.inspection.display_estimator`
in sphinx.
"""
import sys
from docutils.parsers.rst import Directive
from docutils import nodes
from io import StringIO


class ExecuteHTML(Directive):

    has_content = True
    required_arguments = 0
    optional_arguments = 0

    @classmethod
    def execute(cls, code):
        orig_stdout, orig_stderr = sys.stdout, sys.stderr

        output, err = StringIO(), StringIO()

        sys.stdout, sys.stderr = output, err
        exec(code)
        sys.stdout, sys.stderr = orig_stdout, orig_stderr

        return "".join(['<div style="font-size: 1.2em">',
                        output.getvalue(), err.getvalue(), "</div>"])

    def run(self):
        output = []
        code = "\n".join(self.content)
        code_results = self.execute(code)

        input_code = nodes.literal_block(code, code)
        input_code['language'] = 'python'
        output.append(input_code)
        code_results = nodes.raw('', code_results, format='html')
        output.append(code_results)
        return output


def setup(app):
    app.add_directive('display_html', ExecuteHTML)
