from collections import namedtuple
from contextlib import closing
from io import StringIO

from .._config import config_context
from ..base import BaseEstimator
from ..pipeline import Pipeline
from ..pipeline import FeatureUnion
from ..compose import ColumnTransformer
from ..ensemble import VotingClassifier, VotingRegressor


def _estimator_tool_tip(estimator):
    """Replace newlines to allow for css content: attr(...) to properly
    display tooltips.
    """
    return str(estimator).replace('\n', '&#xa;')


def _write_label_html(out, name, tool_tip):
    """Write label to html"""
    out.write('<div class="sk-label" sk-data-tooltip="{}">'
              '{}</div>'.format(tool_tip, name))


_EstHTMLInfo = namedtuple('_EstHTMLInfo',
                          'type, estimators, names, name_tips')


def _type_of_html_estimator(estimator):
    """Generate information about how to display an estimator.
    """
    if isinstance(estimator, str):
        return _EstHTMLInfo('single', estimator, estimator, estimator)

    elif estimator is None:
        return _EstHTMLInfo('single', estimator, 'None', 'None')

    elif isinstance(estimator, Pipeline):
        estimators = [step[1] for step in estimator.steps]
        names = [step[0] for step in estimator.steps]
        name_tips = [_estimator_tool_tip(est) for est in estimators]
        return _EstHTMLInfo('serial', estimators, names, name_tips)

    elif isinstance(estimator, ColumnTransformer):
        estimators = [trans[1] for trans in estimator.transformers]
        names = [trans[0] for trans in estimator.transformers]
        name_tips = [trans[2] for trans in estimator.transformers]
        return _EstHTMLInfo('parallel', estimators, names, name_tips)

    elif isinstance(estimator, FeatureUnion):
        estimators = [trans[1] for trans in estimator.transformer_list]
        names = [trans[0] for trans in estimator.transformer_list]
        name_tips = [_estimator_tool_tip(est) for est in estimators]
        return _EstHTMLInfo('parallel', estimators, names, name_tips)

    elif isinstance(estimator, (VotingClassifier, VotingRegressor)):
        estimators = [est[1] for est in estimator.estimators]
        names = [est[0] for est in estimator.estimators]
        name_tips = [_estimator_tool_tip(est) for est in estimators]
        return _EstHTMLInfo('parallel', estimators, names, name_tips)

    elif isinstance(estimator, BaseEstimator):
        name = estimator.__class__.__name__
        tool_tip = _estimator_tool_tip(estimator)
        return _EstHTMLInfo('single', estimator, name, tool_tip)

    else:
        raise ValueError("Invalid estimator")


def _write_estimator_html(out, estimator, name):
    """Write estimator to html in serial, parallel, or by itself (single).
    """
    est_html_info = _type_of_html_estimator(estimator)

    if est_html_info.type == 'serial':
        out.write('<div class="sk-serial">')
        est_infos = zip(est_html_info.estimators, est_html_info.names,
                        est_html_info.name_tips)
        for est, name, tool_tip in est_infos:
            _write_estimator_html(out, est, name)
        out.write('</div>')  # sk-serial

    elif est_html_info.type == 'parallel':
        out.write('<div class="sk-serial-item sk-dashed-wrapped">')
        if name:
            tool_tip = _estimator_tool_tip(estimator)
            _write_label_html(out, name, tool_tip)
        out.write('<div class="sk-parallel">')

        est_infos = zip(est_html_info.estimators, est_html_info.names,
                        est_html_info.name_tips)
        for est, name, tool_tip in est_infos:
            out.write('<div class="sk-parallel-item">')
            _write_label_html(out, name, tool_tip)
            out.write('<div class="sk-serial">')
            _write_estimator_html(out, est, name)
            out.write('</div></div>')  # sk-parallel-item sk-serial
        out.write('</div></div>')  # sk-parallel sk-serial-item

    elif est_html_info.type == 'single':
        out.write('<div class="sk-serial-item">'
                  '<div class="sk-estimator" sk-data-tooltip="{}">'
                  '{}</div></div>'.format(est_html_info.name_tips,
                                          est_html_info.names))


_STYLE = """
.sk-estimator {
  font-family: monospace;
  background-color: #f0f8ff;
  padding: 0.5em;
  margin: 0.25em 0.25em;
  border: 1px dotted black;
  text-align: center;
}
.sk-parallel-item::after {
  content: "";
  width: 100%;
  border-bottom: 1px solid gray;
  flex-grow: 1;
}
.sk-serial::before {
  content: "";
  position: absolute;
  border-left: 1px solid gray;
  top: 2em;
  bottom: 0;
  left: 50%;
}
.sk-serial {
  display: flex;
  flex-direction: column;
  align-items: center;
  float: left;
  background: white;
}
.sk-parallel {
  display: flex;
  align-items: stretch;
}
.sk-parallel-item {
  display: flex;
  flex-direction: column;
  position: relative;
}
.sk-parallel-item:first-child::after {
  align-self: flex-end;
  width: 50%;
}
.sk-parallel-item:last-child::after {
  align-self: flex-start;
  width: 50%;
}
.sk-final-spacer {
  visibility: hidden;
  font-family: monospace;
  white-space: pre;
}
.sk-dashed-wrapped {
  border: 1px dashed gray;
  padding: 0.25em;
}
.sk-label {
  text-align: center;
  font-family: monospace;
  font-weight: bold;
  margin: 0;
  background: white;
}
.sk-serial-item {
  margin-bottom: 0.25em;
}
.sk-container {
  display: flex;
  flex-direction: column;
  align-items: flex-start;
  position: relative;
  float: left;
}
[sk-data-tooltip] {
  position: relative;
  cursor: pointer;
}
[sk-data-tooltip]:before {
  visibility: hidden;
  opacity: 0;
  pointer-events: none;
  font-weight: 400;
}
[sk-data-tooltip]:before {
  position: absolute;
  top: 0;
  left: 0;
  padding: 0.5em;
  overflow: hidden;
  background-color: #f0f8ff;
  border: 1px solid gray;
  white-space: pre;
  content: attr(sk-data-tooltip);
  text-align: left;
}
[sk-data-tooltip]:hover:before {
  visibility: visible;
  opacity: 1;
  z-index: 2;
}
"""


def export_html(estimator, print_changed_only=True):
    """Build a HTML representation of an estimator

    Parameters
    ----------
    estimator : estimator object
        The estimator to visualize.

    print_changed_only : bool, optional (default=True)
        If True, only the parameters that were set to non-default
        values will be printed when printing an estimator.

    Returns
    -------
    html: str or iPython HTML object
        HTML representation of estimator. When called in jupyter notebook or
        lab, a iPython HTML object is returned.
    """

    with config_context(print_changed_only=print_changed_only), \
            closing(StringIO()) as out:

        if not isinstance(estimator, Pipeline):
            estimator = Pipeline([('', estimator)])

        out.write('<html><head><style>')
        out.write(_STYLE.replace('\n', ''))
        out.write('</style></head><body>')

        out.write('<div class="sk-container">')
        _write_estimator_html(out, estimator, '')
        out.write('</div>')  # sk-container

        # Adds whitespace at the end to allow space for hover info
        out.write('<div class="sk-final-spacer">')
        out.write(_estimator_tool_tip(estimator.steps[-1]))
        out.write('</div>')  # sk-final-spacer
        out.write("</body></html>")

        html_output = out.getvalue()

        # wrap in iPython HTML if in a notebook context
        try:
            cls_name = get_ipython().__class__.__name__
            if cls_name != 'ZMQInteractiveShell':
                return html_output

            from IPython.display import HTML
            return HTML(html_output)
        except (ImportError, NameError):
            return html_output
