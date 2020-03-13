from collections import namedtuple
from contextlib import closing
from io import StringIO
import uuid


def _estimator_details(estimator):
    """Replace newlines to allow for css content: attr(...) to properly
    display estimator details.
    """
    return str(estimator).replace('\n', '&#xa;')


def _write_dropdown_html(out, name, tool_tip, outer_class, inner_class):
    out.write(
        f'<div class="{outer_class}">'
        f'<div class="{inner_class} sk-toggleable">')

    if tool_tip is not None:
        est_id = uuid.uuid4()
        out.write(f'<input class="sk-toggleable__control sk-hidden--visually" '
                  f'id="{est_id}" type="checkbox"/>'
                  f'<label class="sk-toggleable__label" for="{est_id}">'
                  f'{name}</label>'
                  f'<div class="sk-toggleable__content"><pre>{tool_tip}'
                  f'</pre></div>')
    else:
        out.write(f'<label>{name}</label>')
    out.write('</div></div>')  # outer_class inner_class


def _write_label_html(out, name, tool_tip):
    """Write label to html"""
    _write_dropdown_html(out, name, tool_tip, "sk-label-container", "sk-label")


_EstHTMLInfo = namedtuple('_EstHTMLInfo',
                          'type, estimators, names, name_details')


def _type_of_html_estimator(estimator):
    """Generate information about how to display an estimator.
    """
    # import here to avoid circular import from base.py
    from sklearn.pipeline import Pipeline
    from sklearn.pipeline import FeatureUnion
    from sklearn.compose import ColumnTransformer
    from sklearn.ensemble import VotingClassifier, VotingRegressor

    if isinstance(estimator, str):
        return _EstHTMLInfo('single', [estimator], [estimator], [estimator])

    elif estimator is None:
        return _EstHTMLInfo('single', [estimator], ['None'], ['None'])

    elif isinstance(estimator, Pipeline):
        estimators = [step[1] for step in estimator.steps]
        names = [step[0] for step in estimator.steps]
        name_details = [_estimator_details(est) for est in estimators]
        return _EstHTMLInfo('serial', estimators, names, name_details)

    elif isinstance(estimator, ColumnTransformer):
        estimators = [trans[1] for trans in estimator.transformers]
        names = [trans[0] for trans in estimator.transformers]
        name_details = [trans[2] for trans in estimator.transformers]
        return _EstHTMLInfo('parallel', estimators, names, name_details)

    elif isinstance(estimator, FeatureUnion):
        estimators = [trans[1] for trans in estimator.transformer_list]
        names = [trans[0] for trans in estimator.transformer_list]
        name_details = [None] * len(names)
        return _EstHTMLInfo('parallel', estimators, names, name_details)

    elif isinstance(estimator, (VotingClassifier, VotingRegressor)):
        estimators = [est[1] for est in estimator.estimators]
        names = [est[0] for est in estimator.estimators]
        name_details = [None] * len(names)
        return _EstHTMLInfo('parallel', estimators, names, name_details)

    elif hasattr(estimator, "estimator"):
        names = [estimator.__class__.__name__]
        name_details = [_estimator_details(estimator)]
        inner_estimator = estimator.estimator
        return _EstHTMLInfo('single-meta', [inner_estimator], names,
                            name_details)
    # Base estimator
    names = [estimator.__class__.__name__]
    tool_tips = [_estimator_details(estimator)]
    return _EstHTMLInfo('single', [estimator], names, tool_tips)


def _write_estimator_html(out, estimator, name):
    """Write estimator to html in serial, parallel, or by itself (single).
    """
    est_html_info = _type_of_html_estimator(estimator)

    if est_html_info.type == 'serial':
        out.write('<div class="sk-serial">')
        est_infos = zip(est_html_info.estimators, est_html_info.names,
                        est_html_info.name_details)
        for est, name, tool_tip in est_infos:
            _write_estimator_html(out, est, name)
        out.write('</div>')  # sk-serial

    elif est_html_info.type == 'parallel':
        out.write('<div class="sk-serial-item sk-dashed-wrapped">')
        if name:
            tool_tip = _estimator_details(estimator)
            _write_label_html(out, name, tool_tip)
        out.write('<div class="sk-parallel">')

        est_infos = zip(est_html_info.estimators, est_html_info.names,
                        est_html_info.name_details)
        for est, name, tool_tip in est_infos:
            out.write('<div class="sk-parallel-item">')
            _write_label_html(out, name, tool_tip)
            out.write('<div class="sk-serial">')
            _write_estimator_html(out, est, name)
            out.write('</div></div>')  # sk-parallel-item sk-serial
        out.write('</div></div>')  # sk-parallel sk-serial-item

    elif est_html_info.type == 'single-meta':
        out.write('<div class="sk-serial-item sk-dashed-wrapped">')
        _write_label_html(out, est_html_info.names[0],
                          est_html_info.name_details[0])
        _write_estimator_html(out, est_html_info.estimators[0],
                              est_html_info.estimators.__class__.__name__)
        out.write('</div>')  # sk-serial-item # sk-serial

    elif est_html_info.type == 'single':
        _write_dropdown_html(out, est_html_info.names[0],
                             est_html_info.name_details[0],
                             "sk-serial-item", "sk-estimator")


_STYLE = """
div.sk-toggleable {
  background-color: white;
}
label.sk-toggleable__label {
  cursor: pointer;
  display: block;
  width: 100%;
  margin-bottom: 0;
  padding: 0.2em 0.3em;
  box-sizing: border-box;
  text-align: center;
}
div.sk-toggleable__content {
  max-height: 0;
  max-width: 0;
  overflow: hidden;
  text-align: left;
  background-color: #f0f8ff;
}
div.sk-toggleable__content pre {
  margin: 0.2em;
  border-radius: 0.25em;
  background-color: #f0f8ff;
}
input.sk-toggleable__control:checked~div.sk-toggleable__content {
  max-height: 200px;
  max-width: 100%;
  overflow: auto;
}
input.sk-hidden--visually {
  border: 0;
  clip: rect(1px 1px 1px 1px);
  clip: rect(1px, 1px, 1px, 1px);
  height: 1px;
  margin: -1px;
  overflow: hidden;
  padding: 0;
  position: absolute;
  width: 1px;
}
div.sk-estimator {
  font-family: monospace;
  background-color: #f0f8ff;
  margin: 0.25em 0.25em;
  border: 1px dotted black;
  border-radius: 0.25em;
  box-sizing: border-box;
}
div.sk-estimator:hover {
  background-color: #c1e2ff;
}
div.sk-parallel-item::after {
  content: "";
  width: 100%;
  border-bottom: 1px solid gray;
  flex-grow: 1;
}
div.sk-label:hover label.sk-toggleable__label {
  color: #0087fe;
  background-color: rgb(246, 246, 246);
  border-radius: 0.25em;
}
div.sk-serial::before {
  content: "";
  position: absolute;
  border-left: 1px solid gray;
  box-sizing: border-box;
  top: 2em;
  bottom: 0;
  left: 50%;
}
div.sk-serial {
  display: flex;
  flex-direction: column;
  align-items: center;
  float: left;
  background: white;
}
div.sk-serial-item {
  z-index: 1;
}
div.sk-parallel {
  display: flex;
  align-items: stretch;
  justify-content: center;
}
div.sk-parallel-item {
  display: flex;
  flex-direction: column;
  position: relative;
  background: white;
}
div.sk-parallel-item:first-child::after {
  align-self: flex-end;
  width: 50%;
}
div.sk-parallel-item:last-child::after {
  align-self: flex-start;
  width: 50%;
}
div.sk-parallel-item:only-child::after {
  width: 0;
}
div.sk-dashed-wrapped {
  border: 1px dashed gray;
  padding: 0 0.3em 0.3em 0.3em;
  box-sizing: border-box;
}
div.sk-label label {
  font-family: monospace;
  font-weight: bold;
  background: white;
  display: inline-block;
  line-height: 1.4em;
}
div.sk-label-container {
  text-align: center;
  z-index: 1;
}
div.sk-container {
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
  font-weight: 400;
  position: absolute;
  top: 0;
  left: 0;
  padding: 0.5em;
  overflow: hidden;
  background-color: #f0f8ff;
  border: 1px solid gray;
  box-sizing: border-box;
  white-space: pre;
  content: attr(sk-data-tooltip);
  text-align: left;
}
[sk-data-tooltip]:hover:before {
  visibility: visible;
  opacity: 1;
  z-index: 2;
}
""".replace('  ', '').replace('\n', '')


def _estimator_repr_html(estimator, print_changed_only=True):
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
    # import here to avoid circular import from base.py
    from sklearn._config import config_context
    from sklearn.pipeline import Pipeline

    with config_context(print_changed_only=print_changed_only), \
            closing(StringIO()) as out:

        # This forces estimators to always be serial at the first layer
        if not isinstance(estimator, Pipeline):
            estimator = Pipeline([(estimator.__class__.__name__, estimator)])

        out.write('<html><head><style>')
        out.write(_STYLE)
        out.write('</style></head><body>')

        out.write('<div class="sk-top-container"><div class="sk-container">')
        _write_estimator_html(out, estimator, '')
        out.write('</div></div>')  # sk-top-container # sk-container
        out.write('</body></html>')

        html_output = out.getvalue()
        return html_output
