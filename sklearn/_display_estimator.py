from sklearn._config import config_context
from contextlib import closing
from io import StringIO
import uuid


class _VisualBlock:
    """HTML Representation of Estimator

    Parameters
    ----------
    kind : {'serial', 'parallel', 'single'}
        kind of HTML block

    estimators : list of estimators or `_VisualBlock`s or a single estimator
        If kind is in ('parallel', 'serial'), then `estimators` is a list of
        estimators.
        If kind == 'single', then `estimators` is a single estimator.

    names : list of str
        If kind in ('parallel', 'serial'), then `names` corresponds to
        estimators
        If kind is 'single', then `names` is a single string corresponding to
        the single estimator.

    name_details : list of str, str, or None, default=None
        If kind == 'parallel', then `name_details` corresponds to `names`.
        If kind == 'single', then `name_details` is a single string
        corresponding to the single estimator.
        `name_details` is not used when kind == 'single'.

    dash_wrapped : bool, default=True
        If true, wrapped HTML element will be wrapped with a dashed border.
    """
    def __init__(self, kind, estimators, names, name_details=None,
                 dash_wrapped=True):
        self.kind = kind
        self.estimators = estimators
        self.names = names
        self.dash_wrapped = dash_wrapped

        if self.kind == 'parallel' and name_details is None:
            name_details = (None, ) * len(names)

        self.name_details = name_details


def _get_visual_block(estimator):
    """Generate information about how to display an estimator.
    """
    if isinstance(estimator, _VisualBlock):
        return estimator
    elif isinstance(estimator, str):
        return _VisualBlock('single', estimator, estimator, estimator)
    elif estimator is None:
        return _VisualBlock('single', estimator, 'None', 'None')
    # looks like a meta estimator
    elif (hasattr(estimator, 'estimator') and
            hasattr(getattr(estimator, 'estimator'), 'get_params')):
        wrapped_estimator = getattr(estimator, 'estimator')
        wrapped_name = wrapped_estimator.__class__.__name__
        return _VisualBlock('serial', [wrapped_estimator], [wrapped_name])
    return estimator._sk_repr_html()


def _write_label_html(out, name, name_details,
                      outer_class="sk-label-container",
                      inner_class="sk-label",
                      checked=False):
    """Write labeled html with or without a dropdown with named details"""
    out.write(f'<div class="{outer_class}">'
              f'<div class="{inner_class} sk-toggleable">')

    if name_details is not None:
        checked_str = 'checked' if checked else ''
        est_id = uuid.uuid4()
        out.write(f'<input class="sk-toggleable__control sk-hidden--visually" '
                  f'id="{est_id}" type="checkbox" {checked_str}>'
                  f'<label class="sk-toggleable__label" for="{est_id}">'
                  f'{name}</label>'
                  f'<div class="sk-toggleable__content"><pre>{name_details}'
                  f'</pre></div>')
    else:
        out.write(f'<label>{name}</label>')
    out.write('</div></div>')  # outer_class inner_class


def _write_named_label_html(out, estimator, name):
    """Write label with details based on name"""
    if not name or isinstance(estimator, _VisualBlock):
        return
    with config_context(print_changed_only=True):
        name_details = str(estimator)
    _write_label_html(out, name, name_details)


def _write_estimator_html(out, estimator, name, first_call=False):
    """Write estimator to html in serial, parallel, or by itself (single).
    """
    with config_context(print_changed_only=not first_call):
        est_block = _get_visual_block(estimator)

    if est_block.kind == 'serial':
        dashed_wrapped = first_call or est_block.dash_wrapped
        dash_cls = " sk-dashed-wrapped" if dashed_wrapped else ""
        out.write(f'<div class="sk-item{dash_cls}">')

        _write_named_label_html(out, estimator, name)

        out.write('<div class="sk-serial">')
        est_infos = zip(est_block.estimators, est_block.names)
        for est, name in est_infos:
            _write_estimator_html(out, est, name)
        out.write('</div></div>')  # sk-serial sk-item

    elif est_block.kind == 'parallel':
        dashed_wrapped = first_call or est_block.dash_wrapped
        dash_cls = " sk-dashed-wrapped" if dashed_wrapped else ""
        out.write(f'<div class="sk-item{dash_cls}">')

        _write_named_label_html(out, estimator, name)
        out.write('<div class="sk-parallel">')

        est_infos = zip(est_block.estimators, est_block.names,
                        est_block.name_details)
        for est, name, name_details in est_infos:
            out.write('<div class="sk-parallel-item">')
            _write_label_html(out, name, name_details)
            out.write('<div class="sk-serial">')
            _write_estimator_html(out, est, '')
            out.write('</div></div>')  # sk-parallel-item sk-serial
        out.write('</div></div>')  # sk-parallel sk-item

    elif est_block.kind == 'single':
        _write_label_html(out, est_block.names, est_block.name_details,
                          outer_class="sk-item", inner_class="sk-estimator",
                          checked=first_call)


_STYLE = """
div.sk-top-container {
  color: black;
  background-color: white;
}
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
  color: black;
  border-radius: 0.25em;
  background-color: #f0f8ff;
}
input.sk-toggleable__control:checked~div.sk-toggleable__content {
  max-height: 200px;
  max-width: 100%;
  overflow: auto;
}
div.sk-estimator input.sk-toggleable__control:checked~label.sk-toggleable__label {
  background-color: #d4ebff;
}
div.sk-label input.sk-toggleable__control:checked~label.sk-toggleable__label {
  background-color: #d4ebff;
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
  background-color: #d4ebff;
}
div.sk-parallel-item::after {
  content: "";
  width: 100%;
  border-bottom: 1px solid gray;
  flex-grow: 1;
}
div.sk-label:hover label.sk-toggleable__label {
  background-color: #d4ebff;
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
  background-color: white;
}
div.sk-item {
  z-index: 1;
}
div.sk-parallel {
  display: flex;
  align-items: stretch;
  justify-content: center;
  background-color: white;
}
div.sk-parallel-item {
  display: flex;
  flex-direction: column;
  position: relative;
  background-color: white;
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
  margin: 0.2em;
  box-sizing: border-box;
  padding-bottom: 0.1em;
  background-color: white;
  position: relative;
}
div.sk-label label {
  font-family: monospace;
  font-weight: bold;
  background-color: white;
  display: inline-block;
  line-height: 1.2em;
}
div.sk-label-container {
  position: relative;
  z-index: 2;
  text-align: center;
}
div.sk-container {
  display: inline-block;
  position: relative;
}
""".replace('  ', '').replace('\n', '')  # noqa


def _estimator_repr_html(estimator):
    """Build a HTML representation of an estimator

    Parameters
    ----------
    estimator : estimator object
        The estimator to visualize.

    Returns
    -------
    html: str or iPython HTML object
        HTML representation of estimator. When called in jupyter notebook or
        lab, a iPython HTML object is returned.
    """
    with closing(StringIO()) as out:
        out.write(f'<!DOCTYPE html><html lang="en">'
                  f'<head><title>sklearn-viz</title>'
                  f'<style>{_STYLE}</style></head><body>'
                  f'<div class="sk-top-container"><div class="sk-container">')
        _write_estimator_html(out, estimator, estimator.__class__.__name__,
                              first_call=True)
        out.write('</div></div></body></html>')

        html_output = out.getvalue()
        return html_output
