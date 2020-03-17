from sklearn._config import config_context
from contextlib import closing
from io import StringIO
import uuid


class _EstHTMLBlock:
    """HTML Representation of Estimator

    If type == 'single', then the parameters are single items representing the
    single estimator
    if type == 'parallel', then the paramters are list representing the
    parallel estimators
    if type == 'serial', then the parameters are list representing the serial
    estimators
    if type == 'single-meta', then parameters represent the wrapped estimator
    """
    def __init__(self, type, estimators, names, name_details,
                 dash_wrapped=True):
        self.type = type
        self.estimators = estimators
        self.names = names
        self.name_details = name_details
        self.dash_wrapped = dash_wrapped


def _type_of_html_estimator(estimator):
    """Generate information about how to display an estimator.
    """
    if isinstance(estimator, _EstHTMLBlock):
        return estimator
    elif isinstance(estimator, str):
        return _EstHTMLBlock('single', estimator, estimator, estimator)
    elif estimator is None:
        return _EstHTMLBlock('single', estimator, 'None', 'None')
    # looks like a meta estimator
    elif (hasattr(estimator, 'estimator') and
            hasattr(getattr(estimator, 'estimator'), 'get_params')):
        wrapped_estimator = getattr(estimator, 'estimator')
        wrapped_name = wrapped_estimator.__class__.__name__
        return _EstHTMLBlock('single-meta', wrapped_estimator, wrapped_name,
                             None)
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
    if not name or isinstance(estimator, _EstHTMLBlock):
        return
    with config_context(print_changed_only=True):
        name_details = str(estimator)
    _write_label_html(out, name, name_details)


def _write_sk_item(out, dash_wrapped=True):
    dash_cls = " sk-dashed-wrapped" if dash_wrapped else ""
    out.write(f'<div class="sk-item{dash_cls}">')


def _write_estimator_html(out, estimator, name, first_call=False):
    """Write estimator to html in serial, parallel, or by itself (single).
    """
    with config_context(print_changed_only=not first_call):
        est_block = _type_of_html_estimator(estimator)

    if est_block.type == 'serial':
        _write_sk_item(out, dash_wrapped=first_call or est_block.dash_wrapped)
        _write_named_label_html(out, estimator, name)

        out.write('<div class="sk-serial">')
        est_infos = zip(est_block.estimators, est_block.names)
        for est, name in est_infos:
            if name and not isinstance(est, _EstHTMLBlock):
                name = f"{name}: {est.__class__.__name__}"
            _write_estimator_html(out, est, name)
        out.write('</div></div>')  # sk-serial sk-item

    elif est_block.type == 'parallel':
        _write_sk_item(out, dash_wrapped=est_block.dash_wrapped)
        _write_named_label_html(out, estimator, name)
        out.write('<div class="sk-parallel">')

        if est_block.name_details is None:
            name_details = (None,) * len(est_block.estimators)
        else:
            name_details = est_block.name_details

        est_infos = zip(est_block.estimators, est_block.names, name_details)
        for est, name, name_details in est_infos:
            out.write('<div class="sk-parallel-item">')
            _write_label_html(out, name, name_details)
            out.write('<div class="sk-serial">')
            _write_estimator_html(out, est, '')
            out.write('</div></div>')  # sk-parallel-item sk-serial
        out.write('</div></div>')  # sk-parallel sk-item

    elif est_block.type == 'single-meta':
        _write_sk_item(out, dash_wrapped=est_block.dash_wrapped)
        _write_named_label_html(out, estimator, name)
        out.write('<div class="sk-parallel"><div class="sk-parallel-item">')
        _write_estimator_html(out, est_block.estimators, est_block.names)
        # sk-parallel sk-parallel-item sk-item
        out.write('</div></div></div>')

    elif est_block.type == 'single':
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
  display: flex;
  flex-direction: column;
  align-items: flex-start;
  position: relative;
  float: left;
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
