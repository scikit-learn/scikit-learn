from contextlib import closing
from io import StringIO
from inspect import isclass
from string import Template
import html
from pathlib import Path

from .. import config_context


class _IDCounter:
    """Generate sequential ids with a prefix."""

    def __init__(self, prefix):
        self.prefix = prefix
        self.count = 0

    def get_id(self):
        self.count += 1
        return f"{self.prefix}-{self.count}"


_CONTAINER_ID_COUNTER = _IDCounter("sk-container-id")
_ESTIMATOR_ID_COUNTER = _IDCounter("sk-estimator-id")


class _VisualBlock:
    """HTML Representation of Estimator

    Parameters
    ----------
    kind : {'serial', 'parallel', 'single'}
        kind of HTML block

    estimators : list of estimators or `_VisualBlock`s or a single estimator
        If kind != 'single', then `estimators` is a list of
        estimators.
        If kind == 'single', then `estimators` is a single estimator.

    names : list of str, default=None
        If kind != 'single', then `names` corresponds to estimators.
        If kind == 'single', then `names` is a single string corresponding to
        the single estimator.

    name_details : list of str, str, or None, default=None
        If kind != 'single', then `name_details` corresponds to `names`.
        If kind == 'single', then `name_details` is a single string
        corresponding to the single estimator.

    dash_wrapped : bool, default=True
        If true, wrapped HTML element will be wrapped with a dashed border.
        Only active when kind != 'single'.
    """

    def __init__(
        self, kind, estimators, *, names=None, name_details=None, dash_wrapped=True
    ):
        self.kind = kind
        self.estimators = estimators
        self.dash_wrapped = dash_wrapped

        if self.kind in ("parallel", "serial"):
            if names is None:
                names = (None,) * len(estimators)
            if name_details is None:
                name_details = (None,) * len(estimators)

        self.names = names
        self.name_details = name_details

    def _sk_visual_block_(self):
        return self


def _write_label_html(
    out,
    name,
    name_details,
    outer_class="sk-label-container",
    inner_class="sk-label",
    checked=False,
    url_link="",
):
    """Write labeled html with or without a dropdown with named details"""
    out.write(f'<div class="{outer_class}"><div class="{inner_class} sk-toggleable">')
    name = html.escape(name)

    if name_details is not None:
        name_details = html.escape(str(name_details))
        label_class = "sk-toggleable__label sk-toggleable__label-arrow"

        checked_str = "checked" if checked else ""
        est_id = _ESTIMATOR_ID_COUNTER.get_id()

        if url_link:  # if the url_link is valid, use it
            doc_label = "<span>Online documentation</span>"
            if name is not None:
                doc_label = f"<span>Documentation for {name}</span>"
            doc_link = (
                '<a class="sk-estimator-doc-link" target="_blank"'
                f' href="{url_link}">?{doc_label}</a>'
            )
        else:  # no url_link, add no link to the documentation
            doc_link = ""

        fmt_str = f"""<input
        class="sk-toggleable__control sk-hidden--visually" id="{est_id}"
             type="checkbox" {checked_str}><label for="{est_id}"
             class="{label_class}">{name}{doc_link}</label><div
             class="sk-toggleable__content"><pre>{name_details}</pre></div>
            """
        out.write(fmt_str)
    else:
        out.write(f"<label>{name}</label>")
    out.write("</div></div>")  # outer_class inner_class


def _get_visual_block(estimator):
    """Generate information about how to display an estimator."""
    if hasattr(estimator, "_sk_visual_block_"):
        try:
            return estimator._sk_visual_block_()
        except Exception:
            return _VisualBlock(
                "single",
                estimator,
                names=estimator.__class__.__name__,
                name_details=str(estimator),
            )

    if isinstance(estimator, str):
        return _VisualBlock(
            "single", estimator, names=estimator, name_details=estimator
        )
    elif estimator is None:
        return _VisualBlock("single", estimator, names="None", name_details="None")

    # check if estimator looks like a meta estimator wraps estimators
    if hasattr(estimator, "get_params") and not isclass(estimator):
        estimators = [
            (key, est)
            for key, est in estimator.get_params(deep=False).items()
            if hasattr(est, "get_params") and hasattr(est, "fit") and not isclass(est)
        ]
        if estimators:
            return _VisualBlock(
                "parallel",
                [est for _, est in estimators],
                names=[f"{key}: {est.__class__.__name__}" for key, est in estimators],
                name_details=[str(est) for _, est in estimators],
            )

    return _VisualBlock(
        "single",
        estimator,
        names=estimator.__class__.__name__,
        name_details=str(estimator),
    )


def _write_estimator_html(
    out, estimator, estimator_label, estimator_label_details, first_call=False
):
    """Write estimator to html in serial, parallel, or by itself (single)."""
    # Delayed to avoid circular import
    from sklearn.base import BaseEstimator

    if first_call:
        est_block = _get_visual_block(estimator)
    else:
        with config_context(print_changed_only=True):
            est_block = _get_visual_block(estimator)
    # `estimator` can also be an instance of `_VisualBlock`
    if isinstance(estimator, BaseEstimator):
        url_link = estimator._get_url_link()
    else:
        url_link = ""
    if est_block.kind in ("serial", "parallel"):
        dashed_wrapped = first_call or est_block.dash_wrapped
        dash_cls = " sk-dashed-wrapped" if dashed_wrapped else ""
        out.write(f'<div class="sk-item{dash_cls}">')

        if estimator_label:
            _write_label_html(
                out, estimator_label, estimator_label_details, url_link=url_link
            )

        kind = est_block.kind
        out.write(f'<div class="sk-{kind}">')
        est_infos = zip(est_block.estimators, est_block.names, est_block.name_details)

        for est, name, name_details in est_infos:
            if kind == "serial":
                _write_estimator_html(out, est, name, name_details)
            else:  # parallel
                out.write('<div class="sk-parallel-item">')
                # wrap element in a serial visualblock
                serial_block = _VisualBlock("serial", [est], dash_wrapped=False)
                _write_estimator_html(out, serial_block, name, name_details)
                out.write("</div>")  # sk-parallel-item

        out.write("</div></div>")
    elif est_block.kind == "single":
        _write_label_html(
            out,
            est_block.names,
            est_block.name_details,
            outer_class="sk-item",
            inner_class="sk-estimator",
            checked=first_call,
            url_link=url_link,
        )


with open(Path(__file__).with_suffix(".css"), "r") as style_file:
    _STYLE = style_file.read()


def estimator_html_repr(estimator):
    """Build a HTML representation of an estimator.

    Read more in the :ref:`User Guide <visualizing_composite_estimators>`.

    Parameters
    ----------
    estimator : estimator object
        The estimator to visualize.

    Returns
    -------
    html: str
        HTML representation of estimator.
    """
    from sklearn.utils.validation import check_is_fitted
    from sklearn.exceptions import NotFittedError

    try:
        check_is_fitted(estimator)
        bg = "f0f8ff"
        bgh = "d4ebff"
    except NotFittedError:
        bg = "fff5e6"
        bgh = "ffe0b3"

    with closing(StringIO()) as out:
        container_id = _CONTAINER_ID_COUNTER.get_id()
        style_template = Template(_STYLE)

        style_with_id = style_template.substitute(
            id=container_id, background_color=bg, background_color_hover=bgh
        )
        estimator_str = str(estimator)

        # The fallback message is shown by default and loading the CSS sets
        # div.sk-text-repr-fallback to display: none to hide the fallback message.
        #
        # If the notebook is trusted, the CSS is loaded which hides the fallback
        # message. If the notebook is not trusted, then the CSS is not loaded and the
        # fallback message is shown by default.
        #
        # The reverse logic applies to HTML repr div.sk-container.
        # div.sk-container is hidden by default and the loading the CSS displays it.
        fallback_msg = (
            "In a Jupyter environment, please rerun this cell to show the HTML"
            " representation or trust the notebook. <br />On GitHub, the"
            " HTML representation is unable to render, please try loading this page"
            " with nbviewer.org."
        )
        html_template = """<style>{style_with_id}</style>
            <div id="{container_id}" class="sk-top-container">
            <div class="sk-text-repr-fallback">
            <pre>{estimator_str}</pre><b>{fallback_msg}</b>
            </div>
            <div class="sk-container" hidden>
        """.format(
            **{
                "container_id": container_id,
                "style_with_id": style_with_id,
                "estimator_str": html.escape(estimator_str),
                "fallback_msg": fallback_msg,
            }
        )

        out.write(html_template)

        _write_estimator_html(
            out,
            estimator,
            estimator.__class__.__name__,
            estimator_str,
            first_call=True,
        )
        out.write("</div></div>")

        html_output = out.getvalue()
        return html_output
