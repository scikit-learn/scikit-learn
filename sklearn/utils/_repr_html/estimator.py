# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import html
from contextlib import closing
from inspect import isclass
from io import StringIO
from pathlib import Path
from string import Template

from ... import config_context


class _IDCounter:
    """Generate sequential ids with a prefix."""

    def __init__(self, prefix):
        self.prefix = prefix
        self.count = 0

    def get_id(self):
        self.count += 1
        return f"{self.prefix}-{self.count}"


def _get_css_style():
    estimator_css_file = Path(__file__).parent / "estimator.css"
    params_css_file = Path(__file__).parent / "params.css"

    estimator_css = estimator_css_file.read_text(encoding="utf-8")
    params_css = params_css_file.read_text(encoding="utf-8")

    return f"{estimator_css}\n{params_css}"


_CONTAINER_ID_COUNTER = _IDCounter("sk-container-id")
_ESTIMATOR_ID_COUNTER = _IDCounter("sk-estimator-id")
_CSS_STYLE = _get_css_style()


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

    name_caption : str, default=None
        The caption below the name. `None` stands for no caption.
        Only active when kind == 'single'.

    doc_link_label : str, default=None
        The label for the documentation link. If provided, the label would be
        "Documentation for {doc_link_label}". Otherwise it will look for `names`.
        Only active when kind == 'single'.

    dash_wrapped : bool, default=True
        If true, wrapped HTML element will be wrapped with a dashed border.
        Only active when kind != 'single'.
    """

    def __init__(
        self,
        kind,
        estimators,
        *,
        names=None,
        name_details=None,
        name_caption=None,
        doc_link_label=None,
        dash_wrapped=True,
    ):
        self.kind = kind
        self.estimators = estimators
        self.dash_wrapped = dash_wrapped
        self.name_caption = name_caption
        self.doc_link_label = doc_link_label

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
    params,
    name,
    name_details,
    name_caption=None,
    doc_link_label=None,
    outer_class="sk-label-container",
    inner_class="sk-label",
    checked=False,
    doc_link="",
    is_fitted_css_class="",
    is_fitted_icon="",
    param_prefix="",
):
    """Write labeled html with or without a dropdown with named details.

    Parameters
    ----------
    out : file-like object
        The file to write the HTML representation to.
    params: str
        If estimator has `get_params` method, this is the HTML representation
        of the estimator's parameters and their values. When the estimator
        does not have `get_params`, it is an empty string.
    name : str
        The label for the estimator. It corresponds either to the estimator class name
        for a simple estimator or in the case of a `Pipeline` and `ColumnTransformer`,
        it corresponds to the name of the step.
    name_details : str
        The details to show as content in the dropdown part of the toggleable label. It
        can contain information such as non-default parameters or column information for
        `ColumnTransformer`.
    name_caption : str, default=None
        The caption below the name. If `None`, no caption will be created.
    doc_link_label : str, default=None
        The label for the documentation link. If provided, the label would be
        "Documentation for {doc_link_label}". Otherwise it will look for `name`.
    outer_class : {"sk-label-container", "sk-item"}, default="sk-label-container"
        The CSS class for the outer container.
    inner_class : {"sk-label", "sk-estimator"}, default="sk-label"
        The CSS class for the inner container.
    checked : bool, default=False
        Whether the dropdown is folded or not. With a single estimator, we intend to
        unfold the content.
    doc_link : str, default=""
        The link to the documentation for the estimator. If an empty string, no link is
        added to the diagram. This can be generated for an estimator if it uses the
        `_HTMLDocumentationLinkMixin`.
    is_fitted_css_class : {"", "fitted"}
        The CSS class to indicate whether or not the estimator is fitted. The
        empty string means that the estimator is not fitted and "fitted" means that the
        estimator is fitted.
    is_fitted_icon : str, default=""
        The HTML representation to show the fitted information in the diagram. An empty
        string means that no information is shown.
    param_prefix : str, default=""
        The prefix to prepend to parameter names for nested estimators.
    """
    out.write(
        f'<div class="{outer_class}"><div'
        f' class="{inner_class} {is_fitted_css_class} sk-toggleable">'
    )
    name = html.escape(name)
    if name_details is not None:
        name_details = html.escape(str(name_details))
        checked_str = "checked" if checked else ""
        est_id = _ESTIMATOR_ID_COUNTER.get_id()

        if doc_link:
            doc_label = "<span>Online documentation</span>"
            if doc_link_label is not None:
                doc_label = f"<span>Documentation for {doc_link_label}</span>"
            elif name is not None:
                doc_label = f"<span>Documentation for {name}</span>"
            doc_link = (
                f'<a class="sk-estimator-doc-link {is_fitted_css_class}"'
                f' rel="noreferrer" target="_blank" href="{doc_link}">?{doc_label}</a>'
            )

        name_caption_div = (
            ""
            if name_caption is None
            else f'<div class="caption">{html.escape(name_caption)}</div>'
        )
        name_caption_div = f"<div><div>{name}</div>{name_caption_div}</div>"
        links_div = (
            f"<div>{doc_link}{is_fitted_icon}</div>"
            if doc_link or is_fitted_icon
            else ""
        )

        label_html = (
            f'<label for="{est_id}" class="sk-toggleable__label {is_fitted_css_class} '
            f'sk-toggleable__label-arrow">{name_caption_div}{links_div}</label>'
        )

        fmt_str = (
            f'<input class="sk-toggleable__control sk-hidden--visually" id="{est_id}" '
            f'type="checkbox" {checked_str}>{label_html}<div '
            f'class="sk-toggleable__content {is_fitted_css_class}" '
            f'data-param-prefix="{html.escape(param_prefix)}">'
        )

        if params:
            fmt_str = "".join([fmt_str, f"{params}</div>"])
        elif name_details and ("Pipeline" not in name):
            fmt_str = "".join([fmt_str, f"<pre>{name_details}</pre></div>"])

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

    # check if estimator looks like a meta estimator (wraps estimators)
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
    out,
    estimator,
    estimator_label,
    estimator_label_details,
    is_fitted_css_class,
    is_fitted_icon="",
    first_call=False,
    param_prefix="",
):
    """Write estimator to html in serial, parallel, or by itself (single).

    For multiple estimators, this function is called recursively.

    Parameters
    ----------
    out : file-like object
        The file to write the HTML representation to.
    estimator : estimator object
        The estimator to visualize.
    estimator_label : str
        The label for the estimator. It corresponds either to the estimator class name
        for simple estimator or in the case of `Pipeline` and `ColumnTransformer`, it
        corresponds to the name of the step.
    estimator_label_details : str
        The details to show as content in the dropdown part of the toggleable label.
        It can contain information as non-default parameters or column information for
        `ColumnTransformer`.
    is_fitted_css_class : {"", "fitted"}
        The CSS class to indicate whether or not the estimator is fitted or not. The
        empty string means that the estimator is not fitted and "fitted" means that the
        estimator is fitted.
    is_fitted_icon : str, default=""
        The HTML representation to show the fitted information in the diagram. An empty
        string means that no information is shown. If the estimator to be shown is not
        the first estimator (i.e. `first_call=False`), `is_fitted_icon` is always an
        empty string.
    first_call : bool, default=False
        Whether this is the first time this function is called.
    param_prefix : str, default=""
        The prefix to prepend to parameter names for nested estimators.
        For example, in a pipeline this might be "pipeline__stepname__".
    """
    if first_call:
        est_block = _get_visual_block(estimator)
    else:
        is_fitted_icon = ""
        with config_context(print_changed_only=True):
            est_block = _get_visual_block(estimator)
    # `estimator` can also be an instance of `_VisualBlock`
    if hasattr(estimator, "_get_doc_link"):
        doc_link = estimator._get_doc_link()
    else:
        doc_link = ""
    if est_block.kind in ("serial", "parallel"):
        dashed_wrapped = first_call or est_block.dash_wrapped
        dash_cls = " sk-dashed-wrapped" if dashed_wrapped else ""
        out.write(f'<div class="sk-item{dash_cls}">')

        if estimator_label:
            if hasattr(estimator, "get_params") and hasattr(
                estimator, "_get_params_html"
            ):
                params = estimator._get_params_html(deep=False)._repr_html_inner()
            else:
                params = ""

            _write_label_html(
                out,
                params,
                estimator_label,
                estimator_label_details,
                doc_link=doc_link,
                is_fitted_css_class=is_fitted_css_class,
                is_fitted_icon=is_fitted_icon,
                param_prefix=param_prefix,
            )

        kind = est_block.kind
        out.write(f'<div class="sk-{kind}">')
        est_infos = zip(est_block.estimators, est_block.names, est_block.name_details)

        for est, name, name_details in est_infos:
            # Build the parameter prefix for nested estimators

            if param_prefix and hasattr(name, "split"):
                # If we already have a prefix, append the new component
                new_prefix = f"{param_prefix}{name.split(':')[0]}__"
            elif hasattr(name, "split"):
                # If this is the first level, start the prefix
                new_prefix = f"{name.split(':')[0]}__" if name else ""
            else:
                new_prefix = param_prefix

            if kind == "serial":
                _write_estimator_html(
                    out,
                    est,
                    name,
                    name_details,
                    is_fitted_css_class=is_fitted_css_class,
                    param_prefix=new_prefix,
                )
            else:  # parallel
                out.write('<div class="sk-parallel-item">')
                # wrap element in a serial visualblock
                serial_block = _VisualBlock("serial", [est], dash_wrapped=False)
                _write_estimator_html(
                    out,
                    serial_block,
                    name,
                    name_details,
                    is_fitted_css_class=is_fitted_css_class,
                    param_prefix=new_prefix,
                )
                out.write("</div>")  # sk-parallel-item

        out.write("</div></div>")
    elif est_block.kind == "single":
        if hasattr(estimator, "_get_params_html"):
            params = estimator._get_params_html()._repr_html_inner()
        else:
            params = ""

        _write_label_html(
            out,
            params,
            est_block.names,
            est_block.name_details,
            est_block.name_caption,
            est_block.doc_link_label,
            outer_class="sk-item",
            inner_class="sk-estimator",
            checked=first_call,
            doc_link=doc_link,
            is_fitted_css_class=is_fitted_css_class,
            is_fitted_icon=is_fitted_icon,
            param_prefix=param_prefix,
        )


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

    Examples
    --------
    >>> from sklearn.utils._repr_html.estimator import estimator_html_repr
    >>> from sklearn.linear_model import LogisticRegression
    >>> estimator_html_repr(LogisticRegression())
    '<style>#sk-container-id...'
    """
    from sklearn.exceptions import NotFittedError
    from sklearn.utils.validation import check_is_fitted

    if not hasattr(estimator, "fit"):
        status_label = "<span>Not fitted</span>"
        is_fitted_css_class = ""
    else:
        try:
            check_is_fitted(estimator)
            status_label = "<span>Fitted</span>"
            is_fitted_css_class = "fitted"
        except NotFittedError:
            status_label = "<span>Not fitted</span>"
            is_fitted_css_class = ""

    is_fitted_icon = (
        f'<span class="sk-estimator-doc-link {is_fitted_css_class}">'
        f"i{status_label}</span>"
    )
    with closing(StringIO()) as out:
        container_id = _CONTAINER_ID_COUNTER.get_id()
        style_template = Template(_CSS_STYLE)
        style_with_id = style_template.substitute(id=container_id)
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
        html_template = (
            f"<style>{style_with_id}</style>"
            f"<body>"
            f'<div id="{container_id}" class="sk-top-container">'
            '<div class="sk-text-repr-fallback">'
            f"<pre>{html.escape(estimator_str)}</pre><b>{fallback_msg}</b>"
            "</div>"
            '<div class="sk-container" hidden>'
        )

        out.write(html_template)
        _write_estimator_html(
            out,
            estimator,
            estimator.__class__.__name__,
            estimator_str,
            first_call=True,
            is_fitted_css_class=is_fitted_css_class,
            is_fitted_icon=is_fitted_icon,
        )
        with open(str(Path(__file__).parent / "estimator.js"), "r") as f:
            script = f.read()

        html_end = f"</div></div><script>{script}</script></body>"

        out.write(html_end)

        html_output = out.getvalue()
        return html_output
