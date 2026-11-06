# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

"""Human-readable visualisation of metadata routing.

:func:`format_routing` takes a metadata-routing object (as returned by
:func:`~sklearn.utils.metadata_routing.get_routing_for_object`) and returns a
string showing how each metadata parameter flows through an estimator tree,
annotated with per-method status glyphs, plus a per-parameter summary.

:func:`visualise_routing` is a thin wrapper that prints the result.
"""

import os
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from io import StringIO

from sklearn.utils._metadata_requests import (
    COMPOSITE_METHODS,
    SIMPLE_METHODS,
    UNUSED,
    WARN,
    MetadataRequest,
    _routing_repr,
    request_is_alias,
)

_GLYPH = {
    "requested": "✓",
    "ignored": "✗",
    "warns": "⚠",
    "errors": "⛔",
    "unused": "⊘",
}

# Listed most-severe first so the summary's leading glyph and section order
# reflect what is most important to notice.
_CATEGORY_ORDER = ["errors", "warns", "ignored", "requested", "unused"]

# ANSI 8-colour codes. 8-colour rather than 256/truecolor so the output
# adapts to both light and dark terminal themes: basic red/green/yellow/cyan
# render sensibly on any background.
_ANSI = {
    "red": "\033[31m",
    "yellow": "\033[33m",
    "green": "\033[32m",
    "cyan": "\033[36m",
    "reset": "\033[0m",
}

# Ignored/unused are deliberately uncoloured: they represent benign states and
# reserving a colour for each would dilute the red/yellow danger signal.
_CATEGORY_COLOUR = {
    "requested": "green",
    "warns": "yellow",
    "errors": "red",
}


def _colour(text, name, enabled):
    if not enabled or name is None:
        return text
    return f"{_ANSI[name]}{text}{_ANSI['reset']}"


def _glyph(category, colour):
    return _colour(_GLYPH[category], _CATEGORY_COLOUR.get(category), colour)


def _arrow(colour):
    """Coloured arrow used between an alias and its component param name."""
    return _colour("→", "cyan", colour)


def _legend(colour):
    return (
        "Legend: "
        f"{_glyph('requested', colour)} requested   "
        f"{_glyph('ignored', colour)} ignored   "
        f"{_glyph('warns', colour)} warn on use   "
        f"{_glyph('errors', colour)} error on use   "
        f"{_glyph('unused', colour)} unused\n"
        f'Aliased requests are shown as "user_name {_arrow(colour)} component_name".'
    )


def _auto_colour():
    """Best-effort TTY detection with NO_COLOR / FORCE_COLOR honoured."""
    if os.environ.get("NO_COLOR"):
        return False
    if os.environ.get("FORCE_COLOR"):
        return True
    return sys.stdout.isatty()


# --- Status classification ---------------------------------------------------


def _status_category(status):
    """Classify a raw metadata-request value into a summary category."""
    if status is True or request_is_alias(status):
        return "requested"
    if status is False:
        return "ignored"
    if status is None:
        return "errors"
    if status == WARN:
        return "warns"
    if status == UNUSED:
        return "unused"
    return None


def _method_glyph(status, colour):
    # Aliased requests are "requested" by definition — the alias mapping is
    # already visible in the ``name → component`` prefix, so we do not need a
    # distinct glyph to mark them here.
    cat = _status_category(status)
    if cat is None:
        return "?"
    return _glyph(cat, colour)


def _expand(methods):
    """Expand composite method names into their simple parts.

    Composites are kept in the output set alongside their parts so that a
    ``MethodMapping`` edge whose caller is a composite name still matches when
    the composite is reachable.
    """
    out = set()
    for m in methods:
        out.add(m)
        if m in COMPOSITE_METHODS:
            out.update(COMPOSITE_METHODS[m])
    return out


# --- Node tree built from a MetadataRequest / MetadataRouter ----------------


@dataclass
class _Node:
    owner_repr: str
    step_name: object  # str on children; None at root
    path: str  # "root_repr/mapping_name/..."
    statuses: dict = field(default_factory=dict)  # param -> {method: status}
    reachable: set = field(default_factory=set)  # simple methods here
    children: list = field(default_factory=list)
    # (param, method) -> set of root simple methods that can invoke this
    # consumption. Populated in the same pass as statuses.
    origins: dict = field(default_factory=dict)


def _build(obj, path, step_name, incoming_by_root):
    """Recursively build a :class:`_Node` from a routing object.

    ``incoming_by_root`` maps each root simple-method to the set of callee
    method names that reach this node under that root. Keeping it keyed by
    root lets us populate ``origins`` without a second traversal.
    """
    all_reachable = set()
    for callers in incoming_by_root.values():
        all_reachable |= _expand(callers)

    statuses = defaultdict(dict)
    origins = {}

    def record(consumer):
        # Iterate simple methods only: composite methods (``fit_transform``)
        # are synthesised views over their simple parts on ``MetadataRequest``
        # and would produce duplicate rows in the summary.
        for method in SIMPLE_METHODS:
            if method not in all_reachable or not hasattr(consumer, method):
                continue
            for param, status in getattr(consumer, method).requests.items():
                statuses[param][method] = status
                for root, callers in incoming_by_root.items():
                    if method in _expand(callers):
                        origins.setdefault((param, method), set()).add(root)

    if isinstance(obj, MetadataRequest):
        record(obj)
        return _Node(
            owner_repr=_routing_repr(obj.owner),
            step_name=step_name,
            path=path,
            statuses=dict(statuses),
            reachable=all_reachable,
            origins=origins,
        )

    if obj._self_request is not None:
        record(obj._self_request)

    children = []
    for name, pair in obj._route_mappings.items():
        child_incoming = {}
        for root, callers in incoming_by_root.items():
            expanded = _expand(callers)
            callees = {p.callee for p in pair.mapping if p.caller in expanded}
            if callees:
                child_incoming[root] = callees
        if not child_incoming:
            continue
        children.append(_build(pair.router, f"{path}/{name}", name, child_incoming))

    return _Node(
        owner_repr=_routing_repr(obj.owner),
        step_name=step_name,
        path=path,
        statuses=dict(statuses),
        reachable=all_reachable,
        children=children,
        origins=origins,
    )


def _build_tree(routing):
    # Seed every simple method as its own root so that every path in the tree
    # knows which invocations can reach it.
    initial = {m: {m} for m in SIMPLE_METHODS}
    return _build(routing, _routing_repr(routing.owner), None, initial)


# --- Filter handling --------------------------------------------------------


def _as_set(value):
    if value is None:
        return None
    if isinstance(value, str):
        return {value}
    return set(value)


def _coerce_method_filter(value):
    """Normalise the ``method=`` filter.

    Composite names in the filter are expanded to their simple parts so that
    ``method="fit_transform"`` matches the ``fit`` and ``transform`` roots.
    """
    as_set = _as_set(value)
    if as_set is None:
        return None
    expanded = set()
    for m in as_set:
        if m in COMPOSITE_METHODS:
            expanded.update(COMPOSITE_METHODS[m])
        else:
            expanded.add(m)
    return expanded


def _user_param(comp_param, status):
    """User-facing parameter name — the alias string if ``status`` is one."""
    if request_is_alias(status):
        return status
    return comp_param


def _keep_pair(comp_param, method, status, origins, method_filter, param_filter):
    if method_filter is not None:
        roots = origins.get((comp_param, method), set())
        if not roots & method_filter:
            return False
    if param_filter is not None and _user_param(comp_param, status) not in param_filter:
        return False
    return True


# --- Tree rendering ---------------------------------------------------------


def _format_param(comp_param, method_statuses, colour):
    """Render one parameter's per-method glyphs as an inline string.

    When some methods alias and others do not, the output is split into two
    comma-separated segments so the ``alias → param`` prefix only annotates
    the aliasing methods.
    """
    aliased = []
    direct = []
    alias_name = None
    for method in SIMPLE_METHODS:
        if method not in method_statuses:
            continue
        status = method_statuses[method]
        token = f"{method}{_method_glyph(status, colour)}"
        if request_is_alias(status):
            alias_name = status
            aliased.append(token)
        else:
            direct.append(token)

    parts = []
    if aliased:
        parts.append(f"{alias_name} {_arrow(colour)} {comp_param}[{','.join(aliased)}]")
    if direct:
        parts.append(f"{comp_param}[{','.join(direct)}]")
    return ", ".join(parts)


def _write_params(node, out, indent, method_filter, param_filter, colour):
    for param in sorted(node.statuses):
        visible = {
            m: s
            for m, s in node.statuses[param].items()
            if _keep_pair(param, m, s, node.origins, method_filter, param_filter)
        }
        if visible:
            out.write(f"{indent}    ➤ {_format_param(param, visible, colour)}\n")


def _write_children(node, out, prefix, method_filter, param_filter, colour):
    for i, child in enumerate(node.children):
        is_last = i == len(node.children) - 1
        connector = "└── " if is_last else "├── "
        label = (
            f"{child.step_name} ({child.owner_repr})"
            if child.step_name is not None
            else child.owner_repr
        )
        out.write(f"{prefix}{connector}{label}\n")
        child_indent = prefix + ("    " if is_last else "│   ")
        _write_params(child, out, child_indent, method_filter, param_filter, colour)
        _write_children(child, out, child_indent, method_filter, param_filter, colour)


def _render_tree(root, out, method_filter, param_filter, colour):
    out.write(f"{root.owner_repr}\n")
    _write_params(root, out, "", method_filter, param_filter, colour)
    _write_children(root, out, "", method_filter, param_filter, colour)


# --- Summary rendering ------------------------------------------------------


def _collect_summary(node, summary):
    """Flatten a node tree into ``{root_method: {user_param: {cat: [paths]}}}``."""
    for comp_param, m_to_s in node.statuses.items():
        for method, status in m_to_s.items():
            cat = _status_category(status)
            if cat is None:
                continue
            uparam = _user_param(comp_param, status)
            location = f"{node.path}.{method}"
            # origins may be empty if this method is listed on a node that no
            # root actually reaches — shouldn't happen given the traversal, but
            # fall back to the method itself to avoid silently dropping rows.
            roots = node.origins.get((comp_param, method)) or {method}
            for root in roots:
                summary[root][uparam][cat].append(location)
    for child in node.children:
        _collect_summary(child, summary)


def _leading_category(cats):
    if "errors" in cats:
        return "errors"
    if "warns" in cats:
        return "warns"
    if set(cats) == {"ignored"}:
        return "ignored"
    return "requested"


def _render_summary(root, out, method_filter, param_filter, colour):
    summary = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    _collect_summary(root, summary)

    filtered = {}
    for root_method, by_param in summary.items():
        if method_filter is not None and root_method not in method_filter:
            continue
        params = {
            p: cats
            for p, cats in by_param.items()
            if param_filter is None or p in param_filter
        }
        if params:
            filtered[root_method] = params

    if not filtered:
        return

    out.write("\nParameter summary:\n")
    for root_method in SIMPLE_METHODS:
        if root_method not in filtered:
            continue
        out.write(f"{root_method}\n")
        for param in sorted(filtered[root_method]):
            cats = filtered[root_method][param]
            out.write(f" ├─ {_glyph(_leading_category(cats), colour)} {param}\n")
            for cat in _CATEGORY_ORDER:
                if cat not in cats:
                    continue
                out.write(f" │   • {_glyph(cat, colour)} {cat}:\n")
                for location in cats[cat]:
                    out.write(f" │       - {location}\n")
        out.write("\n")


# --- Public API -------------------------------------------------------------


def format_routing(
    routing,
    *,
    tree=True,
    summary=True,
    method=None,
    param=None,
    legend=True,
    colour=False,
):
    """Return a human-readable string describing metadata routing.

    Parameters
    ----------
    routing : MetadataRequest or MetadataRouter
        Routing object, typically obtained from
        :func:`~sklearn.utils.metadata_routing.get_routing_for_object`.

    tree : bool, default=True
        Include the hierarchical tree view.

    summary : bool, default=True
        Include the per-parameter summary grouped by root method.

    method : str or iterable of str, default=None
        If given, restrict the output to parameters reachable from these root
        methods. Composite names (``"fit_transform"``, ``"fit_predict"``) are
        expanded to their simple parts.

    param : str or iterable of str, default=None
        If given, restrict the output to these user-facing parameter names.
        Aliases are resolved first, so passing the alias matches.

    legend : bool, default=True
        Prepend a glyph legend.

    colour : bool, default=False
        If True, wrap status glyphs and the alias arrow in ANSI escape codes:
        red for errors, yellow for warns, green for requested, cyan for the
        alias arrow. Defaults to False so the returned string stays plain text
        when captured to logs or tested; :func:`visualise_routing` enables it
        automatically when stdout is a terminal.

    Returns
    -------
    text : str
        The formatted visualisation.
    """
    method_filter = _coerce_method_filter(method)
    param_filter = _as_set(param)

    root = _build_tree(routing)

    out = StringIO()
    if legend:
        out.write(_legend(colour) + "\n\n")
    if tree:
        _render_tree(root, out, method_filter, param_filter, colour)
    if summary:
        _render_summary(root, out, method_filter, param_filter, colour)
    return out.getvalue()


def visualise_routing(routing, **kwargs):
    """Print a human-readable metadata routing diagram.

    See :func:`format_routing` for keyword arguments. ``colour`` defaults to
    True when stdout is a terminal (and ``NO_COLOR`` is unset), False
    otherwise.
    """
    kwargs.setdefault("colour", _auto_colour())
    print(format_routing(routing, **kwargs))
