"""sklearn.utils._metadata_routing_visualise
================================================
Utilities to *inspect* and *visualise* scikit-learn's experimental
*metadata-routing* configuration.  The module is **purely diagnostic** – it
never mutates routing information – and can be imported inside notebooks or
unit-tests to print human-readable ASCII diagrams that explain how
``sample_weight``, ``groups`` … flow through complex estimator pipelines.

Overview of sections
--------------------
1.  **Public entry-points**
    • ``visualise_routing`` – prints a hierarchical tree with parameters
      annotated inline (the primary public helper).

2.  **Core data gathering**
    ``_collect_routing_info`` walks a (potentially nested) structure of
    :class:`~sklearn.utils.metadata_routing.MetadataRequest` and
    :class:`~sklearn.utils.metadata_routing.MetadataRouter` objects and returns a
    *rich* mapping that contains, for **every** component path:

    ::

        {
            "params"          : set[str],             # names understood by component
            "methods"         : {param -> {methods}}, # caller methods that *request*
            "aliases"         : {param -> alias},     # component-param → user-param
            "statuses"        : {param -> {method -> status}},
            "existing_methods": set[str],             # methods implemented on object
        }

    The heavy lifting of reachability analysis (only consider methods that can
    actually be invoked) lives in ``_compute_reachable_methods``.

3.  **Rendering helpers**
    Once the mapping is built we only need pretty-printers.  Each view resides
    in its own function (``_display_tree_new``, ``_build_flow_diagram`` …) so
    that UI changes remain isolated from the routing logic.

4.  **Test fixtures** (bottom of file)
    The last section sets up an end-to-end *example* pipeline that is executed
    when running the file directly (or when doctests need a quick demo).  The
    code lives at module level on purpose so our unit-tests can import the
    objects without duplicating boilerplate.

Conventions & design notes
--------------------------
* A *path* is a ``/``-separated string that encodes the position of a component
  inside the estimator tree (e.g. ``"Pipeline/StandardScaler"``).  Paths are
  unique keys in the mapping produced by ``_collect_routing_info``.
* Unicode glyphs (``✓ ✗ ↗ ⊘ ⚠``) are used to keep the output compact yet
  expressive.  They are generated via ``_get_status_indicator``.
* The module relies *only* on the public helpers defined in
  ``sklearn.utils._metadata_requests`` and keeps a strict separation between
  *data collection* and *rendering* so that new visualisation styles can be
  added without touching the traversal logic.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from collections import defaultdict

from sklearn.utils._metadata_requests import (
    COMPOSITE_METHODS,
    SIMPLE_METHODS,
    UNUSED,
    WARN,
    MetadataRequest,
    MetadataRouter,
    _routing_repr,
    request_is_alias,
)

# -----------------------------------------------------------------------------
# General-purpose helpers (shared across collectors / renderers)
# -----------------------------------------------------------------------------


def _expand_methods(methods):
    """Expand *composite* methods (like ``fit_transform``) into their constituent
    *simple* methods (``fit`` and ``transform``).

    Parameters
    ----------
    methods : Iterable[str]
        Iterable of method names which may include composite entries present in
        ``sklearn.utils._metadata_requests.COMPOSITE_METHODS``.

    Returns
    -------
    set[str]
        A set containing every *simple* method after expansion.
    """
    expanded: set[str] = set()
    for m in methods:
        if m in COMPOSITE_METHODS:
            expanded.update(COMPOSITE_METHODS[m])
        expanded.add(m)
    return expanded


# -----------------------------------------------------------------------------
# Status classification helpers & glyph mapping
# -----------------------------------------------------------------------------

# Mapping of *high-level* status categories to the glyph shown in the UI.  Note
# that both "errors" (raise) and "warns" (warn-only) share the same ⚠ symbol –
# the distinction is made textual when required.
_STATUS_GLYPH: dict[str, str] = {
    "requested": "✓",  # param requested (True) or alias
    "ignored": "✗",  # param explicitly *not* requested (False)
    "warns": "⚠",  # warn if passed (WARN sentinel)
    "errors": "⛔",  # error if passed (None)
    "unused": "⊘",  # UNUSED sentinel
}

# Ordered list so that the parameter summary is printed in a predictable and
# meaningful order (most severe first).
_CATEGORY_ORDER: list[tuple[str, str]] = [
    ("errors", _STATUS_GLYPH["errors"]),
    ("warns", _STATUS_GLYPH["warns"]),
    ("ignored", _STATUS_GLYPH["ignored"]),
    ("requested", _STATUS_GLYPH["requested"]),
    ("unused", _STATUS_GLYPH["unused"]),
]


def _status_category(status):
    """Return the *category* name for a raw request *status* value."""
    # `True` → requested -------------------------------------------------------
    if status is True or request_is_alias(status):
        return "requested"

    # `False` → explicitly ignored -------------------------------------------
    if status is False:
        return "ignored"

    # `None`   → raise error if metadata provided ----------------------------
    if status is None:
        return "errors"

    # WARN sentinel → warn only ----------------------------------------------
    if status == WARN:
        return "warns"

    # UNUSED sentinel → parameter is accepted but never used -----------------
    if status == UNUSED:
        return "unused"

    # Fallback (should not happen) -------------------------------------------
    return "unknown"


def visualise_routing(routing_info):
    """
    Visualize metadata routing information.

    This diagnostic always displays *all* possible metadata parameters with
    their status indicators in a unified, comprehensive view.
    """
    print("\n=== METADATA ROUTING TREE ===")

    # Get all routing information
    routing_map = _collect_routing_info(routing_info)

    # Display tree structure without duplicate root entry
    _display_tree(
        routing_info,
        routing_map,
        prefix="",
        is_last=True,
        parent_path="",
        root=True,
    )

    # ------------------------------------------------------------------
    # New style: root-level method blocks
    # ------------------------------------------------------------------

    summary_by_method = _summarise_params_by_method(routing_info, routing_map)

    if summary_by_method:
        print("\nParameter summary:")

        glyph_for = dict(_CATEGORY_ORDER)

        for root_method in SIMPLE_METHODS:
            if root_method not in summary_by_method:
                continue

            print(f"{root_method}")

            params_for_method = summary_by_method[root_method]

            for param in sorted(params_for_method):
                cats = params_for_method[param]

                # Leading glyph precedence (error > warn > ignored > requested)
                if "errors" in cats:
                    leading = glyph_for["errors"]
                elif "warns" in cats:
                    leading = glyph_for["warns"]
                elif set(cats.keys()) == {"ignored"}:
                    leading = glyph_for["ignored"]
                else:
                    leading = glyph_for["requested"]

                print(f" ├─ {leading} {param}")

                for cat, glyph in _CATEGORY_ORDER:
                    if cat in cats:
                        print(f" │   • {glyph} {cat}:")
                        for p in cats[cat]:
                            print(f" │       - {_shorten_path(p)}")

            print()  # blank line between methods

    else:
        print("\nNo parameter summary.")


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================


def _compute_method_origins(router):
    """Return mapping {"path.method": set(root_methods)}."""

    origins = defaultdict(set)

    def _walk(obj, incoming_methods: set[str], root_method: str, path=""):
        current_path = (
            f"{path}/{_routing_repr(obj.owner)}" if path else _routing_repr(obj.owner)
        )

        inc_expanded = _expand_methods(incoming_methods)

        for m in inc_expanded:
            origins[f"{current_path}.{m}"].add(root_method)

        if isinstance(obj, MetadataRequest):
            return  # leaf

        # handle self request
        if obj._self_request is not None:
            _walk(obj._self_request, inc_expanded, root_method, current_path)

        # children via mappings
        for name, mpair in obj._route_mappings.items():
            child_obj = mpair.router
            child_methods_raw = {
                p.callee for p in mpair.mapping if p.caller in inc_expanded
            }
            child_methods = _expand_methods(child_methods_raw)
            if child_methods:
                _walk(child_obj, child_methods, root_method, f"{current_path}/{name}")

    for root in SIMPLE_METHODS:
        _walk(router, {root}, root)

    return origins


def _collect_routing_info(router):
    """Return a *rich* mapping describing every parameter/method/alias path.

    The returned ``defaultdict`` has *paths* ( ``"Pipeline/StandardScaler"``)
    as keys and for each path stores::

        {
            "params"          : set[str],             # every param encountered
            "methods"         : {param -> {methods}}, # methods that actively request
            "aliases"         : {param -> alias},     # component-param → user-param
            "statuses"        : {param -> {method -> status}},
            "existing_methods": set[str],             # methods implemented on object
        }

    The heavy lifting is delegated to two small helpers so that the control-flow
    is easier to follow:

    * :func:`_collect_request_info` handles *simple* consumers
      (:class:`~sklearn.utils.metadata_routing.MetadataRequest`).
    * :func:`_collect_router_info` handles routers which can also include a
      *self-request* consumer.
    """

    info: dict[str, dict] = defaultdict(
        lambda: {
            "params": set(),
            "methods": defaultdict(set),
            "aliases": {},
            "statuses": defaultdict(dict),
            "existing_methods": set(),
        }
    )

    # ------------------------------------------------------------------
    # Helper functions – defined *inside* to keep them private to algorithm
    # ------------------------------------------------------------------

    def _record_status(current_path: str, param: str, method: str, status):
        """Record *status* for (*path*, *param*, *method*) in all auxiliary
        structures.

        The helper is on the hot-path of the DFS traversal.  We therefore cache
        expensive predicate checks (`request_is_alias`) so each call evaluates
        them only once.
        """

        # Cache predicate evaluations once ------------------------------------
        is_alias = request_is_alias(status)
        is_requested = (status is True) or is_alias

        # 1. Raw status mapping (used later for glyph generation)
        info[current_path]["statuses"][param][method] = status

        # 2. Which methods *actively* request the parameter -------------------
        if is_requested:
            info[current_path]["methods"][param].add(method)

        # 3. Alias bookkeeping (component_param → user_param) -----------------
        if is_alias and status != param:
            info[current_path]["aliases"][param] = status

    # ---------------------------------------
    # Simple consumer (MetadataRequest) branch
    # ---------------------------------------

    def _collect_request_info(obj, current_path: str, reachable_here: set[str]):
        """Handle the *consumer* case.

        Parameters
        ----------
        obj : MetadataRequest
        current_path : str
            Absolute path in the routing tree.
        reachable_here : set[str]
            Methods that can *actually* be invoked given parent→child mappings.
        """
        for method in SIMPLE_METHODS:
            if method not in reachable_here or not hasattr(obj, method):
                continue

            info[current_path]["existing_methods"].add(method)
            method_req = getattr(obj, method)

            # Consider only the parameters that *this specific method* can
            # handle so that unrelated parameters for other methods never show
            # up as "not requested" here.
            param_source = method_req.requests.keys()

            for param in param_source:
                info[current_path]["params"].add(param)
                alias = method_req.requests.get(param, False)
                _record_status(current_path, param, method, alias)

    # ---------------------------------------
    # Router branch (can include a *self* request)
    # ---------------------------------------

    def _collect_router_info(obj, current_path: str, reachable_here: set[str]):
        """Handle :class:`MetadataRouter` instances (including *self* request)."""
        if obj._self_request is not None:
            _collect_request_info(obj._self_request, current_path, reachable_here)

        # Recurse into children -------------------------------------------------
        for name, mapping in obj._route_mappings.items():
            child_methods_raw = {
                pair.callee for pair in mapping.mapping if pair.caller in reachable_here
            }
            child_methods = _expand_methods(child_methods_raw)
            if child_methods:
                _collect(
                    mapping.router,
                    f"{current_path}/{name}",  # becomes *path* for child
                    child_methods,
                )

    # ------------------------------------------------------------------
    # Main DFS traversal
    # ------------------------------------------------------------------

    def _collect(
        obj, path: "str" = "", incoming_methods: set[str] = set(SIMPLE_METHODS)
    ):
        current_path = (
            f"{path}/{_routing_repr(obj.owner)}" if path else _routing_repr(obj.owner)
        )
        reachable_here = _expand_methods(incoming_methods)

        if isinstance(obj, MetadataRequest):
            _collect_request_info(obj, current_path, reachable_here)
        elif isinstance(obj, MetadataRouter):
            _collect_router_info(obj, current_path, reachable_here)

    _collect(router)
    return info


def _get_status_indicator(status):
    """Return the glyph corresponding to *status* (✓, ✗ …)."""
    if request_is_alias(status):
        # Alias itself gets the alias arrow glyph.
        return "↗"

    cat = _status_category(status)
    return _STATUS_GLYPH.get(cat, "?")


def _format_param_with_status(param, statuses, aliases):
    """Return a compact textual representation of *param* annotated with
    method-level status indicators.

    The implementation groups methods by *alias* vs *direct* usage in a single
    pass, avoiding the duplication present in the earlier version while
    preserving the exact output format.  When some methods use an alias and
    others do not, two comma-separated segments are produced, e.g.::

        inner_weights→sample_weight[fit↗], sample_weight[partial_fit⚠]
    """

    # `param_statuses` maps method → raw status (True/False/alias/…)
    param_statuses = statuses.get(param, {})
    if not param_statuses:
        # No status information – fall back to simple param/alias name.
        return f"{aliases.get(param, param)}→{param}" if param in aliases else param

    # Group methods by whether they rely on an alias (str) or not.
    grouped: dict[bool, list[str]] = {True: [], False: []}
    for method in SIMPLE_METHODS:
        if method not in param_statuses:
            continue
        status = param_statuses[method]
        indicator = _get_status_indicator(status)
        grouped[request_is_alias(status)].append(f"{method}{indicator}")

    parts: list[str] = []
    # Iterate in stable order so alias segment (if any) comes first.
    for uses_alias in (True, False):
        if not grouped[uses_alias]:
            continue

        # Decide *display name* for this segment.
        if uses_alias:
            # All alias statuses share the same alias string by construction.
            alias_name = next(s for s in param_statuses.values() if request_is_alias(s))
            name = f"{alias_name}→{param}"
        else:
            name = param

        parts.append(f"{name}[{','.join(grouped[uses_alias])}]")

    return ", ".join(parts)


# -----------------------------------------------------------------------------
# Tree rendering helpers (branch / prefix calculation)
# -----------------------------------------------------------------------------


def _branch_connector(root: bool, prefix: str, is_last: bool) -> str:
    """Return the unicode branch connector for the current node."""
    if root and prefix == "":
        return ""
    return "└── " if is_last else "├── "


def _param_prefix(prefix: str, root: bool, is_last: bool) -> str:
    """Indentation for parameter lines that hang under a node."""
    if root and prefix == "":
        return "    "
    return prefix + ("    " if is_last else "│   ") + "    "


def _child_prefix(prefix: str, root: bool, is_last: bool) -> str:
    """Prefix to pass when recursing into child nodes."""
    if root and prefix == "":
        return prefix  # keep root prefix empty for the very top-level call
    return prefix + ("    " if is_last else "│   ")


def _display_tree(
    router,
    routing_map,
    prefix="",
    is_last=True,
    step_name=None,
    parent_path="",
    root=False,
):
    """Display the routing tree with proper formatting and inline parameters."""
    # Get current path
    if parent_path:
        current_path = f"{parent_path}/{_routing_repr(router.owner)}"
    else:
        current_path = _routing_repr(router.owner)

    # Determine if this node has parameters
    node_info = routing_map.get(current_path, {})
    has_params = bool(node_info.get("params"))

    # Build the display line
    connector = _branch_connector(root, prefix, is_last)

    display_parts = []
    if step_name:
        display_parts.append(f"{step_name} ({_routing_repr(router.owner)})")
    else:
        display_parts.append(_routing_repr(router.owner))

    # Collect parameters for separate printing (one per line)
    param_strs = []
    if has_params:
        for param in sorted(node_info["params"]):
            param_str = _format_param_with_status(
                param,
                node_info["statuses"],
                node_info.get("aliases", {}),
            )
            param_strs.append(param_str)

    display_line = "".join(display_parts)
    print(f"{prefix}{connector}{display_line}")

    # Print each parameter on its own indented line
    if param_strs:
        param_prefix = _param_prefix(prefix, root, is_last)
        for p in param_strs:
            print(f"{param_prefix}➤ {p}")

    # Process children
    if isinstance(router, MetadataRouter):
        children = list(router._route_mappings.items())
        new_prefix = _child_prefix(prefix, root, is_last)

        for i, (name, mapping) in enumerate(children):
            is_last_child = i == len(children) - 1

            # Pass the *mapping name* on the path so that it matches the key
            # format used in `_collect_routing_info`.
            _display_tree(
                mapping.router,
                routing_map,
                new_prefix,
                is_last_child,
                name,
                f"{current_path}/{name}",  # parent_path for child
                root=False,
            )


def _summarise_params_by_method(router, routing_map):
    """Return a mapping ``{root_method -> user_param -> category -> [paths]}``.

    A three-level ``defaultdict`` removes the tedious nested ``setdefault``
    chains while keeping the public return type identical (converted back to a
    plain ``dict`` at the end).
    """

    from collections import defaultdict

    # summary[root][param][cat] → list[str]
    summary = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    origins = _compute_method_origins(router)

    for path, info in routing_map.items():
        for comp_param in info["params"]:
            statuses = info["statuses"].get(comp_param, {})

            for method, status in statuses.items():
                # Determine user-facing name: only *real* aliases count.
                user_param = status if request_is_alias(status) else comp_param

                cat = _status_category(status)
                if cat == "unknown":
                    continue

                roots = origins.get(f"{path}.{method}", {method.split(".")[0]})
                for root_method in roots:
                    summary[root_method][user_param][cat].append(f"{path}.{method}")

    # Convert back to plain dicts for stable downstream expectations.
    return {
        r: {p: dict(cats) for p, cats in params.items()}
        for r, params in summary.items()
    }


# -----------------------------------------------------------------------------
# Formatting helpers
# -----------------------------------------------------------------------------


def _shorten_path(path_method: str) -> str:
    """Remove estimator class tokens from a path for display purposes."""
    if "/" not in path_method:
        return path_method

    try:
        path, method = path_method.rsplit(".", 1)
    except ValueError:
        # no method part
        return path_method

    tokens = path.split("/")
    if len(tokens) <= 2:
        # nothing to shorten
        return path_method

    root = tokens[0]
    rest = tokens[1:]
    # Keep every other token starting from index0 of *rest* (mapping names).
    filtered = [root] + rest[::2]
    short_path = "/".join(filtered)
    return f"{short_path}.{method}"
