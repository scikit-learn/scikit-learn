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
    request_is_alias,
)


def visualise_routing(routing_info):
    """
    Visualize metadata routing information.

    This diagnostic always displays *all* possible metadata parameters with
    their status indicators. The previous ``show_all_metadata`` toggle has been
    removed – comprehensive output is now the default (and only) mode.
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

    # Get user-facing parameters (including aliases)
    user_params = sorted(_get_user_facing_params(routing_map))

    if user_params:
        print("\nParameters:")
        for param in user_params:
            consumers = []

            # Find all consumers for this user-facing parameter
            for comp_path, info in routing_map.items():
                comp_name = comp_path.split("/")[-1]

                # Check if this component consumes the param directly (no alias)
                if param in info["params"] and param not in info.get("aliases", {}):
                    consumers.append(comp_name)

                # Check if this component consumes the param via alias
                elif param in info.get("aliases", {}).values():
                    # Find the component parameter that maps to this user param
                    for comp_param, alias in info["aliases"].items():
                        if alias == param:
                            consumers.append(f"{comp_name} (as {comp_param})")

            if consumers:
                print(f"  • {param} → {', '.join(consumers)}")
    else:
        print("\nNo parameters are being routed.")

    # Print parameter summary
    summary = _summarise_params(routing_map)
    if summary:
        print("\nParameter summary:")

        glyph_for = dict(_CATEGORY_ORDER)

        for param in sorted(summary):
            cats = summary[param]

            # Choose leading glyph based on first category present in order.
            leading_glyph = next(glyph for cat, glyph in _CATEGORY_ORDER if cat in cats)

            print(f"{leading_glyph} {param}")

            for cat, glyph in _CATEGORY_ORDER:
                if cat in cats:
                    paths = ", ".join(cats[cat])
                    print(f"    • {glyph} {cat}: {paths}")

    else:
        print("\nNo parameter summary.")


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================


def _compute_reachable_methods(router):
    """Return a mapping ``{component_path: {methods}}`` of *reachable* methods.

    A method on a component is considered reachable if there exists a chain of
    ``MethodMapping`` objects starting from *any* root‐level simple method that
    results in that method being invoked.  This respects the caller → callee
    relationships stored in each :class:`~sklearn.utils.metadata_routing.MethodMapping`.
    """

    def _expand_methods(methods):
        """Return a set with composite methods broken into their simple parts."""
        expanded = set()
        for m in methods:
            if m in COMPOSITE_METHODS:
                expanded.update(COMPOSITE_METHODS[m])
            expanded.add(m)
        return expanded

    reachable = defaultdict(set)

    def _walk(obj, incoming_methods, path=""):
        # Build path in the exact same way as in _collect_routing_info so that
        # look-ups can be shared.
        current_path = f"{path}/{obj.owner}" if path else obj.owner

        incoming_methods = _expand_methods(incoming_methods)

        if isinstance(obj, MetadataRequest):
            for method in incoming_methods:
                if hasattr(obj, method):
                    reachable[current_path].add(method)
            # No children – stop.
            return

        # obj is a MetadataRouter -------------------------------------------
        # 1. Handle potential self-request (router acts as a consumer).
        if obj._self_request:
            for method in incoming_methods:
                if hasattr(obj._self_request, method):
                    reachable[current_path].add(method)

        # 2. Propagate reachability to children via method mappings.
        for name, mapping_pair in obj._route_mappings.items():
            child_obj = mapping_pair.router
            # Collect child methods reachable through *any* caller that is itself
            # reachable on *this* router.
            child_methods_raw = {
                pair.callee
                for pair in mapping_pair.mapping
                if pair.caller in incoming_methods
            }
            child_methods = _expand_methods(child_methods_raw)
            if child_methods:
                # Propagate the *mapping name* into the path so that it stays
                # in-sync with the scheme used in `_collect_routing_info`.
                _walk(child_obj, child_methods, f"{current_path}/{name}")

    _walk(router, set(SIMPLE_METHODS))
    return reachable


def _collect_routing_info(router, top_router=None):
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
    if top_router is None:
        top_router = router

    reachable_map = _compute_reachable_methods(top_router)

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

    def _record_status(current_path: str, param: str, method: str, alias):
        """Populate *info* dictionaries in a single place."""

        # Record the raw status so we can later derive visual indicators.
        info[current_path]["statuses"][param][method] = alias

        # Only count the method as *requesting* the parameter when the status
        # truly indicates a request (True) or a *user* alias. Special markers
        # such as WARN/UNUSED (*do not* constitute a request) and must be
        # skipped, otherwise they would incorrectly show up as "requested".
        if alias is True or request_is_alias(alias):
            info[current_path]["methods"][param].add(method)

        # Maintain the mapping {component_param -> user_alias}. Again, ignore
        # special placeholders such as WARN/UNUSED which are *not* aliases.
        if request_is_alias(alias) and alias != param:
            info[current_path]["aliases"][param] = alias

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
            # Include the *mapping name* in the path so that identical child
            # objects appearing in different locations (e.g. several
            # ``StandardScaler`` instances) do not share the same key in the
            # routing map.  This fixes cases where their individual requests
            # get merged and displayed incorrectly.
            _collect(
                mapping.router,
                f"{current_path}/{name}",  # becomes *path* for child
            )

    # ------------------------------------------------------------------
    # Main DFS traversal
    # ------------------------------------------------------------------

    def _collect(obj, path=""):
        current_path = f"{path}/{obj.owner}" if path else obj.owner
        reachable_here = reachable_map.get(current_path, set())

        if isinstance(obj, MetadataRequest):
            _collect_request_info(obj, current_path, reachable_here)
        elif isinstance(obj, MetadataRouter):
            _collect_router_info(obj, current_path, reachable_here)

    _collect(router)
    return info


def _get_status_indicator(alias):
    """Get a visual indicator for parameter status."""
    if alias is True:
        return "✓"  # requested
    elif alias is False:
        return "✗"  # not requested
    elif alias is None:
        return "⛔"  # error if passed (unrequested)
    elif alias == WARN:
        return "⚠"  # warn status
    elif alias == UNUSED:
        return "⊘"  # unused status
    elif isinstance(alias, str):
        return "↗"  # alias
    else:
        return "?"  # unknown status


def _format_param_with_status(param, methods, statuses, aliases):
    """Format a parameter with its status indicators and methods."""
    if param in aliases:
        alias = aliases[param]
        param_display = f"{alias}→{param}"
    else:
        param_display = param

    # Split methods by whether they use an alias (str) or not. This allows us
    # to display alias relationships only for the methods that actually use
    # them. We need to generate *two* entries when a parameter is aliased in
    # some methods but not others, e.g. ``inner_weights→sample_weight[fit↗]``
    # and ``sample_weight[partial_fit⚠]``.

    # Collect the status mapping once for convenience
    param_statuses = statuses.get(param, {})

    alias_methods = {m: s for m, s in param_statuses.items() if request_is_alias(s)}
    non_alias_methods = {
        m: s for m, s in param_statuses.items() if not request_is_alias(s)
    }

    parts = []

    # Helper to build the method part respecting the `show_all_metadata` flag.
    def _build_method_part(method_map):
        method_parts = []
        # Always in full-metadata mode now.
        for method in SIMPLE_METHODS:
            if method in method_map:
                status = method_map[method]
                indicator = _get_status_indicator(status)
                method_parts.append(f"{method}{indicator}")
        return method_parts

    # 1. Alias part (where status is a str)
    if alias_methods:
        alias = next(iter(alias_methods.values()))  # They should all be the same
        method_parts = _build_method_part(alias_methods)
        alias_display = f"{alias}→{param}"
        if method_parts:
            parts.append(f"{alias_display}[{','.join(method_parts)}]")
        else:
            parts.append(alias_display)

    # 2. Non-alias part
    if non_alias_methods or (not alias_methods):
        method_parts = _build_method_part(non_alias_methods)
        if method_parts:
            parts.append(f"{param}[{','.join(method_parts)}]")
        else:
            parts.append(param)

    # Join the parts with ", " to mimic the desired output style.
    return ", ".join(parts)


def _get_user_facing_params(routing_map):
    """Get all user-facing parameter names (including aliases)."""
    user_params = set()

    for comp_path, info in routing_map.items():
        # Add all component parameters that don't have aliases (they are user-facing)
        for param in info["params"]:
            if param not in info.get("aliases", {}):
                user_params.add(param)

        # Add all alias names (these are the user-facing names for aliased params)
        for alias in info.get("aliases", {}).values():
            user_params.add(alias)

    return user_params


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
        current_path = f"{parent_path}/{router.owner}"
    else:
        current_path = router.owner

    # Determine if this node has parameters
    node_info = routing_map.get(current_path, {})
    has_params = bool(node_info.get("params"))

    # Build the display line
    # For the root node we don't want branch connectors
    if root and prefix == "":
        connector = ""
    else:
        connector = "└── " if is_last else "├── "

    display_parts = []
    if step_name:
        display_parts.append(f"{step_name} ({router.owner})")
    else:
        display_parts.append(router.owner)

    # Collect parameters for separate printing (one per line)
    param_strs = []
    if has_params:
        for param in sorted(node_info["params"]):
            param_str = _format_param_with_status(
                param,
                node_info["methods"][param],
                node_info["statuses"],
                node_info.get("aliases", {}),
            )
            param_strs.append(param_str)

    display_line = "".join(display_parts)
    print(f"{prefix}{connector}{display_line}")

    # Print each parameter on its own indented line
    if param_strs:
        # Special-case root to avoid surplus indentation
        if root and prefix == "":
            param_prefix = "    "
        else:
            param_prefix = prefix + ("    " if is_last else "│   ") + "    "
        for p in param_strs:
            print(f"{param_prefix}➤ {p}")

    # Process children
    if isinstance(router, MetadataRouter):
        children = list(router._route_mappings.items())
        extension = "    " if (is_last or root and prefix == "") else "│   "
        new_prefix = prefix + ("" if root and prefix == "" else extension)

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


_CATEGORY_ORDER = [
    ("requested", "✓"),
    ("ignored", "✗"),
    ("warns", "⚠"),
    ("errors", "⛔"),
    ("unused", "⊘"),
]


def _summarise_params(routing_map):
    """Aggregate parameter statuses across the whole tree.

    Returns
    -------
    dict
        { user_param -> {category -> [path.method, ...]} }
    """
    summary = {}

    for path, info in routing_map.items():
        for comp_param in info["params"]:
            statuses = info["statuses"].get(comp_param, {})

            for method, status in statuses.items():
                # Determine user-facing name
                user_param = status if isinstance(status, str) else comp_param

                # Classify category
                if status is True or isinstance(status, str):
                    cat = "requested"
                elif status is False:
                    cat = "ignored"
                elif status == WARN:
                    cat = "warns"
                elif status is None:
                    cat = "errors"
                elif status == UNUSED:
                    cat = "unused"
                else:
                    continue

                summary.setdefault(user_param, {}).setdefault(cat, []).append(
                    f"{path}.{method}"
                )

    return summary
