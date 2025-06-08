# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from collections import defaultdict

from sklearn import set_config
from sklearn.compose import ColumnTransformer
from sklearn.feature_selection import SelectPercentile, chi2
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import get_scorer
from sklearn.model_selection import GroupKFold, RandomizedSearchCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.utils._metadata_requests import (
    COMPOSITE_METHODS,
    SIMPLE_METHODS,
    UNUSED,
    WARN,
    MetadataRequest,
    MetadataRouter,
)
from sklearn.utils.metadata_routing import get_routing_for_object

# Enable metadata routing
set_config(enable_metadata_routing=True)


def _param_names(router):
    res = set()
    for method in SIMPLE_METHODS:
        res = res.union(
            router._get_param_names(
                method=method, return_alias=False, ignore_self_request=False
            )
        )
    return res


# ============================================================================
# NEW VISUALIZATION FUNCTIONS
# ============================================================================


def visualise_routing(routing_info, show_method_mappings=False, show_all_metadata=True):
    """
    Visualize metadata routing information.

    Parameters
    ----------
    routing_info : MetadataRouter or MetadataRequest
        The routing information to visualize
    show_method_mappings : bool, default=False
        Whether to show method mappings (e.g., fit→fit_transform) in the output.
    show_all_metadata : bool, default=True
        Whether to show all possible metadata parameters with status indicators,
        or only the parameters that are explicitly requested.
    """
    visualise_tree(routing_info, show_method_mappings, show_all_metadata)


def visualise_tree(routing_info, show_method_mappings=False, show_all_metadata=True):
    """Show a consolidated tree view with parameters annotated inline."""
    print("\n=== METADATA ROUTING TREE ===")

    # Get all routing information
    routing_map = _collect_routing_info(
        routing_info, show_all_metadata=show_all_metadata
    )

    # Display tree structure
    print(f"Root: {routing_info.owner}")
    print("│")
    _display_tree_new(
        routing_info,
        routing_map,
        show_method_mappings=show_method_mappings,
        show_all_metadata=show_all_metadata,
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
                _walk(child_obj, child_methods, current_path)

    _walk(router, set(SIMPLE_METHODS))
    return reachable


def _collect_routing_info(router, top_router=None, show_all_metadata=True):
    """Collect all routing information into a structured format."""
    if top_router is None:
        top_router = router

    # Compute reachability map so we can ignore unreachable methods.
    reachable_map = _compute_reachable_methods(top_router)

    info = defaultdict(
        lambda: {
            "params": set(),
            "methods": defaultdict(set),
            "aliases": {},
            "statuses": defaultdict(dict),
            # Track which methods actually exist for this component
            "existing_methods": set(),
        }
    )

    def _collect(obj, path=""):
        # Build the current path correctly
        if path:
            current_path = f"{path}/{obj.owner}"
        else:
            current_path = obj.owner

        # Fetch reachable methods for *this* component (may be empty)
        reachable_here = reachable_map.get(current_path, set())

        if isinstance(obj, MetadataRequest):
            # Get detailed info about each method's requests
            for method in SIMPLE_METHODS:
                # Skip if method not reachable or not implemented
                if method not in reachable_here or not hasattr(obj, method):
                    continue

                info[current_path]["existing_methods"].add(method)
                method_req = getattr(obj, method)

                if show_all_metadata:
                    # Get all possible metadata parameters for this estimator
                    all_possible_params = _param_names(obj)
                    for param in all_possible_params:
                        info[current_path]["params"].add(param)

                        # Determine alias/status for this method
                        if param in method_req.requests:
                            alias = method_req.requests[param]
                        else:
                            alias = False  # reachable but not requested

                        info[current_path]["statuses"][param][method] = alias

                        # Only add to requested-methods set if actually requested
                        if (
                            alias is not False
                            and alias != WARN
                            and alias is not None
                            and alias != UNUSED
                        ):
                            info[current_path]["methods"][param].add(method)

                        # Track alias relationships
                        if isinstance(alias, str) and alias != param:
                            info[current_path]["aliases"][param] = alias
                else:
                    # Original behavior: only show explicitly set parameters
                    for param, alias in method_req.requests.items():
                        # Include ALL parameters regardless of status
                        info[current_path]["params"].add(param)
                        # Record status only if this method actually defines the
                        # parameter
                        if param in method_req.requests:
                            alias = method_req.requests[param]
                            info[current_path]["statuses"][param][method] = alias

                            # Only add to requested-methods set if it is actively
                            # requested
                            if (
                                alias is not False
                                and alias != WARN
                                and alias is not None
                                and alias != UNUSED
                            ):
                                info[current_path]["methods"][param].add(method)

                            # Track alias relationships
                            if isinstance(alias, str) and alias != param:
                                info[current_path]["aliases"][param] = alias

        elif isinstance(obj, MetadataRouter):
            if obj._self_request:
                for method in SIMPLE_METHODS:
                    # Skip if method not reachable or not implemented
                    if method not in reachable_here or not hasattr(
                        obj._self_request, method
                    ):
                        continue

                    info[current_path]["existing_methods"].add(method)
                    method_req = getattr(obj._self_request, method)

                    if show_all_metadata:
                        # Get all possible metadata parameters for this estimator
                        all_possible_params = _param_names(obj._self_request)
                        for param in all_possible_params:
                            info[current_path]["params"].add(param)

                            if param in method_req.requests:
                                alias = method_req.requests[param]
                            else:
                                alias = False

                            info[current_path]["statuses"][param][method] = alias

                            if (
                                alias is not False
                                and alias != WARN
                                and alias is not None
                                and alias != UNUSED
                            ):
                                info[current_path]["methods"][param].add(method)

                            if isinstance(alias, str) and alias != param:
                                info[current_path]["aliases"][param] = alias
                    else:
                        # Original behavior: only show explicitly set parameters
                        for param, alias in method_req.requests.items():
                            # Include ALL parameters regardless of status
                            info[current_path]["params"].add(param)
                            # Only record status for methods that define the parameter
                            if param in method_req.requests:
                                alias = method_req.requests[param]
                                info[current_path]["statuses"][param][method] = alias

                                if (
                                    alias is not False
                                    and alias != WARN
                                    and alias is not None
                                    and alias != UNUSED
                                ):
                                    info[current_path]["methods"][param].add(method)

                                if isinstance(alias, str) and alias != param:
                                    # Store the mapping from component param to user
                                    # param for display
                                    info[current_path]["aliases"][param] = alias

            for name, mapping in obj._route_mappings.items():
                _collect(mapping.router, current_path)

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


def _format_param_with_status(
    param, methods, statuses, aliases, show_all_metadata=True
):
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

    alias_methods = {m: s for m, s in param_statuses.items() if isinstance(s, str)}
    non_alias_methods = {
        m: s for m, s in param_statuses.items() if not isinstance(s, str)
    }

    parts = []

    # Helper to build the method part respecting the `show_all_metadata` flag.
    def _build_method_part(method_map):
        method_parts = []
        if show_all_metadata:
            # Show all methods but only include indicators for those present in
            # the supplied map.
            for method in SIMPLE_METHODS:
                if method in method_map:
                    status = method_map[method]
                    indicator = _get_status_indicator(status)
                    method_parts.append(f"{method}{indicator}")
        else:
            for method, status in method_map.items():
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
    if non_alias_methods or (show_all_metadata and not alias_methods):
        method_parts = _build_method_part(non_alias_methods)
        if method_parts:
            parts.append(f"{param}[{','.join(method_parts)}]")
        else:
            parts.append(param)

    # Join the parts with ", " to mimic the desired output style.
    return ", ".join(parts)


def _get_original_params(router):
    """Get all original parameter names before aliasing."""
    params = set()

    def _search(obj):
        if isinstance(obj, MetadataRequest):
            for method in SIMPLE_METHODS:
                method_request = getattr(obj, method)
                for param in method_request.requests:
                    params.add(param)
        elif isinstance(obj, MetadataRouter):
            if obj._self_request:
                for method in SIMPLE_METHODS:
                    method_request = getattr(obj._self_request, method)
                    for param in method_request.requests:
                        params.add(param)
            for mapping in obj._route_mappings.values():
                _search(mapping.router)

    _search(router)
    return params


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


def _display_tree_new(
    router,
    routing_map,
    prefix="",
    is_last=True,
    show_method_mappings=False,
    show_all_metadata=True,
    step_name=None,
    parent_path="",
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
                show_all_metadata=show_all_metadata,
            )
            param_strs.append(param_str)

    display_line = "".join(display_parts)
    print(f"{prefix}{connector}{display_line}")

    # Print each parameter on its own indented line
    if param_strs:
        # Use spaces after the branch connector to avoid an unnecessary vertical bar
        param_prefix = prefix + ("    " if is_last else "│   ") + "    "
        for p in param_strs:
            print(f"{param_prefix}➤ {p}")

    # Show method mappings if requested
    if show_method_mappings and isinstance(router, MetadataRouter):
        extension = "    " if is_last else "│   "
        new_prefix = prefix + extension

        # Collect all method mappings
        all_mappings = []
        for mapping in router._route_mappings.values():
            for method_pair in mapping.mapping:
                all_mappings.append(f"{method_pair.caller} → {method_pair.callee}")

        if all_mappings:
            print(f"{new_prefix}├─ Method mappings:")
            for i, mapping_str in enumerate(sorted(all_mappings)):
                mapping_connector = "└─" if i == len(all_mappings) - 1 else "├─"
                print(f"{new_prefix}│ {mapping_connector} {mapping_str}")

    # Process children
    if isinstance(router, MetadataRouter):
        children = list(router._route_mappings.items())
        extension = "    " if is_last else "│   "
        new_prefix = prefix + extension

        for i, (name, mapping) in enumerate(children):
            is_last_child = i == len(children) - 1
            full_path = f"{current_path}/{mapping.router.owner}"

            # Display child with step name
            _display_tree_new(
                mapping.router,
                routing_map,
                new_prefix,
                is_last_child,
                show_method_mappings,
                show_all_metadata,
                name,
                current_path,
            )


def _collect_all_components(router, prefix=""):
    """Collect all component names in the routing hierarchy."""
    components = []

    if isinstance(router, MetadataRequest):
        if router.owner:
            components.append(router.owner)
    elif isinstance(router, MetadataRouter):
        if router.owner:
            components.append(router.owner)

        for name, mapping in router._route_mappings.items():
            sub_comps = _collect_all_components(mapping.router, f"{prefix}{name}.")
            components.extend(sub_comps)

    return components


def _find_consumers(router, param):
    """Find all components that consume a parameter."""
    consumers = []

    def _search(obj):
        if isinstance(obj, MetadataRequest):
            for method in SIMPLE_METHODS:
                if obj.consumes(method, [param]):
                    consumers.append(obj.owner)
                    break

        elif isinstance(obj, MetadataRouter):
            if obj._self_request:
                for method in SIMPLE_METHODS:
                    if obj._self_request.consumes(method, [param]):
                        consumers.append(obj.owner)
                        break

            for mapping in obj._route_mappings.values():
                _search(mapping.router)

    _search(router)
    return consumers


def _find_consumers_with_methods(router, param):
    """Find all components that consume a parameter with their methods."""
    consumers = defaultdict(set)

    def _search(obj):
        if isinstance(obj, MetadataRequest):
            for method in SIMPLE_METHODS:
                if obj.consumes(method, [param]):
                    consumers[obj.owner].add(method)

        elif isinstance(obj, MetadataRouter):
            if obj._self_request:
                for method in SIMPLE_METHODS:
                    if obj._self_request.consumes(method, [param]):
                        consumers[obj.owner].add(method)

            for mapping in obj._route_mappings.values():
                _search(mapping.router)

    _search(router)
    return consumers


def _build_compact_structure(
    router, level=0, top_router=None, parent_path="", show_all_metadata=True
):
    """Build a compact structure representation."""
    if top_router is None:
        top_router = router

    # Build current path
    if parent_path:
        current_path = f"{parent_path}/{router.owner}"
    else:
        current_path = router.owner

    structure = {
        "name": router.owner,
        "type": type(router).__name__,
        "level": level,
        "path": current_path,
        "children": [],
    }

    # Add routing info
    routing_map = _collect_routing_info(top_router, show_all_metadata=show_all_metadata)
    if current_path in routing_map:
        structure["params"] = routing_map[current_path]["params"]
        structure["methods"] = routing_map[current_path]["methods"]

    if isinstance(router, MetadataRouter):
        for name, mapping in router._route_mappings.items():
            child_structure = _build_compact_structure(
                mapping.router, level + 1, top_router, current_path, show_all_metadata
            )
            child_structure["mapping_name"] = name

            # Collect method mappings if any
            mappings = []
            for method_pair in mapping.mapping:
                if method_pair.caller != method_pair.callee:
                    mappings.append(f"{method_pair.caller} → {method_pair.callee}")
            if mappings:
                child_structure["mappings"] = mappings

            structure["children"].append(child_structure)

    return structure


def _print_compact_structure(
    structure,
    indent=0,
    show_method_mappings=False,
    show_all_metadata=True,
    mapping_info=None,
    routing_map=None,
):
    """Print the compact structure."""
    prefix = "  " * indent

    # For child nodes with mappings, combine on one line
    if mapping_info and indent > 0:
        # This is a mapped component - show mapping name and component together
        name_part = f"↳ via {mapping_info['name']}: {structure['name']}"

        # Prepare parameter strings if any
        param_strs = []
        if structure.get("params"):
            for param in sorted(structure["params"]):
                if routing_map and structure["path"] in routing_map:
                    param_str = _format_param_with_status(
                        param,
                        routing_map[structure["path"]]["methods"][param],
                        routing_map[structure["path"]]["statuses"],
                        routing_map[structure["path"]].get("aliases", {}),
                        show_all_metadata=show_all_metadata,
                    )
                else:
                    # Fallback if routing_map not available
                    param_str = param
                param_strs.append(param_str)

        # First, print the component line
        print(f"{prefix}{name_part}")

        # Then, print each parameter on its own indented line
        if param_strs:
            # Use spaces after the branch connector to avoid an unnecessary vertical bar
            param_prefix = prefix + "    "
            for ps in param_strs:
                print(f"{param_prefix}◆ {ps}")

        # Show method mappings if requested
        if show_method_mappings and mapping_info.get("mappings"):
            for mapping in sorted(mapping_info["mappings"]):
                print(f"{prefix}    └─ {mapping}")
    else:
        # Root node or no mapping
        print(f"{prefix}▸ {structure['name']}")

        # Show parameters for root
        if structure.get("params") and indent == 0:
            param_strs = []
            for param in sorted(structure["params"]):
                if routing_map and structure["path"] in routing_map:
                    param_str = _format_param_with_status(
                        param,
                        routing_map[structure["path"]]["methods"][param],
                        routing_map[structure["path"]]["statuses"],
                        routing_map[structure["path"]].get("aliases", {}),
                        show_all_metadata=show_all_metadata,
                    )
                else:
                    # Fallback if routing_map not available
                    param_str = param
                param_strs.append(param_str)
            print(f"{prefix}  Parameters: {', '.join(param_strs)}")

    # Process children
    for child in structure.get("children", []):
        mapping_info = {
            "name": child["mapping_name"],
            "mappings": child.get("mappings", []),
        }
        _print_compact_structure(
            child,
            indent + 2,
            show_method_mappings,
            show_all_metadata,
            mapping_info,
            routing_map,
        )


def _find_routing_paths(router, param):
    """Find all routing paths for a parameter."""
    paths = []

    def _trace_path(obj, current_path):
        if isinstance(obj, MetadataRequest):
            for method in SIMPLE_METHODS:
                if obj.consumes(method, [param]):
                    paths.append(current_path + [f"{obj.owner}.{method}"])

        elif isinstance(obj, MetadataRouter):
            current_path.append(obj.owner)

            if obj._self_request:
                for method in SIMPLE_METHODS:
                    if obj._self_request.consumes(method, [param]):
                        paths.append(current_path + [f"self.{method}"])

            for name, mapping in obj._route_mappings.items():
                for m in mapping.mapping:
                    if mapping.router.consumes(m.callee, [param]):
                        path_copy = current_path.copy()
                        path_copy.append(f"{m.caller} → {name}.{m.callee}")
                        _trace_path(mapping.router, path_copy)

    _trace_path(router, [])
    return paths


def _build_flow_diagram(router, show_method_mappings=False, show_all_metadata=True):
    """Build an ASCII flow diagram."""
    lines = []

    # Header
    lines.append("┌" + "─" * 50 + "┐")
    lines.append(f"│ {router.owner:^48} │")
    lines.append("└" + "─" * 25 + "┬" + "─" * 25 + "┘")

    # Collect all routing information
    routing_map = _collect_routing_info(router)
    original_params = sorted(_get_original_params(router))

    if original_params:
        lines.append(" " * 26 + "│")
        lines.append(" " * 20 + f"Parameters: {', '.join(original_params)}")
        lines.append(" " * 26 + "│")

    # Build the flow structure recursively
    def _add_flow_node(obj, indent=0, prefix="", is_last=True, parent_path=""):
        # Build current path
        if parent_path:
            current_path = f"{parent_path}/{obj.owner}"
        else:
            current_path = obj.owner

        if isinstance(obj, MetadataRouter):
            # Add connector for non-root nodes
            if indent > 0:
                connector = "└──" if is_last else "├──"
                lines.append(f"{prefix}{connector}─┐")
                lines.append(f"{prefix}{'   ' if is_last else '│  '} │")

            # Add component box
            comp_prefix = prefix + ("   " if is_last else "│  ") if indent > 0 else ""

            # Get consumed parameters for this component from routing map
            param_info = ""
            if current_path in routing_map:
                info = routing_map[current_path]
                if info["params"]:
                    param_details = []
                    for param in sorted(info["params"]):
                        param_str = _format_param_with_status(
                            param,
                            info["methods"][param],
                            info["statuses"],
                            info.get("aliases", {}),
                            show_all_metadata=show_all_metadata,
                        )
                        param_details.append(param_str)
                    param_info = f" ◆ {', '.join(param_details)}"

            # Draw component box
            box_content = obj.owner + param_info
            box_width = max(40, len(box_content) + 4)
            lines.append(f"{comp_prefix}┌" + "─" * box_width + "┐")
            lines.append(f"{comp_prefix}│ {box_content:<{box_width - 2}} │")
            lines.append(f"{comp_prefix}└" + "─" * box_width + "┘")

            # Process children
            children = list(obj._route_mappings.items())
            if children:
                lines.append(f"{comp_prefix}" + " " * (box_width // 2) + "│")

                for i, (name, mapping) in enumerate(children):
                    is_last_child = i == len(children) - 1

                    # Show method mapping
                    method_info = ""
                    if show_method_mappings and getattr(mapping, "mapping", []):
                        unique_maps = set()
                        for m in mapping.mapping:
                            if m.caller != m.callee:
                                unique_maps.add(f"{m.caller}→{m.callee}")
                        if unique_maps:
                            method_info = f" [{', '.join(sorted(unique_maps))}]"

                    # Add branch
                    branch_char = "└" if is_last_child else "├"
                    lines.append(
                        f"{comp_prefix}"
                        + " " * (box_width // 2)
                        + f"{branch_char}── {name}{method_info}"
                    )

                    # Recurse
                    new_prefix = (
                        comp_prefix
                        + " " * (box_width // 2)
                        + ("    " if is_last_child else "│   ")
                    )
                    _add_flow_node(
                        mapping.router, indent + 1, new_prefix, True, current_path
                    )

        elif isinstance(obj, MetadataRequest):
            # Terminal node - show consumed parameters from routing map
            param_str = ""
            if current_path in routing_map:
                info = routing_map[current_path]
                if info["params"]:
                    param_details = []
                    for param in sorted(info["params"]):
                        param_str = _format_param_with_status(
                            param,
                            info["methods"][param],
                            info["statuses"],
                            info.get("aliases", {}),
                            show_all_metadata=show_all_metadata,
                        )
                        param_details.append(param_str)
                    param_str = f" ◆ {', '.join(param_details)}"

            # Show component with or without parameters
            lines.append(f"{prefix}    ◆ {obj.owner}{param_str}")

    # Start building the flow
    _add_flow_node(router)

    # Add parameter flow summary at the bottom with aliases
    if original_params:
        lines.append("")
        lines.append("─" * 60)
        lines.append("Parameter Flow Summary:")
        for param in original_params:
            consumers = []
            # Find all consumers and their aliases
            for comp_path, info in routing_map.items():
                if param in info["params"]:
                    comp_name = comp_path.split("/")[-1]
                    if param in info.get("aliases", {}):
                        alias = info["aliases"][param]
                        consumers.append(f"{comp_name} (as {alias})")
                    else:
                        consumers.append(comp_name)

            if consumers:
                lines.append(f"  • {param} → {', '.join(consumers)}")

    return lines


# ============================================================================
# TEST CASES
# ============================================================================

numeric_features = ["age", "fare"]
numeric_transformer = Pipeline(
    steps=[
        ("imputer", SimpleImputer(strategy="median")),
        (
            "scaler",
            StandardScaler()
            .set_fit_request(sample_weight="inner_weights")
            .set_transform_request(copy=True),
        ),
    ]
)

categorical_features = ["embarked", "sex", "pclass"]
categorical_transformer = Pipeline(
    steps=[
        ("encoder", OneHotEncoder(handle_unknown="ignore")),
        ("selector", SelectPercentile(chi2, percentile=50)),
    ]
)
preprocessor = ColumnTransformer(
    transformers=[
        ("num", numeric_transformer, numeric_features),
        ("cat", categorical_transformer, categorical_features),
    ]
)

# %%
# Append classifier to preprocessing pipeline.
# Now we have a full prediction pipeline.
clf = Pipeline(
    steps=[
        ("preprocessor", preprocessor),
        ("classifier", LogisticRegression().set_fit_request(sample_weight=True)),
    ]
)

param_grid = {
    "preprocessor__num__imputer__strategy": ["mean", "median"],
    "preprocessor__cat__selector__percentile": [10, 30, 50, 70],
    "classifier__C": [0.1, 1.0, 10, 100],
}

scorer = get_scorer("accuracy").set_score_request(sample_weight=True)

search_cv = RandomizedSearchCV(
    clf, param_grid, cv=GroupKFold(), scoring=scorer, random_state=0
)


# Get the routing information
test = get_routing_for_object(search_cv)

visualise_routing(test, show_method_mappings=False, show_all_metadata=True)
