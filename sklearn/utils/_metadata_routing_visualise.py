# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from collections import defaultdict
from io import StringIO

from sklearn import set_config
from sklearn.compose import ColumnTransformer
from sklearn.feature_selection import SelectPercentile, chi2
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold, RandomizedSearchCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.utils._metadata_requests import (
    SIMPLE_METHODS,
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
                method=method, return_alias=True, ignore_self_request=False
            )
        )
    return res


# ============================================================================
# NEW VISUALIZATION FUNCTIONS
# ============================================================================


def visualise_routing(routing_info, view="tree", show_method_mappings=False):
    """
    Visualize metadata routing information.

    Parameters
    ----------
    routing_info : MetadataRouter or MetadataRequest
        The routing information to visualize
    view : str, default='tree'
        The visualization style. Options: 'tree', 'matrix', 'compact', 'flow'
    show_method_mappings : bool, default=False
        Whether to show method mappings (e.g., fit→fit_transform) in the output.
    """
    if view == "tree":
        visualise_tree(routing_info, show_method_mappings)
    elif view == "matrix":
        visualise_matrix(routing_info)
    elif view == "compact":
        visualise_compact(routing_info, show_method_mappings)
    elif view == "flow":
        visualise_flow(routing_info, show_method_mappings)
    else:
        # Default to old behavior for backward compatibility
        visualise_routing_old(routing_info)


def visualise_tree(routing_info, show_method_mappings=False):
    """Show a consolidated tree view with parameters annotated inline."""
    print("\n=== METADATA ROUTING TREE ===")

    # Get all routing information
    routing_map = _collect_routing_info(routing_info)

    # Display tree structure
    print(f"Root: {routing_info.owner}")
    print("│")
    _display_tree_new(
        routing_info, routing_map, show_method_mappings=show_method_mappings
    )

    # Get original parameters (before aliasing)
    original_params = sorted(_get_original_params(routing_info))

    if original_params:
        print("\nParameters:")
        for param in original_params:
            consumers = []
            aliases_info = []

            # Find all consumers and their aliases
            for comp_path, info in routing_map.items():
                if param in info["params"]:
                    comp_name = comp_path.split("/")[-1]
                    if param in info.get("aliases", {}):
                        alias = info["aliases"][param]
                        consumers.append(f"{comp_name} (as {alias})")
                        aliases_info.append((comp_name, alias))
                    else:
                        consumers.append(comp_name)

            if consumers:
                print(f"  • {param} → {', '.join(consumers)}")
    else:
        print("\nNo parameters are being routed.")


def visualise_matrix(routing_info):
    """Show a matrix view of components vs parameters."""
    print("\n=== METADATA ROUTING MATRIX ===")

    # Collect all components and original parameters
    components = _collect_all_components(routing_info)
    params = sorted(_get_original_params(routing_info))

    if not params:
        print("No parameters are being routed.")
        return

    # Build the consumption matrix with alias info
    routing_map = _collect_routing_info(routing_info)
    matrix = defaultdict(lambda: defaultdict(lambda: {"methods": set(), "alias": None}))

    for comp_path, info in routing_map.items():
        comp_name = comp_path.split("/")[-1]  # Get the component name from path
        for param in info["params"]:
            matrix[comp_name][param]["methods"] = info["methods"][param]
            if param in info.get("aliases", {}):
                matrix[comp_name][param]["alias"] = info["aliases"][param]

    # Calculate column widths - need to check all method strings that will be displayed
    comp_width = max(len(comp) for comp in components) + 2
    param_widths = {}
    for param in params:
        # Start with parameter name length
        max_width = len(param)

        # Check all method strings that will be displayed for this parameter
        for comp in components:
            cell_info = matrix[comp][param]
            if cell_info["methods"]:
                methods = sorted(cell_info["methods"])
                method_str = ",".join(methods)
                if cell_info["alias"]:
                    display_str = f"{method_str}→{cell_info['alias']}"
                else:
                    display_str = method_str
                max_width = max(max_width, len(display_str))

        # Add padding
        param_widths[param] = max_width + 2

    # Print header
    print(f"{'Component':<{comp_width}}", end="")
    for param in params:
        print(f"{param:^{param_widths[param]}}", end="")
    print()

    # Print separator
    print("=" * comp_width, end="")
    for param in params:
        print("=" * param_widths[param], end="")
    print()

    # Print rows
    for comp in components:
        print(f"{comp:<{comp_width}}", end="")
        for param in params:
            cell_info = matrix[comp][param]
            if cell_info["methods"]:
                methods = sorted(cell_info["methods"])
                method_str = ",".join(methods)
                if cell_info["alias"]:
                    display_str = f"{method_str}→{cell_info['alias']}"
                else:
                    display_str = method_str
                print(f"{display_str:^{param_widths[param]}}", end="")
            else:
                print(f"{'-':^{param_widths[param]}}", end="")
        print()


def visualise_compact(routing_info, show_method_mappings=False):
    """Show a compact hierarchical view."""
    print("\n=== COMPACT METADATA ROUTING ===")

    # Build compact representation
    structure = _build_compact_structure(routing_info)
    _print_compact_structure(
        structure, indent=0, show_method_mappings=show_method_mappings
    )

    # Show routing paths
    original_params = sorted(_get_original_params(routing_info))
    if original_params:
        print("\nRouting Paths:")
        for param in original_params:
            paths = _find_routing_paths(routing_info, param)
            if paths:
                print(f"\n{param}:")
                for path in paths:
                    print(f"  → {' → '.join(path)}")


def visualise_flow(routing_info, show_method_mappings=False):
    """Show an ASCII flow diagram."""
    print("\n=== METADATA ROUTING FLOW ===")

    # Build flow representation
    flow = _build_flow_diagram(routing_info, show_method_mappings)

    # Display the flow
    for line in flow:
        print(line)


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================


def _collect_routing_info(router, top_router=None):
    """Collect all routing information into a structured format."""
    if top_router is None:
        top_router = router

    info = defaultdict(
        lambda: {"params": set(), "methods": defaultdict(set), "aliases": {}}
    )

    def _collect(obj, path=""):
        # Build the current path correctly
        if path:
            current_path = f"{path}/{obj.owner}"
        else:
            current_path = obj.owner

        if isinstance(obj, MetadataRequest):
            # Get detailed info about each method's requests
            for method in SIMPLE_METHODS:
                method_req = getattr(obj, method)
                for param, alias in method_req.requests.items():
                    if alias is False or alias == WARN:
                        continue
                    # Track the original parameter name
                    info[current_path]["params"].add(param)
                    info[current_path]["methods"][param].add(method)
                    # If it's an alias (string), track the mapping
                    if isinstance(alias, str) and alias != param:
                        info[current_path]["aliases"][param] = alias

        elif isinstance(obj, MetadataRouter):
            if obj._self_request:
                for method in SIMPLE_METHODS:
                    method_req = getattr(obj._self_request, method)
                    for param, alias in method_req.requests.items():
                        if alias is False or alias == WARN:
                            continue
                        info[current_path]["params"].add(param)
                        info[current_path]["methods"][param].add(method)
                        if isinstance(alias, str) and alias != param:
                            info[current_path]["aliases"][param] = alias

            for name, mapping in obj._route_mappings.items():
                _collect(mapping.router, current_path)

    _collect(router)
    return info


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


def _display_tree_new(
    router,
    routing_map,
    prefix="",
    is_last=True,
    show_method_mappings=False,
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

    # Add parameters if any
    if has_params:
        param_strs = []
        for param in sorted(node_info["params"]):
            methods = sorted(node_info["methods"][param])
            method_str = ",".join(methods)
            # Check if this parameter has an alias
            if param in node_info.get("aliases", {}):
                alias = node_info["aliases"][param]
                param_strs.append(f"{param}→{alias}[{method_str}]")
            else:
                param_strs.append(f"{param}[{method_str}]")
        display_parts.append(f"  ➤ {', '.join(param_strs)}")

    display_line = "".join(display_parts)
    print(f"{prefix}{connector}{display_line}")

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


def _build_compact_structure(router, level=0, top_router=None, parent_path=""):
    """Build a compact structure representation."""
    if top_router is None:
        top_router = router

    # Get current path
    if parent_path:
        current_path = f"{parent_path}/{router.owner}"
    else:
        current_path = router.owner

    structure = {
        "name": router.owner,
        "type": type(router).__name__,
        "params": defaultdict(set),
        "children": [],
        "path": current_path,
    }

    # Get routing info with aliases
    routing_map = _collect_routing_info(top_router)
    structure["routing_map"] = routing_map

    # Get params for this component from routing map
    if current_path in routing_map:
        info = routing_map[current_path]
        for param in info["params"]:
            structure["params"][param] = info["methods"][param]

    elif isinstance(router, MetadataRouter):
        for name, mapping in router._route_mappings.items():
            child = _build_compact_structure(
                mapping.router, level + 1, top_router, current_path
            )
            child["mapping_name"] = name

            # Collect method mappings if any
            mappings = []
            for method_pair in mapping.mapping:  # Fix attribute access
                if method_pair.caller != method_pair.callee:
                    mappings.append(f"{method_pair.caller} → {method_pair.callee}")
            if mappings:
                child["mappings"] = mappings

            structure["children"].append(child)

    return structure


def _print_compact_structure(
    structure, indent=0, show_method_mappings=False, mapping_info=None
):
    """Print the compact structure."""
    prefix = "  " * indent

    # Get routing info for aliases
    routing_map = structure.get("routing_map", {})

    # For child nodes with mappings, combine on one line
    if mapping_info and indent > 0:
        # This is a mapped component - show mapping name and component together
        name_part = f"↳ via {mapping_info['name']}: {structure['name']}"

        # Add parameters if any
        if structure["params"]:
            param_strs = []
            for param, methods in sorted(structure["params"].items()):
                method_str = ",".join(sorted(methods))
                # Check for alias in routing map
                comp_path = structure.get("path", structure["name"])
                if comp_path in routing_map and param in routing_map[comp_path].get(
                    "aliases", {}
                ):
                    alias = routing_map[comp_path]["aliases"][param]
                    param_strs.append(f"{param}→{alias}[{method_str}]")
                else:
                    param_strs.append(f"{param}[{method_str}]")
            print(f"{prefix}{name_part} ◆ {', '.join(param_strs)}")
        else:
            print(f"{prefix}{name_part}")

        # Show method mappings if requested
        if show_method_mappings and mapping_info.get("mappings"):
            for mapping in sorted(mapping_info["mappings"]):
                print(f"{prefix}    └─ {mapping}")
    else:
        # Root node or no mapping
        print(f"{prefix}▸ {structure['name']}")

        # Show parameters for root
        if structure["params"] and indent == 0:
            param_strs = []
            for param, methods in sorted(structure["params"].items()):
                method_str = ",".join(sorted(methods))
                comp_path = structure.get("path", structure["name"])
                if comp_path in routing_map and param in routing_map[comp_path].get(
                    "aliases", {}
                ):
                    alias = routing_map[comp_path]["aliases"][param]
                    param_strs.append(f"{param}→{alias}[{method_str}]")
                else:
                    param_strs.append(f"{param}[{method_str}]")
            print(f"{prefix}  Parameters: {', '.join(param_strs)}")

    # Process children
    for child in structure["children"]:
        # Check if this child has a mapping
        mapping_info = None
        if "mapping_name" in child:
            mapping_info = {
                "name": child["mapping_name"],
                "mappings": child.get("mappings", []),
            }
        _print_compact_structure(child, indent + 2, show_method_mappings, mapping_info)


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


def _build_flow_diagram(router, show_method_mappings=False):
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
                        methods = ",".join(sorted(info["methods"][param]))
                        # Show alias if exists
                        if param in info.get("aliases", {}):
                            alias = info["aliases"][param]
                            param_details.append(f"{param}→{alias}[{methods}]")
                        else:
                            param_details.append(f"{param}[{methods}]")
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
                    if show_method_mappings and getattr(
                        mapping, "mapping", []
                    ):  # Fix attribute access
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
                        methods = ",".join(sorted(info["methods"][param]))
                        # Show alias if exists
                        if param in info.get("aliases", {}):
                            alias = info["aliases"][param]
                            param_details.append(f"{param}→{alias}[{methods}]")
                        else:
                            param_details.append(f"{param}[{methods}]")
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
# OLD VISUALIZATION (for backward compatibility)
# ============================================================================


def visualise_routing_old(routing_info):
    params = _param_names(routing_info)
    for param in params:
        output = StringIO()
        print(f"Visualising routing for {param}", file=output)
        routed = _visualise_param(param, routing_info, indent=0, output=output)
        if routed:
            print(output.getvalue())


def _visualise_param(param, routing_info, indent, output):
    if isinstance(routing_info, MetadataRouter):
        return _visualise_metadata_router(param, routing_info, indent, output)
    elif isinstance(routing_info, MetadataRequest):
        return _visualise_param_in_metadata_request(param, routing_info, indent, output)
    else:
        raise ValueError(f"Unknown type {type(routing_info)}")


def _visualise_param_in_metadata_request(param, request, indent, output):
    print(" ." * indent, request.owner, file=output)
    routed = False
    for method in SIMPLE_METHODS:
        if request.consumes(method, [param]):
            print(" ." * (indent + 1), f"{method} consumes {param}", file=output)
            routed = True
    return routed


def _visualise_metadata_router(param, router, indent, output):
    routed = False
    print(" ." * indent, router.owner, file=output)
    indent += 1
    if router._self_request:
        sub_output = StringIO()
        sub_routed = _visualise_param_in_metadata_request(
            param, router._self_request, indent=indent + 1, output=sub_output
        )
        if sub_routed:
            print(" ." * indent, "Self request", file=output)
            print(sub_output.getvalue(), end="", file=output)
            routed = True

    for name, router_mapping_pair in router._route_mappings.items():
        inner_router = router_mapping_pair.router
        for mappings in router_mapping_pair.mapping:
            caller, callee = mappings.caller, mappings.callee
            if inner_router.consumes(method=callee, params=[param]):
                sub_output = StringIO()
                sub_routed = _visualise_param(
                    param, inner_router, indent=indent + 1, output=sub_output
                )
                if sub_routed:
                    print(" ." * indent, f"{caller} -> {name}.{callee}", file=output)
                    print(sub_output.getvalue(), end="", file=output)
                    routed = True

    return routed


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

search_cv = RandomizedSearchCV(clf, param_grid, cv=GroupKFold(), random_state=0)


# Get the routing information
test = get_routing_for_object(search_cv)

print("=" * 80)
print("OLD VERBOSE OUTPUT:")
print("=" * 80)
visualise_routing_old(test)

print("\n" * 3)
print("=" * 80)
print("NEW VISUALIZATION OPTIONS:")
print("=" * 80)

# Show all new visualization styles
for view_type in ["tree", "matrix", "compact", "flow"]:
    visualise_routing(test, view=view_type, show_method_mappings=False)
    print("\n")
