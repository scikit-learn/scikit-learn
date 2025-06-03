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
from sklearn.utils._metadata_requests import SIMPLE_METHODS
from sklearn.utils.metadata_routing import (
    MetadataRequest,
    MetadataRouter,
    get_routing_for_object,
)

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


def visualise_routing(routing_info, view="tree"):
    """
    Visualize metadata routing information.

    Parameters
    ----------
    routing_info : MetadataRouter or MetadataRequest
        The routing information to visualize
    view : str, default='tree'
        The visualization style. Options: 'tree', 'matrix', 'compact', 'flow'
    """
    if view == "tree":
        visualise_tree(routing_info)
    elif view == "matrix":
        visualise_matrix(routing_info)
    elif view == "compact":
        visualise_compact(routing_info)
    elif view == "flow":
        visualise_flow(routing_info)
    else:
        # Default to old behavior for backward compatibility
        visualise_routing_old(routing_info)


def visualise_tree(routing_info):
    """Show a consolidated tree view with parameters annotated inline."""
    print("\n=== METADATA ROUTING TREE ===")
    print(f"Root: {routing_info.owner}")

    # Collect all routing information first
    routing_map = _collect_routing_info(routing_info)

    # Display the tree
    _display_tree(routing_info, routing_map, indent=0)

    # Show parameter summary
    params = _param_names(routing_info)
    if params:
        print("\nParameters tracked:")
        for param in sorted(params):
            consumers = _find_consumers(routing_info, param)
            print(f"  • {param}: consumed by {len(consumers)} component(s)")


def visualise_matrix(routing_info):
    """Show a matrix view of components vs parameters."""
    print("\n=== METADATA ROUTING MATRIX ===")

    # Collect all components and parameters
    components = _collect_all_components(routing_info)
    params = sorted(_param_names(routing_info))

    if not params:
        print("No parameters are being routed.")
        return

    # Build the consumption matrix
    matrix = defaultdict(lambda: defaultdict(set))
    for param in params:
        consumers = _find_consumers_with_methods(routing_info, param)
        for comp_name, methods in consumers.items():
            matrix[comp_name][param] = methods

    # Calculate column widths
    comp_width = max(len(comp) for comp in components) + 2
    param_widths = {param: max(len(param), 8) + 2 for param in params}

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
            methods = matrix[comp][param]
            if methods:
                method_str = ",".join(sorted(methods))
                if len(method_str) > param_widths[param] - 2:
                    method_str = method_str[: param_widths[param] - 5] + "..."
                print(f"{method_str:^{param_widths[param]}}", end="")
            else:
                print(f"{'-':^{param_widths[param]}}", end="")
        print()


def visualise_compact(routing_info):
    """Show a compact hierarchical view."""
    print("\n=== COMPACT METADATA ROUTING ===")

    # Build compact representation
    structure = _build_compact_structure(routing_info)
    _print_compact_structure(structure, indent=0)

    # Show routing paths
    params = _param_names(routing_info)
    if params:
        print("\nRouting paths:")
        for param in sorted(params):
            paths = _find_routing_paths(routing_info, param)
            print(f"\n{param}:")
            for path in paths:
                print(f"  → {' → '.join(path)}")


def visualise_flow(routing_info):
    """Show an ASCII flow diagram."""
    print("\n=== METADATA ROUTING FLOW ===")

    # Build flow representation
    flow = _build_flow_diagram(routing_info)

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

    info = defaultdict(lambda: {"params": set(), "methods": defaultdict(set)})

    def _collect(obj, path=""):
        current_path = f"{path}/{obj.owner}" if path else obj.owner

        if isinstance(obj, MetadataRequest):
            # Get all parameters from the top-level router
            all_params = _param_names(top_router)
            for method in SIMPLE_METHODS:
                for param in all_params:
                    if obj.consumes(method, [param]):
                        info[current_path]["params"].add(param)
                        info[current_path]["methods"][param].add(method)

        elif isinstance(obj, MetadataRouter):
            if obj._self_request:
                all_params = _param_names(top_router)
                for method in SIMPLE_METHODS:
                    for param in all_params:
                        if obj._self_request.consumes(method, [param]):
                            info[current_path]["params"].add(param)
                            info[current_path]["methods"][param].add(method)

            for name, mapping in obj._route_mappings.items():
                _collect(mapping.router, current_path)

    _collect(router)
    return info


def _display_tree(router, routing_map, indent=0, prefix="", is_last=True):
    """Display the routing tree with inline parameter annotations."""
    connector = "└── " if is_last else "├── "

    # Get info for current node
    info = routing_map.get(router.owner, {"params": set(), "methods": defaultdict(set)})
    param_str = ""
    if info["params"]:
        param_details = []
        for param in sorted(info["params"]):
            methods = ",".join(sorted(info["methods"][param]))
            param_details.append(f"{param}[{methods}]")
        param_str = f" ({', '.join(param_details)})"

    if indent == 0:
        print("│")
    else:
        print(f"{prefix}{connector}{router.owner}{param_str}")

    # Process children
    if isinstance(router, MetadataRouter):
        children = list(router._route_mappings.items())
        for i, (name, mapping) in enumerate(children):
            is_last_child = i == len(children) - 1
            new_prefix = prefix + ("    " if is_last else "│   ")

            # Show method mapping
            if mapping.mapping:
                method_maps = []
                for m in mapping.mapping:
                    if m.caller != m.callee:
                        method_maps.append(f"{m.caller}→{m.callee}")
                    else:
                        method_maps.append(m.caller)
                method_str = f" [{', '.join(method_maps)}]"
            else:
                method_str = ""

            print(f"{new_prefix}{connector}{name}{method_str}")

            # Recurse
            newer_prefix = new_prefix + ("    " if is_last_child else "│   ")
            _display_tree_inner(mapping.router, routing_map, newer_prefix)


def _display_tree_inner(router, routing_map, prefix):
    """Helper for displaying inner tree nodes."""
    if isinstance(router, MetadataRequest):
        # Get parameter info from routing map
        info = routing_map.get(
            router.owner, {"params": set(), "methods": defaultdict(set)}
        )

        if info["params"]:
            params = []
            for param in sorted(info["params"]):
                methods = ",".join(sorted(info["methods"][param]))
                params.append(f"{param}[{methods}]")
            param_str = f" consumes: {', '.join(params)}"
            print(f"{prefix}    {router.owner}{param_str}")

    elif isinstance(router, MetadataRouter):
        # Display the router
        _display_tree(router, routing_map, len(prefix) // 4, prefix, True)


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


def _build_compact_structure(router, level=0, top_router=None):
    """Build a compact structure representation."""
    if top_router is None:
        top_router = router

    structure = {
        "name": router.owner,
        "type": type(router).__name__,
        "params": defaultdict(set),
        "children": [],
    }

    if isinstance(router, MetadataRequest):
        all_params = _param_names(top_router)
        for method in SIMPLE_METHODS:
            for param in all_params:
                if router.consumes(method, [param]):
                    structure["params"][param].add(method)

    elif isinstance(router, MetadataRouter):
        if router._self_request:
            all_params = _param_names(top_router)
            for method in SIMPLE_METHODS:
                for param in all_params:
                    if router._self_request.consumes(method, [param]):
                        structure["params"][param].add(method)

        for name, mapping in router._route_mappings.items():
            child = _build_compact_structure(mapping.router, level + 1, top_router)
            child["mapping_name"] = name
            child["mappings"] = [(m.caller, m.callee) for m in mapping.mapping]
            structure["children"].append(child)

    return structure


def _print_compact_structure(structure, indent=0):
    """Print the compact structure."""
    prefix = "  " * indent

    # Component name with type
    comp_str = structure["name"]
    if structure["params"]:
        param_strs = []
        for param, methods in structure["params"].items():
            param_strs.append(f"{param}[{','.join(sorted(methods))}]")
        comp_str += f" ◆ {', '.join(param_strs)}"

    print(f"{prefix}▸ {comp_str}")

    # Children
    for child in structure["children"]:
        mapping_str = ""
        if "mapping_name" in child:
            mapping_str = f" via {child['mapping_name']}"
            if child["mappings"]:
                unique_mappings = set()
                for caller, callee in child["mappings"]:
                    if caller != callee:
                        unique_mappings.add(f"{caller}→{callee}")
                if unique_mappings:
                    mapping_str += f" ({', '.join(unique_mappings)})"

        if mapping_str:
            print(f"{prefix}  ↳{mapping_str}")

        _print_compact_structure(child, indent + 2)


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


def _build_flow_diagram(router):
    """Build an ASCII flow diagram."""
    lines = []

    # Header
    lines.append("┌" + "─" * 50 + "┐")
    lines.append(f"│ {router.owner:^48} │")
    lines.append("└" + "─" * 25 + "┬" + "─" * 25 + "┘")

    # Collect all routing information
    routing_map = _collect_routing_info(router)
    params = sorted(_param_names(router))

    # Also create a map of which parameters each component consumes
    consumers_map = defaultdict(lambda: {"params": set(), "methods": defaultdict(set)})
    for param in params:
        consumers = _find_consumers(router, param)
        for consumer in consumers:
            consumers_map[consumer]["params"].add(param)

            # Find which methods consume this param for this consumer
            def _find_methods(obj, target_owner, param):
                if isinstance(obj, MetadataRequest) and obj.owner == target_owner:
                    for method in SIMPLE_METHODS:
                        if obj.consumes(method, [param]):
                            consumers_map[target_owner]["methods"][param].add(method)
                elif isinstance(obj, MetadataRouter):
                    if obj._self_request and obj.owner == target_owner:
                        for method in SIMPLE_METHODS:
                            if obj._self_request.consumes(method, [param]):
                                consumers_map[target_owner]["methods"][param].add(
                                    method
                                )
                    for mapping in obj._route_mappings.values():
                        _find_methods(mapping.router, target_owner, param)

            _find_methods(router, consumer, param)

    if params:
        lines.append(" " * 26 + "│")
        lines.append(" " * 20 + f"Parameters: {', '.join(params)}")
        lines.append(" " * 26 + "│")

    # Build the flow structure recursively
    def _add_flow_node(obj, indent=0, prefix="", is_last=True):
        if isinstance(obj, MetadataRouter):
            # Add connector for non-root nodes
            if indent > 0:
                connector = "└──" if is_last else "├──"
                lines.append(f"{prefix}{connector}─┐")
                lines.append(f"{prefix}{'   ' if is_last else '│  '} │")

            # Add component box
            comp_prefix = prefix + ("   " if is_last else "│  ") if indent > 0 else ""

            # Get consumed parameters for this component
            comp_info = consumers_map.get(
                obj.owner, {"params": set(), "methods": defaultdict(set)}
            )
            param_info = ""
            if comp_info["params"]:
                param_details = []
                for param in sorted(comp_info["params"]):
                    methods = ",".join(sorted(comp_info["methods"][param]))
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
                    if mapping.mapping:
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
                    _add_flow_node(mapping.router, indent + 1, new_prefix, True)

        elif isinstance(obj, MetadataRequest):
            # Terminal node - show consumed parameters
            comp_info = consumers_map.get(
                obj.owner, {"params": set(), "methods": defaultdict(set)}
            )
            param_str = ""
            if comp_info["params"]:
                param_details = []
                for param in sorted(comp_info["params"]):
                    methods = ",".join(sorted(comp_info["methods"][param]))
                    param_details.append(f"{param}[{methods}]")
                param_str = f" ◆ {', '.join(param_details)}"

            # Show component with or without parameters
            lines.append(f"{prefix}    ◆ {obj.owner}{param_str}")

    # Start building the flow
    _add_flow_node(router)

    # Add parameter flow summary at the bottom
    if params:
        lines.append("")
        lines.append("─" * 60)
        lines.append("Parameter Flow Summary:")
        for param in params:
            consumers = _find_consumers(router, param)
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
            .set_fit_request(sample_weight=True)
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
    visualise_routing(test, view=view_type)
    print("\n")
