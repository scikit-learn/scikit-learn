"""
Metadata Routing Tree Visualization

This module provides functions to visualize metadata routing information
in a hierarchical tree format, making it easier to understand the routing
relationships between components.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from io import StringIO

from sklearn.utils._metadata_requests import SIMPLE_METHODS
from sklearn.utils.metadata_routing import (
    MetadataRequest,
    MetadataRouter,
)


class MetadataTreeNode:
    """A node in the metadata routing tree."""

    def __init__(self, name, node_type):
        self.name = name
        self.node_type = node_type  # "router" or "consumer"
        self.children = {}  # name -> MetadataTreeNode
        self.parameters = {}  # param -> {method -> caller/callee}
        self.method_mappings = {}  # caller -> [(callee, child_name)]


def _build_tree(router, node_name, parent=None):
    """Build a tree representation of the routing information."""
    if isinstance(router, MetadataRouter):
        node = MetadataTreeNode(node_name, "router")

        # Handle self request if router is also a consumer
        if router._self_request:
            self_node = _build_tree(router._self_request, "self_request", node)
            node.children["self"] = self_node

        # Process each child router/consumer
        for name, router_mapping_pair in router._route_mappings.items():
            inner_router = router_mapping_pair.router

            # Store method mappings
            for mapping in router_mapping_pair.mapping:
                caller, callee = mapping.caller, mapping.callee
                if caller not in node.method_mappings:
                    node.method_mappings[caller] = []
                node.method_mappings[caller].append((callee, name))

            # Recursively build child node
            child_node = _build_tree(inner_router, name, node)
            node.children[name] = child_node

        return node

    elif isinstance(router, MetadataRequest):
        node = MetadataTreeNode(node_name, "consumer")

        # Process parameters for each method
        for method in SIMPLE_METHODS:
            method_request = getattr(router, method)
            for param, alias in method_request.requests.items():
                if param not in node.parameters:
                    node.parameters[param] = {}
                node.parameters[param][method] = alias

        return node

    else:
        raise ValueError(f"Unknown routing info type: {type(router)}")


def _param_names(router):
    """Get all parameter names from a router."""
    res = set()
    for method in SIMPLE_METHODS:
        res = res.union(
            router._get_param_names(
                method=method, return_alias=True, ignore_self_request=False
            )
        )
    return res


def _render_node(node, param_filter=None, indent=0, output=None):
    """Render a node in the tree."""
    if output is None:
        output = StringIO()

    # Box drawing characters for better visualization
    indent_char = "  "
    corner = "└─"
    tee = "├─"
    vertical = "│ "
    horizontal = "──"

    # Print node name with appropriate type indicator
    if indent == 0:
        print(f"┌─ {node.name}", file=output)
    else:
        prefix = indent_char * (indent - 1)
        if node.name == "self":
            print(f"{prefix}{corner} Self request", file=output)
        else:
            print(f"{prefix}{corner} {node.name}", file=output)

    # Print parameters consumed by this node
    if node.parameters and (
        not param_filter or any(p in param_filter for p in node.parameters)
    ):
        prefix = indent_char * indent
        print(f"{prefix}{tee} Parameters:", file=output)

        param_list = sorted(node.parameters.keys())
        if param_filter:
            param_list = [p for p in param_list if p in param_filter]

        for i, param in enumerate(param_list):
            methods = sorted(node.parameters[param].keys())
            is_last = i == len(param_list) - 1

            param_prefix = f"{prefix}{vertical if not is_last else indent_char}"
            method_str = ", ".join(f"{m}" for m in methods)
            print(
                f"{param_prefix}{corner if is_last else tee} [{param}] → {method_str}",
                file=output,
            )

    # Print method mappings
    if node.method_mappings:
        prefix = indent_char * indent
        if node.parameters:
            print(f"{prefix}{tee} Method mappings:", file=output)
        else:
            print(f"{prefix}{tee} Method mappings:", file=output)

        mappings = []
        for caller, callees in node.method_mappings.items():
            for callee, child_name in callees:
                mappings.append((caller, callee, child_name))

        # Sort by caller, then child_name
        mappings.sort(key=lambda x: (x[0], x[2]))

        for i, (caller, callee, child_name) in enumerate(mappings):
            is_last = i == len(mappings) - 1 and not node.children
            map_prefix = f"{prefix}{vertical if not is_last else indent_char}"
            print(
                f"{map_prefix}{corner if is_last else tee} {caller} → "
                f"{child_name}.{callee}",
                file=output,
            )

    # Print children
    if node.children:
        child_names = sorted(node.children.keys())
        if "self" in child_names:  # Always process "self" first
            child_names.remove("self")
            child_names = ["self"] + child_names

        for i, child_name in enumerate(child_names):
            child = node.children[child_name]
            _render_node(child, param_filter, indent + 1, output)

    return output


def visualise_routing_tree(routing_info, param=None):
    """Visualize routing information as a hierarchical tree.

    This function creates a tree-based visualization of the metadata routing
    information, showing how parameters are routed through different components.

    Parameters
    ----------
    routing_info : MetadataRouter or MetadataRequest
        The routing information to visualize.
    param : str, optional
        If provided, only show routing for this specific parameter.

    Returns
    -------
    str
        The formatted tree representation as a string.
    """
    # Build the tree
    root = _build_tree(routing_info, routing_info.owner)

    # Get parameters to filter by
    params_to_show = None
    if param:
        params_to_show = {param}

    # Render the tree
    output = _render_node(root, params_to_show)
    return output.getvalue()


def visualise_routing(routing_info):
    """Visualize routing information for all parameters.

    Parameters
    ----------
    routing_info : MetadataRouter or MetadataRequest
        The routing information to visualize.
    """
    params = _param_names(routing_info)
    if not params:
        print(f"No metadata routing found for {routing_info.owner}")
        return

    print(visualise_routing_tree(routing_info))
