"""Utilities to route metadata within scikit-learn estimators."""

# This module is not a separate sub-folder since that would result in a circular
# import issue.
#
# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from ._metadata_requests import (  # noqa: F401
    UNCHANGED,
    UNUSED,
    WARN,
    MetadataRequest,
    MetadataRouter,
    MethodMapping,
    _MetadataRequester,
    _raise_for_params,
    _raise_for_unsupported_routing,
    _routing_enabled,
    _RoutingNotSupportedMixin,
    get_routing_for_object,
    process_routing,
)
