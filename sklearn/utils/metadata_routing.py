"""Utilities to route metadata within scikit-learn estimators."""

# This module is not a separate sub-folder since that would result in a circular
# import issue.
#
# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from ._metadata_requests import WARN, UNUSED, UNCHANGED  # noqa
from ._metadata_requests import get_routing_for_object
from ._metadata_requests import MetadataRouter
from ._metadata_requests import MetadataRequest
from ._metadata_requests import MethodMapping
from ._metadata_requests import process_routing
from ._metadata_requests import _MetadataRequester  # noqa
from ._metadata_requests import _routing_enabled  # noqa
from ._metadata_requests import _raise_for_params  # noqa
from ._metadata_requests import _RoutingNotSupportedMixin  # noqa
from ._metadata_requests import _raise_for_unsupported_routing  # noqa


__all__ = [
    "MetadataRequest",
    "MetadataRouter",
    "MethodMapping",
    "get_routing_for_object",
    "process_routing",
]


def __dir__():
    return __all__
