"""Utilities to route metadata within scikit-learn estimators."""

# This module is not a separate sub-folder since that would result in a circular
# import issue.
#
# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from ._metadata_requests import WARN, UNUSED, UNCHANGED  # noqa: F401
from ._metadata_requests import get_routing_for_object  # noqa: F401
from ._metadata_requests import MetadataRouter  # noqa: F401
from ._metadata_requests import MetadataRequest  # noqa: F401
from ._metadata_requests import MethodMapping  # noqa: F401
from ._metadata_requests import process_routing  # noqa: F401
from ._metadata_requests import _MetadataRequester  # noqa: F401
from ._metadata_requests import _routing_enabled  # noqa: F401
from ._metadata_requests import _raise_for_params  # noqa: F401
from ._metadata_requests import _RoutingNotSupportedMixin  # noqa: F401
from ._metadata_requests import _raise_for_unsupported_routing  # noqa: F401
