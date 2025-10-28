from __future__ import annotations

from narwhals._compliant.typing import CompliantNamespaceT_co
from narwhals._namespace import Namespace as NwNamespace
from narwhals._utils import Version

__all__ = ["Namespace"]


class Namespace(NwNamespace[CompliantNamespaceT_co], version=Version.V2): ...
