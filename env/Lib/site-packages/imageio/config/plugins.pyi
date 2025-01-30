from typing import Any, Dict, Optional
from ..core.v3_plugin_api import PluginV3

class PluginConfig:
    name: str
    class_name: str
    module_name: str
    is_legacy: bool
    package_name: Optional[str] = None
    install_name: Optional[str] = None
    legacy_args: Optional[dict] = None
    @property
    def format(self) -> Any: ...
    @property
    def plugin_class(self) -> PluginV3: ...
    def __init__(
        self,
        name: str,
        class_name: str,
        module_name: str,
        *,
        is_legacy: bool = False,
        package_name: str = None,
        install_name: str = None,
        legacy_args: dict = None,
    ) -> None: ...

known_plugins: Dict[str, PluginConfig]
