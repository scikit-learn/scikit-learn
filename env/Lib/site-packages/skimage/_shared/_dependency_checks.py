from .version_requirements import is_installed
import sys
import platform

has_mpl = is_installed("matplotlib", ">=3.3")

is_wasm = (sys.platform == "emscripten") or (platform.machine() in ["wasm32", "wasm64"])
