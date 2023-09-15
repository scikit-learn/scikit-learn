#!/bin/bash

set -e

git clone https://github.com/emscripten-core/emsdk.git
cd emsdk
./emsdk install $EMSCRIPTEN_VERSION
./emsdk activate $EMSCRIPTEN_VERSION
source emsdk_env.sh
cd -

pip install pyodide-build==$PYODIDE_VERSION pyodide-cli

# TODO: temporary fix for pywasmcross.py taken from https://github.com/pyodide/pyodide/pull/4136
PYODIDE_BUILD_DIR=$(python -c 'import os; import pyodide_build; print(os.path.dirname(pyodide_build.__file__))')
patch ${PYODIDE_BUILD_DIR}/pywasmcross.py << EOF
diff --git a/pyodide-build/pyodide_build/pywasmcross.py b/pyodide-build/pyodide_build/pywasmcross.py
index 16632e1f..92ca04ef 100755
--- a/pyodide-build/pyodide_build/pywasmcross.py
+++ b/pyodide-build/pyodide_build/pywasmcross.py
@@ -226,14 +226,19 @@ def replay_genargs_handle_dashI(arg: str, target_install_dir: str) -> str | None
         The new argument, or None to delete the argument.
     """
     assert arg.startswith("-I")
-    if (
-        str(Path(arg[2:]).resolve()).startswith(sys.prefix + "/include/python")
-        and "site-packages" not in arg
-    ):
-        return arg.replace("-I" + sys.prefix, "-I" + target_install_dir)
+
     # Don't include any system directories
     if arg[2:].startswith("/usr"):
         return None
+
+    # Replace local Python include paths with the cross compiled ones
+    include_path = str(Path(arg[2:]).resolve())
+    if include_path.startswith(sys.prefix + "/include/python"):
+        return arg.replace("-I" + sys.prefix, "-I" + target_install_dir)
+
+    if include_path.startswith(sys.base_prefix + "/include/python"):
+        return arg.replace("-I" + sys.base_prefix, "-I" + target_install_dir)
+
     return arg
EOF

pyodide build

ls -ltrh dist
