import os

PLOTLY_DIR = os.environ.get(
    "PLOTLY_DIR", os.path.join(os.path.expanduser("~"), ".plotly")
)
TEST_FILE = os.path.join(PLOTLY_DIR, ".permission_test")


def _permissions():
    try:
        if not os.path.exists(PLOTLY_DIR):
            try:
                os.mkdir(PLOTLY_DIR)
            except Exception:
                # in case of race
                if not os.path.isdir(PLOTLY_DIR):
                    raise
        with open(TEST_FILE, "w") as f:
            f.write("testing\n")
        try:
            os.remove(TEST_FILE)
        except Exception:
            pass
        return True
    except Exception:  # Do not trap KeyboardInterrupt.
        return False


_file_permissions = None


def ensure_writable_plotly_dir():
    # Cache permissions status
    global _file_permissions
    if _file_permissions is None:
        _file_permissions = _permissions()
    return _file_permissions
