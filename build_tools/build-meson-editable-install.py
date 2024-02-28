import re
import shlex
import subprocess
from pathlib import Path


def main():
    pyproject_path = Path("pyproject.toml")

    if not pyproject_path.exists():
        raise SystemExit(
            "Can not find pyproject.toml. You should run this script from the"
            " scikit-learn root folder."
        )

    old_pyproject_content = pyproject_path.read_text(encoding="utf-8")
    if 'build-backend = "mesonpy"' not in old_pyproject_content:
        new_pyproject_content = re.sub(
            r"\[build-system\]",
            r'[build-system]\nbuild-backend = "mesonpy"',
            old_pyproject_content,
        )
        pyproject_path.write_text(new_pyproject_content, encoding="utf-8")

    command = shlex.split(
        "pip install --editable .  --verbose --no-build-isolation "
        "--config-settings editable-verbose=true"
    )

    exception = None
    try:
        subprocess.check_call(command)
    except Exception as e:
        exception = e
    finally:
        pyproject_path.write_text(old_pyproject_content, encoding="utf-8")

    if exception is not None:
        raise RuntimeError(
            "There was some error when running the script"
        ) from exception


if __name__ == "__main__":
    main()
