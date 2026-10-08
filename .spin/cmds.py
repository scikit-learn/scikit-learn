import shutil
import sys

import click
from spin.cmds import util


@click.command()
def clean():
    """🪥 Clean build folder.

    Very rarely needed since meson-python recompiles as needed when sklearn is
    imported.

    One known use case where "spin clean" is useful: avoid compilation errors
    when switching from numpy<2 to numpy>=2 in the same conda environment or
    virtualenv.
    """
    util.run([sys.executable, "-m", "pip", "uninstall", "scikit-learn", "-y"])
    default_meson_build_dir = (
        f"build/cp{sys.version_info.major}{sys.version_info.minor}"
    )
    click.secho(
        f"removing default Meson build dir: {default_meson_build_dir}",
        bold=True,
        fg="bright_blue",
    )

    shutil.rmtree(default_meson_build_dir, ignore_errors=True)
