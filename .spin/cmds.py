import shutil
import sys

import click


@click.command()
def clean():
    """🔧 Clean Meson build.

    Very rarely needed since meson-python recompile as needed when sklearn is
    imported.

    One known use case where "spin clean" is useful: to avoid compilation
    errors when switching from numpy<2 to numpy>=2 in the same conda
    environment or virtualenv.
    """
    default_meson_build_dir = (
        f"build/cp{sys.version_info.major}{sys.version_info.minor}"
    )
    click.secho(
        f"removing default Meson build dir: {default_meson_build_dir}",
        bold=True,
        fg="bright_blue",
    )

    shutil.rmtree(default_meson_build_dir, ignore_errors=True)
