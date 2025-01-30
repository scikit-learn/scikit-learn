"""Record metadata for the build process."""

from __future__ import annotations

from typing import TYPE_CHECKING

from sphinx.locale import __
from sphinx.util._serialise import stable_hash

if TYPE_CHECKING:
    from collections.abc import Set
    from pathlib import Path

    from sphinx.config import Config, _ConfigRebuild
    from sphinx.util.tags import Tags


class BuildInfo:
    """buildinfo file manipulator.

    HTMLBuilder and its family are storing their own envdata to ``.buildinfo``.
    This class is a manipulator for the file.
    """

    @classmethod
    def load(cls: type[BuildInfo], filename: Path, /) -> BuildInfo:
        content = filename.read_text(encoding='utf-8')
        lines = content.splitlines()

        version = lines[0].rstrip()
        if version != '# Sphinx build info version 1':
            msg = __('failed to read broken build info file (unknown version)')
            raise ValueError(msg)

        if not lines[2].startswith('config: '):
            msg = __('failed to read broken build info file (missing config entry)')
            raise ValueError(msg)
        if not lines[3].startswith('tags: '):
            msg = __('failed to read broken build info file (missing tags entry)')
            raise ValueError(msg)

        build_info = BuildInfo()
        build_info.config_hash = lines[2].removeprefix('config: ').strip()
        build_info.tags_hash = lines[3].removeprefix('tags: ').strip()
        return build_info

    def __init__(
        self,
        config: Config | None = None,
        tags: Tags | None = None,
        config_categories: Set[_ConfigRebuild] = frozenset(),
    ) -> None:
        self.config_hash = ''
        self.tags_hash = ''

        if config:
            values = {c.name: c.value for c in config.filter(config_categories)}
            self.config_hash = stable_hash(values)

        if tags:
            self.tags_hash = stable_hash(sorted(tags))

    def __eq__(self, other: BuildInfo) -> bool:  # type: ignore[override]
        return (
            self.config_hash == other.config_hash and self.tags_hash == other.tags_hash
        )

    def dump(self, filename: Path, /) -> None:
        build_info = (
            '# Sphinx build info version 1\n'
            '# This file records the configuration used when building these files. '
            'When it is not found, a full rebuild will be done.\n'
            f'config: {self.config_hash}\n'
            f'tags: {self.tags_hash}\n'
        )
        filename.write_text(build_info, encoding='utf-8')
