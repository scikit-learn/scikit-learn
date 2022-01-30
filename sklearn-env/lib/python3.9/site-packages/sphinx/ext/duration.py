"""
    sphinx.ext.duration
    ~~~~~~~~~~~~~~~~~~~

    Measure durations of Sphinx processing.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

from datetime import datetime, timedelta
from itertools import islice
from operator import itemgetter
from typing import Any, Dict, List, cast

from docutils import nodes

from sphinx.application import Sphinx
from sphinx.domains import Domain
from sphinx.locale import __
from sphinx.util import logging

logger = logging.getLogger(__name__)


class DurationDomain(Domain):
    """A domain for durations of Sphinx processing."""
    name = 'duration'

    @property
    def reading_durations(self) -> Dict[str, timedelta]:
        return self.data.setdefault('reading_durations', {})

    def note_reading_duration(self, duration: timedelta) -> None:
        self.reading_durations[self.env.docname] = duration

    def clear(self) -> None:
        self.reading_durations.clear()

    def clear_doc(self, docname: str) -> None:
        self.reading_durations.pop(docname, None)

    def merge_domaindata(self, docnames: List[str], otherdata: Dict[str, timedelta]) -> None:
        for docname, duration in otherdata.items():
            if docname in docnames:
                self.reading_durations[docname] = duration


def on_builder_inited(app: Sphinx) -> None:
    """Initialize DurationDomain on bootstrap.

    This clears results of last build.
    """
    domain = cast(DurationDomain, app.env.get_domain('duration'))
    domain.clear()


def on_source_read(app: Sphinx, docname: str, content: List[str]) -> None:
    """Start to measure reading duration."""
    app.env.temp_data['started_at'] = datetime.now()


def on_doctree_read(app: Sphinx, doctree: nodes.document) -> None:
    """Record a reading duration."""
    started_at = app.env.temp_data.get('started_at')
    duration = datetime.now() - started_at
    domain = cast(DurationDomain, app.env.get_domain('duration'))
    domain.note_reading_duration(duration)


def on_build_finished(app: Sphinx, error: Exception) -> None:
    """Display duration ranking on current build."""
    domain = cast(DurationDomain, app.env.get_domain('duration'))
    durations = sorted(domain.reading_durations.items(), key=itemgetter(1), reverse=True)
    if not durations:
        return

    logger.info('')
    logger.info(__('====================== slowest reading durations ======================='))
    for docname, d in islice(durations, 5):
        logger.info('%d.%03d %s', d.seconds, d.microseconds / 1000, docname)


def setup(app: Sphinx) -> Dict[str, Any]:
    app.add_domain(DurationDomain)
    app.connect('builder-inited', on_builder_inited)
    app.connect('source-read', on_source_read)
    app.connect('doctree-read', on_doctree_read)
    app.connect('build-finished', on_build_finished)

    return {
        'version': 'builtin',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
