from __future__ import absolute_import

import collections
import logging
import os
import re

from pip._vendor import pkg_resources, six
from pip._vendor.packaging.utils import canonicalize_name
from pip._vendor.pkg_resources import RequirementParseError

from pip._internal.exceptions import InstallationError
from pip._internal.req import InstallRequirement
from pip._internal.req.req_file import COMMENT_RE
from pip._internal.utils.deprecation import deprecated
from pip._internal.utils.misc import (
    dist_is_editable, get_installed_distributions,
)

logger = logging.getLogger(__name__)


def freeze(
        requirement=None,
        find_links=None, local_only=None, user_only=None, skip_regex=None,
        isolated=False,
        wheel_cache=None,
        exclude_editable=False,
        skip=()):
    find_links = find_links or []
    skip_match = None

    if skip_regex:
        skip_match = re.compile(skip_regex).search

    dependency_links = []

    for dist in pkg_resources.working_set:
        if dist.has_metadata('dependency_links.txt'):
            dependency_links.extend(
                dist.get_metadata_lines('dependency_links.txt')
            )
    for link in find_links:
        if '#egg=' in link:
            dependency_links.append(link)
    for link in find_links:
        yield '-f %s' % link
    installations = {}
    for dist in get_installed_distributions(local_only=local_only,
                                            skip=(),
                                            user_only=user_only):
        try:
            req = FrozenRequirement.from_dist(
                dist,
                dependency_links
            )
        except RequirementParseError:
            logger.warning(
                "Could not parse requirement: %s",
                dist.project_name
            )
            continue
        if exclude_editable and req.editable:
            continue
        installations[req.name] = req

    if requirement:
        # the options that don't get turned into an InstallRequirement
        # should only be emitted once, even if the same option is in multiple
        # requirements files, so we need to keep track of what has been emitted
        # so that we don't emit it again if it's seen again
        emitted_options = set()
        # keep track of which files a requirement is in so that we can
        # give an accurate warning if a requirement appears multiple times.
        req_files = collections.defaultdict(list)
        for req_file_path in requirement:
            with open(req_file_path) as req_file:
                for line in req_file:
                    if (not line.strip() or
                            line.strip().startswith('#') or
                            (skip_match and skip_match(line)) or
                            line.startswith((
                                '-r', '--requirement',
                                '-Z', '--always-unzip',
                                '-f', '--find-links',
                                '-i', '--index-url',
                                '--pre',
                                '--trusted-host',
                                '--process-dependency-links',
                                '--extra-index-url'))):
                        line = line.rstrip()
                        if line not in emitted_options:
                            emitted_options.add(line)
                            yield line
                        continue

                    if line.startswith('-e') or line.startswith('--editable'):
                        if line.startswith('-e'):
                            line = line[2:].strip()
                        else:
                            line = line[len('--editable'):].strip().lstrip('=')
                        line_req = InstallRequirement.from_editable(
                            line,
                            isolated=isolated,
                            wheel_cache=wheel_cache,
                        )
                    else:
                        line_req = InstallRequirement.from_line(
                            COMMENT_RE.sub('', line).strip(),
                            isolated=isolated,
                            wheel_cache=wheel_cache,
                        )

                    if not line_req.name:
                        logger.info(
                            "Skipping line in requirement file [%s] because "
                            "it's not clear what it would install: %s",
                            req_file_path, line.strip(),
                        )
                        logger.info(
                            "  (add #egg=PackageName to the URL to avoid"
                            " this warning)"
                        )
                    elif line_req.name not in installations:
                        # either it's not installed, or it is installed
                        # but has been processed already
                        if not req_files[line_req.name]:
                            logger.warning(
                                "Requirement file [%s] contains %s, but that "
                                "package is not installed",
                                req_file_path,
                                COMMENT_RE.sub('', line).strip(),
                            )
                        else:
                            req_files[line_req.name].append(req_file_path)
                    else:
                        yield str(installations[line_req.name]).rstrip()
                        del installations[line_req.name]
                        req_files[line_req.name].append(req_file_path)

        # Warn about requirements that were included multiple times (in a
        # single requirements file or in different requirements files).
        for name, files in six.iteritems(req_files):
            if len(files) > 1:
                logger.warning("Requirement %s included multiple times [%s]",
                               name, ', '.join(sorted(set(files))))

        yield(
            '## The following requirements were added by '
            'pip freeze:'
        )
    for installation in sorted(
            installations.values(), key=lambda x: x.name.lower()):
        if canonicalize_name(installation.name) not in skip:
            yield str(installation).rstrip()


class FrozenRequirement(object):
    def __init__(self, name, req, editable, comments=()):
        self.name = name
        self.req = req
        self.editable = editable
        self.comments = comments

    _rev_re = re.compile(r'-r(\d+)$')
    _date_re = re.compile(r'-(20\d\d\d\d\d\d)$')

    @classmethod
    def from_dist(cls, dist, dependency_links):
        location = os.path.normcase(os.path.abspath(dist.location))
        comments = []
        from pip._internal.vcs import vcs, get_src_requirement
        if dist_is_editable(dist) and vcs.get_backend_name(location):
            editable = True
            try:
                req = get_src_requirement(dist, location)
            except InstallationError as exc:
                logger.warning(
                    "Error when trying to get requirement for VCS system %s, "
                    "falling back to uneditable format", exc
                )
                req = None
            if req is None:
                logger.warning(
                    'Could not determine repository location of %s', location
                )
                comments.append(
                    '## !! Could not determine repository location'
                )
                req = dist.as_requirement()
                editable = False
        else:
            editable = False
            req = dist.as_requirement()
            specs = req.specs
            assert len(specs) == 1 and specs[0][0] in ["==", "==="], \
                'Expected 1 spec with == or ===; specs = %r; dist = %r' % \
                (specs, dist)
            version = specs[0][1]
            ver_match = cls._rev_re.search(version)
            date_match = cls._date_re.search(version)
            if ver_match or date_match:
                svn_backend = vcs.get_backend('svn')
                if svn_backend:
                    svn_location = svn_backend().get_location(
                        dist,
                        dependency_links,
                    )
                if not svn_location:
                    logger.warning(
                        'Warning: cannot find svn location for %s', req,
                    )
                    comments.append(
                        '## FIXME: could not find svn URL in dependency_links '
                        'for this package:'
                    )
                else:
                    deprecated(
                        "SVN editable detection based on dependency links "
                        "will be dropped in the future.",
                        replacement=None,
                        gone_in="18.2",
                        issue=4187,
                    )
                    comments.append(
                        '# Installing as editable to satisfy requirement %s:' %
                        req
                    )
                    if ver_match:
                        rev = ver_match.group(1)
                    else:
                        rev = '{%s}' % date_match.group(1)
                    editable = True
                    req = '%s@%s#egg=%s' % (
                        svn_location,
                        rev,
                        cls.egg_name(dist)
                    )
        return cls(dist.project_name, req, editable, comments)

    @staticmethod
    def egg_name(dist):
        name = dist.egg_name()
        match = re.search(r'-py\d\.\d$', name)
        if match:
            name = name[:match.start()]
        return name

    def __str__(self):
        req = self.req
        if self.editable:
            req = '-e %s' % req
        return '\n'.join(list(self.comments) + [str(req)]) + '\n'
