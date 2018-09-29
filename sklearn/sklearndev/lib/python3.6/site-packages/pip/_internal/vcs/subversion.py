from __future__ import absolute_import

import logging
import os
import re

from pip._vendor.six.moves.urllib import parse as urllib_parse

from pip._internal.index import Link
from pip._internal.utils.logging import indent_log
from pip._internal.utils.misc import display_path, remove_auth_from_url, rmtree
from pip._internal.vcs import VersionControl, vcs

_svn_xml_url_re = re.compile('url="([^"]+)"')
_svn_rev_re = re.compile(r'committed-rev="(\d+)"')
_svn_url_re = re.compile(r'URL: (.+)')
_svn_revision_re = re.compile(r'Revision: (.+)')
_svn_info_xml_rev_re = re.compile(r'\s*revision="(\d+)"')
_svn_info_xml_url_re = re.compile(r'<url>(.*)</url>')


logger = logging.getLogger(__name__)


class Subversion(VersionControl):
    name = 'svn'
    dirname = '.svn'
    repo_name = 'checkout'
    schemes = ('svn', 'svn+ssh', 'svn+http', 'svn+https', 'svn+svn')

    def get_base_rev_args(self, rev):
        return ['-r', rev]

    def get_info(self, location):
        """Returns (url, revision), where both are strings"""
        assert not location.rstrip('/').endswith(self.dirname), \
            'Bad directory: %s' % location
        output = self.run_command(
            ['info', location],
            show_stdout=False,
            extra_environ={'LANG': 'C'},
        )
        match = _svn_url_re.search(output)
        if not match:
            logger.warning(
                'Cannot determine URL of svn checkout %s',
                display_path(location),
            )
            logger.debug('Output that cannot be parsed: \n%s', output)
            return None, None
        url = match.group(1).strip()
        match = _svn_revision_re.search(output)
        if not match:
            logger.warning(
                'Cannot determine revision of svn checkout %s',
                display_path(location),
            )
            logger.debug('Output that cannot be parsed: \n%s', output)
            return url, None
        return url, match.group(1)

    def export(self, location):
        """Export the svn repository at the url to the destination location"""
        url, rev_options = self.get_url_rev_options(self.url)

        logger.info('Exporting svn repository %s to %s', url, location)
        with indent_log():
            if os.path.exists(location):
                # Subversion doesn't like to check out over an existing
                # directory --force fixes this, but was only added in svn 1.5
                rmtree(location)
            cmd_args = ['export'] + rev_options.to_args() + [url, location]
            self.run_command(cmd_args, show_stdout=False)

    def fetch_new(self, dest, url, rev_options):
        rev_display = rev_options.to_display()
        logger.info(
            'Checking out %s%s to %s',
            url,
            rev_display,
            display_path(dest),
        )
        cmd_args = ['checkout', '-q'] + rev_options.to_args() + [url, dest]
        self.run_command(cmd_args)

    def switch(self, dest, url, rev_options):
        cmd_args = ['switch'] + rev_options.to_args() + [url, dest]
        self.run_command(cmd_args)

    def update(self, dest, rev_options):
        cmd_args = ['update'] + rev_options.to_args() + [dest]
        self.run_command(cmd_args)

    def get_location(self, dist, dependency_links):
        for url in dependency_links:
            egg_fragment = Link(url).egg_fragment
            if not egg_fragment:
                continue
            if '-' in egg_fragment:
                # FIXME: will this work when a package has - in the name?
                key = '-'.join(egg_fragment.split('-')[:-1]).lower()
            else:
                key = egg_fragment
            if key == dist.key:
                return url.split('#', 1)[0]
        return None

    def get_revision(self, location):
        """
        Return the maximum revision for all files under a given location
        """
        # Note: taken from setuptools.command.egg_info
        revision = 0

        for base, dirs, files in os.walk(location):
            if self.dirname not in dirs:
                dirs[:] = []
                continue    # no sense walking uncontrolled subdirs
            dirs.remove(self.dirname)
            entries_fn = os.path.join(base, self.dirname, 'entries')
            if not os.path.exists(entries_fn):
                # FIXME: should we warn?
                continue

            dirurl, localrev = self._get_svn_url_rev(base)

            if base == location:
                base = dirurl + '/'   # save the root url
            elif not dirurl or not dirurl.startswith(base):
                dirs[:] = []
                continue    # not part of the same svn tree, skip it
            revision = max(revision, localrev)
        return revision

    def get_url_rev(self, url):
        # hotfix the URL scheme after removing svn+ from svn+ssh:// readd it
        url, rev = super(Subversion, self).get_url_rev(url)
        if url.startswith('ssh://'):
            url = 'svn+' + url
        return url, rev

    def get_url_rev_args(self, url):
        extra_args = get_rev_options_args(url)
        url = remove_auth_from_url(url)

        return url, extra_args

    def get_url(self, location):
        # In cases where the source is in a subdirectory, not alongside
        # setup.py we have to look up in the location until we find a real
        # setup.py
        orig_location = location
        while not os.path.exists(os.path.join(location, 'setup.py')):
            last_location = location
            location = os.path.dirname(location)
            if location == last_location:
                # We've traversed up to the root of the filesystem without
                # finding setup.py
                logger.warning(
                    "Could not find setup.py for directory %s (tried all "
                    "parent directories)",
                    orig_location,
                )
                return None

        return self._get_svn_url_rev(location)[0]

    def _get_svn_url_rev(self, location):
        from pip._internal.exceptions import InstallationError

        entries_path = os.path.join(location, self.dirname, 'entries')
        if os.path.exists(entries_path):
            with open(entries_path) as f:
                data = f.read()
        else:  # subversion >= 1.7 does not have the 'entries' file
            data = ''

        if (data.startswith('8') or
                data.startswith('9') or
                data.startswith('10')):
            data = list(map(str.splitlines, data.split('\n\x0c\n')))
            del data[0][0]  # get rid of the '8'
            url = data[0][3]
            revs = [int(d[9]) for d in data if len(d) > 9 and d[9]] + [0]
        elif data.startswith('<?xml'):
            match = _svn_xml_url_re.search(data)
            if not match:
                raise ValueError('Badly formatted data: %r' % data)
            url = match.group(1)    # get repository URL
            revs = [int(m.group(1)) for m in _svn_rev_re.finditer(data)] + [0]
        else:
            try:
                # subversion >= 1.7
                xml = self.run_command(
                    ['info', '--xml', location],
                    show_stdout=False,
                )
                url = _svn_info_xml_url_re.search(xml).group(1)
                revs = [
                    int(m.group(1)) for m in _svn_info_xml_rev_re.finditer(xml)
                ]
            except InstallationError:
                url, revs = None, []

        if revs:
            rev = max(revs)
        else:
            rev = 0

        return url, rev

    def get_src_requirement(self, dist, location):
        repo = self.get_url(location)
        if repo is None:
            return None
        # FIXME: why not project name?
        egg_project_name = dist.egg_name().split('-', 1)[0]
        rev = self.get_revision(location)
        return 'svn+%s@%s#egg=%s' % (repo, rev, egg_project_name)

    def is_commit_id_equal(self, dest, name):
        """Always assume the versions don't match"""
        return False


def get_rev_options_args(url):
    """
    Return the extra arguments to pass to RevOptions.
    """
    r = urllib_parse.urlsplit(url)
    if hasattr(r, 'username'):
        # >= Python-2.5
        username, password = r.username, r.password
    else:
        netloc = r[1]
        if '@' in netloc:
            auth = netloc.split('@')[0]
            if ':' in auth:
                username, password = auth.split(':', 1)
            else:
                username, password = auth, None
        else:
            username, password = None, None

    extra_args = []
    if username:
        extra_args += ['--username', username]
    if password:
        extra_args += ['--password', password]

    return extra_args


vcs.register(Subversion)
