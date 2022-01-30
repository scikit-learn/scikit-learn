"""Coverage controllers for use by pytest-cov and nose-cov."""
import contextlib
import copy
import functools
import os
import random
import socket
import sys

import coverage
from coverage.data import CoverageData

from .compat import StringIO
from .embed import cleanup


class _NullFile:
    @staticmethod
    def write(v):
        pass


@contextlib.contextmanager
def _backup(obj, attr):
    backup = getattr(obj, attr)
    try:
        setattr(obj, attr, copy.copy(backup))
        yield
    finally:
        setattr(obj, attr, backup)


def _ensure_topdir(meth):
    @functools.wraps(meth)
    def ensure_topdir_wrapper(self, *args, **kwargs):
        try:
            original_cwd = os.getcwd()
        except OSError:
            # Looks like it's gone, this is non-ideal because a side-effect will
            # be introduced in the tests here but we can't do anything about it.
            original_cwd = None
        os.chdir(self.topdir)
        try:
            return meth(self, *args, **kwargs)
        finally:
            if original_cwd is not None:
                os.chdir(original_cwd)

    return ensure_topdir_wrapper


class CovController:
    """Base class for different plugin implementations."""

    def __init__(self, cov_source, cov_report, cov_config, cov_append, cov_branch, config=None, nodeid=None):
        """Get some common config used by multiple derived classes."""
        self.cov_source = cov_source
        self.cov_report = cov_report
        self.cov_config = cov_config
        self.cov_append = cov_append
        self.cov_branch = cov_branch
        self.config = config
        self.nodeid = nodeid

        self.cov = None
        self.combining_cov = None
        self.data_file = None
        self.node_descs = set()
        self.failed_workers = []
        self.topdir = os.getcwd()
        self.is_collocated = None

    @contextlib.contextmanager
    def ensure_topdir(self):
        original_cwd = os.getcwd()
        os.chdir(self.topdir)
        yield
        os.chdir(original_cwd)

    @_ensure_topdir
    def pause(self):
        self.cov.stop()
        self.unset_env()

    @_ensure_topdir
    def resume(self):
        self.cov.start()
        self.set_env()

    @_ensure_topdir
    def set_env(self):
        """Put info about coverage into the env so that subprocesses can activate coverage."""
        if self.cov_source is None:
            os.environ['COV_CORE_SOURCE'] = os.pathsep
        else:
            os.environ['COV_CORE_SOURCE'] = os.pathsep.join(self.cov_source)
        config_file = os.path.abspath(self.cov_config)
        if os.path.exists(config_file):
            os.environ['COV_CORE_CONFIG'] = config_file
        else:
            os.environ['COV_CORE_CONFIG'] = os.pathsep
        os.environ['COV_CORE_DATAFILE'] = os.path.abspath(self.cov.config.data_file)
        if self.cov_branch:
            os.environ['COV_CORE_BRANCH'] = 'enabled'

    @staticmethod
    def unset_env():
        """Remove coverage info from env."""
        os.environ.pop('COV_CORE_SOURCE', None)
        os.environ.pop('COV_CORE_CONFIG', None)
        os.environ.pop('COV_CORE_DATAFILE', None)
        os.environ.pop('COV_CORE_BRANCH', None)
        os.environ.pop('COV_CORE_CONTEXT', None)

    @staticmethod
    def get_node_desc(platform, version_info):
        """Return a description of this node."""

        return 'platform {}, python {}'.format(platform, '%s.%s.%s-%s-%s' % version_info[:5])

    @staticmethod
    def sep(stream, s, txt):
        if hasattr(stream, 'sep'):
            stream.sep(s, txt)
        else:
            sep_total = max((70 - 2 - len(txt)), 2)
            sep_len = sep_total // 2
            sep_extra = sep_total % 2
            out = f'{s * sep_len} {txt} {s * (sep_len + sep_extra)}\n'
            stream.write(out)

    @_ensure_topdir
    def summary(self, stream):
        """Produce coverage reports."""
        total = None

        if not self.cov_report:
            with _backup(self.cov, "config"):
                return self.cov.report(show_missing=True, ignore_errors=True, file=_NullFile)

        # Output coverage section header.
        if len(self.node_descs) == 1:
            self.sep(stream, '-', 'coverage: %s' % ''.join(self.node_descs))
        else:
            self.sep(stream, '-', 'coverage')
            for node_desc in sorted(self.node_descs):
                self.sep(stream, ' ', '%s' % node_desc)

        # Report on any failed workers.
        if self.failed_workers:
            self.sep(stream, '-', 'coverage: failed workers')
            stream.write('The following workers failed to return coverage data, '
                         'ensure that pytest-cov is installed on these workers.\n')
            for node in self.failed_workers:
                stream.write('%s\n' % node.gateway.id)

        # Produce terminal report if wanted.
        if any(x in self.cov_report for x in ['term', 'term-missing']):
            options = {
                'show_missing': ('term-missing' in self.cov_report) or None,
                'ignore_errors': True,
                'file': stream,
            }
            skip_covered = isinstance(self.cov_report, dict) and 'skip-covered' in self.cov_report.values()
            options.update({'skip_covered': skip_covered or None})
            with _backup(self.cov, "config"):
                total = self.cov.report(**options)

        # Produce annotated source code report if wanted.
        if 'annotate' in self.cov_report:
            annotate_dir = self.cov_report['annotate']

            with _backup(self.cov, "config"):
                self.cov.annotate(ignore_errors=True, directory=annotate_dir)
            # We need to call Coverage.report here, just to get the total
            # Coverage.annotate don't return any total and we need it for --cov-fail-under.

            with _backup(self.cov, "config"):
                total = self.cov.report(ignore_errors=True, file=_NullFile)
            if annotate_dir:
                stream.write('Coverage annotated source written to dir %s\n' % annotate_dir)
            else:
                stream.write('Coverage annotated source written next to source\n')

        # Produce html report if wanted.
        if 'html' in self.cov_report:
            output = self.cov_report['html']
            with _backup(self.cov, "config"):
                total = self.cov.html_report(ignore_errors=True, directory=output)
            stream.write('Coverage HTML written to dir %s\n' % (self.cov.config.html_dir if output is None else output))

        # Produce xml report if wanted.
        if 'xml' in self.cov_report:
            output = self.cov_report['xml']
            with _backup(self.cov, "config"):
                total = self.cov.xml_report(ignore_errors=True, outfile=output)
            stream.write('Coverage XML written to file %s\n' % (self.cov.config.xml_output if output is None else output))

        return total


class Central(CovController):
    """Implementation for centralised operation."""

    @_ensure_topdir
    def start(self):
        cleanup()

        self.cov = coverage.Coverage(source=self.cov_source,
                                     branch=self.cov_branch,
                                     data_suffix=True,
                                     config_file=self.cov_config)
        self.combining_cov = coverage.Coverage(source=self.cov_source,
                                               branch=self.cov_branch,
                                               data_suffix=True,
                                               data_file=os.path.abspath(self.cov.config.data_file),
                                               config_file=self.cov_config)

        # Erase or load any previous coverage data and start coverage.
        if not self.cov_append:
            self.cov.erase()
        self.cov.start()
        self.set_env()

    @_ensure_topdir
    def finish(self):
        """Stop coverage, save data to file and set the list of coverage objects to report on."""

        self.unset_env()
        self.cov.stop()
        self.cov.save()

        self.cov = self.combining_cov
        self.cov.load()
        self.cov.combine()
        self.cov.save()

        node_desc = self.get_node_desc(sys.platform, sys.version_info)
        self.node_descs.add(node_desc)


class DistMaster(CovController):
    """Implementation for distributed master."""

    @_ensure_topdir
    def start(self):
        cleanup()

        # Ensure coverage rc file rsynced if appropriate.
        if self.cov_config and os.path.exists(self.cov_config):
            self.config.option.rsyncdir.append(self.cov_config)

        self.cov = coverage.Coverage(source=self.cov_source,
                                     branch=self.cov_branch,
                                     data_suffix=True,
                                     config_file=self.cov_config)
        self.cov._warn_no_data = False
        self.cov._warn_unimported_source = False
        self.cov._warn_preimported_source = False
        self.combining_cov = coverage.Coverage(source=self.cov_source,
                                               branch=self.cov_branch,
                                               data_suffix=True,
                                               data_file=os.path.abspath(self.cov.config.data_file),
                                               config_file=self.cov_config)
        if not self.cov_append:
            self.cov.erase()
        self.cov.start()
        self.cov.config.paths['source'] = [self.topdir]

    def configure_node(self, node):
        """Workers need to know if they are collocated and what files have moved."""

        node.workerinput.update({
            'cov_master_host': socket.gethostname(),
            'cov_master_topdir': self.topdir,
            'cov_master_rsync_roots': [str(root) for root in node.nodemanager.roots],
        })

    def testnodedown(self, node, error):
        """Collect data file name from worker."""

        # If worker doesn't return any data then it is likely that this
        # plugin didn't get activated on the worker side.
        output = getattr(node, 'workeroutput', {})
        if 'cov_worker_node_id' not in output:
            self.failed_workers.append(node)
            return

        # If worker is not collocated then we must save the data file
        # that it returns to us.
        if 'cov_worker_data' in output:
            data_suffix = '%s.%s.%06d.%s' % (
                socket.gethostname(), os.getpid(),
                random.randint(0, 999999),
                output['cov_worker_node_id']
            )

            cov = coverage.Coverage(source=self.cov_source,
                                    branch=self.cov_branch,
                                    data_suffix=data_suffix,
                                    config_file=self.cov_config)
            cov.start()
            if coverage.version_info < (5, 0):
                data = CoverageData()
                data.read_fileobj(StringIO(output['cov_worker_data']))
                cov.data.update(data)
            else:
                data = CoverageData(no_disk=True)
                data.loads(output['cov_worker_data'])
                cov.get_data().update(data)
            cov.stop()
            cov.save()
            path = output['cov_worker_path']
            self.cov.config.paths['source'].append(path)

        # Record the worker types that contribute to the data file.
        rinfo = node.gateway._rinfo()
        node_desc = self.get_node_desc(rinfo.platform, rinfo.version_info)
        self.node_descs.add(node_desc)

    @_ensure_topdir
    def finish(self):
        """Combines coverage data and sets the list of coverage objects to report on."""

        # Combine all the suffix files into the data file.
        self.cov.stop()
        self.cov.save()
        self.cov = self.combining_cov
        self.cov.load()
        self.cov.combine()
        self.cov.save()


class DistWorker(CovController):
    """Implementation for distributed workers."""

    @_ensure_topdir
    def start(self):

        cleanup()

        # Determine whether we are collocated with master.
        self.is_collocated = (socket.gethostname() == self.config.workerinput['cov_master_host'] and
                              self.topdir == self.config.workerinput['cov_master_topdir'])

        # If we are not collocated then rewrite master paths to worker paths.
        if not self.is_collocated:
            master_topdir = self.config.workerinput['cov_master_topdir']
            worker_topdir = self.topdir
            if self.cov_source is not None:
                self.cov_source = [source.replace(master_topdir, worker_topdir)
                                   for source in self.cov_source]
            self.cov_config = self.cov_config.replace(master_topdir, worker_topdir)

        # Erase any previous data and start coverage.
        self.cov = coverage.Coverage(source=self.cov_source,
                                     branch=self.cov_branch,
                                     data_suffix=True,
                                     config_file=self.cov_config)
        self.cov.start()
        self.set_env()

    @_ensure_topdir
    def finish(self):
        """Stop coverage and send relevant info back to the master."""
        self.unset_env()
        self.cov.stop()

        if self.is_collocated:
            # We don't combine data if we're collocated - we can get
            # race conditions in the .combine() call (it's not atomic)
            # The data is going to be combined in the master.
            self.cov.save()

            # If we are collocated then just inform the master of our
            # data file to indicate that we have finished.
            self.config.workeroutput['cov_worker_node_id'] = self.nodeid
        else:
            self.cov.combine()
            self.cov.save()
            # If we are not collocated then add the current path
            # and coverage data to the output so we can combine
            # it on the master node.

            # Send all the data to the master over the channel.
            if coverage.version_info < (5, 0):
                buff = StringIO()
                self.cov.data.write_fileobj(buff)
                data = buff.getvalue()
            else:
                data = self.cov.get_data().dumps()

            self.config.workeroutput.update({
                'cov_worker_path': self.topdir,
                'cov_worker_node_id': self.nodeid,
                'cov_worker_data': data,
            })

    def summary(self, stream):
        """Only the master reports so do nothing."""

        pass
