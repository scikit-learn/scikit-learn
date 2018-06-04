"""Coverage controllers for use by pytest-cov and nose-cov."""

import os
import random
import socket
import sys

import coverage
from coverage.data import CoverageData

from .compat import StringIO


class CovController(object):
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
        self.node_descs = set()
        self.failed_slaves = []
        self.topdir = os.getcwd()

    def set_env(self):
        """Put info about coverage into the env so that subprocesses can activate coverage."""
        if self.cov_source is None:
            os.environ['COV_CORE_SOURCE'] = ''
        else:
            os.environ['COV_CORE_SOURCE'] = os.pathsep.join(self.cov_source)
        config_file = os.path.abspath(self.cov_config)
        if os.path.exists(config_file):
            os.environ['COV_CORE_CONFIG'] = config_file
        else:
            os.environ['COV_CORE_CONFIG'] = ''
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

    @staticmethod
    def get_node_desc(platform, version_info):
        """Return a description of this node."""

        return 'platform %s, python %s' % (platform, '%s.%s.%s-%s-%s' % version_info[:5])

    @staticmethod
    def sep(stream, s, txt):
        if hasattr(stream, 'sep'):
            stream.sep(s, txt)
        else:
            sep_total = max((70 - 2 - len(txt)), 2)
            sep_len = sep_total // 2
            sep_extra = sep_total % 2
            out = '%s %s %s\n' % (s * sep_len, txt, s * (sep_len + sep_extra))
            stream.write(out)

    def summary(self, stream):
        """Produce coverage reports."""
        total = 0

        if not self.cov_report:
            with open(os.devnull, 'w') as null:
                total = self.cov.report(show_missing=True, ignore_errors=True, file=null)
                return total

        # Output coverage section header.
        if len(self.node_descs) == 1:
            self.sep(stream, '-', 'coverage: %s' % ''.join(self.node_descs))
        else:
            self.sep(stream, '-', 'coverage')
            for node_desc in sorted(self.node_descs):
                self.sep(stream, ' ', '%s' % node_desc)

        # Produce terminal report if wanted.
        if any(x in self.cov_report for x in ['term', 'term-missing']):
            options = {
                'show_missing': ('term-missing' in self.cov_report) or None,
                'ignore_errors': True,
                'file': stream,
            }
            skip_covered = isinstance(self.cov_report, dict) and 'skip-covered' in self.cov_report.values()
            if hasattr(coverage, 'version_info') and coverage.version_info[0] >= 4:
                options.update({'skip_covered': skip_covered or None})
            total = self.cov.report(**options)

        # Produce annotated source code report if wanted.
        if 'annotate' in self.cov_report:
            annotate_dir = self.cov_report['annotate']
            self.cov.annotate(ignore_errors=True, directory=annotate_dir)
            # We need to call Coverage.report here, just to get the total
            # Coverage.annotate don't return any total and we need it for --cov-fail-under.
            total = self.cov.report(ignore_errors=True, file=StringIO())
            if annotate_dir:
                stream.write('Coverage annotated source written to dir %s\n' % annotate_dir)
            else:
                stream.write('Coverage annotated source written next to source\n')

        # Produce html report if wanted.
        if 'html' in self.cov_report:
            total = self.cov.html_report(ignore_errors=True, directory=self.cov_report['html'])
            stream.write('Coverage HTML written to dir %s\n' % self.cov.config.html_dir)

        # Produce xml report if wanted.
        if 'xml' in self.cov_report:
            total = self.cov.xml_report(ignore_errors=True, outfile=self.cov_report['xml'])
            stream.write('Coverage XML written to file %s\n' % self.cov.config.xml_output)

        # Report on any failed slaves.
        if self.failed_slaves:
            self.sep(stream, '-', 'coverage: failed slaves')
            stream.write('The following slaves failed to return coverage data, '
                         'ensure that pytest-cov is installed on these slaves.\n')
            for node in self.failed_slaves:
                stream.write('%s\n' % node.gateway.id)

        return total


class Central(CovController):
    """Implementation for centralised operation."""

    def start(self):
        """Erase any previous coverage data and start coverage."""

        self.cov = coverage.coverage(source=self.cov_source,
                                     branch=self.cov_branch,
                                     config_file=self.cov_config)
        if self.cov_append:
            self.cov.load()
        else:
            self.cov.erase()
        self.cov.start()
        self.set_env()

    def finish(self):
        """Stop coverage, save data to file and set the list of coverage objects to report on."""

        self.unset_env()
        self.cov.stop()
        self.cov.combine()
        self.cov.save()
        node_desc = self.get_node_desc(sys.platform, sys.version_info)
        self.node_descs.add(node_desc)


class DistMaster(CovController):
    """Implementation for distributed master."""

    def start(self):
        """Ensure coverage rc file rsynced if appropriate."""

        if self.cov_config and os.path.exists(self.cov_config):
            self.config.option.rsyncdir.append(self.cov_config)

        self.cov = coverage.coverage(source=self.cov_source,
                                     branch=self.cov_branch,
                                     config_file=self.cov_config)
        if self.cov_append:
            self.cov.load()
        else:
            self.cov.erase()
        self.cov.start()
        self.cov.config.paths['source'] = [self.topdir]

    def configure_node(self, node):
        """Slaves need to know if they are collocated and what files have moved."""

        node.slaveinput['cov_master_host'] = socket.gethostname()
        node.slaveinput['cov_master_topdir'] = self.topdir
        node.slaveinput['cov_master_rsync_roots'] = [str(root) for root in node.nodemanager.roots]

    def testnodedown(self, node, error):
        """Collect data file name from slave."""

        # If slave doesn't return any data then it is likely that this
        # plugin didn't get activated on the slave side.
        if not (hasattr(node, 'slaveoutput') and 'cov_slave_node_id' in node.slaveoutput):
            self.failed_slaves.append(node)
            return

        # If slave is not collocated then we must save the data file
        # that it returns to us.
        if 'cov_slave_data' in node.slaveoutput:
            data_suffix = '%s.%s.%06d.%s' % (
                socket.gethostname(), os.getpid(),
                random.randint(0, 999999),
                node.slaveoutput['cov_slave_node_id']
                )

            cov = coverage.coverage(source=self.cov_source,
                                    branch=self.cov_branch,
                                    data_suffix=data_suffix,
                                    config_file=self.cov_config)
            cov.start()
            if hasattr(self.cov.data, 'read_fileobj'):  # for coverage 4.0
                data = CoverageData()
                data.read_fileobj(StringIO(node.slaveoutput['cov_slave_data']))
                cov.data.update(data)
            else:
                cov.data.lines, cov.data.arcs = node.slaveoutput['cov_slave_data']
            cov.stop()
            cov.save()
            path = node.slaveoutput['cov_slave_path']
            self.cov.config.paths['source'].append(path)

        # Record the slave types that contribute to the data file.
        rinfo = node.gateway._rinfo()
        node_desc = self.get_node_desc(rinfo.platform, rinfo.version_info)
        self.node_descs.add(node_desc)

    def finish(self):
        """Combines coverage data and sets the list of coverage objects to report on."""

        # Combine all the suffix files into the data file.
        self.cov.stop()
        self.cov.combine()
        self.cov.save()


class DistSlave(CovController):
    """Implementation for distributed slaves."""

    def start(self):
        """Determine what data file and suffix to contribute to and start coverage."""

        # Determine whether we are collocated with master.
        self.is_collocated = (socket.gethostname() == self.config.slaveinput['cov_master_host'] and
                              self.topdir == self.config.slaveinput['cov_master_topdir'])

        # If we are not collocated then rewrite master paths to slave paths.
        if not self.is_collocated:
            master_topdir = self.config.slaveinput['cov_master_topdir']
            slave_topdir = self.topdir
            self.cov_source = [source.replace(master_topdir, slave_topdir)
                               for source in self.cov_source]
            self.cov_config = self.cov_config.replace(master_topdir, slave_topdir)

        # Erase any previous data and start coverage.
        self.cov = coverage.coverage(source=self.cov_source,
                                     branch=self.cov_branch,
                                     data_suffix=True,
                                     config_file=self.cov_config)
        if self.cov_append:
            self.cov.load()
        else:
            self.cov.erase()
        self.cov.start()
        self.set_env()

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
            self.config.slaveoutput['cov_slave_node_id'] = self.nodeid
        else:
            self.cov.combine()
            self.cov.save()
            # If we are not collocated then add the current path
            # and coverage data to the output so we can combine
            # it on the master node.

            # Send all the data to the master over the channel.
            self.config.slaveoutput['cov_slave_path'] = self.topdir
            self.config.slaveoutput['cov_slave_node_id'] = self.nodeid
            if hasattr(self.cov.data, 'write_fileobj'):  # for coverage 4.0
                buff = StringIO()
                self.cov.data.write_fileobj(buff)
                self.config.slaveoutput['cov_slave_data'] = buff.getvalue()
            else:
                self.config.slaveoutput['cov_slave_data'] = self.cov.data.lines, self.cov.data.arcs

    def summary(self, stream):
        """Only the master reports so do nothing."""

        pass
