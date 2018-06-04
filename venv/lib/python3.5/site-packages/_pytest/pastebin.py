""" submit failure or test session information to a pastebin service. """
from __future__ import absolute_import, division, print_function

import pytest
import six
import sys
import tempfile


def pytest_addoption(parser):
    group = parser.getgroup("terminal reporting")
    group._addoption('--pastebin', metavar="mode",
                     action='store', dest="pastebin", default=None,
                     choices=['failed', 'all'],
                     help="send failed|all info to bpaste.net pastebin service.")


@pytest.hookimpl(trylast=True)
def pytest_configure(config):
    if config.option.pastebin == "all":
        tr = config.pluginmanager.getplugin('terminalreporter')
        # if no terminal reporter plugin is present, nothing we can do here;
        # this can happen when this function executes in a slave node
        # when using pytest-xdist, for example
        if tr is not None:
            # pastebin file will be utf-8 encoded binary file
            config._pastebinfile = tempfile.TemporaryFile('w+b')
            oldwrite = tr._tw.write

            def tee_write(s, **kwargs):
                oldwrite(s, **kwargs)
                if isinstance(s, six.text_type):
                    s = s.encode('utf-8')
                config._pastebinfile.write(s)

            tr._tw.write = tee_write


def pytest_unconfigure(config):
    if hasattr(config, '_pastebinfile'):
        # get terminal contents and delete file
        config._pastebinfile.seek(0)
        sessionlog = config._pastebinfile.read()
        config._pastebinfile.close()
        del config._pastebinfile
        # undo our patching in the terminal reporter
        tr = config.pluginmanager.getplugin('terminalreporter')
        del tr._tw.__dict__['write']
        # write summary
        tr.write_sep("=", "Sending information to Paste Service")
        pastebinurl = create_new_paste(sessionlog)
        tr.write_line("pastebin session-log: %s\n" % pastebinurl)


def create_new_paste(contents):
    """
    Creates a new paste using bpaste.net service.

    :contents: paste contents as utf-8 encoded bytes
    :returns: url to the pasted contents
    """
    import re
    if sys.version_info < (3, 0):
        from urllib import urlopen, urlencode
    else:
        from urllib.request import urlopen
        from urllib.parse import urlencode

    params = {
        'code': contents,
        'lexer': 'python3' if sys.version_info[0] == 3 else 'python',
        'expiry': '1week',
    }
    url = 'https://bpaste.net'
    response = urlopen(url, data=urlencode(params).encode('ascii')).read()
    m = re.search(r'href="/raw/(\w+)"', response.decode('utf-8'))
    if m:
        return '%s/show/%s' % (url, m.group(1))
    else:
        return 'bad response: ' + response


def pytest_terminal_summary(terminalreporter):
    import _pytest.config
    if terminalreporter.config.option.pastebin != "failed":
        return
    tr = terminalreporter
    if 'failed' in tr.stats:
        terminalreporter.write_sep("=", "Sending information to Paste Service")
        for rep in terminalreporter.stats.get('failed'):
            try:
                msg = rep.longrepr.reprtraceback.reprentries[-1].reprfileloc
            except AttributeError:
                msg = tr._getfailureheadline(rep)
            tw = _pytest.config.create_terminal_writer(terminalreporter.config, stringio=True)
            rep.toterminal(tw)
            s = tw.stringio.getvalue()
            assert len(s)
            pastebinurl = create_new_paste(s)
            tr.write_line("%s --> %s" % (msg, pastebinurl))
