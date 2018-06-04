# Copyright 2017 Palantir Technologies, Inc.
import logging
import re
from pyls import hookimpl, _utils

log = logging.getLogger(__name__)

SPHINX = re.compile(r"\s*:param\s+(?P<param>\w+):\s*(?P<doc>[^\n]+)")
EPYDOC = re.compile(r"\s*@param\s+(?P<param>\w+):\s*(?P<doc>[^\n]+)")
GOOGLE = re.compile(r"\s*(?P<param>\w+).*:\s*(?P<doc>[^\n]+)")

DOC_REGEX = [SPHINX, EPYDOC, GOOGLE]


@hookimpl
def pyls_signature_help(document, position):
    signatures = document.jedi_script(position).call_signatures()

    if not signatures:
        return {'signatures': []}

    s = signatures[0]
    sig = {
        'label': s.docstring().splitlines()[0],
        'documentation': _utils.format_docstring(s.docstring(raw=True))
    }

    # If there are params, add those
    if s.params:
        sig['parameters'] = [{
            'label': p.name,
            'documentation': _param_docs(s.docstring(), p.name)
        } for p in s.params]

    # We only return a single signature because Python doesn't allow overloading
    sig_info = {'signatures': [sig], 'activeSignature': 0}

    if s.index is not None and s.params:
        # Then we know which parameter we're looking at
        sig_info['activeParameter'] = s.index

    return sig_info


def _param_docs(docstring, param_name):
    for line in docstring.splitlines():
        for regex in DOC_REGEX:
            m = regex.match(line)
            if not m:
                continue
            if m.group('param') != param_name:
                continue
            return m.group('doc') or ""
