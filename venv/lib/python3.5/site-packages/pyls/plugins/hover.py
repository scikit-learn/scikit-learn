# Copyright 2017 Palantir Technologies, Inc.
import logging
from pyls import hookimpl, _utils

log = logging.getLogger(__name__)


@hookimpl
def pyls_hover(document, position):
    definitions = document.jedi_script(position).goto_definitions()
    word = document.word_at_position(position)

    # Find an exact match for a completion
    definitions = [d for d in definitions if d.name == word]

    if not definitions:
        # :(
        return {'contents': ''}

    return {'contents': _utils.format_docstring(definitions[0].docstring()) or ""}
