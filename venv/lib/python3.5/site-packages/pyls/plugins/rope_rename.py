# Copyright 2017 Palantir Technologies, Inc.
import logging
import os

from rope.base import libutils
from rope.refactor.rename import Rename

from pyls import hookimpl, uris

log = logging.getLogger(__name__)


@hookimpl
def pyls_rename(config, workspace, document, position, new_name):
    rope_config = config.settings(document_path=document.path).get('rope', {})
    rope_project = workspace._rope_project_builder(rope_config)

    rename = Rename(
        rope_project,
        libutils.path_to_resource(rope_project, document.path),
        document.offset_at_position(position)
    )

    log.debug("Executing rename of %s to %s", document.word_at_position(position), new_name)
    changeset = rename.get_changes(new_name, in_hierarchy=True, docs=True)
    log.debug("Finished rename: %s", changeset.changes)
    return {
        'documentChanges': [{
            'textDocument': {
                'uri': uris.uri_with(
                    document.uri, path=os.path.join(workspace.root_path, change.resource.path)
                ),
                'version': workspace.get_document(document.uri).version
            },
            'edits': [{
                'range': {
                    'start': {'line': 0, 'character': 0},
                    'end': {'line': _num_lines(change.resource), 'character': 0},
                },
                'newText': change.new_contents
            }]
        } for change in changeset.changes]
    }


def _num_lines(resource):
    "Count the number of lines in a `File` resource."
    return len(resource.read().splitlines())
