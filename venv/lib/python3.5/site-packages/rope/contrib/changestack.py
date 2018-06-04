"""For performing many refactorings as a single command

`changestack` module can be used to perform many refactorings on top
of each other as one bigger command.  It can be used like::

  stack = ChangeStack(project, 'my big command')

  #..
  stack.push(refactoring1.get_changes())
  #..
  stack.push(refactoring2.get_changes())
  #..
  stack.push(refactoringX.get_changes())

  stack.pop_all()
  changes = stack.merged()

Now `changes` can be previewed or performed as before.
"""

from rope.base import change


class ChangeStack(object):

    def __init__(self, project, description='merged changes'):
        self.project = project
        self.description = description
        self.stack = []

    def push(self, changes):
        self.stack.append(changes)
        self.project.do(changes)

    def pop_all(self):
        for i in range(len(self.stack)):
            self.project.history.undo(drop=True)

    def merged(self):
        result = change.ChangeSet(self.description)
        for changes in self.stack:
            for c in self._basic_changes(changes):
                result.add_change(c)
        return result

    def _basic_changes(self, changes):
        if isinstance(changes, change.ChangeSet):
            for child in changes.changes:
                for atom in self._basic_changes(child):
                    yield atom
        else:
            yield changes
