import os


class ResourceObserver(object):
    """Provides the interface for observing resources

    `ResourceObserver`\s can be registered using `Project.
    add_observer()`.  But most of the time `FilteredResourceObserver`
    should be used.  `ResourceObserver`\s report all changes passed
    to them and they don't report changes to all resources.  For
    example if a folder is removed, it only calls `removed()` for that
    folder and not its contents.  You can use
    `FilteredResourceObserver` if you are interested in changes only
    to a list of resources.  And you want changes to be reported on
    individual resources.

    """

    def __init__(self, changed=None, moved=None, created=None,
                 removed=None, validate=None):
        self.changed = changed
        self.moved = moved
        self.created = created
        self.removed = removed
        self._validate = validate

    def resource_changed(self, resource):
        """It is called when the resource changes"""
        if self.changed is not None:
            self.changed(resource)

    def resource_moved(self, resource, new_resource):
        """It is called when a resource is moved"""
        if self.moved is not None:
            self.moved(resource, new_resource)

    def resource_created(self, resource):
        """Is called when a new resource is created"""
        if self.created is not None:
            self.created(resource)

    def resource_removed(self, resource):
        """Is called when a new resource is removed"""
        if self.removed is not None:
            self.removed(resource)

    def validate(self, resource):
        """Validate the existence of this resource and its children.

        This function is called when rope need to update its resource
        cache about the files that might have been changed or removed
        by other processes.

        """
        if self._validate is not None:
            self._validate(resource)


class FilteredResourceObserver(object):
    """A useful decorator for `ResourceObserver`

    Most resource observers have a list of resources and are
    interested only in changes to those files.  This class satisfies
    this need.  It dispatches resource changed and removed messages.
    It performs these tasks:

    * Changes to files and folders are analyzed to check whether any
      of the interesting resources are changed or not.  If they are,
      it reports these changes to `resource_observer` passed to the
      constructor.
    * When a resource is removed it checks whether any of the
      interesting resources are contained in that folder and reports
      them to `resource_observer`.
    * When validating a folder it validates all of the interesting
      files in that folder.

    Since most resource observers are interested in a list of
    resources that change over time, `add_resource` and
    `remove_resource` might be useful.

    """

    def __init__(self, resource_observer, initial_resources=None,
                 timekeeper=None):
        self.observer = resource_observer
        self.resources = {}
        if timekeeper is not None:
            self.timekeeper = timekeeper
        else:
            self.timekeeper = ChangeIndicator()
        if initial_resources is not None:
            for resource in initial_resources:
                self.add_resource(resource)

    def add_resource(self, resource):
        """Add a resource to the list of interesting resources"""
        if resource.exists():
            self.resources[resource] = self.timekeeper.get_indicator(resource)
        else:
            self.resources[resource] = None

    def remove_resource(self, resource):
        """Add a resource to the list of interesting resources"""
        if resource in self.resources:
            del self.resources[resource]

    def clear_resources(self):
        """Removes all registered resources"""
        self.resources.clear()

    def resource_changed(self, resource):
        changes = _Changes()
        self._update_changes_caused_by_changed(changes, resource)
        self._perform_changes(changes)

    def _update_changes_caused_by_changed(self, changes, changed):
        if changed in self.resources:
            changes.add_changed(changed)
        if self._is_parent_changed(changed):
            changes.add_changed(changed.parent)

    def _update_changes_caused_by_moved(self, changes, resource,
                                        new_resource=None):
        if resource in self.resources:
            changes.add_removed(resource, new_resource)
        if new_resource in self.resources:
            changes.add_created(new_resource)
        if resource.is_folder():
            for file in list(self.resources):
                if resource.contains(file):
                    new_file = self._calculate_new_resource(
                        resource, new_resource, file)
                    changes.add_removed(file, new_file)
        if self._is_parent_changed(resource):
            changes.add_changed(resource.parent)
        if new_resource is not None:
            if self._is_parent_changed(new_resource):
                changes.add_changed(new_resource.parent)

    def _is_parent_changed(self, child):
        return child.parent in self.resources

    def resource_moved(self, resource, new_resource):
        changes = _Changes()
        self._update_changes_caused_by_moved(changes, resource, new_resource)
        self._perform_changes(changes)

    def resource_created(self, resource):
        changes = _Changes()
        self._update_changes_caused_by_created(changes, resource)
        self._perform_changes(changes)

    def _update_changes_caused_by_created(self, changes, resource):
        if resource in self.resources:
            changes.add_created(resource)
        if self._is_parent_changed(resource):
            changes.add_changed(resource.parent)

    def resource_removed(self, resource):
        changes = _Changes()
        self._update_changes_caused_by_moved(changes, resource)
        self._perform_changes(changes)

    def _perform_changes(self, changes):
        for resource in changes.changes:
            self.observer.resource_changed(resource)
            self.resources[resource] = self.timekeeper.get_indicator(resource)
        for resource, new_resource in changes.moves.items():
            self.resources[resource] = None
            if new_resource is not None:
                self.observer.resource_moved(resource, new_resource)
            else:
                self.observer.resource_removed(resource)
        for resource in changes.creations:
            self.observer.resource_created(resource)
            self.resources[resource] = self.timekeeper.get_indicator(resource)

    def validate(self, resource):
        changes = _Changes()
        for file in self._search_resource_moves(resource):
            if file in self.resources:
                self._update_changes_caused_by_moved(changes, file)
        for file in self._search_resource_changes(resource):
            if file in self.resources:
                self._update_changes_caused_by_changed(changes, file)
        for file in self._search_resource_creations(resource):
            if file in self.resources:
                changes.add_created(file)
        self._perform_changes(changes)

    def _search_resource_creations(self, resource):
        creations = set()
        if resource in self.resources and resource.exists() and \
           self.resources[resource] is None:
            creations.add(resource)
        if resource.is_folder():
            for file in self.resources:
                if file.exists() and resource.contains(file) and \
                   self.resources[file] is None:
                    creations.add(file)
        return creations

    def _search_resource_moves(self, resource):
        all_moved = set()
        if resource in self.resources and not resource.exists():
            all_moved.add(resource)
        if resource.is_folder():
            for file in self.resources:
                if resource.contains(file):
                    if not file.exists():
                        all_moved.add(file)
        moved = set(all_moved)
        for folder in [file for file in all_moved if file.is_folder()]:
            if folder in moved:
                for file in list(moved):
                    if folder.contains(file):
                        moved.remove(file)
        return moved

    def _search_resource_changes(self, resource):
        changed = set()
        if resource in self.resources and self._is_changed(resource):
            changed.add(resource)
        if resource.is_folder():
            for file in self.resources:
                if file.exists() and resource.contains(file):
                    if self._is_changed(file):
                        changed.add(file)
        return changed

    def _is_changed(self, resource):
        if self.resources[resource] is None:
            return False
        return self.resources[resource] != \
            self.timekeeper.get_indicator(resource)

    def _calculate_new_resource(self, main, new_main, resource):
        if new_main is None:
            return None
        diff = resource.path[len(main.path):]
        return resource.project.get_resource(new_main.path + diff)


class ChangeIndicator(object):

    def get_indicator(self, resource):
        """Return the modification time and size of a `Resource`."""
        path = resource.real_path
        # on dos, mtime does not change for a folder when files are added
        if os.name != 'posix' and os.path.isdir(path):
            return (os.path.getmtime(path),
                    len(os.listdir(path)),
                    os.path.getsize(path))
        return (os.path.getmtime(path),
                os.path.getsize(path))


class _Changes(object):

    def __init__(self):
        self.changes = set()
        self.creations = set()
        self.moves = {}

    def add_changed(self, resource):
        self.changes.add(resource)

    def add_removed(self, resource, new_resource=None):
        self.moves[resource] = new_resource

    def add_created(self, resource):
        self.creations.add(resource)
