import rope.base.codeanalyze
import rope.base.evaluate
import rope.base.pyobjects
from rope.base import taskhandle, exceptions, worder
from rope.contrib import fixsyntax
from rope.refactor import occurrences


def find_occurrences(project, resource, offset, unsure=False, resources=None,
                     in_hierarchy=False,
                     task_handle=taskhandle.NullTaskHandle()):
    """Return a list of `Location`\s

    If `unsure` is `True`, possible matches are returned, too.  You
    can use `Location.unsure` to see which are unsure occurrences.
    `resources` can be a list of `rope.base.resource.File`\s that
    should be searched for occurrences; if `None` all python files
    in the project are searched.

    """
    name = worder.get_name_at(resource, offset)
    this_pymodule = project.get_pymodule(resource)
    primary, pyname = rope.base.evaluate.eval_location2(
        this_pymodule, offset)

    def is_match(occurrence):
        return unsure
    finder = occurrences.create_finder(
        project, name, pyname, unsure=is_match,
        in_hierarchy=in_hierarchy, instance=primary)
    if resources is None:
        resources = project.get_python_files()
    job_set = task_handle.create_jobset('Finding Occurrences',
                                        count=len(resources))
    return _find_locations(finder, resources, job_set)


def find_implementations(project, resource, offset, resources=None,
                         task_handle=taskhandle.NullTaskHandle()):
    """Find the places a given method is overridden.

    Finds the places a method is implemented.  Returns a list of
    `Location`\s.
    """
    name = worder.get_name_at(resource, offset)
    this_pymodule = project.get_pymodule(resource)
    pyname = rope.base.evaluate.eval_location(this_pymodule, offset)
    if pyname is not None:
        pyobject = pyname.get_object()
        if not isinstance(pyobject, rope.base.pyobjects.PyFunction) or \
           pyobject.get_kind() != 'method':
            raise exceptions.BadIdentifierError('Not a method!')
    else:
        raise exceptions.BadIdentifierError('Cannot resolve the identifier!')

    def is_defined(occurrence):
        if not occurrence.is_defined():
            return False

    def not_self(occurrence):
        if occurrence.get_pyname().get_object() == pyname.get_object():
            return False
    filters = [is_defined, not_self,
               occurrences.InHierarchyFilter(pyname, True)]
    finder = occurrences.Finder(project, name, filters=filters)
    if resources is None:
        resources = project.get_python_files()
    job_set = task_handle.create_jobset('Finding Implementations',
                                        count=len(resources))
    return _find_locations(finder, resources, job_set)


def find_definition(project, code, offset, resource=None, maxfixes=1):
    """Return the definition location of the python name at `offset`

    A `Location` object is returned if the definition location can be
    determined, otherwise ``None`` is returned.
    """
    fixer = fixsyntax.FixSyntax(project, code, resource, maxfixes)
    pyname = fixer.pyname_at(offset)
    if pyname is not None:
        module, lineno = pyname.get_definition_location()
        name = rope.base.worder.Worder(code).get_word_at(offset)
        if lineno is not None:
            start = module.lines.get_line_start(lineno)

            def check_offset(occurrence):
                if occurrence.offset < start:
                    return False
            pyname_filter = occurrences.PyNameFilter(pyname)
            finder = occurrences.Finder(project, name,
                                        [check_offset, pyname_filter])
            for occurrence in finder.find_occurrences(pymodule=module):
                return Location(occurrence)


class Location(object):

    def __init__(self, occurrence):
        self.resource = occurrence.resource
        self.region = occurrence.get_word_range()
        self.offset = self.region[0]
        self.unsure = occurrence.is_unsure()
        self.lineno = occurrence.lineno


def _find_locations(finder, resources, job_set):
    result = []
    for resource in resources:
        job_set.started_job(resource.path)
        for occurrence in finder.find_occurrences(resource):
            result.append(Location(occurrence))
        job_set.finished_job()
    return result
