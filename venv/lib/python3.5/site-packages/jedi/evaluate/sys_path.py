import os

from jedi._compatibility import unicode, force_unicode, all_suffixes
from jedi.evaluate.cache import evaluator_method_cache
from jedi.evaluate.base_context import ContextualizedNode
from jedi.evaluate.helpers import is_string
from jedi.common.utils import traverse_parents
from jedi.parser_utils import get_cached_code_lines
from jedi import settings
from jedi import debug


def _abs_path(module_context, path):
    if os.path.isabs(path):
        return path

    module_path = module_context.py__file__()
    if module_path is None:
        # In this case we have no idea where we actually are in the file
        # system.
        return None

    base_dir = os.path.dirname(module_path)
    path = force_unicode(path)
    return os.path.abspath(os.path.join(base_dir, path))


def _paths_from_assignment(module_context, expr_stmt):
    """
    Extracts the assigned strings from an assignment that looks as follows::

        sys.path[0:0] = ['module/path', 'another/module/path']

    This function is in general pretty tolerant (and therefore 'buggy').
    However, it's not a big issue usually to add more paths to Jedi's sys_path,
    because it will only affect Jedi in very random situations and by adding
    more paths than necessary, it usually benefits the general user.
    """
    for assignee, operator in zip(expr_stmt.children[::2], expr_stmt.children[1::2]):
        try:
            assert operator in ['=', '+=']
            assert assignee.type in ('power', 'atom_expr') and \
                len(assignee.children) > 1
            c = assignee.children
            assert c[0].type == 'name' and c[0].value == 'sys'
            trailer = c[1]
            assert trailer.children[0] == '.' and trailer.children[1].value == 'path'
            # TODO Essentially we're not checking details on sys.path
            # manipulation. Both assigment of the sys.path and changing/adding
            # parts of the sys.path are the same: They get added to the end of
            # the current sys.path.
            """
            execution = c[2]
            assert execution.children[0] == '['
            subscript = execution.children[1]
            assert subscript.type == 'subscript'
            assert ':' in subscript.children
            """
        except AssertionError:
            continue

        cn = ContextualizedNode(module_context.create_context(expr_stmt), expr_stmt)
        for lazy_context in cn.infer().iterate(cn):
            for context in lazy_context.infer():
                if is_string(context):
                    abs_path = _abs_path(module_context, context.get_safe_value())
                    if abs_path is not None:
                        yield abs_path


def _paths_from_list_modifications(module_context, trailer1, trailer2):
    """ extract the path from either "sys.path.append" or "sys.path.insert" """
    # Guarantee that both are trailers, the first one a name and the second one
    # a function execution with at least one param.
    if not (trailer1.type == 'trailer' and trailer1.children[0] == '.'
            and trailer2.type == 'trailer' and trailer2.children[0] == '('
            and len(trailer2.children) == 3):
        return

    name = trailer1.children[1].value
    if name not in ['insert', 'append']:
        return
    arg = trailer2.children[1]
    if name == 'insert' and len(arg.children) in (3, 4):  # Possible trailing comma.
        arg = arg.children[2]

    for context in module_context.create_context(arg).eval_node(arg):
        if is_string(context):
            abs_path = _abs_path(module_context, context.get_safe_value())
            if abs_path is not None:
                yield abs_path


@evaluator_method_cache(default=[])
def check_sys_path_modifications(module_context):
    """
    Detect sys.path modifications within module.
    """
    def get_sys_path_powers(names):
        for name in names:
            power = name.parent.parent
            if power.type in ('power', 'atom_expr'):
                c = power.children
                if c[0].type == 'name' and c[0].value == 'sys' \
                        and c[1].type == 'trailer':
                    n = c[1].children[1]
                    if n.type == 'name' and n.value == 'path':
                        yield name, power

    if module_context.tree_node is None:
        return []

    added = []
    try:
        possible_names = module_context.tree_node.get_used_names()['path']
    except KeyError:
        pass
    else:
        for name, power in get_sys_path_powers(possible_names):
            expr_stmt = power.parent
            if len(power.children) >= 4:
                added.extend(
                    _paths_from_list_modifications(
                        module_context, *power.children[2:4]
                    )
                )
            elif expr_stmt is not None and expr_stmt.type == 'expr_stmt':
                added.extend(_paths_from_assignment(module_context, expr_stmt))
    return added


def discover_buildout_paths(evaluator, script_path):
    buildout_script_paths = set()

    for buildout_script_path in _get_buildout_script_paths(script_path):
        for path in _get_paths_from_buildout_script(evaluator, buildout_script_path):
            buildout_script_paths.add(path)

    return buildout_script_paths


def _get_paths_from_buildout_script(evaluator, buildout_script_path):
    try:
        module_node = evaluator.parse(
            path=buildout_script_path,
            cache=True,
            cache_path=settings.cache_directory
        )
    except IOError:
        debug.warning('Error trying to read buildout_script: %s', buildout_script_path)
        return

    from jedi.evaluate.context import ModuleContext
    module = ModuleContext(
        evaluator, module_node, buildout_script_path,
        code_lines=get_cached_code_lines(evaluator.grammar, buildout_script_path),
    )
    for path in check_sys_path_modifications(module):
        yield path


def _get_parent_dir_with_file(path, filename):
    for parent in traverse_parents(path):
        if os.path.isfile(os.path.join(parent, filename)):
            return parent
    return None


def _get_buildout_script_paths(search_path):
    """
    if there is a 'buildout.cfg' file in one of the parent directories of the
    given module it will return a list of all files in the buildout bin
    directory that look like python files.

    :param search_path: absolute path to the module.
    :type search_path: str
    """
    project_root = _get_parent_dir_with_file(search_path, 'buildout.cfg')
    if not project_root:
        return
    bin_path = os.path.join(project_root, 'bin')
    if not os.path.exists(bin_path):
        return

    for filename in os.listdir(bin_path):
        try:
            filepath = os.path.join(bin_path, filename)
            with open(filepath, 'r') as f:
                firstline = f.readline()
                if firstline.startswith('#!') and 'python' in firstline:
                    yield filepath
        except (UnicodeDecodeError, IOError) as e:
            # Probably a binary file; permission error or race cond. because
            # file got deleted. Ignore it.
            debug.warning(unicode(e))
            continue


def dotted_path_in_sys_path(sys_path, module_path):
    """
    Returns the dotted path inside a sys.path.
    """
    # First remove the suffix.
    for suffix in all_suffixes():
        if module_path.endswith(suffix):
            module_path = module_path[:-len(suffix)]
        break
    else:
        # There should always be a suffix in a valid Python file on the path.
        return None

    if module_path.startswith(os.path.sep):
        # The paths in sys.path most of the times don't end with a slash.
        module_path = module_path[1:]

    for p in sys_path:
        if module_path.startswith(p):
            rest = module_path[len(p):]
            if rest:
                split = rest.split(os.path.sep)
                for string in split:
                    if not string or '.' in string:
                        return None
                return '.'.join(split)

    return None
