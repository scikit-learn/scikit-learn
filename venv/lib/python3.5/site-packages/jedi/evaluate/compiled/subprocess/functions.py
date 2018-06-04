import sys
import os

from jedi._compatibility import find_module, cast_path, force_unicode, \
    iter_modules, all_suffixes, print_to_stderr
from jedi.evaluate.compiled import access
from jedi import parser_utils


def get_sys_path():
    return list(map(cast_path, sys.path))


def load_module(evaluator, **kwargs):
    return access.load_module(evaluator, **kwargs)


def get_compiled_method_return(evaluator, id, attribute, *args, **kwargs):
    handle = evaluator.compiled_subprocess.get_access_handle(id)
    return getattr(handle.access, attribute)(*args, **kwargs)


def get_special_object(evaluator, identifier):
    return access.get_special_object(evaluator, identifier)


def create_simple_object(evaluator, obj):
    return access.create_access_path(evaluator, obj)


def get_module_info(evaluator, sys_path=None, full_name=None, **kwargs):
    if sys_path is not None:
        sys.path, temp = sys_path, sys.path
    try:
        module_file, module_path, is_pkg = find_module(full_name=full_name, **kwargs)
    except ImportError:
        return None, None, None
    finally:
        if sys_path is not None:
            sys.path = temp

    code = None
    if is_pkg:
        # In this case, we don't have a file yet. Search for the
        # __init__ file.
        if module_path.endswith(('.zip', '.egg')):
            code = module_file.loader.get_source(full_name)
        else:
            module_path = _get_init_path(module_path)
    elif module_file:
        if module_path.endswith(('.zip', '.egg')):
            # Unfortunately we are reading unicode here already, not byes.
            # It seems however hard to get bytes, because the zip importer
            # logic just unpacks the zip file and returns a file descriptor
            # that we cannot as easily access. Therefore we just read it as
            # a string.
            code = module_file.read()
        else:
            # Read the code with a binary file, because the binary file
            # might not be proper unicode. This is handled by the parser
            # wrapper.
            with open(module_path, 'rb') as f:
                code = f.read()

        module_file.close()

    return code, cast_path(module_path), is_pkg


def list_module_names(evaluator, search_path):
    return [
        name
        for module_loader, name, is_pkg in iter_modules(search_path)
    ]


def get_builtin_module_names(evaluator):
    return list(map(force_unicode, sys.builtin_module_names))


def _test_raise_error(evaluator, exception_type):
    """
    Raise an error to simulate certain problems for unit tests.
    """
    raise exception_type


def _test_print(evaluator, stderr=None, stdout=None):
    """
    Force some prints in the subprocesses. This exists for unit tests.
    """
    if stderr is not None:
        print_to_stderr(stderr)
        sys.stderr.flush()
    if stdout is not None:
        print(stdout)
        sys.stdout.flush()


def _get_init_path(directory_path):
    """
    The __init__ file can be searched in a directory. If found return it, else
    None.
    """
    for suffix in all_suffixes():
        path = os.path.join(directory_path, '__init__' + suffix)
        if os.path.exists(path):
            return path
    return None


def safe_literal_eval(evaluator, value):
    return parser_utils.safe_literal_eval(value)
