import copy

from .. import Options


def backup_Options():
    backup = {}
    for name, value in vars(Options).items():
        # we need a deep copy of _directive_defaults, because they can be changed
        if name == '_directive_defaults':
            value = copy.deepcopy(value)
        backup[name] = value
    return backup


def restore_Options(backup):
    no_value = object()
    for name, orig_value in backup.items():
        if getattr(Options, name, no_value) != orig_value:
            setattr(Options, name, orig_value)
    # strip Options from new keys that might have been added:
    for name in vars(Options).keys():
        if name not in backup:
            delattr(Options, name)


def check_global_options(expected_options, white_list=[]):
    """
    returns error message of "" if check Ok
    """
    no_value = object()
    for name, orig_value in expected_options.items():
        if name not in white_list:
            if getattr(Options, name, no_value) != orig_value:
                return "error in option " + name
    return ""
