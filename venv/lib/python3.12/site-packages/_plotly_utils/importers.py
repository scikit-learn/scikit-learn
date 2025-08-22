import importlib


def relative_import(parent_name, rel_modules=(), rel_classes=()):
    """
    Helper function to import submodules lazily in Python 3.7+

    Parameters
    ----------
    rel_modules: list of str
        list of submodules to import, of the form .submodule
    rel_classes: list of str
        list of submodule classes/variables to import, of the form ._submodule.Foo

    Returns
    -------
    tuple
        Tuple that should be assigned to __all__, __getattr__ in the caller
    """
    module_names = {rel_module.split(".")[-1]: rel_module for rel_module in rel_modules}
    class_names = {rel_path.split(".")[-1]: rel_path for rel_path in rel_classes}

    def __getattr__(import_name):
        # In Python 3.7+, lazy import submodules

        # Check for submodule
        if import_name in module_names:
            rel_import = module_names[import_name]
            return importlib.import_module(rel_import, parent_name)

        # Check for submodule class
        if import_name in class_names:
            rel_path_parts = class_names[import_name].split(".")
            rel_module = ".".join(rel_path_parts[:-1])
            class_name = import_name
            class_module = importlib.import_module(rel_module, parent_name)
            return getattr(class_module, class_name)

        raise AttributeError(
            "module {__name__!r} has no attribute {name!r}".format(
                name=import_name, __name__=parent_name
            )
        )

    __all__ = list(module_names) + list(class_names)

    def __dir__():
        return __all__

    return __all__, __getattr__, __dir__
