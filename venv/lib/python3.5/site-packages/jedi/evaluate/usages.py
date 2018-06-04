from jedi.evaluate import imports
from jedi.evaluate.filters import TreeNameDefinition
from jedi.evaluate.context import ModuleContext


def _resolve_names(definition_names, avoid_names=()):
    for name in definition_names:
        if name in avoid_names:
            # Avoiding recursions here, because goto on a module name lands
            # on the same module.
            continue

        if not isinstance(name, imports.SubModuleName):
            # SubModuleNames are not actually existing names but created
            # names when importing something like `import foo.bar.baz`.
            yield name

        if name.api_type == 'module':
            for name in _resolve_names(name.goto(), definition_names):
                yield name


def _dictionarize(names):
    return dict(
        (n if n.tree_name is None else n.tree_name, n)
        for n in names
    )


def _find_names(module_context, tree_name):
    context = module_context.create_context(tree_name)
    name = TreeNameDefinition(context, tree_name)
    found_names = set(name.goto())
    found_names.add(name)
    return _dictionarize(_resolve_names(found_names))


def usages(module_context, tree_name):
    search_name = tree_name.value
    found_names = _find_names(module_context, tree_name)
    modules = set(d.get_root_context() for d in found_names.values())
    modules = set(m for m in modules if isinstance(m, ModuleContext))

    non_matching_usage_maps = {}
    for m in imports.get_modules_containing_name(module_context.evaluator, modules, search_name):
        for name_leaf in m.tree_node.get_used_names().get(search_name, []):
            new = _find_names(m, name_leaf)
            if any(tree_name in found_names for tree_name in new):
                found_names.update(new)
                for tree_name in new:
                    for dct in non_matching_usage_maps.get(tree_name, []):
                        # A usage that was previously searched for matches with
                        # a now found name. Merge.
                        found_names.update(dct)
                    try:
                        del non_matching_usage_maps[tree_name]
                    except KeyError:
                        pass
            else:
                for name in new:
                    non_matching_usage_maps.setdefault(name, []).append(new)
    return found_names.values()
