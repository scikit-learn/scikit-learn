import re
import os

from parso import python_bytes_to_unicode

from jedi.evaluate.cache import evaluator_method_cache
from jedi._compatibility import iter_modules, all_suffixes
from jedi.evaluate.filters import GlobalNameFilter, ContextNameMixin, \
    AbstractNameDefinition, ParserTreeFilter, DictFilter, MergedFilter
from jedi.evaluate import compiled
from jedi.evaluate.base_context import TreeContext
from jedi.evaluate.imports import SubModuleName, infer_import


class _ModuleAttributeName(AbstractNameDefinition):
    """
    For module attributes like __file__, __str__ and so on.
    """
    api_type = u'instance'

    def __init__(self, parent_module, string_name):
        self.parent_context = parent_module
        self.string_name = string_name

    def infer(self):
        return compiled.get_string_context_set(self.parent_context.evaluator)


class ModuleName(ContextNameMixin, AbstractNameDefinition):
    start_pos = 1, 0

    def __init__(self, context, name):
        self._context = context
        self._name = name

    @property
    def string_name(self):
        return self._name


class ModuleContext(TreeContext):
    api_type = u'module'
    parent_context = None

    def __init__(self, evaluator, module_node, path, code_lines):
        super(ModuleContext, self).__init__(evaluator, parent_context=None)
        self.tree_node = module_node
        self._path = path
        self.code_lines = code_lines

    def get_filters(self, search_global, until_position=None, origin_scope=None):
        yield MergedFilter(
            ParserTreeFilter(
                self.evaluator,
                context=self,
                until_position=until_position,
                origin_scope=origin_scope
            ),
            GlobalNameFilter(self, self.tree_node),
        )
        yield DictFilter(self._sub_modules_dict())
        yield DictFilter(self._module_attributes_dict())
        for star_module in self.star_imports():
            yield next(star_module.get_filters(search_global))

    # I'm not sure if the star import cache is really that effective anymore
    # with all the other really fast import caches. Recheck. Also we would need
    # to push the star imports into Evaluator.module_cache, if we reenable this.
    @evaluator_method_cache([])
    def star_imports(self):
        modules = []
        for i in self.tree_node.iter_imports():
            if i.is_star_import():
                name = i.get_paths()[-1][-1]
                new = infer_import(self, name)
                for module in new:
                    if isinstance(module, ModuleContext):
                        modules += module.star_imports()
                modules += new
        return modules

    @evaluator_method_cache()
    def _module_attributes_dict(self):
        names = ['__file__', '__package__', '__doc__', '__name__']
        # All the additional module attributes are strings.
        return dict((n, _ModuleAttributeName(self, n)) for n in names)

    @property
    def _string_name(self):
        """ This is used for the goto functions. """
        if self._path is None:
            return ''  # no path -> empty name
        else:
            sep = (re.escape(os.path.sep),) * 2
            r = re.search(r'([^%s]*?)(%s__init__)?(\.py|\.so)?$' % sep, self._path)
            # Remove PEP 3149 names
            return re.sub(r'\.[a-z]+-\d{2}[mud]{0,3}$', '', r.group(1))

    @property
    @evaluator_method_cache()
    def name(self):
        return ModuleName(self, self._string_name)

    def _get_init_directory(self):
        """
        :return: The path to the directory of a package. None in case it's not
                 a package.
        """
        for suffix in all_suffixes():
            ending = '__init__' + suffix
            py__file__ = self.py__file__()
            if py__file__ is not None and py__file__.endswith(ending):
                # Remove the ending, including the separator.
                return self.py__file__()[:-len(ending) - 1]
        return None

    def py__name__(self):
        for name, module in self.evaluator.module_cache.iterate_modules_with_names():
            if module == self and name != '':
                return name

        return '__main__'

    def py__file__(self):
        """
        In contrast to Python's __file__ can be None.
        """
        if self._path is None:
            return None

        return os.path.abspath(self._path)

    def py__package__(self):
        if self._get_init_directory() is None:
            return re.sub(r'\.?[^.]+$', '', self.py__name__())
        else:
            return self.py__name__()

    def _py__path__(self):
        search_path = self.evaluator.get_sys_path()
        init_path = self.py__file__()
        if os.path.basename(init_path) == '__init__.py':
            with open(init_path, 'rb') as f:
                content = python_bytes_to_unicode(f.read(), errors='replace')
                # these are strings that need to be used for namespace packages,
                # the first one is ``pkgutil``, the second ``pkg_resources``.
                options = ('declare_namespace(__name__)', 'extend_path(__path__')
                if options[0] in content or options[1] in content:
                    # It is a namespace, now try to find the rest of the
                    # modules on sys_path or whatever the search_path is.
                    paths = set()
                    for s in search_path:
                        other = os.path.join(s, self.name.string_name)
                        if os.path.isdir(other):
                            paths.add(other)
                    if paths:
                        return list(paths)
                    # TODO I'm not sure if this is how nested namespace
                    # packages work. The tests are not really good enough to
                    # show that.
        # Default to this.
        return [self._get_init_directory()]

    @property
    def py__path__(self):
        """
        Not seen here, since it's a property. The callback actually uses a
        variable, so use it like::

            foo.py__path__(sys_path)

        In case of a package, this returns Python's __path__ attribute, which
        is a list of paths (strings).
        Raises an AttributeError if the module is not a package.
        """
        path = self._get_init_directory()

        if path is None:
            raise AttributeError('Only packages have __path__ attributes.')
        else:
            return self._py__path__

    @evaluator_method_cache()
    def _sub_modules_dict(self):
        """
        Lists modules in the directory of this module (if this module is a
        package).
        """
        path = self._path
        names = {}
        if path is not None and path.endswith(os.path.sep + '__init__.py'):
            mods = iter_modules([os.path.dirname(path)])
            for module_loader, name, is_pkg in mods:
                # It's obviously a relative import to the current module.
                names[name] = SubModuleName(self, name)

        # TODO add something like this in the future, its cleaner than the
        #   import hacks.
        # ``os.path`` is a hardcoded exception, because it's a
        # ``sys.modules`` modification.
        # if str(self.name) == 'os':
        #     names.append(Name('path', parent_context=self))

        return names

    def py__class__(self):
        return compiled.get_special_object(self.evaluator, u'MODULE_CLASS')

    def __repr__(self):
        return "<%s: %s@%s-%s>" % (
            self.__class__.__name__, self._string_name,
            self.tree_node.start_pos[0], self.tree_node.end_pos[0])
