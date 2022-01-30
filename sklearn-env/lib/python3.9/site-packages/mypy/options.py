from mypy.backports import OrderedDict
import re
import pprint
import sys

from typing_extensions import Final, TYPE_CHECKING
from typing import Dict, List, Mapping, Optional, Pattern, Set, Tuple, Callable, Any

from mypy import defaults
from mypy.util import get_class_descriptors, replace_object_state

if TYPE_CHECKING:
    from mypy.errors import ErrorCode


class BuildType:
    STANDARD: Final = 0
    MODULE: Final = 1
    PROGRAM_TEXT: Final = 2


PER_MODULE_OPTIONS: Final = {
    # Please keep this list sorted
    "allow_redefinition",
    "allow_untyped_globals",
    "always_false",
    "always_true",
    "check_untyped_defs",
    "debug_cache",
    "disallow_any_decorated",
    "disallow_any_explicit",
    "disallow_any_expr",
    "disallow_any_generics",
    "disallow_any_unimported",
    "disallow_incomplete_defs",
    "disallow_subclassing_any",
    "disallow_untyped_calls",
    "disallow_untyped_decorators",
    "disallow_untyped_defs",
    "follow_imports",
    "follow_imports_for_stubs",
    "ignore_errors",
    "ignore_missing_imports",
    "implicit_reexport",
    "local_partial_types",
    "mypyc",
    "no_implicit_optional",
    "show_none_errors",
    "strict_equality",
    "strict_optional",
    "strict_optional_whitelist",
    "warn_no_return",
    "warn_return_any",
    "warn_unreachable",
    "warn_unused_ignores",
}

OPTIONS_AFFECTING_CACHE: Final = (PER_MODULE_OPTIONS | {"platform", "bazel", "plugins"}) - {
    "debug_cache"
}


class Options:
    """Options collected from flags."""

    def __init__(self) -> None:
        # Cache for clone_for_module()
        self._per_module_cache: Optional[Dict[str, Options]] = None

        # -- build options --
        self.build_type = BuildType.STANDARD
        self.python_version: Tuple[int, int] = sys.version_info[:2]
        # The executable used to search for PEP 561 packages. If this is None,
        # then mypy does not search for PEP 561 packages.
        self.python_executable: Optional[str] = sys.executable
        self.platform = sys.platform
        self.custom_typing_module: Optional[str] = None
        self.custom_typeshed_dir: Optional[str] = None
        self.mypy_path: List[str] = []
        self.report_dirs: Dict[str, str] = {}
        # Show errors in PEP 561 packages/site-packages modules
        self.no_silence_site_packages = False
        self.no_site_packages = False
        self.ignore_missing_imports = False
        # Is ignore_missing_imports set in a per-module section
        self.ignore_missing_imports_per_module = False
        self.follow_imports = 'normal'  # normal|silent|skip|error
        # Whether to respect the follow_imports setting even for stub files.
        # Intended to be used for disabling specific stubs.
        self.follow_imports_for_stubs = False
        # PEP 420 namespace packages
        # This allows definitions of packages without __init__.py and allows packages to span
        # multiple directories. This flag affects both import discovery and the association of
        # input files/modules/packages to the relevant file and fully qualified module name.
        self.namespace_packages = False
        # Use current directory and MYPYPATH to determine fully qualified module names of files
        # passed by automatically considering their subdirectories as packages. This is only
        # relevant if namespace packages are enabled, since otherwise examining __init__.py's is
        # sufficient to determine module names for files. As a possible alternative, add a single
        # top-level __init__.py to your packages.
        self.explicit_package_bases = False
        # File names, directory names or subpaths to avoid checking
        self.exclude: List[str] = []

        # disallow_any options
        self.disallow_any_generics = False
        self.disallow_any_unimported = False
        self.disallow_any_expr = False
        self.disallow_any_decorated = False
        self.disallow_any_explicit = False

        # Disallow calling untyped functions from typed ones
        self.disallow_untyped_calls = False

        # Disallow defining untyped (or incompletely typed) functions
        self.disallow_untyped_defs = False

        # Disallow defining incompletely typed functions
        self.disallow_incomplete_defs = False

        # Type check unannotated functions
        self.check_untyped_defs = False

        # Disallow decorating typed functions with untyped decorators
        self.disallow_untyped_decorators = False

        # Disallow subclassing values of type 'Any'
        self.disallow_subclassing_any = False

        # Also check typeshed for missing annotations
        self.warn_incomplete_stub = False

        # Warn about casting an expression to its inferred type
        self.warn_redundant_casts = False

        # Warn about falling off the end of a function returning non-None
        self.warn_no_return = True

        # Warn about returning objects of type Any when the function is
        # declared with a precise type
        self.warn_return_any = False

        # Warn about unused '# type: ignore' comments
        self.warn_unused_ignores = False

        # Warn about unused '[mypy-<pattern>]'  or '[[tool.mypy.overrides]]' config sections
        self.warn_unused_configs = False

        # Files in which to ignore all non-fatal errors
        self.ignore_errors = False

        # Apply strict None checking
        self.strict_optional = True

        # Show "note: In function "foo":" messages.
        self.show_error_context = False

        # Use nicer output (when possible).
        self.color_output = True
        self.error_summary = True

        # Files in which to allow strict-Optional related errors
        # TODO: Kill this in favor of show_none_errors
        self.strict_optional_whitelist: Optional[List[str]] = None

        # Alternate way to show/hide strict-None-checking related errors
        self.show_none_errors = True

        # Don't assume arguments with default values of None are Optional
        self.no_implicit_optional = False

        # Don't re-export names unless they are imported with `from ... as ...`
        self.implicit_reexport = True

        # Suppress toplevel errors caused by missing annotations
        self.allow_untyped_globals = False

        # Allow variable to be redefined with an arbitrary type in the same block
        # and the same nesting level as the initialization
        self.allow_redefinition = False

        # Prohibit equality, identity, and container checks for non-overlapping types.
        # This makes 1 == '1', 1 in ['1'], and 1 is '1' errors.
        self.strict_equality = False

        # Report an error for any branches inferred to be unreachable as a result of
        # type analysis.
        self.warn_unreachable = False

        # Variable names considered True
        self.always_true: List[str] = []

        # Variable names considered False
        self.always_false: List[str] = []

        # Error codes to disable
        self.disable_error_code: List[str] = []
        self.disabled_error_codes: Set[ErrorCode] = set()

        # Error codes to enable
        self.enable_error_code: List[str] = []
        self.enabled_error_codes: Set[ErrorCode] = set()

        # Use script name instead of __main__
        self.scripts_are_modules = False

        # Config file name
        self.config_file: Optional[str] = None

        # A filename containing a JSON mapping from filenames to
        # mtime/size/hash arrays, used to avoid having to recalculate
        # source hashes as often.
        self.quickstart_file: Optional[str] = None

        # A comma-separated list of files/directories for mypy to type check;
        # supports globbing
        self.files: Optional[List[str]] = None

        # Write junit.xml to given file
        self.junit_xml: Optional[str] = None

        # Caching and incremental checking options
        self.incremental = True
        self.cache_dir = defaults.CACHE_DIR
        self.sqlite_cache = False
        self.debug_cache = False
        self.skip_version_check = False
        self.skip_cache_mtime_checks = False
        self.fine_grained_incremental = False
        # Include fine-grained dependencies in written cache files
        self.cache_fine_grained = False
        # Read cache files in fine-grained incremental mode (cache must include dependencies)
        self.use_fine_grained_cache = False

        # Tune certain behaviors when being used as a front-end to mypyc. Set per-module
        # in modules being compiled. Not in the config file or command line.
        self.mypyc = False

        # Disable the memory optimization of freeing ASTs when
        # possible. This isn't exposed as a command line option
        # because it is intended for software integrating with
        # mypy. (Like mypyc.)
        self.preserve_asts = False

        # Paths of user plugins
        self.plugins: List[str] = []

        # Per-module options (raw)
        self.per_module_options: OrderedDict[str, Dict[str, object]] = OrderedDict()
        self._glob_options: List[Tuple[str, Pattern[str]]] = []
        self.unused_configs: Set[str] = set()

        # -- development options --
        self.verbosity = 0  # More verbose messages (for troubleshooting)
        self.pdb = False
        self.show_traceback = False
        self.raise_exceptions = False
        self.dump_type_stats = False
        self.dump_inference_stats = False
        self.dump_build_stats = False

        # -- test options --
        # Stop after the semantic analysis phase
        self.semantic_analysis_only = False

        # Use stub builtins fixtures to speed up tests
        self.use_builtins_fixtures = False

        # -- experimental options --
        self.shadow_file: Optional[List[List[str]]] = None
        self.show_column_numbers: bool = False
        self.show_error_codes = False
        # Use soft word wrap and show trimmed source snippets with error location markers.
        self.pretty = False
        self.dump_graph = False
        self.dump_deps = False
        self.logical_deps = False
        # If True, partial types can't span a module top level and a function
        self.local_partial_types = False
        # Some behaviors are changed when using Bazel (https://bazel.build).
        self.bazel = False
        # If True, export inferred types for all expressions as BuildResult.types
        self.export_types = False
        # List of package roots -- directories under these are packages even
        # if they don't have __init__.py.
        self.package_root: List[str] = []
        self.cache_map: Dict[str, Tuple[str, str]] = {}
        # Don't properly free objects on exit, just kill the current process.
        self.fast_exit = True
        # Used to transform source code before parsing if not None
        # TODO: Make the type precise (AnyStr -> AnyStr)
        self.transform_source: Optional[Callable[[Any], Any]] = None
        # Print full path to each file in the report.
        self.show_absolute_path: bool = False
        # Install missing stub packages if True
        self.install_types = False
        # Install missing stub packages in non-interactive mode (don't prompt for
        # confirmation, and don't show any errors)
        self.non_interactive = False
        # When we encounter errors that may cause many additional errors,
        # skip most errors after this many messages have been reported.
        # -1 means unlimited.
        self.many_errors_threshold = defaults.MANY_ERRORS_THRESHOLD

    # To avoid breaking plugin compatibility, keep providing new_semantic_analyzer
    @property
    def new_semantic_analyzer(self) -> bool:
        return True

    def snapshot(self) -> object:
        """Produce a comparable snapshot of this Option"""
        # Under mypyc, we don't have a __dict__, so we need to do worse things.
        d = dict(getattr(self, '__dict__', ()))
        for k in get_class_descriptors(Options):
            if hasattr(self, k) and k != "new_semantic_analyzer":
                d[k] = getattr(self, k)
        # Remove private attributes from snapshot
        d = {k: v for k, v in d.items() if not k.startswith('_')}
        return d

    def __repr__(self) -> str:
        return 'Options({})'.format(pprint.pformat(self.snapshot()))

    def apply_changes(self, changes: Dict[str, object]) -> 'Options':
        new_options = Options()
        # Under mypyc, we don't have a __dict__, so we need to do worse things.
        replace_object_state(new_options, self, copy_dict=True)
        for key, value in changes.items():
            setattr(new_options, key, value)
        if changes.get("ignore_missing_imports"):
            # This is the only option for which a per-module and a global
            # option sometimes beheave differently.
            new_options.ignore_missing_imports_per_module = True
        return new_options

    def build_per_module_cache(self) -> None:
        self._per_module_cache = {}

        # Config precedence is as follows:
        #  1. Concrete section names: foo.bar.baz
        #  2. "Unstructured" glob patterns: foo.*.baz, in the order
        #     they appear in the file (last wins)
        #  3. "Well-structured" wildcard patterns: foo.bar.*, in specificity order.

        # Since structured configs inherit from structured configs above them in the hierarchy,
        # we need to process per-module configs in a careful order.
        # We have to process foo.* before foo.bar.* before foo.bar,
        # and we need to apply *.bar to foo.bar but not to foo.bar.*.
        # To do this, process all well-structured glob configs before non-glob configs and
        # exploit the fact that foo.* sorts earlier ASCIIbetically (unicodebetically?)
        # than foo.bar.*.
        # (A section being "processed last" results in its config "winning".)
        # Unstructured glob configs are stored and are all checked for each module.
        unstructured_glob_keys = [k for k in self.per_module_options.keys()
                                  if '*' in k[:-1]]
        structured_keys = [k for k in self.per_module_options.keys()
                           if '*' not in k[:-1]]
        wildcards = sorted(k for k in structured_keys if k.endswith('.*'))
        concrete = [k for k in structured_keys if not k.endswith('.*')]

        for glob in unstructured_glob_keys:
            self._glob_options.append((glob, self.compile_glob(glob)))

        # We (for ease of implementation) treat unstructured glob
        # sections as used if any real modules use them or if any
        # concrete config sections use them. This means we need to
        # track which get used while constructing.
        self.unused_configs = set(unstructured_glob_keys)

        for key in wildcards + concrete:
            # Find what the options for this key would be, just based
            # on inheriting from parent configs.
            options = self.clone_for_module(key)
            # And then update it with its per-module options.
            self._per_module_cache[key] = options.apply_changes(self.per_module_options[key])

        # Add the more structured sections into unused configs, since
        # they only count as used if actually used by a real module.
        self.unused_configs.update(structured_keys)

    def clone_for_module(self, module: str) -> 'Options':
        """Create an Options object that incorporates per-module options.

        NOTE: Once this method is called all Options objects should be
        considered read-only, else the caching might be incorrect.
        """
        if self._per_module_cache is None:
            self.build_per_module_cache()
        assert self._per_module_cache is not None

        # If the module just directly has a config entry, use it.
        if module in self._per_module_cache:
            self.unused_configs.discard(module)
            return self._per_module_cache[module]

        # If not, search for glob paths at all the parents. So if we are looking for
        # options for foo.bar.baz, we search foo.bar.baz.*, foo.bar.*, foo.*,
        # in that order, looking for an entry.
        # This is technically quadratic in the length of the path, but module paths
        # don't actually get all that long.
        options = self
        path = module.split('.')
        for i in range(len(path), 0, -1):
            key = '.'.join(path[:i] + ['*'])
            if key in self._per_module_cache:
                self.unused_configs.discard(key)
                options = self._per_module_cache[key]
                break

        # OK and *now* we need to look for unstructured glob matches.
        # We only do this for concrete modules, not structured wildcards.
        if not module.endswith('.*'):
            for key, pattern in self._glob_options:
                if pattern.match(module):
                    self.unused_configs.discard(key)
                    options = options.apply_changes(self.per_module_options[key])

        # We could update the cache to directly point to modules once
        # they have been looked up, but in testing this made things
        # slower and not faster, so we don't bother.

        return options

    def compile_glob(self, s: str) -> Pattern[str]:
        # Compile one of the glob patterns to a regex so that '.*' can
        # match *zero or more* module sections. This means we compile
        # '.*' into '(\..*)?'.
        parts = s.split('.')
        expr = re.escape(parts[0]) if parts[0] != '*' else '.*'
        for part in parts[1:]:
            expr += re.escape('.' + part) if part != '*' else r'(\..*)?'
        return re.compile(expr + '\\Z')

    def select_options_affecting_cache(self) -> Mapping[str, object]:
        return {opt: getattr(self, opt) for opt in OPTIONS_AFFECTING_CACHE}
