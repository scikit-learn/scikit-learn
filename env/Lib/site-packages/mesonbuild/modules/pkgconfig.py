# SPDX-License-Identifier: Apache-2.0
# Copyright 2015-2022 The Meson development team

from __future__ import annotations
from collections import defaultdict
from dataclasses import dataclass
from pathlib import PurePath
import os
import typing as T

from . import NewExtensionModule, ModuleInfo
from . import ModuleReturnValue
from .. import build
from .. import dependencies
from .. import mesonlib
from ..options import OptionKey
from .. import mlog
from ..options import BUILTIN_DIR_OPTIONS
from ..dependencies.pkgconfig import PkgConfigDependency, PkgConfigInterface
from ..interpreter.type_checking import D_MODULE_VERSIONS_KW, INSTALL_DIR_KW, VARIABLES_KW, NoneType
from ..interpreterbase import FeatureNew, FeatureDeprecated
from ..interpreterbase.decorators import ContainerTypeInfo, KwargInfo, typed_kwargs, typed_pos_args

if T.TYPE_CHECKING:
    from typing_extensions import TypedDict

    from . import ModuleState
    from .. import mparser
    from ..interpreter import Interpreter

    ANY_DEP = T.Union[dependencies.Dependency, build.BuildTargetTypes, str]
    LIBS = T.Union[build.LibTypes, str]

    class GenerateKw(TypedDict):

        version: T.Optional[str]
        name: T.Optional[str]
        filebase: T.Optional[str]
        description: T.Optional[str]
        url: str
        subdirs: T.List[str]
        conflicts: T.List[str]
        dataonly: bool
        libraries: T.List[ANY_DEP]
        libraries_private: T.List[ANY_DEP]
        requires: T.List[T.Union[str, build.StaticLibrary, build.SharedLibrary, dependencies.Dependency]]
        requires_private: T.List[T.Union[str, build.StaticLibrary, build.SharedLibrary, dependencies.Dependency]]
        install_dir: T.Optional[str]
        d_module_versions: T.List[T.Union[str, int]]
        extra_cflags: T.List[str]
        variables: T.Dict[str, str]
        uninstalled_variables: T.Dict[str, str]
        unescaped_variables: T.Dict[str, str]
        unescaped_uninstalled_variables: T.Dict[str, str]


_PKG_LIBRARIES: KwargInfo[T.List[T.Union[str, dependencies.Dependency, build.SharedLibrary, build.StaticLibrary, build.CustomTarget, build.CustomTargetIndex]]] = KwargInfo(
    'libraries',
    ContainerTypeInfo(list, (str, dependencies.Dependency,
                             build.SharedLibrary, build.StaticLibrary,
                             build.CustomTarget, build.CustomTargetIndex)),
    default=[],
    listify=True,
)

_PKG_REQUIRES: KwargInfo[T.List[T.Union[str, build.SharedLibrary, build.StaticLibrary, dependencies.Dependency]]] = KwargInfo(
    'requires',
    ContainerTypeInfo(list, (str, build.SharedLibrary, build.StaticLibrary, dependencies.Dependency)),
    default=[],
    listify=True,
)


def _as_str(obj: object) -> str:
    assert isinstance(obj, str)
    return obj


@dataclass
class MetaData:

    filebase: str
    display_name: str
    location: mparser.BaseNode
    warned: bool = False


class DependenciesHelper:
    def __init__(self, state: ModuleState, name: str, metadata: T.Dict[str, MetaData]) -> None:
        self.state = state
        self.name = name
        self.metadata = metadata
        self.pub_libs: T.List[LIBS] = []
        self.pub_reqs: T.List[str] = []
        self.priv_libs: T.List[LIBS] = []
        self.priv_reqs: T.List[str] = []
        self.cflags: T.List[str] = []
        self.version_reqs: T.DefaultDict[str, T.Set[str]] = defaultdict(set)
        self.link_whole_targets: T.List[T.Union[build.CustomTarget, build.CustomTargetIndex, build.StaticLibrary]] = []
        self.uninstalled_incdirs: mesonlib.OrderedSet[str] = mesonlib.OrderedSet()

    def add_pub_libs(self, libs: T.List[ANY_DEP]) -> None:
        p_libs, reqs, cflags = self._process_libs(libs, True)
        self.pub_libs = p_libs + self.pub_libs # prepend to preserve dependencies
        self.pub_reqs += reqs
        self.cflags += cflags

    def add_priv_libs(self, libs: T.List[ANY_DEP]) -> None:
        p_libs, reqs, _ = self._process_libs(libs, False)
        self.priv_libs = p_libs + self.priv_libs
        self.priv_reqs += reqs

    def add_pub_reqs(self, reqs: T.List[T.Union[str, build.StaticLibrary, build.SharedLibrary, dependencies.Dependency]]) -> None:
        self.pub_reqs += self._process_reqs(reqs)

    def add_priv_reqs(self, reqs: T.List[T.Union[str, build.StaticLibrary, build.SharedLibrary, dependencies.Dependency]]) -> None:
        self.priv_reqs += self._process_reqs(reqs)

    def _check_generated_pc_deprecation(self, obj: T.Union[build.CustomTarget, build.CustomTargetIndex, build.StaticLibrary, build.SharedLibrary]) -> None:
        if obj.get_id() in self.metadata:
            return
        data = self.metadata[obj.get_id()]
        if data.warned:
            return
        mlog.deprecation('Library', mlog.bold(obj.name), 'was passed to the '
                         '"libraries" keyword argument of a previous call '
                         'to generate() method instead of first positional '
                         'argument.', 'Adding', mlog.bold(data.display_name),
                         'to "Requires" field, but this is a deprecated '
                         'behaviour that will change in version 2.0 '
                         'of Meson. Please report the issue if this '
                         'warning cannot be avoided in your case.',
                         location=data.location)
        data.warned = True

    def _process_reqs(self, reqs: T.Sequence[T.Union[str, build.StaticLibrary, build.SharedLibrary, dependencies.Dependency]]) -> T.List[str]:
        '''Returns string names of requirements'''
        processed_reqs: T.List[str] = []
        for obj in mesonlib.listify(reqs):
            if not isinstance(obj, str):
                FeatureNew.single_use('pkgconfig.generate requirement from non-string object', '0.46.0', self.state.subproject)
            if (isinstance(obj, (build.CustomTarget, build.CustomTargetIndex, build.SharedLibrary, build.StaticLibrary))
                    and obj.get_id() in self.metadata):
                self._check_generated_pc_deprecation(obj)
                processed_reqs.append(self.metadata[obj.get_id()].filebase)
            elif isinstance(obj, PkgConfigDependency):
                if obj.found():
                    processed_reqs.append(obj.name)
                    self.add_version_reqs(obj.name, obj.version_reqs)
            elif isinstance(obj, str):
                name, version_req = self.split_version_req(obj)
                processed_reqs.append(name)
                self.add_version_reqs(name, [version_req] if version_req is not None else None)
            elif isinstance(obj, dependencies.Dependency) and not obj.found():
                pass
            elif isinstance(obj, dependencies.ExternalDependency) and obj.name == 'threads':
                pass
            else:
                raise mesonlib.MesonException('requires argument not a string, '
                                              'library with pkgconfig-generated file '
                                              f'or pkgconfig-dependency object, got {obj!r}')
        return processed_reqs

    def add_cflags(self, cflags: T.List[str]) -> None:
        self.cflags += mesonlib.stringlistify(cflags)

    def _add_uninstalled_incdirs(self, incdirs: T.List[build.IncludeDirs], subdir: T.Optional[str] = None) -> None:
        for i in incdirs:
            curdir = i.get_curdir()
            for d in i.get_incdirs():
                path = os.path.join(curdir, d)
                self.uninstalled_incdirs.add(path)
        if subdir is not None:
            self.uninstalled_incdirs.add(subdir)

    def _process_libs(
            self, libs: T.List[ANY_DEP], public: bool
            ) -> T.Tuple[T.List[T.Union[str, build.SharedLibrary, build.StaticLibrary, build.CustomTarget, build.CustomTargetIndex]], T.List[str], T.List[str]]:
        libs = mesonlib.listify(libs)
        processed_libs: T.List[T.Union[str, build.SharedLibrary, build.StaticLibrary, build.CustomTarget, build.CustomTargetIndex]] = []
        processed_reqs: T.List[str] = []
        processed_cflags: T.List[str] = []
        for obj in libs:
            if (isinstance(obj, (build.CustomTarget, build.CustomTargetIndex, build.SharedLibrary, build.StaticLibrary))
                    and obj.get_id() in self.metadata):
                self._check_generated_pc_deprecation(obj)
                processed_reqs.append(self.metadata[obj.get_id()].filebase)
            elif isinstance(obj, dependencies.ExternalDependency) and obj.name == 'valgrind':
                pass
            elif isinstance(obj, PkgConfigDependency):
                if obj.found():
                    processed_reqs.append(obj.name)
                    self.add_version_reqs(obj.name, obj.version_reqs)
            elif isinstance(obj, dependencies.InternalDependency):
                if obj.found():
                    if obj.objects:
                        raise mesonlib.MesonException('.pc file cannot refer to individual object files.')

                    # Ensure BothLibraries are resolved:
                    if self.pub_libs and isinstance(self.pub_libs[0], build.StaticLibrary):
                        obj = obj.get_as_static(recursive=True)
                    else:
                        obj = obj.get_as_shared(recursive=True)

                    processed_libs += obj.get_link_args()
                    processed_cflags += obj.get_compile_args()
                    self._add_lib_dependencies(obj.libraries, obj.whole_libraries, obj.ext_deps, public, private_external_deps=True)
                    self._add_uninstalled_incdirs(obj.get_include_dirs())
            elif isinstance(obj, dependencies.Dependency):
                if obj.found():
                    processed_libs += obj.get_link_args()
                    processed_cflags += obj.get_compile_args()
            elif isinstance(obj, build.SharedLibrary) and obj.shared_library_only:
                # Do not pull dependencies for shared libraries because they are
                # only required for static linking. Adding private requires has
                # the side effect of exposing their cflags, which is the
                # intended behaviour of pkg-config but force Debian to add more
                # than needed build deps.
                # See https://bugs.freedesktop.org/show_bug.cgi?id=105572
                processed_libs.append(obj)
                self._add_uninstalled_incdirs(obj.get_include_dirs(), obj.get_subdir())
            elif isinstance(obj, (build.SharedLibrary, build.StaticLibrary)):
                processed_libs.append(obj)
                self._add_uninstalled_incdirs(obj.get_include_dirs(), obj.get_subdir())
                # If there is a static library in `Libs:` all its deps must be
                # public too, otherwise the generated pc file will never be
                # usable without --static.
                self._add_lib_dependencies(obj.link_targets,
                                           obj.link_whole_targets,
                                           obj.external_deps,
                                           isinstance(obj, build.StaticLibrary) and public)
            elif isinstance(obj, (build.CustomTarget, build.CustomTargetIndex)):
                if not obj.is_linkable_target():
                    raise mesonlib.MesonException('library argument contains a not linkable custom_target.')
                FeatureNew.single_use('custom_target in pkgconfig.generate libraries', '0.58.0', self.state.subproject)
                processed_libs.append(obj)
            elif isinstance(obj, str):
                processed_libs.append(obj)
            else:
                raise mesonlib.MesonException(f'library argument of type {type(obj).__name__} not a string, library or dependency object.')

        return processed_libs, processed_reqs, processed_cflags

    def _add_lib_dependencies(
            self, link_targets: T.Sequence[build.BuildTargetTypes],
            link_whole_targets: T.Sequence[T.Union[build.StaticLibrary, build.CustomTarget, build.CustomTargetIndex]],
            external_deps: T.List[dependencies.Dependency],
            public: bool,
            private_external_deps: bool = False) -> None:
        add_libs = self.add_pub_libs if public else self.add_priv_libs
        # Recursively add all linked libraries
        for t in link_targets:
            # Internal libraries (uninstalled static library) will be promoted
            # to link_whole, treat them as such here.
            if t.is_internal():
                # `is_internal` shouldn't return True for anything but a
                # StaticLibrary, or a CustomTarget that is a StaticLibrary
                assert isinstance(t, (build.StaticLibrary, build.CustomTarget, build.CustomTargetIndex)), 'for mypy'
                self._add_link_whole(t, public)
            else:
                add_libs([t])
        for t in link_whole_targets:
            self._add_link_whole(t, public)
        # And finally its external dependencies
        if private_external_deps:
            self.add_priv_libs(T.cast('T.List[ANY_DEP]', external_deps))
        else:
            add_libs(T.cast('T.List[ANY_DEP]', external_deps))

    def _add_link_whole(self, t: T.Union[build.CustomTarget, build.CustomTargetIndex, build.StaticLibrary], public: bool) -> None:
        # Don't include static libraries that we link_whole. But we still need to
        # include their dependencies: a static library we link_whole
        # could itself link to a shared library or an installed static library.
        # Keep track of link_whole_targets so we can remove them from our
        # lists in case a library is link_with and link_whole at the same time.
        # See remove_dups() below.
        self.link_whole_targets.append(t)
        if isinstance(t, build.BuildTarget):
            self._add_lib_dependencies(t.link_targets, t.link_whole_targets, t.external_deps, public)

    def add_version_reqs(self, name: str, version_reqs: T.Optional[T.List[str]]) -> None:
        if version_reqs:
            # Note that pkg-config is picky about whitespace.
            # 'foo > 1.2' is ok but 'foo>1.2' is not.
            # foo, bar' is ok, but 'foo,bar' is not.
            self.version_reqs[name].update(version_reqs)

    def split_version_req(self, s: str) -> T.Tuple[str, T.Optional[str]]:
        for op in ['>=', '<=', '!=', '==', '=', '>', '<']:
            pos = s.find(op)
            if pos > 0:
                return s[0:pos].strip(), s[pos:].strip()
        return s, None

    def format_vreq(self, vreq: str) -> str:
        # vreq are '>=1.0' and pkgconfig wants '>= 1.0'
        for op in ['>=', '<=', '!=', '==', '=', '>', '<']:
            if vreq.startswith(op):
                return op + ' ' + vreq[len(op):]
        return vreq

    def format_reqs(self, reqs: T.List[str]) -> str:
        result: T.List[str] = []
        for name in reqs:
            vreqs = self.version_reqs.get(name, None)
            if vreqs:
                result += [name + ' ' + self.format_vreq(vreq) for vreq in sorted(vreqs)]
            else:
                result += [name]
        return ', '.join(result)

    def remove_dups(self) -> None:
        # Set of ids that have already been handled and should not be added any more
        exclude: T.Set[str] = set()

        # We can't just check if 'x' is excluded because we could have copies of
        # the same SharedLibrary object for example.
        def _ids(x: T.Union[str, build.CustomTarget, build.CustomTargetIndex, build.StaticLibrary, build.SharedLibrary]) -> T.Iterable[str]:
            if isinstance(x, str):
                yield x
            else:
                if x.get_id() in self.metadata:
                    yield self.metadata[x.get_id()].display_name
                yield x.get_id()

        # Exclude 'x' in all its forms and return if it was already excluded
        def _add_exclude(x: T.Union[str, build.CustomTarget, build.CustomTargetIndex, build.StaticLibrary, build.SharedLibrary]) -> bool:
            was_excluded = False
            for i in _ids(x):
                if i in exclude:
                    was_excluded = True
                else:
                    exclude.add(i)
            return was_excluded

        # link_whole targets are already part of other targets, exclude them all.
        for t in self.link_whole_targets:
            _add_exclude(t)

        # Mypy thinks these overlap, but since List is invariant they don't,
        # `List[str]`` is not a valid input to `List[str | BuildTarget]`.
        # pylance/pyright gets this right, but for mypy we have to ignore the
        # error
        @T.overload
        def _fn(xs: T.List[str], libs: bool = False) -> T.List[str]: ...  # type: ignore

        @T.overload
        def _fn(xs: T.List[LIBS], libs: bool = False) -> T.List[LIBS]: ...

        def _fn(xs: T.Union[T.List[str], T.List[LIBS]], libs: bool = False) -> T.Union[T.List[str], T.List[LIBS]]:
            # Remove duplicates whilst preserving original order
            result = []
            for x in xs:
                # Don't de-dup unknown strings to avoid messing up arguments like:
                # ['-framework', 'CoreAudio', '-framework', 'CoreMedia']
                known_flags = ['-pthread']
                cannot_dedup = libs and isinstance(x, str) and \
                    not x.startswith(('-l', '-L')) and \
                    x not in known_flags
                if not cannot_dedup and _add_exclude(x):
                    continue
                result.append(x)
            return result

        # Handle lists in priority order: public items can be excluded from
        # private and Requires can excluded from Libs.
        self.pub_reqs = _fn(self.pub_reqs)
        self.pub_libs = _fn(self.pub_libs, True)
        self.priv_reqs = _fn(self.priv_reqs)
        self.priv_libs = _fn(self.priv_libs, True)
        # Reset exclude list just in case some values can be both cflags and libs.
        exclude = set()
        self.cflags = _fn(self.cflags)

class PkgConfigModule(NewExtensionModule):

    INFO = ModuleInfo('pkgconfig')

    # Track already generated pkg-config files This is stored as a class
    # variable so that multiple `import()`s share metadata
    devenv: T.Optional[mesonlib.EnvironmentVariables] = None
    _metadata: T.ClassVar[T.Dict[str, MetaData]] = {}

    def __init__(self) -> None:
        super().__init__()
        self.methods.update({
            'generate': self.generate,
        })

    def postconf_hook(self, b: build.Build) -> None:
        if self.devenv is not None:
            b.devenv.append(self.devenv)

    def _get_lname(self, l: T.Union[build.SharedLibrary, build.StaticLibrary, build.CustomTarget, build.CustomTargetIndex],
                   msg: str, pcfile: str) -> str:
        if isinstance(l, (build.CustomTargetIndex, build.CustomTarget)):
            basename = os.path.basename(l.get_filename())
            name = os.path.splitext(basename)[0]
            if name.startswith('lib'):
                name = name[3:]
            return name
        # Nothing special
        if not l.name_prefix_set:
            return l.name
        # Sometimes people want the library to start with 'lib' everywhere,
        # which is achieved by setting name_prefix to '' and the target name to
        # 'libfoo'. In that case, try to get the pkg-config '-lfoo' arg correct.
        if l.prefix == '' and l.name.startswith('lib'):
            return l.name[3:]
        # If the library is imported via an import library which is always
        # named after the target name, '-lfoo' is correct.
        if isinstance(l, build.SharedLibrary) and l.import_filename:
            return l.name
        # In other cases, we can't guarantee that the compiler will be able to
        # find the library via '-lfoo', so tell the user that.
        mlog.warning(msg.format(l.name, 'name_prefix', l.name, pcfile))
        return l.name

    def _escape(self, value: T.Union[str, PurePath]) -> str:
        '''
        We cannot use quote_arg because it quotes with ' and " which does not
        work with pkg-config and pkgconf at all.
        '''
        # We should always write out paths with / because pkg-config requires
        # spaces to be quoted with \ and that messes up on Windows:
        # https://bugs.freedesktop.org/show_bug.cgi?id=103203
        if isinstance(value, PurePath):
            value = value.as_posix()
        return value.replace(' ', r'\ ')

    def _make_relative(self, prefix: T.Union[PurePath, str], subdir: T.Union[PurePath, str]) -> str:
        prefix = PurePath(prefix)
        subdir = PurePath(subdir)
        try:
            libdir = subdir.relative_to(prefix)
        except ValueError:
            libdir = subdir
        # pathlib joining makes sure absolute libdir is not appended to '${prefix}'
        return ('${prefix}' / libdir).as_posix()

    def _generate_pkgconfig_file(self, state: ModuleState, deps: DependenciesHelper,
                                 subdirs: T.List[str], name: str,
                                 description: str, url: str, version: str,
                                 pcfile: str, conflicts: T.List[str],
                                 variables: T.List[T.Tuple[str, str]],
                                 unescaped_variables: T.List[T.Tuple[str, str]],
                                 uninstalled: bool = False, dataonly: bool = False,
                                 pkgroot: T.Optional[str] = None) -> None:
        coredata = state.environment.get_coredata()
        referenced_vars = set()
        optnames = [x.name for x in BUILTIN_DIR_OPTIONS.keys()]

        if not dataonly:
            # includedir is always implied, although libdir may not be
            # needed for header-only libraries
            referenced_vars |= {'prefix', 'includedir'}
            if deps.pub_libs or deps.priv_libs:
                referenced_vars |= {'libdir'}
        # also automatically infer variables referenced in other variables
        implicit_vars_warning = False
        redundant_vars_warning = False
        varnames = set()
        varstrings = set()
        for k, v in variables + unescaped_variables:
            varnames |= {k}
            varstrings |= {v}
        for optname in optnames:
            optvar = f'${{{optname}}}'
            if any(x.startswith(optvar) for x in varstrings):
                if optname in varnames:
                    redundant_vars_warning = True
                else:
                    # these 3 vars were always "implicit"
                    if dataonly or optname not in {'prefix', 'includedir', 'libdir'}:
                        implicit_vars_warning = True
                    referenced_vars |= {'prefix', optname}
        if redundant_vars_warning:
            FeatureDeprecated.single_use('pkgconfig.generate variable for builtin directories', '0.62.0',
                                         state.subproject, 'They will be automatically included when referenced',
                                         state.current_node)
        if implicit_vars_warning:
            FeatureNew.single_use('pkgconfig.generate implicit variable for builtin directories', '0.62.0',
                                  state.subproject, location=state.current_node)

        if uninstalled:
            outdir = os.path.join(state.environment.build_dir, 'meson-uninstalled')
            if not os.path.exists(outdir):
                os.mkdir(outdir)
            prefix = PurePath(state.environment.get_build_dir())
            srcdir = PurePath(state.environment.get_source_dir())
        else:
            outdir = state.environment.scratch_dir
            prefix = PurePath(_as_str(coredata.get_option(OptionKey('prefix'))))
            if pkgroot:
                pkgroot_ = PurePath(pkgroot)
                if not pkgroot_.is_absolute():
                    pkgroot_ = prefix / pkgroot
                elif prefix not in pkgroot_.parents:
                    raise mesonlib.MesonException('Pkgconfig prefix cannot be outside of the prefix '
                                                  'when pkgconfig.relocatable=true. '
                                                  f'Pkgconfig prefix is {pkgroot_.as_posix()}.')
                prefix = PurePath('${pcfiledir}', os.path.relpath(prefix, pkgroot_))
        fname = os.path.join(outdir, pcfile)
        with open(fname, 'w', encoding='utf-8') as ofile:
            for optname in optnames:
                if optname in referenced_vars - varnames:
                    if optname == 'prefix':
                        ofile.write('prefix={}\n'.format(self._escape(prefix)))
                    else:
                        dirpath = PurePath(_as_str(coredata.get_option(OptionKey(optname))))
                        ofile.write('{}={}\n'.format(optname, self._escape('${prefix}' / dirpath)))
            if uninstalled and not dataonly:
                ofile.write('srcdir={}\n'.format(self._escape(srcdir)))
            if variables or unescaped_variables:
                ofile.write('\n')
            for k, v in variables:
                ofile.write('{}={}\n'.format(k, self._escape(v)))
            for k, v in unescaped_variables:
                ofile.write(f'{k}={v}\n')
            ofile.write('\n')
            ofile.write(f'Name: {name}\n')
            if len(description) > 0:
                ofile.write(f'Description: {description}\n')
            if len(url) > 0:
                ofile.write(f'URL: {url}\n')
            ofile.write(f'Version: {version}\n')
            reqs_str = deps.format_reqs(deps.pub_reqs)
            if len(reqs_str) > 0:
                ofile.write(f'Requires: {reqs_str}\n')
            reqs_str = deps.format_reqs(deps.priv_reqs)
            if len(reqs_str) > 0:
                ofile.write(f'Requires.private: {reqs_str}\n')
            if len(conflicts) > 0:
                ofile.write('Conflicts: {}\n'.format(' '.join(conflicts)))

            def generate_libs_flags(libs: T.List[LIBS]) -> T.Iterable[str]:
                msg = 'Library target {0!r} has {1!r} set. Compilers ' \
                      'may not find it from its \'-l{2}\' linker flag in the ' \
                      '{3!r} pkg-config file.'
                Lflags = []
                for l in libs:
                    if isinstance(l, str):
                        yield l
                    else:
                        install_dir: T.Union[str, bool]
                        if uninstalled:
                            install_dir = os.path.dirname(state.backend.get_target_filename_abs(l))
                        else:
                            _i = l.get_custom_install_dir()
                            install_dir = _i[0] if _i else None
                        if install_dir is False:
                            continue
                        if isinstance(l, build.BuildTarget) and 'cs' in l.compilers:
                            if isinstance(install_dir, str):
                                Lflag = '-r{}/{}'.format(self._escape(self._make_relative(prefix, install_dir)), l.filename)
                            else:  # install_dir is True
                                Lflag = '-r${libdir}/%s' % l.filename
                        else:
                            if isinstance(install_dir, str):
                                Lflag = '-L{}'.format(self._escape(self._make_relative(prefix, install_dir)))
                            else:  # install_dir is True
                                Lflag = '-L${libdir}'
                        if Lflag not in Lflags:
                            Lflags.append(Lflag)
                            yield Lflag
                        lname = self._get_lname(l, msg, pcfile)
                        # If using a custom suffix, the compiler may not be able to
                        # find the library
                        if isinstance(l, build.BuildTarget) and l.name_suffix_set:
                            mlog.warning(msg.format(l.name, 'name_suffix', lname, pcfile))
                        if isinstance(l, (build.CustomTarget, build.CustomTargetIndex)) or 'cs' not in l.compilers:
                            yield f'-l{lname}'

            if len(deps.pub_libs) > 0:
                ofile.write('Libs: {}\n'.format(' '.join(generate_libs_flags(deps.pub_libs))))
            if len(deps.priv_libs) > 0:
                ofile.write('Libs.private: {}\n'.format(' '.join(generate_libs_flags(deps.priv_libs))))

            cflags: T.List[str] = []
            if uninstalled:
                for d in deps.uninstalled_incdirs:
                    for basedir in ['${prefix}', '${srcdir}']:
                        path = self._escape(PurePath(basedir, d).as_posix())
                        cflags.append(f'-I{path}')
            else:
                for d in subdirs:
                    if d == '.':
                        cflags.append('-I${includedir}')
                    else:
                        cflags.append(self._escape(PurePath('-I${includedir}') / d))
            cflags += [self._escape(f) for f in deps.cflags]
            if cflags and not dataonly:
                ofile.write('Cflags: {}\n'.format(' '.join(cflags)))

    @typed_pos_args('pkgconfig.generate', optargs=[(build.SharedLibrary, build.StaticLibrary)])
    @typed_kwargs(
        'pkgconfig.generate',
        D_MODULE_VERSIONS_KW.evolve(since='0.43.0'),
        INSTALL_DIR_KW,
        KwargInfo('conflicts', ContainerTypeInfo(list, str), default=[], listify=True),
        KwargInfo('dataonly', bool, default=False, since='0.54.0'),
        KwargInfo('description', (str, NoneType)),
        KwargInfo('extra_cflags', ContainerTypeInfo(list, str), default=[], listify=True, since='0.42.0'),
        KwargInfo('filebase', (str, NoneType), validator=lambda x: 'must not be an empty string' if x == '' else None),
        KwargInfo('name', (str, NoneType), validator=lambda x: 'must not be an empty string' if x == '' else None),
        KwargInfo('subdirs', ContainerTypeInfo(list, str), default=[], listify=True),
        KwargInfo('url', str, default=''),
        KwargInfo('version', (str, NoneType)),
        VARIABLES_KW.evolve(name="unescaped_uninstalled_variables", since='0.59.0'),
        VARIABLES_KW.evolve(name="unescaped_variables", since='0.59.0'),
        VARIABLES_KW.evolve(name="uninstalled_variables", since='0.54.0', since_values={dict: '0.56.0'}),
        VARIABLES_KW.evolve(since='0.41.0', since_values={dict: '0.56.0'}),
        _PKG_LIBRARIES,
        _PKG_LIBRARIES.evolve(name='libraries_private'),
        _PKG_REQUIRES,
        _PKG_REQUIRES.evolve(name='requires_private'),
    )
    def generate(self, state: ModuleState,
                 args: T.Tuple[T.Optional[T.Union[build.SharedLibrary, build.StaticLibrary]]],
                 kwargs: GenerateKw) -> ModuleReturnValue:
        default_version = state.project_version
        default_install_dir: T.Optional[str] = None
        default_description: T.Optional[str] = None
        default_name: T.Optional[str] = None
        mainlib: T.Optional[T.Union[build.SharedLibrary, build.StaticLibrary]] = None
        default_subdirs = ['.']
        if args[0]:
            FeatureNew.single_use('pkgconfig.generate optional positional argument', '0.46.0', state.subproject)
            mainlib = args[0]
            default_name = mainlib.name
            default_description = state.project_name + ': ' + mainlib.name
            install_dir = mainlib.get_custom_install_dir()
            if install_dir and isinstance(install_dir[0], str):
                default_install_dir = os.path.join(install_dir[0], 'pkgconfig')
        else:
            if kwargs['version'] is None:
                FeatureNew.single_use('pkgconfig.generate implicit version keyword', '0.46.0', state.subproject)
            msg = ('pkgconfig.generate: if a library is not passed as a '
                   'positional argument, the {!r} keyword argument is '
                   'required.')
            if kwargs['name'] is None:
                raise build.InvalidArguments(msg.format('name'))
            if kwargs['description'] is None:
                raise build.InvalidArguments(msg.format('description'))

        dataonly = kwargs['dataonly']
        if dataonly:
            default_subdirs = []
            blocked_vars = ['libraries', 'libraries_private', 'requires_private', 'extra_cflags', 'subdirs']
            # Mypy can't figure out that this TypedDict index is correct, without repeating T.Literal for the entire list
            if any(kwargs[k] for k in blocked_vars):  # type: ignore
                raise mesonlib.MesonException(f'Cannot combine dataonly with any of {blocked_vars}')
            default_install_dir = os.path.join(state.environment.get_datadir(), 'pkgconfig')

        subdirs = kwargs['subdirs'] or default_subdirs
        version = kwargs['version'] if kwargs['version'] is not None else default_version
        name = kwargs['name'] if kwargs['name'] is not None else default_name
        assert isinstance(name, str), 'for mypy'
        filebase = kwargs['filebase'] if kwargs['filebase'] is not None else name
        description = kwargs['description'] if kwargs['description'] is not None else default_description
        url = kwargs['url']
        conflicts = kwargs['conflicts']

        # Prepend the main library to public libraries list. This is required
        # so dep.add_pub_libs() can handle dependency ordering correctly and put
        # extra libraries after the main library.
        libraries = kwargs['libraries'].copy()
        if mainlib:
            libraries.insert(0, mainlib)

        deps = DependenciesHelper(state, filebase, self._metadata)
        deps.add_pub_libs(libraries)
        deps.add_priv_libs(kwargs['libraries_private'])
        deps.add_pub_reqs(kwargs['requires'])
        deps.add_priv_reqs(kwargs['requires_private'])
        deps.add_cflags(kwargs['extra_cflags'])

        dversions = kwargs['d_module_versions']
        if dversions:
            compiler = state.environment.coredata.compilers.host.get('d')
            if compiler:
                deps.add_cflags(compiler.get_feature_args(
                    {'versions': dversions, 'import_dirs': [], 'debug': [], 'unittest': False}, None))

        deps.remove_dups()

        def parse_variable_list(vardict: T.Dict[str, str]) -> T.List[T.Tuple[str, str]]:
            reserved = ['prefix', 'libdir', 'includedir']
            variables = []
            for name, value in vardict.items():
                if not value:
                    FeatureNew.single_use('empty variable value in pkg.generate', '1.4.0', state.subproject, location=state.current_node)
                if not dataonly and name in reserved:
                    raise mesonlib.MesonException(f'Variable "{name}" is reserved')
                variables.append((name, value))
            return variables

        variables = parse_variable_list(kwargs['variables'])
        unescaped_variables = parse_variable_list(kwargs['unescaped_variables'])

        pcfile = filebase + '.pc'
        pkgroot = pkgroot_name = kwargs['install_dir'] or default_install_dir
        if pkgroot is None:
            if mesonlib.is_freebsd():
                pkgroot = os.path.join(_as_str(state.environment.coredata.get_option(OptionKey('prefix'))), 'libdata', 'pkgconfig')
                pkgroot_name = os.path.join('{prefix}', 'libdata', 'pkgconfig')
            elif mesonlib.is_haiku():
                pkgroot = os.path.join(_as_str(state.environment.coredata.get_option(OptionKey('prefix'))), 'develop', 'lib', 'pkgconfig')
                pkgroot_name = os.path.join('{prefix}', 'develop', 'lib', 'pkgconfig')
            else:
                pkgroot = os.path.join(_as_str(state.environment.coredata.get_option(OptionKey('libdir'))), 'pkgconfig')
                pkgroot_name = os.path.join('{libdir}', 'pkgconfig')
        relocatable = state.get_option('pkgconfig.relocatable')
        self._generate_pkgconfig_file(state, deps, subdirs, name, description, url,
                                      version, pcfile, conflicts, variables,
                                      unescaped_variables, False, dataonly,
                                      pkgroot=pkgroot if relocatable else None)
        res = build.Data([mesonlib.File(True, state.environment.get_scratch_dir(), pcfile)], pkgroot, pkgroot_name, None, state.subproject, install_tag='devel')
        variables = parse_variable_list(kwargs['uninstalled_variables'])
        unescaped_variables = parse_variable_list(kwargs['unescaped_uninstalled_variables'])

        pcfile = filebase + '-uninstalled.pc'
        self._generate_pkgconfig_file(state, deps, subdirs, name, description, url,
                                      version, pcfile, conflicts, variables,
                                      unescaped_variables, uninstalled=True, dataonly=dataonly)
        # Associate the main library with this generated pc file. If the library
        # is used in any subsequent call to the generated, it will generate a
        # 'Requires:' or 'Requires.private:'.
        # Backward compatibility: We used to set 'generated_pc' on all public
        # libraries instead of just the main one. Keep doing that but warn if
        # anyone is relying on that deprecated behaviour.
        if mainlib:
            if mainlib.get_id() not in self._metadata:
                self._metadata[mainlib.get_id()] = MetaData(
                    filebase, name, state.current_node)
            else:
                mlog.warning('Already generated a pkg-config file for', mlog.bold(mainlib.name))
        else:
            for lib in deps.pub_libs:
                if not isinstance(lib, str) and lib.get_id() not in self._metadata:
                    self._metadata[lib.get_id()] = MetaData(
                        filebase, name, state.current_node)
        if self.devenv is None:
            self.devenv = PkgConfigInterface.get_env(state.environment, mesonlib.MachineChoice.HOST, uninstalled=True)
        return ModuleReturnValue(res, [res])


def initialize(interp: Interpreter) -> PkgConfigModule:
    return PkgConfigModule()
