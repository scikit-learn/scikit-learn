# SPDX-License-Identifier: Apache-2.0
# Copyright 2014-2021 The Meson development team

from __future__ import annotations

import functools, uuid, os, operator, re
import typing as T

from . import backends
from .. import build
from .. import mesonlib
from .. import mlog
from ..arglist import CompilerArgs
from ..mesonlib import MesonBugException, MesonException
from ..options import OptionKey

if T.TYPE_CHECKING:
    from ..build import BuildTarget
    from ..compilers import Compiler
    from ..interpreter import Interpreter

INDENT = '\t'
XCODETYPEMAP = {'c': 'sourcecode.c.c',
                'a': 'archive.ar',
                'cc': 'sourcecode.cpp.cpp',
                'cxx': 'sourcecode.cpp.cpp',
                'cpp': 'sourcecode.cpp.cpp',
                'c++': 'sourcecode.cpp.cpp',
                'm': 'sourcecode.c.objc',
                'mm': 'sourcecode.cpp.objcpp',
                'h': 'sourcecode.c.h',
                'hpp': 'sourcecode.cpp.h',
                'hxx': 'sourcecode.cpp.h',
                'hh': 'sourcecode.cpp.hh',
                'inc': 'sourcecode.c.h',
                'swift': 'sourcecode.swift',
                'dylib': 'compiled.mach-o.dylib',
                'o': 'compiled.mach-o.objfile',
                's': 'sourcecode.asm',
                'asm': 'sourcecode.nasm',
                'metal': 'sourcecode.metal',
                'glsl': 'sourcecode.glsl',
                }
NEEDS_CUSTOM_RULES = {'nasm': 'sourcecode.nasm',
                      }
LANGNAMEMAP = {'c': 'C',
               'cpp': 'CPLUSPLUS',
               'objc': 'OBJC',
               'objcpp': 'OBJCPLUSPLUS',
               'swift': 'SWIFT_'
               }
OPT2XCODEOPT = {'plain': None,
                '0': '0',
                'g': '0',
                '1': '1',
                '2': '2',
                '3': '3',
                's': 's',
                }
BOOL2XCODEBOOL = {True: 'YES', False: 'NO'}
LINKABLE_EXTENSIONS = {'.o', '.a', '.obj', '.so', '.dylib'}
XCODEVERSIONS = {'1500': ('Xcode 15.0', 60),
                 '1400': ('Xcode 14.0', 56),
                 '1300': ('Xcode 13.0', 55),
                 '1200': ('Xcode 12.0', 54),
                 '1140': ('Xcode 11.4', 53),
                 '1100': ('Xcode 11.0', 52),
                 '1000': ('Xcode 10.0', 51),
                 '930': ('Xcode 9.3', 50),
                 '800': ('Xcode 8.0', 48),
                 '630': ('Xcode 6.3', 47),
                 '320': ('Xcode 3.2', 46),
                 '310': ('Xcode 3.1', 45)
                 }

def autodetect_xcode_version() -> T.Tuple[str, int]:
    try:
        pc, stdout, stderr = mesonlib.Popen_safe(['xcodebuild', '-version'])
    except FileNotFoundError:
        raise MesonException('Could not detect Xcode. Please install it if you wish to use the Xcode backend.')
    if pc.returncode != 0:
        raise MesonException(f'An error occurred while detecting Xcode: {stderr}')
    version = int(''.join(re.search(r'\d*\.\d*\.*\d*', stdout).group(0).split('.')))
    # If the version number does not have two decimal points, pretend it does.
    if stdout.count('.') < 2:
        version *= 10
    for v, r in XCODEVERSIONS.items():
        if int(v) <= version:
            return r
    raise MesonException('Your Xcode installation is too old and is not supported.')

class FileTreeEntry:

    def __init__(self) -> None:
        self.subdirs: T.Dict[str, FileTreeEntry] = {}
        self.targets: T.List[build.BuildTarget] = []

class PbxArray:
    def __init__(self) -> None:
        self.items: T.List[PbxArrayItem] = []

    def add_item(self, item: T.Union[PbxArrayItem, str], comment: str = '') -> None:
        if isinstance(item, PbxArrayItem):
            self.items.append(item)
        else:
            self.items.append(PbxArrayItem(item, comment))

    def write(self, ofile: T.TextIO, indent_level: int) -> None:
        ofile.write('(\n')
        indent_level += 1
        for i in self.items:
            if i.comment:
                ofile.write(indent_level*INDENT + f'{i.value} {i.comment},\n')
            else:
                ofile.write(indent_level*INDENT + f'{i.value},\n')
        indent_level -= 1
        ofile.write(indent_level*INDENT + ');\n')

class PbxArrayItem:
    def __init__(self, value: str, comment: str = ''):
        self.value = value
        if comment:
            if '/*' in comment:
                self.comment = comment
            else:
                self.comment = f'/* {comment} */'
        else:
            self.comment = comment

class PbxComment:
    def __init__(self, text: str):
        assert isinstance(text, str)
        assert '/*' not in text
        self.text = f'/* {text} */'

    def write(self, ofile: T.TextIO, indent_level: int) -> None:
        ofile.write(f'\n{self.text}\n')

class PbxDictItem:
    def __init__(self, key: str, value: T.Union[PbxArray, PbxDict, str, int], comment: str = ''):
        self.key = key
        if isinstance(value, str):
            self.value = self.quote_value(value)
        else:
            self.value = value

        if comment:
            if '/*' in comment:
                self.comment = comment
            else:
                self.comment = f'/* {comment} */'
        else:
            self.comment = comment

    def quote_value(self, value: str) -> str:
        quoted = f'"{value}"'

        if not value:
            return quoted

        if set(' +@$<>/').isdisjoint(value) or value[0] == '"':
            return value

        return quoted

class PbxDict:
    def __init__(self) -> None:
        # This class is a bit weird, because we want to write PBX dicts in
        # defined order _and_ we want to write intermediate comments also in order.
        self.keys: T.Set[str] = set()
        self.items: T.List[T.Union[PbxDictItem, PbxComment]] = []

    def add_item(self, key: str, value: T.Union[PbxArray, PbxDict, str, int], comment: str = '') -> None:
        assert key not in self.keys
        item = PbxDictItem(key, value, comment)
        self.keys.add(key)
        self.items.append(item)

    def get_item(self, key: str) -> PbxDictItem:
        assert key in self.keys
        for item in self.items:
            if not isinstance(item, PbxDictItem):
                continue
            if item.key == key:
                return item
        return None

    def has_item(self, key: str) -> bool:
        return key in self.keys

    def add_comment(self, comment: PbxComment) -> None:
        assert isinstance(comment, PbxComment)
        self.items.append(comment)

    def write(self, ofile: T.TextIO, indent_level: int) -> None:
        ofile.write('{\n')
        indent_level += 1
        for i in self.items:
            if isinstance(i, PbxComment):
                i.write(ofile, indent_level)
            elif isinstance(i, PbxDictItem):
                if isinstance(i.value, (str, int)):
                    if i.comment:
                        ofile.write(indent_level*INDENT + f'{i.key} = {i.value} {i.comment};\n')
                    else:
                        ofile.write(indent_level*INDENT + f'{i.key} = {i.value};\n')
                elif isinstance(i.value, PbxDict):
                    if i.comment:
                        ofile.write(indent_level*INDENT + f'{i.key} {i.comment} = ')
                    else:
                        ofile.write(indent_level*INDENT + f'{i.key} = ')
                    i.value.write(ofile, indent_level)
                elif isinstance(i.value, PbxArray):
                    if i.comment:
                        ofile.write(indent_level*INDENT + f'{i.key} {i.comment} = ')
                    else:
                        ofile.write(indent_level*INDENT + f'{i.key} = ')
                    i.value.write(ofile, indent_level)
                else:
                    print(i)
                    print(i.key)
                    print(i.value)
                    raise RuntimeError('missing code')
            else:
                print(i)
                raise RuntimeError('missing code2')

        indent_level -= 1
        ofile.write(indent_level*INDENT + '}')
        if indent_level == 0:
            ofile.write('\n')
        else:
            ofile.write(';\n')

class XCodeBackend(backends.Backend):

    name = 'xcode'

    def __init__(self, build: T.Optional[build.Build], interpreter: T.Optional[Interpreter]):
        super().__init__(build, interpreter)
        self.project_uid = self.environment.coredata.lang_guids['default'].replace('-', '')[:24]
        self.buildtype = T.cast('str', self.environment.coredata.optstore.get_value_for(OptionKey('buildtype')))
        self.project_conflist = self.gen_id()
        self.maingroup_id = self.gen_id()
        self.all_id = self.gen_id()
        self.all_buildconf_id = self.gen_id()
        self.buildtypes = [self.buildtype]
        self.test_id = self.gen_id()
        self.test_command_id = self.gen_id()
        self.test_buildconf_id = self.gen_id()
        self.regen_id = self.gen_id()
        self.regen_command_id = self.gen_id()
        self.regen_buildconf_id = self.gen_id()
        self.regen_dependency_id = self.gen_id()
        self.top_level_dict = PbxDict()
        self.generator_outputs = {}
        self.arch = self.build.environment.machines.host.cpu
        if self.arch == 'aarch64':
            self.arch = 'arm64'
        self.xcodeversion, self.objversion = autodetect_xcode_version()
        # In Xcode files are not accessed via their file names, but rather every one of them
        # gets an unique id. More precisely they get one unique id per target they are used
        # in. If you generate only one id per file and use them, compilation will work but the
        # UI will only show the file in one target but not the others. Thus they key is
        # a tuple containing the target and filename.
        self.buildfile_ids = {}
        # That is not enough, though. Each target/file combination also gets a unique id
        # in the file reference section. Because why not. This means that a source file
        # that is used in two targets gets a total of four unique ID numbers.
        self.fileref_ids = {}

    def write_pbxfile(self, top_level_dict, ofilename) -> None:
        tmpname = ofilename + '.tmp'
        with open(tmpname, 'w', encoding='utf-8') as ofile:
            ofile.write('// !$*UTF8*$!\n')
            top_level_dict.write(ofile, 0)
        os.replace(tmpname, ofilename)

    def gen_id(self) -> str:
        return str(uuid.uuid4()).upper().replace('-', '')[:24]

    @functools.lru_cache(maxsize=None)
    def get_target_dir(self, target: T.Union[build.Target, build.CustomTargetIndex]) -> str:
        dirname = os.path.join(target.get_subdir(), T.cast('str', self.environment.coredata.optstore.get_value_for(OptionKey('buildtype'))))
        #os.makedirs(os.path.join(self.environment.get_build_dir(), dirname), exist_ok=True)
        return dirname

    def get_custom_target_output_dir(self, target: T.Union[build.Target, build.CustomTargetIndex]) -> str:
        dirname = target.get_subdir()
        os.makedirs(os.path.join(self.environment.get_build_dir(), dirname), exist_ok=True)
        return dirname

    def object_filename_from_source(self, target: build.BuildTarget, source: mesonlib.FileOrString, targetdir: T.Optional[str] = None) -> str:
        # Xcode has the following naming scheme:
        # projectname.build/debug/prog@exe.build/Objects-normal/x86_64/func.o
        project = self.build.project_name
        buildtype = self.buildtype
        tname = target.get_id()
        if isinstance(source, mesonlib.File):
            source = source.fname
        stem = os.path.splitext(os.path.basename(source))[0]
        # Append "build" before the actual object path to match OBJROOT
        obj_path = f'build/{project}.build/{buildtype}/{tname}.build/Objects-normal/{self.arch}/{stem}.o'
        return obj_path

    def determine_swift_dep_dirs(self, target: build.BuildTarget) -> T.List[str]:
        result: T.List[str] = []
        for l in target.link_targets:
            # Xcode does not recognize our private directories, so we have to use its build directories instead.
            result.append(os.path.join(self.environment.get_build_dir(), self.get_target_dir(l)))
        return result

    def generate(self, capture: bool = False, vslite_ctx: T.Optional[T.Dict] = None) -> None:
        # Check for (currently) unexpected capture arg use cases -
        if capture:
            raise MesonBugException('We do not expect the xcode backend to generate with \'capture = True\'')
        if vslite_ctx:
            raise MesonBugException('We do not expect the xcode backend to be given a valid \'vslite_ctx\'')
        self.serialize_tests()
        # Cache the result as the method rebuilds the array every time it is called.
        self.build_targets = self.build.get_build_targets()
        self.custom_targets = self.build.get_custom_targets()
        self.generate_filemap()
        if self.objversion < 50:
            self.generate_buildstylemap()
        self.generate_build_phase_map()
        self.generate_build_configuration_map()
        self.generate_build_configurationlist_map()
        self.generate_project_configurations_map()
        self.generate_buildall_configurations_map()
        self.generate_test_configurations_map()
        self.generate_native_target_map()
        self.generate_native_frameworks_map()
        self.generate_custom_target_map()
        self.generate_native_target_build_rules_map()
        self.generate_generator_target_map()
        self.generate_source_phase_map()
        self.generate_target_dependency_map()
        self.generate_pbxdep_map()
        self.generate_containerproxy_map()
        self.generate_target_file_maps()
        self.generate_build_file_maps()
        self.proj_dir = os.path.join(self.environment.get_build_dir(), self.build.project_name + '.xcodeproj')
        os.makedirs(self.proj_dir, exist_ok=True)
        self.proj_file = os.path.join(self.proj_dir, 'project.pbxproj')
        objects_dict = self.generate_prefix(self.top_level_dict)
        objects_dict.add_comment(PbxComment('Begin PBXAggregateTarget section'))
        self.generate_pbx_aggregate_target(objects_dict)
        objects_dict.add_comment(PbxComment('End PBXAggregateTarget section'))
        objects_dict.add_comment(PbxComment('Begin PBXBuildFile section'))
        self.generate_pbx_build_file(objects_dict)
        objects_dict.add_comment(PbxComment('End PBXBuildFile section'))
        objects_dict.add_comment(PbxComment('Begin PBXBuildRule section'))
        self.generate_pbx_build_rule(objects_dict)
        objects_dict.add_comment(PbxComment('End PBXBuildRule section'))
        objects_dict.add_comment(PbxComment('Begin PBXBuildStyle section'))
        if self.objversion < 50:
            self.generate_pbx_build_style(objects_dict)
            objects_dict.add_comment(PbxComment('End PBXBuildStyle section'))
            objects_dict.add_comment(PbxComment('Begin PBXContainerItemProxy section'))
        self.generate_pbx_container_item_proxy(objects_dict)
        objects_dict.add_comment(PbxComment('End PBXContainerItemProxy section'))
        objects_dict.add_comment(PbxComment('Begin PBXFileReference section'))
        self.generate_pbx_file_reference(objects_dict)
        objects_dict.add_comment(PbxComment('End PBXFileReference section'))
        objects_dict.add_comment(PbxComment('Begin PBXFrameworksBuildPhase section'))
        self.generate_pbx_frameworks_buildphase(objects_dict)
        objects_dict.add_comment(PbxComment('End PBXFrameworksBuildPhase section'))
        objects_dict.add_comment(PbxComment('Begin PBXGroup section'))
        self.generate_pbx_group(objects_dict)
        objects_dict.add_comment(PbxComment('End PBXGroup section'))
        objects_dict.add_comment(PbxComment('Begin PBXNativeTarget section'))
        self.generate_pbx_native_target(objects_dict)
        objects_dict.add_comment(PbxComment('End PBXNativeTarget section'))
        objects_dict.add_comment(PbxComment('Begin PBXProject section'))
        self.generate_pbx_project(objects_dict)
        objects_dict.add_comment(PbxComment('End PBXProject section'))
        objects_dict.add_comment(PbxComment('Begin PBXShellScriptBuildPhase section'))
        self.generate_pbx_shell_build_phase(objects_dict)
        objects_dict.add_comment(PbxComment('End PBXShellScriptBuildPhase section'))
        objects_dict.add_comment(PbxComment('Begin PBXSourcesBuildPhase section'))
        self.generate_pbx_sources_build_phase(objects_dict)
        objects_dict.add_comment(PbxComment('End PBXSourcesBuildPhase section'))
        objects_dict.add_comment(PbxComment('Begin PBXTargetDependency section'))
        self.generate_pbx_target_dependency(objects_dict)
        objects_dict.add_comment(PbxComment('End PBXTargetDependency section'))
        objects_dict.add_comment(PbxComment('Begin XCBuildConfiguration section'))
        self.generate_xc_build_configuration(objects_dict)
        objects_dict.add_comment(PbxComment('End XCBuildConfiguration section'))
        objects_dict.add_comment(PbxComment('Begin XCConfigurationList section'))
        self.generate_xc_configurationList(objects_dict)
        objects_dict.add_comment(PbxComment('End XCConfigurationList section'))
        self.generate_suffix(self.top_level_dict)
        self.write_pbxfile(self.top_level_dict, self.proj_file)
        self.generate_regen_info()

    def get_xcodetype(self, fname: str) -> str:
        extension = fname.split('.')[-1]
        if extension == 'C':
            extension = 'cpp'
        xcodetype = XCODETYPEMAP.get(extension.lower())
        if not xcodetype:
            xcodetype = 'sourcecode.unknown'
        return xcodetype

    def generate_filemap(self) -> None:
        self.filemap = {} # Key is source file relative to src root.
        self.foldermap = {}
        self.target_filemap = {}
        for name, t in self.build_targets.items():
            for s in t.sources:
                if isinstance(s, mesonlib.File):
                    if '/' in s.fname:
                        # From the top level down, add the folders containing the source file.
                        folder = os.path.split(os.path.dirname(s.fname))
                        while folder:
                            fpath = os.path.join(*folder)
                            # Multiple targets might use the same folders, so store their targets with them.
                            # Otherwise, folders and their source files will appear in the wrong places in Xcode.
                            if (fpath, t) not in self.foldermap:
                                self.foldermap[(fpath, t)] = self.gen_id()
                            else:
                                break
                            folder = folder[:-1]
                    s = os.path.join(s.subdir, s.fname)
                    self.filemap[s] = self.gen_id()
            for o in t.objects:
                if isinstance(o, str):
                    o = os.path.join(t.subdir, o)
                    self.filemap[o] = self.gen_id()
            for e in t.extra_files:
                if isinstance(e, mesonlib.File):
                    e = os.path.join(e.subdir, e.fname)
                    self.filemap[e] = self.gen_id()
                else:
                    e = os.path.join(t.subdir, e)
                    self.filemap[e] = self.gen_id()
            self.target_filemap[name] = self.gen_id()

    def generate_buildstylemap(self) -> None:
        self.buildstylemap = {self.buildtype: self.gen_id()}

    def generate_build_phase_map(self) -> None:
        for tname, t in self.build_targets.items():
            # generate id for our own target-name
            t.buildphasemap = {}
            t.buildphasemap[tname] = self.gen_id()
            # each target can have its own Frameworks/Sources/..., generate id's for those
            t.buildphasemap['Frameworks'] = self.gen_id()
            t.buildphasemap['Resources'] = self.gen_id()
            t.buildphasemap['Sources'] = self.gen_id()

    def generate_build_configuration_map(self) -> None:
        self.buildconfmap = {}
        for t in self.build_targets:
            bconfs = {self.buildtype: self.gen_id()}
            self.buildconfmap[t] = bconfs
        for t in self.custom_targets:
            bconfs = {self.buildtype: self.gen_id()}
            self.buildconfmap[t] = bconfs

    def generate_project_configurations_map(self) -> None:
        self.project_configurations = {self.buildtype: self.gen_id()}

    def generate_buildall_configurations_map(self) -> None:
        self.buildall_configurations = {self.buildtype: self.gen_id()}

    def generate_test_configurations_map(self) -> None:
        self.test_configurations = {self.buildtype: self.gen_id()}

    def generate_build_configurationlist_map(self) -> None:
        self.buildconflistmap = {}
        for t in self.build_targets:
            self.buildconflistmap[t] = self.gen_id()
        for t in self.custom_targets:
            self.buildconflistmap[t] = self.gen_id()

    def generate_native_target_map(self) -> None:
        self.native_targets = {}
        for t in self.build_targets:
            self.native_targets[t] = self.gen_id()

    def generate_native_target_build_rules_map(self) -> None:
        self.build_rules = {}
        for name, target in self.build_targets.items():
            languages = {}
            for language in target.compilers:
                if language not in NEEDS_CUSTOM_RULES:
                    continue
                languages[language] = self.gen_id()
            self.build_rules[name] = languages

    def generate_custom_target_map(self) -> None:
        self.shell_targets = {}
        self.custom_target_output_buildfile = {}
        self.custom_target_output_fileref = {}
        for tname, t in self.custom_targets.items():
            self.shell_targets[tname] = self.gen_id()
            if not isinstance(t, build.CustomTarget):
                continue
            (srcs, ofilenames, cmd) = self.eval_custom_target_command(t)
            for o in ofilenames:
                self.custom_target_output_buildfile[o] = self.gen_id()
                self.custom_target_output_fileref[o] = self.gen_id()

    def generate_generator_target_map(self) -> None:
        # Generator objects do not have natural unique ids
        # so use a counter.
        self.generator_fileref_ids = {}
        self.generator_buildfile_ids = {}
        for tname, t in self.build_targets.items():
            generator_id = 0
            for genlist in t.generated:
                if not isinstance(genlist, build.GeneratedList):
                    continue
                self.gen_single_target_map(genlist, tname, t, generator_id)
                generator_id += 1
        # FIXME add outputs.
        for tname, t in self.custom_targets.items():
            generator_id = 0
            for genlist in t.sources:
                if not isinstance(genlist, build.GeneratedList):
                    continue
                self.gen_single_target_map(genlist, tname, t, generator_id)
                generator_id += 1

    def gen_single_target_map(self, genlist, tname, t, generator_id) -> None:
        k = (tname, generator_id)
        assert k not in self.shell_targets
        self.shell_targets[k] = self.gen_id()
        ofile_abs = []
        for i in genlist.get_inputs():
            for o_base in genlist.get_outputs_for(i):
                o = os.path.join(self.get_target_private_dir(t), o_base)
                ofile_abs.append(os.path.join(self.environment.get_build_dir(), o))
        assert k not in self.generator_outputs
        self.generator_outputs[k] = ofile_abs
        buildfile_ids = []
        fileref_ids = []
        for i in range(len(ofile_abs)):
            buildfile_ids.append(self.gen_id())
            fileref_ids.append(self.gen_id())
        self.generator_buildfile_ids[k] = buildfile_ids
        self.generator_fileref_ids[k] = fileref_ids

    def generate_native_frameworks_map(self) -> None:
        self.native_frameworks = {}
        self.native_frameworks_fileref = {}
        for t in self.build_targets.values():
            for dep in t.get_external_deps():
                if dep.name == 'appleframeworks' and dep.found():
                    for f in dep.frameworks:
                        self.native_frameworks[f] = self.gen_id()
                        self.native_frameworks_fileref[f] = self.gen_id()

    def generate_target_dependency_map(self) -> None:
        self.target_dependency_map = {}
        for tname, t in self.build_targets.items():
            for target in t.link_targets:
                if isinstance(target, build.CustomTargetIndex):
                    k = (tname, target.target.get_basename())
                    if k in self.target_dependency_map:
                        continue
                else:
                    k = (tname, target.get_basename())
                    assert k not in self.target_dependency_map
                self.target_dependency_map[k] = self.gen_id()
        for tname, t in self.custom_targets.items():
            k = tname
            assert k not in self.target_dependency_map
            self.target_dependency_map[k] = self.gen_id()

    def generate_pbxdep_map(self) -> None:
        self.pbx_dep_map = {}
        self.pbx_custom_dep_map = {}
        for t in self.build_targets:
            self.pbx_dep_map[t] = self.gen_id()
        for t in self.custom_targets:
            self.pbx_custom_dep_map[t] = self.gen_id()

    def generate_containerproxy_map(self) -> None:
        self.containerproxy_map = {}
        for t in self.build_targets:
            self.containerproxy_map[t] = self.gen_id()

    def generate_target_file_maps(self) -> None:
        self.generate_target_file_maps_impl(self.build_targets)
        self.generate_target_file_maps_impl(self.custom_targets)

    def generate_target_file_maps_impl(self, targets) -> None:
        for tname, t in targets.items():
            for s in t.sources:
                if isinstance(s, mesonlib.File):
                    s = os.path.join(s.subdir, s.fname)
                if not isinstance(s, str):
                    continue
                k = (tname, s)
                assert k not in self.buildfile_ids
                self.buildfile_ids[k] = self.gen_id()
                assert k not in self.fileref_ids
                self.fileref_ids[k] = self.gen_id()
            if not hasattr(t, 'objects'):
                continue
            for o in t.objects:
                if isinstance(o, build.ExtractedObjects):
                    # Extracted objects do not live in "the Xcode world".
                    continue
                if isinstance(o, mesonlib.File):
                    o = os.path.join(o.subdir, o.fname)
                if isinstance(o, str):
                    o = os.path.join(t.subdir, o)
                    k = (tname, o)
                    assert k not in self.buildfile_ids
                    self.buildfile_ids[k] = self.gen_id()
                    assert k not in self.fileref_ids
                    self.fileref_ids[k] = self.gen_id()
                else:
                    raise RuntimeError('Unknown input type ' + str(o))
            for e in t.extra_files:
                if isinstance(e, mesonlib.File):
                    e = os.path.join(e.subdir, e.fname)
                if isinstance(e, str):
                    e = os.path.join(t.subdir, e)
                    k = (tname, e)
                    assert k not in self.buildfile_ids
                    self.buildfile_ids[k] = self.gen_id()
                    assert k not in self.fileref_ids
                    self.fileref_ids[k] = self.gen_id()

    def generate_build_file_maps(self) -> None:
        for buildfile in self.interpreter.get_build_def_files():
            assert isinstance(buildfile, str)
            self.buildfile_ids[buildfile] = self.gen_id()
            self.fileref_ids[buildfile] = self.gen_id()

    def generate_source_phase_map(self) -> None:
        self.source_phase = {}
        for t in self.build_targets:
            self.source_phase[t] = self.gen_id()

    def generate_pbx_aggregate_target(self, objects_dict: PbxDict) -> None:
        self.custom_aggregate_targets = {}
        self.build_all_tdep_id = self.gen_id()
        target_dependencies = []
        custom_target_dependencies = []
        for tname, t in self.get_build_by_default_targets().items():
            if isinstance(t, build.CustomTarget):
                custom_target_dependencies.append(self.pbx_custom_dep_map[t.get_id()])
            elif isinstance(t, build.BuildTarget):
                target_dependencies.append(self.pbx_dep_map[t.get_id()])
        aggregated_targets = []
        aggregated_targets.append((self.all_id,
                                   'ALL_BUILD',
                                   self.all_buildconf_id,
                                   [],
                                   [self.regen_dependency_id] + target_dependencies + custom_target_dependencies))
        aggregated_targets.append((self.test_id,
                                   'RUN_TESTS',
                                   self.test_buildconf_id,
                                   [self.test_command_id],
                                   [self.regen_dependency_id, self.build_all_tdep_id]))
        aggregated_targets.append((self.regen_id,
                                   'REGENERATE',
                                   self.regen_buildconf_id,
                                   [self.regen_command_id],
                                   []))
        for tname, t in self.build.get_custom_targets().items():
            ct_id = self.gen_id()
            self.custom_aggregate_targets[tname] = ct_id
            build_phases = []
            dependencies = [self.regen_dependency_id]
            generator_id = 0
            for d in t.dependencies:
                if isinstance(d, build.CustomTarget):
                    dependencies.append(self.pbx_custom_dep_map[d.get_id()])
                elif isinstance(d, build.BuildTarget):
                    dependencies.append(self.pbx_dep_map[d.get_id()])
            for s in t.sources:
                if isinstance(s, build.GeneratedList):
                    build_phases.append(self.shell_targets[(tname, generator_id)])
                    for d in s.depends:
                        dependencies.append(self.pbx_custom_dep_map[d.get_id()])
                    generator_id += 1
                elif isinstance(s, build.ExtractedObjects):
                    source_target_id = self.pbx_dep_map[s.target.get_id()]
                    if source_target_id not in dependencies:
                        dependencies.append(source_target_id)
            build_phases.append(self.shell_targets[tname])
            aggregated_targets.append((ct_id, tname, self.buildconflistmap[tname], build_phases, dependencies))

        # Sort objects by ID before writing
        sorted_aggregated_targets = sorted(aggregated_targets, key=operator.itemgetter(0))
        for t in sorted_aggregated_targets:
            agt_dict = PbxDict()
            name = t[1]
            buildconf_id = t[2]
            build_phases = t[3]
            dependencies = t[4]
            agt_dict.add_item('isa', 'PBXAggregateTarget')
            agt_dict.add_item('buildConfigurationList', buildconf_id, f'Build configuration list for PBXAggregateTarget "{name}"')
            bp_arr = PbxArray()
            agt_dict.add_item('buildPhases', bp_arr)
            for bp in build_phases:
                bp_arr.add_item(bp, 'ShellScript')
            dep_arr = PbxArray()
            agt_dict.add_item('dependencies', dep_arr)
            for td in dependencies:
                dep_arr.add_item(td, 'PBXTargetDependency')
            agt_dict.add_item('name', name)
            agt_dict.add_item('productName', name)
            objects_dict.add_item(t[0], agt_dict, name)

    def generate_pbx_build_file(self, objects_dict: PbxDict) -> None:
        for tname, t in self.build_targets.items():
            for dep in t.get_external_deps():
                if dep.name == 'appleframeworks' and dep.found():
                    for f in dep.frameworks:
                        fw_dict = PbxDict()
                        fwkey = self.native_frameworks[f]
                        if fwkey not in objects_dict.keys:
                            objects_dict.add_item(fwkey, fw_dict, f'{f}.framework in Frameworks')
                        fw_dict.add_item('isa', 'PBXBuildFile')
                        fw_dict.add_item('fileRef', self.native_frameworks_fileref[f], f)

            for s in t.sources:
                in_build_dir = False
                if isinstance(s, mesonlib.File):
                    if s.is_built:
                        in_build_dir = True
                    s = os.path.join(s.subdir, s.fname)

                if not isinstance(s, str):
                    continue
                sdict = PbxDict()
                k = (tname, s)
                idval = self.buildfile_ids[k]
                fileref = self.fileref_ids[k]
                if in_build_dir:
                    fullpath = os.path.join(self.environment.get_build_dir(), s)
                else:
                    fullpath = os.path.join(self.environment.get_source_dir(), s)
                sdict.add_item('isa', 'PBXBuildFile')
                sdict.add_item('fileRef', fileref, fullpath)
                objects_dict.add_item(idval, sdict)

            for o in t.objects:
                if isinstance(o, build.ExtractedObjects):
                    # Object files are not source files as such. We add them
                    # by hand in linker flags. It is also not particularly
                    # clear how to define build files in Xcode's file format.
                    continue
                if isinstance(o, mesonlib.File):
                    o = os.path.join(o.subdir, o.fname)
                elif isinstance(o, str):
                    o = os.path.join(t.subdir, o)
                idval = self.buildfile_ids[(tname, o)]
                k = (tname, o)
                fileref = self.fileref_ids[k]
                assert o not in self.filemap
                self.filemap[o] = idval
                fullpath = os.path.join(self.environment.get_source_dir(), o)
                fullpath2 = fullpath
                o_dict = PbxDict()
                objects_dict.add_item(idval, o_dict, fullpath)
                o_dict.add_item('isa', 'PBXBuildFile')
                o_dict.add_item('fileRef', fileref, fullpath2)

            generator_id = 0
            for g in t.generated:
                if not isinstance(g, build.GeneratedList):
                    continue
                self.create_generator_shellphase(objects_dict, tname, generator_id)
                generator_id += 1

        # Custom targets are shell build phases in Xcode terminology.
        for tname, t in self.custom_targets.items():
            if not isinstance(t, build.CustomTarget):
                continue
            (srcs, ofilenames, cmd) = self.eval_custom_target_command(t)
            for o in ofilenames:
                custom_dict = PbxDict()
                objects_dict.add_item(self.custom_target_output_buildfile[o], custom_dict, f'/* {o} */')
                custom_dict.add_item('isa', 'PBXBuildFile')
                custom_dict.add_item('fileRef', self.custom_target_output_fileref[o])
            generator_id = 0
            for g in t.sources:
                if not isinstance(g, build.GeneratedList):
                    continue
                self.create_generator_shellphase(objects_dict, tname, generator_id)
                generator_id += 1

    def create_generator_shellphase(self, objects_dict, tname, generator_id) -> None:
        file_ids = self.generator_buildfile_ids[(tname, generator_id)]
        ref_ids = self.generator_fileref_ids[(tname, generator_id)]
        assert len(ref_ids) == len(file_ids)
        for file_o, ref_id in zip(file_ids, ref_ids):
            odict = PbxDict()
            objects_dict.add_item(file_o, odict)
            odict.add_item('isa', 'PBXBuildFile')
            odict.add_item('fileRef', ref_id)

    # This is skipped if Xcode 9 or above is installed, as PBXBuildStyle was removed on that version.
    def generate_pbx_build_style(self, objects_dict: PbxDict) -> None:
        for name, idval in self.buildstylemap.items():
            styledict = PbxDict()
            objects_dict.add_item(idval, styledict, name)
            styledict.add_item('isa', 'PBXBuildStyle')
            settings_dict = PbxDict()
            styledict.add_item('buildSettings', settings_dict)
            settings_dict.add_item('COPY_PHASE_STRIP', 'NO')
            styledict.add_item('name', name)

    def to_shell_script(self, args: CompilerArgs) -> str:
        quoted_cmd = []
        for c in args:
            quoted_cmd.append(c.replace('"', chr(92) + '"'))
        cmd = ' '.join(quoted_cmd)
        return f"\"#!/bin/sh\\n{cmd}\\n\""

    def generate_pbx_build_rule(self, objects_dict: PbxDict) -> None:
        for name, languages in self.build_rules.items():
            target: BuildTarget = self.build_targets[name]
            for language, idval in languages.items():
                compiler: Compiler = target.compilers[language]
                buildrule = PbxDict()
                buildrule.add_item('isa', 'PBXBuildRule')
                buildrule.add_item('compilerSpec', 'com.apple.compilers.proxy.script')
                if compiler.get_id() != 'yasm':
                    # Yasm doesn't generate escaped build rules
                    buildrule.add_item('dependencyFile', '$(DERIVED_FILE_DIR)/$(INPUT_FILE_BASE).d')
                buildrule.add_item('fileType', NEEDS_CUSTOM_RULES[language])
                inputfiles = PbxArray()
                buildrule.add_item('inputFiles', inputfiles)
                buildrule.add_item('isEditable', '0')
                outputfiles = PbxArray()
                outputfiles.add_item('$(DERIVED_FILE_DIR)/$(INPUT_FILE_BASE).o')
                buildrule.add_item('outputFiles', outputfiles)
                # Do NOT use this parameter. Xcode will accept it from the UI,
                # but the parser will break down inconsistently upon next
                # opening. rdar://FB12144055
                # outputargs = PbxArray()
                # args = self.generate_basic_compiler_args(target, compiler)
                # outputargs.add_item(self.to_shell_script(args))
                # buildrule.add_item('outputFilesCompilerFlags', outputargs)
                commands = CompilerArgs(compiler)
                commands += compiler.get_exelist()
                if compiler.get_id() == 'yasm':
                    # Yasm doesn't generate escaped build rules
                    commands += self.compiler_to_generator_args(target, compiler, output='"$SCRIPT_OUTPUT_FILE_0"', input='"$SCRIPT_INPUT_FILE"', depfile=None)
                else:
                    commands += self.compiler_to_generator_args(target,
                                                                compiler,
                                                                output='"$SCRIPT_OUTPUT_FILE_0"',
                                                                input='"$SCRIPT_INPUT_FILE"',
                                                                depfile='"$(dirname "$SCRIPT_OUTPUT_FILE_0")/$(basename "$SCRIPT_OUTPUT_FILE_0" .o).d"',
                                                                extras=['$OTHER_INPUT_FILE_FLAGS'])
                buildrule.add_item('script', self.to_shell_script(commands))
                objects_dict.add_item(idval, buildrule, 'PBXBuildRule')

    def generate_pbx_container_item_proxy(self, objects_dict: PbxDict) -> None:
        for t in self.build_targets:
            proxy_dict = PbxDict()
            objects_dict.add_item(self.containerproxy_map[t], proxy_dict, 'PBXContainerItemProxy')
            proxy_dict.add_item('isa', 'PBXContainerItemProxy')
            proxy_dict.add_item('containerPortal', self.project_uid, 'Project object')
            proxy_dict.add_item('proxyType', '1')
            proxy_dict.add_item('remoteGlobalIDString', self.native_targets[t])
            proxy_dict.add_item('remoteInfo', t)

    def generate_pbx_file_reference(self, objects_dict: PbxDict) -> None:
        for tname, t in self.build_targets.items():
            for dep in t.get_external_deps():
                if dep.name == 'appleframeworks' and dep.found():
                    for f in dep.frameworks:
                        fw_dict = PbxDict()
                        framework_fileref = self.native_frameworks_fileref[f]
                        if objects_dict.has_item(framework_fileref):
                            continue
                        objects_dict.add_item(framework_fileref, fw_dict, f)
                        fw_dict.add_item('isa', 'PBXFileReference')
                        fw_dict.add_item('lastKnownFileType', 'wrapper.framework')
                        fw_dict.add_item('name', f'{f}.framework')
                        fw_dict.add_item('path', f'System/Library/Frameworks/{f}.framework')
                        fw_dict.add_item('sourceTree', 'SDKROOT')
            for s in t.sources:
                in_build_dir = False
                if isinstance(s, mesonlib.File):
                    if s.is_built:
                        in_build_dir = True
                    s = os.path.join(s.subdir, s.fname)
                if not isinstance(s, str):
                    continue
                idval = self.fileref_ids[(tname, s)]
                fullpath = os.path.join(self.environment.get_source_dir(), s)
                src_dict = PbxDict()
                xcodetype = self.get_xcodetype(s)
                name = os.path.basename(s)
                path = s
                objects_dict.add_item(idval, src_dict, fullpath)
                src_dict.add_item('isa', 'PBXFileReference')
                src_dict.add_item('explicitFileType', xcodetype)
                src_dict.add_item('fileEncoding', '4')
                if in_build_dir:
                    src_dict.add_item('name', name)
                    # This makes no sense. This should say path instead of name
                    # but then the path gets added twice.
                    src_dict.add_item('path', name)
                    src_dict.add_item('sourceTree', 'BUILD_ROOT')
                else:
                    src_dict.add_item('name', name)
                    src_dict.add_item('path', path)
                    src_dict.add_item('sourceTree', 'SOURCE_ROOT')

            generator_id = 0
            for g in t.generated:
                if not isinstance(g, build.GeneratedList):
                    continue
                outputs = self.generator_outputs[(tname, generator_id)]
                ref_ids = self.generator_fileref_ids[tname, generator_id]
                assert len(ref_ids) == len(outputs)
                for o, ref_id in zip(outputs, ref_ids):
                    odict = PbxDict()
                    name = os.path.basename(o)
                    objects_dict.add_item(ref_id, odict, o)
                    xcodetype = self.get_xcodetype(o)
                    rel_name = mesonlib.relpath(o, self.environment.get_source_dir())
                    odict.add_item('isa', 'PBXFileReference')
                    odict.add_item('explicitFileType', xcodetype)
                    odict.add_item('fileEncoding', '4')
                    odict.add_item('name', name)
                    odict.add_item('path', rel_name)
                    odict.add_item('sourceTree', 'SOURCE_ROOT')

                generator_id += 1

            for o in t.objects:
                if isinstance(o, build.ExtractedObjects):
                    # Same as with pbxbuildfile.
                    continue
                if isinstance(o, mesonlib.File):
                    fullpath = o.absolute_path(self.environment.get_source_dir(), self.environment.get_build_dir())
                    o = os.path.join(o.subdir, o.fname)
                else:
                    o = os.path.join(t.subdir, o)
                    fullpath = os.path.join(self.environment.get_source_dir(), o)
                idval = self.fileref_ids[(tname, o)]
                rel_name = mesonlib.relpath(fullpath, self.environment.get_source_dir())
                o_dict = PbxDict()
                name = os.path.basename(o)
                objects_dict.add_item(idval, o_dict, fullpath)
                o_dict.add_item('isa', 'PBXFileReference')
                o_dict.add_item('explicitFileType', self.get_xcodetype(o))
                o_dict.add_item('fileEncoding', '4')
                o_dict.add_item('name', name)
                o_dict.add_item('path', rel_name)
                o_dict.add_item('sourceTree', 'SOURCE_ROOT')

            for e in t.extra_files:
                if isinstance(e, mesonlib.File):
                    e = os.path.join(e.subdir, e.fname)
                else:
                    e = os.path.join(t.subdir, e)
                idval = self.fileref_ids[(tname, e)]
                fullpath = os.path.join(self.environment.get_source_dir(), e)
                e_dict = PbxDict()
                xcodetype = self.get_xcodetype(e)
                name = os.path.basename(e)
                path = e
                objects_dict.add_item(idval, e_dict, fullpath)
                e_dict.add_item('isa', 'PBXFileReference')
                e_dict.add_item('explicitFileType', xcodetype)
                e_dict.add_item('name', name)
                e_dict.add_item('path', path)
                e_dict.add_item('sourceTree', 'SOURCE_ROOT')
        for tname, idval in self.target_filemap.items():
            target_dict = PbxDict()
            objects_dict.add_item(idval, target_dict, tname)
            t = self.build_targets[tname]
            fname = t.get_filename()
            reftype = 0
            if isinstance(t, build.Executable):
                typestr = 'compiled.mach-o.executable'
                path = fname
            elif isinstance(t, build.SharedLibrary):
                typestr = self.get_xcodetype('dummy.dylib')
                path = fname
            else:
                typestr = self.get_xcodetype(fname)
                path = '"%s"' % t.get_filename()
            target_dict.add_item('isa', 'PBXFileReference')
            target_dict.add_item('explicitFileType', typestr)
            target_dict.add_item('path', path)
            target_dict.add_item('refType', reftype)
            target_dict.add_item('sourceTree', 'BUILT_PRODUCTS_DIR')

        for tname, t in self.custom_targets.items():
            if not isinstance(t, build.CustomTarget):
                continue
            (srcs, ofilenames, cmd) = self.eval_custom_target_command(t)
            for s in t.sources:
                if isinstance(s, mesonlib.File):
                    s = os.path.join(s.subdir, s.fname)
                elif isinstance(s, str):
                    s = os.path.join(t.subdir, s)
                else:
                    continue
                custom_dict = PbxDict()
                typestr = self.get_xcodetype(s)
                custom_dict.add_item('isa', 'PBXFileReference')
                custom_dict.add_item('explicitFileType', typestr)
                custom_dict.add_item('name', s)
                custom_dict.add_item('path', s)
                custom_dict.add_item('refType', 0)
                custom_dict.add_item('sourceTree', 'SOURCE_ROOT')
                objects_dict.add_item(self.fileref_ids[(tname, s)], custom_dict)
            for o in ofilenames:
                custom_dict = PbxDict()
                typestr = self.get_xcodetype(o)
                custom_dict.add_item('isa', 'PBXFileReference')
                custom_dict.add_item('explicitFileType', typestr)
                custom_dict.add_item('name', o)
                custom_dict.add_item('path', os.path.join(self.src_to_build, o))
                custom_dict.add_item('refType', 0)
                custom_dict.add_item('sourceTree', 'SOURCE_ROOT')
                objects_dict.add_item(self.custom_target_output_fileref[o], custom_dict)

        for buildfile in self.interpreter.get_build_def_files():
            basename = os.path.split(buildfile)[1]
            buildfile_dict = PbxDict()
            typestr = self.get_xcodetype(buildfile)
            buildfile_dict.add_item('isa', 'PBXFileReference')
            buildfile_dict.add_item('explicitFileType', typestr)
            buildfile_dict.add_item('name', basename)
            buildfile_dict.add_item('path', buildfile)
            buildfile_dict.add_item('refType', 0)
            buildfile_dict.add_item('sourceTree', 'SOURCE_ROOT')
            objects_dict.add_item(self.fileref_ids[buildfile], buildfile_dict)

    def generate_pbx_frameworks_buildphase(self, objects_dict: PbxDict) -> None:
        for t in self.build_targets.values():
            bt_dict = PbxDict()
            objects_dict.add_item(t.buildphasemap['Frameworks'], bt_dict, 'Frameworks')
            bt_dict.add_item('isa', 'PBXFrameworksBuildPhase')
            bt_dict.add_item('buildActionMask', 2147483647)
            file_list = PbxArray()
            bt_dict.add_item('files', file_list)
            for dep in t.get_external_deps():
                if dep.name == 'appleframeworks' and dep.found():
                    for f in dep.frameworks:
                        file_list.add_item(self.native_frameworks[f], f'{f}.framework in Frameworks')
            bt_dict.add_item('runOnlyForDeploymentPostprocessing', 0)

    def generate_pbx_group(self, objects_dict: PbxDict) -> None:
        groupmap = {}
        target_src_map = {}
        for t in self.build_targets:
            groupmap[t] = self.gen_id()
            target_src_map[t] = self.gen_id()
        for t in self.custom_targets:
            groupmap[t] = self.gen_id()
            target_src_map[t] = self.gen_id()
        projecttree_id = self.gen_id()
        resources_id = self.gen_id()
        products_id = self.gen_id()
        frameworks_id = self.gen_id()
        main_dict = PbxDict()
        objects_dict.add_item(self.maingroup_id, main_dict)
        main_dict.add_item('isa', 'PBXGroup')
        main_children = PbxArray()
        main_dict.add_item('children', main_children)
        main_children.add_item(projecttree_id, 'Project tree')
        main_children.add_item(resources_id, 'Resources')
        main_children.add_item(products_id, 'Products')
        main_children.add_item(frameworks_id, 'Frameworks')
        main_dict.add_item('sourceTree', '<group>')

        # Define each folder as a group in Xcode. That way, it can build the file tree correctly.
        # This must be done before the project tree group is generated, as source files are added during that phase.
        for (path, target), id in self.foldermap.items():
            folder_dict = PbxDict()
            objects_dict.add_item(id, folder_dict, path)
            folder_dict.add_item('isa', 'PBXGroup')
            folder_children = PbxArray()
            folder_dict.add_item('children', folder_children)
            folder_dict.add_item('name', '"{}"'.format(path.rsplit('/', 1)[-1]))
            folder_dict.add_item('path', f'"{path}"')
            folder_dict.add_item('sourceTree', 'SOURCE_ROOT')

            # Add any detected subdirectories (not declared as subdir()) here, but only one level higher.
            # Example: In "root", add "root/sub", but not "root/sub/subtwo".
            for path_dep, target_dep in self.foldermap:
                if path_dep.startswith(path) and path_dep.split('/', 1)[0] == path.split('/', 1)[0] and path_dep != path and path_dep.count('/') == path.count('/') + 1 and target == target_dep:
                    folder_children.add_item(self.foldermap[(path_dep, target)], path_dep)

        self.add_projecttree(objects_dict, projecttree_id)

        resource_dict = PbxDict()
        objects_dict.add_item(resources_id, resource_dict, 'Resources')
        resource_dict.add_item('isa', 'PBXGroup')
        resource_children = PbxArray()
        resource_dict.add_item('children', resource_children)
        resource_dict.add_item('name', 'Resources')
        resource_dict.add_item('sourceTree', '<group>')

        frameworks_dict = PbxDict()
        objects_dict.add_item(frameworks_id, frameworks_dict, 'Frameworks')
        frameworks_dict.add_item('isa', 'PBXGroup')
        frameworks_children = PbxArray()
        frameworks_dict.add_item('children', frameworks_children)
        # write frameworks

        for t in self.build_targets.values():
            for dep in t.get_external_deps():
                if dep.name == 'appleframeworks' and dep.found():
                    for f in dep.frameworks:
                        frameworks_children.add_item(self.native_frameworks_fileref[f], f)

        frameworks_dict.add_item('name', 'Frameworks')
        frameworks_dict.add_item('sourceTree', '<group>')

        for tname, t in self.custom_targets.items():
            target_dict = PbxDict()
            objects_dict.add_item(groupmap[tname], target_dict, tname)
            target_dict.add_item('isa', 'PBXGroup')
            target_children = PbxArray()
            target_dict.add_item('children', target_children)
            target_children.add_item(target_src_map[tname], 'Source files')
            if t.subproject:
                target_dict.add_item('name', f'{t.subproject} • {t.name}')
            else:
                target_dict.add_item('name', t.name)
            target_dict.add_item('sourceTree', '<group>')
            source_files_dict = PbxDict()
            objects_dict.add_item(target_src_map[tname], source_files_dict, 'Source files')
            source_files_dict.add_item('isa', 'PBXGroup')
            source_file_children = PbxArray()
            source_files_dict.add_item('children', source_file_children)
            for s in t.sources:
                if isinstance(s, mesonlib.File):
                    s = os.path.join(s.subdir, s.fname)
                elif isinstance(s, str):
                    s = os.path.join(t.subdir, s)
                else:
                    continue
                source_file_children.add_item(self.fileref_ids[(tname, s)], s)
            source_files_dict.add_item('name', 'Source files')
            source_files_dict.add_item('sourceTree', '<group>')

        # And finally products
        product_dict = PbxDict()
        objects_dict.add_item(products_id, product_dict, 'Products')
        product_dict.add_item('isa', 'PBXGroup')
        product_children = PbxArray()
        product_dict.add_item('children', product_children)
        for t in self.build_targets:
            product_children.add_item(self.target_filemap[t], t)
        product_dict.add_item('name', 'Products')
        product_dict.add_item('sourceTree', '<group>')

    def write_group_target_entry(self, objects_dict, t):
        tid = t.get_id()
        group_id = self.gen_id()
        target_dict = PbxDict()
        folder_ids = set()
        objects_dict.add_item(group_id, target_dict, tid)
        target_dict.add_item('isa', 'PBXGroup')
        target_children = PbxArray()
        target_dict.add_item('children', target_children)
        target_dict.add_item('name', f'{t} · target')
        target_dict.add_item('sourceTree', '<group>')
        source_files_dict = PbxDict()
        for s in t.sources:
            if isinstance(s, mesonlib.File):
                # If the file is in a folder, add it to the group representing that folder.
                if '/' in s.fname:
                    folder = '/'.join(s.fname.split('/')[:-1])
                    folder_dict = objects_dict.get_item(self.foldermap[(folder, t)]).value.get_item('children').value
                    temp = os.path.join(s.subdir, s.fname)
                    folder_dict.add_item(self.fileref_ids[(tid, temp)], temp)
                    if self.foldermap[(folder, t)] in folder_ids:
                        continue
                    if len(folder.split('/')) == 1:
                        target_children.add_item(self.foldermap[(folder, t)], folder)
                        folder_ids.add(self.foldermap[(folder, t)])
                    continue
                s = os.path.join(s.subdir, s.fname)
            elif isinstance(s, str):
                s = os.path.join(t.subdir, s)
            else:
                continue
            target_children.add_item(self.fileref_ids[(tid, s)], s)
        for o in t.objects:
            if isinstance(o, build.ExtractedObjects):
                # Do not show built object files in the project tree.
                continue
            if isinstance(o, mesonlib.File):
                o = os.path.join(o.subdir, o.fname)
            else:
                o = os.path.join(t.subdir, o)
            target_children.add_item(self.fileref_ids[(tid, o)], o)
        for e in t.extra_files:
            if isinstance(e, mesonlib.File):
                e = os.path.join(e.subdir, e.fname)
            elif isinstance(e, str):
                e = os.path.join(t.subdir, e)
            else:
                continue
            target_children.add_item(self.fileref_ids[(tid, e)], e)
        source_files_dict.add_item('name', 'Source files')
        source_files_dict.add_item('sourceTree', '<group>')
        return group_id

    def add_projecttree(self, objects_dict, projecttree_id) -> None:
        root_dict = PbxDict()
        objects_dict.add_item(projecttree_id, root_dict, "Root of project tree")
        root_dict.add_item('isa', 'PBXGroup')
        target_children = PbxArray()
        root_dict.add_item('children', target_children)
        root_dict.add_item('name', 'Project root')
        root_dict.add_item('sourceTree', '<group>')

        project_tree = self.generate_project_tree()
        self.write_tree(objects_dict, project_tree, target_children, '')

    def write_tree(self, objects_dict, tree_node, children_array, current_subdir) -> None:
        for subdir_name, subdir_node in tree_node.subdirs.items():
            subdir_dict = PbxDict()
            subdir_children = PbxArray()
            subdir_id = self.gen_id()
            objects_dict.add_item(subdir_id, subdir_dict)
            children_array.add_item(subdir_id)
            subdir_dict.add_item('isa', 'PBXGroup')
            subdir_dict.add_item('children', subdir_children)
            subdir_dict.add_item('name', subdir_name)
            subdir_dict.add_item('sourceTree', '<group>')
            self.write_tree(objects_dict, subdir_node, subdir_children, os.path.join(current_subdir, subdir_name))
        for target in tree_node.targets:
            group_id = self.write_group_target_entry(objects_dict, target)
            children_array.add_item(group_id)
        potentials = [os.path.join(current_subdir, 'meson.build'),
                      os.path.join(current_subdir, 'meson.options'),
                      os.path.join(current_subdir, 'meson_options.txt')]
        for bf in potentials:
            i = self.fileref_ids.get(bf, None)
            if i:
                children_array.add_item(i)

    def generate_project_tree(self) -> FileTreeEntry:
        tree_info = FileTreeEntry()
        for tname, t in self.build_targets.items():
            self.add_target_to_tree(tree_info, t)
        return tree_info

    def add_target_to_tree(self, tree_root: FileTreeEntry, t: build.BuildTarget) -> None:
        current_node = tree_root
        path_segments = t.subdir.split('/')
        for s in path_segments:
            if not s:
                continue
            if s not in current_node.subdirs:
                current_node.subdirs[s] = FileTreeEntry()
            current_node = current_node.subdirs[s]
        current_node.targets.append(t)

    def generate_pbx_native_target(self, objects_dict: PbxDict) -> None:
        for tname, idval in self.native_targets.items():
            ntarget_dict = PbxDict()
            t = self.build_targets[tname]
            objects_dict.add_item(idval, ntarget_dict, tname)
            ntarget_dict.add_item('isa', 'PBXNativeTarget')
            ntarget_dict.add_item('buildConfigurationList', self.buildconflistmap[tname], f'Build configuration list for PBXNativeTarget "{tname}"')
            buildphases_array = PbxArray()
            ntarget_dict.add_item('buildPhases', buildphases_array)
            generator_id = 0
            for g in t.generated:
                # Custom target are handled via inter-target dependencies.
                # Generators are built as a shellscriptbuildphase.
                if isinstance(g, build.GeneratedList):
                    buildphases_array.add_item(self.shell_targets[(tname, generator_id)], f'Generator {generator_id}/{tname}')
                    generator_id += 1
            for bpname, bpval in t.buildphasemap.items():
                buildphases_array.add_item(bpval, f'{bpname} yyy')
            build_rules = PbxArray()
            for language, build_rule_idval in self.build_rules[tname].items():
                build_rules.add_item(build_rule_idval, f'{language}')
            ntarget_dict.add_item('buildRules', build_rules)
            dep_array = PbxArray()
            ntarget_dict.add_item('dependencies', dep_array)
            dep_array.add_item(self.regen_dependency_id)
            # These dependencies only tell Xcode that the deps must be built
            # before this one. They don't set up linkage or anything
            # like that. Those are set up in the XCBuildConfiguration.
            for lt in self.build_targets[tname].link_targets:
                # NOT DOCUMENTED, may need to make different links
                # to same target have different targetdependency item.
                if isinstance(lt, build.CustomTarget):
                    dep_array.add_item(self.pbx_custom_dep_map[lt.get_id()], lt.name)
                elif isinstance(lt, build.CustomTargetIndex):
                    dep_array.add_item(self.pbx_custom_dep_map[lt.target.get_id()], lt.target.name)
                else:
                    idval = self.pbx_dep_map[lt.get_id()]
                    dep_array.add_item(idval, 'PBXTargetDependency')
            for o in t.objects:
                if isinstance(o, build.ExtractedObjects):
                    source_target_id = o.target.get_id()
                    idval = self.pbx_dep_map[source_target_id]
                    dep_array.add_item(idval, 'PBXTargetDependency')
            generator_id = 0
            for o in t.generated:
                if isinstance(o, build.CustomTarget):
                    dep_array.add_item(self.pbx_custom_dep_map[o.get_id()], o.name)
                elif isinstance(o, build.CustomTargetIndex):
                    dep_array.add_item(self.pbx_custom_dep_map[o.target.get_id()], o.target.name)

                generator_id += 1

            ntarget_dict.add_item('name', tname)
            ntarget_dict.add_item('productName', tname)
            ntarget_dict.add_item('productReference', self.target_filemap[tname], tname)
            if isinstance(t, build.Executable):
                typestr = 'com.apple.product-type.tool'
            elif isinstance(t, build.StaticLibrary):
                typestr = 'com.apple.product-type.library.static'
            elif isinstance(t, build.SharedLibrary):
                typestr = 'com.apple.product-type.library.dynamic'
            else:
                raise MesonException('Unknown target type for %s' % tname)
            ntarget_dict.add_item('productType', f'"{typestr}"')

    def generate_pbx_project(self, objects_dict: PbxDict) -> None:
        project_dict = PbxDict()
        objects_dict.add_item(self.project_uid, project_dict, 'Project object')
        project_dict.add_item('isa', 'PBXProject')
        attr_dict = PbxDict()
        project_dict.add_item('attributes', attr_dict)
        attr_dict.add_item('BuildIndependentTargetsInParallel', 'YES')
        project_dict.add_item('buildConfigurationList', self.project_conflist, f'Build configuration list for PBXProject "{self.build.project_name}"')
        project_dict.add_item('buildSettings', PbxDict())
        if self.objversion < 50:
            style_arr = PbxArray()
            project_dict.add_item('buildStyles', style_arr)
            for name, idval in self.buildstylemap.items():
                style_arr.add_item(idval, name)
        project_dict.add_item('compatibilityVersion', self.xcodeversion)
        project_dict.add_item('hasScannedForEncodings', 0)
        project_dict.add_item('mainGroup', self.maingroup_id)
        project_dict.add_item('projectDirPath', self.environment.get_source_dir())
        project_dict.add_item('projectRoot', '')
        targets_arr = PbxArray()
        project_dict.add_item('targets', targets_arr)
        targets_arr.add_item(self.all_id, 'ALL_BUILD')
        targets_arr.add_item(self.test_id, 'RUN_TESTS')
        targets_arr.add_item(self.regen_id, 'REGENERATE')
        for t in self.build_targets:
            targets_arr.add_item(self.native_targets[t], t)
        for t in self.custom_targets:
            targets_arr.add_item(self.custom_aggregate_targets[t], t)

    def generate_pbx_shell_build_phase(self, objects_dict: PbxDict) -> None:
        self.generate_test_shell_build_phase(objects_dict)
        self.generate_regen_shell_build_phase(objects_dict)
        self.generate_custom_target_shell_build_phases(objects_dict)
        self.generate_generator_target_shell_build_phases(objects_dict)

    def generate_test_shell_build_phase(self, objects_dict: PbxDict) -> None:
        shell_dict = PbxDict()
        objects_dict.add_item(self.test_command_id, shell_dict, 'ShellScript')
        shell_dict.add_item('isa', 'PBXShellScriptBuildPhase')
        shell_dict.add_item('buildActionMask', 2147483647)
        shell_dict.add_item('files', PbxArray())
        shell_dict.add_item('inputPaths', PbxArray())
        shell_dict.add_item('outputPaths', PbxArray())
        shell_dict.add_item('runOnlyForDeploymentPostprocessing', 0)
        shell_dict.add_item('shellPath', '/bin/sh')
        cmd = mesonlib.get_meson_command() + ['test', '--no-rebuild', '-C', self.environment.get_build_dir()]
        cmdstr = ' '.join(["'%s'" % i for i in cmd])
        shell_dict.add_item('shellScript', cmdstr)
        shell_dict.add_item('showEnvVarsInLog', 0)

    def generate_regen_shell_build_phase(self, objects_dict: PbxDict) -> None:
        shell_dict = PbxDict()
        objects_dict.add_item(self.regen_command_id, shell_dict, 'ShellScript')
        shell_dict.add_item('isa', 'PBXShellScriptBuildPhase')
        shell_dict.add_item('buildActionMask', 2147483647)
        shell_dict.add_item('files', PbxArray())
        shell_dict.add_item('inputPaths', PbxArray())
        shell_dict.add_item('outputPaths', PbxArray())
        shell_dict.add_item('runOnlyForDeploymentPostprocessing', 0)
        shell_dict.add_item('shellPath', '/bin/sh')
        cmd = mesonlib.get_meson_command() + ['--internal', 'regencheck', os.path.join(self.environment.get_build_dir(), 'meson-private')]
        cmdstr = ' '.join(["'%s'" % i for i in cmd])
        shell_dict.add_item('shellScript', cmdstr)
        shell_dict.add_item('showEnvVarsInLog', 0)

    def generate_custom_target_shell_build_phases(self, objects_dict: PbxDict) -> None:
        # Custom targets are shell build phases in Xcode terminology.
        for tname, t in self.custom_targets.items():
            if not isinstance(t, build.CustomTarget):
                continue
            (srcs, ofilenames, cmd) = self.eval_custom_target_command(t, absolute_outputs=True)
            fixed_cmd, _ = self.as_meson_exe_cmdline(cmd[0],
                                                     cmd[1:],
                                                     capture=ofilenames[0] if t.capture else None,
                                                     feed=srcs[0] if t.feed else None,
                                                     env=t.env)
            custom_dict = PbxDict()
            objects_dict.add_item(self.shell_targets[tname], custom_dict, f'/* Custom target {tname} */')
            custom_dict.add_item('isa', 'PBXShellScriptBuildPhase')
            custom_dict.add_item('buildActionMask', 2147483647)
            custom_dict.add_item('files', PbxArray())
            custom_dict.add_item('inputPaths', PbxArray())
            outarray = PbxArray()
            custom_dict.add_item('name', '"Generate {}."'.format(ofilenames[0]))
            custom_dict.add_item('outputPaths', outarray)
            for o in ofilenames:
                outarray.add_item(os.path.join(self.environment.get_build_dir(), o))
            custom_dict.add_item('runOnlyForDeploymentPostprocessing', 0)
            custom_dict.add_item('shellPath', '/bin/sh')
            workdir = self.environment.get_build_dir()
            quoted_cmd = []
            for c in fixed_cmd:
                quoted_cmd.append(c.replace('"', chr(92) + '"'))
            cmdstr = ' '.join([f"\\'{x}\\'" for x in quoted_cmd])
            custom_dict.add_item('shellScript', f'cd \'{workdir}\'; {cmdstr}')
            custom_dict.add_item('showEnvVarsInLog', 0)

    def generate_generator_target_shell_build_phases(self, objects_dict: PbxDict) -> None:
        for tname, t in self.build_targets.items():
            generator_id = 0
            for genlist in t.generated:
                if isinstance(genlist, build.GeneratedList):
                    self.generate_single_generator_phase(tname, t, genlist, generator_id, objects_dict)
                    generator_id += 1
        for tname, t in self.custom_targets.items():
            generator_id = 0
            for genlist in t.sources:
                if isinstance(genlist, build.GeneratedList):
                    self.generate_single_generator_phase(tname, t, genlist, generator_id, objects_dict)
                    generator_id += 1

    def generate_single_generator_phase(self, tname, t, genlist, generator_id, objects_dict) -> None:
        # TODO: this should be rewritten to use the meson wrapper, like the other generators do
        # Currently it doesn't handle a host binary that requires an exe wrapper correctly.
        generator = genlist.get_generator()
        exe = generator.get_exe()
        exe_arr = self.build_target_to_cmd_array(exe)
        workdir = self.environment.get_build_dir()
        target_private_dir = self.relpath(self.get_target_private_dir(t), self.get_target_dir(t))
        gen_dict = PbxDict()
        objects_dict.add_item(self.shell_targets[(tname, generator_id)], gen_dict, f'Generator {generator_id}/{tname}')
        infilelist = genlist.get_inputs()
        outfilelist = genlist.get_outputs()
        gen_dict.add_item('isa', 'PBXShellScriptBuildPhase')
        gen_dict.add_item('buildActionMask', 2147483647)
        gen_dict.add_item('files', PbxArray())
        gen_dict.add_item('inputPaths', PbxArray())
        gen_dict.add_item('name', f'Generator {generator_id}/{tname}')
        commands = [["cd", workdir]] # Array of arrays, each one a single command, will get concatenated below.
        k = (tname, generator_id)
        ofile_abs = self.generator_outputs[k]
        outarray = PbxArray()
        gen_dict.add_item('outputPaths', outarray)
        for of in ofile_abs:
            outarray.add_item(of)
        for i in infilelist:
            # This might be needed to be added to inputPaths. It's not done yet as it is
            # unclear whether it is necessary, what actually happens when it is defined
            # and currently the build works without it.
            #infile_abs = i.absolute_path(self.environment.get_source_dir(), self.environment.get_build_dir())
            infilename = i.rel_to_builddir(self.build_to_src, target_private_dir)
            base_args = generator.get_arglist(infilename)
            for o_base in genlist.get_outputs_for(i):
                o = os.path.join(self.get_target_private_dir(t), o_base)
                args = []
                for arg in base_args:
                    arg = arg.replace("@INPUT@", infilename)
                    arg = arg.replace('@OUTPUT@', o).replace('@BUILD_DIR@', self.get_target_private_dir(t))
                    arg = arg.replace("@CURRENT_SOURCE_DIR@", os.path.join(self.build_to_src, t.subdir))
                    args.append(arg)
                args = self.replace_outputs(args, self.get_target_private_dir(t), outfilelist)
                args = self.replace_extra_args(args, genlist)
                if generator.capture:
                    # When capturing, stdout is the output. Forward it with the shell.
                    full_command = ['('] + exe_arr + args + ['>', o, ')']
                else:
                    full_command = exe_arr + args
                commands.append(full_command)
        gen_dict.add_item('runOnlyForDeploymentPostprocessing', 0)
        gen_dict.add_item('shellPath', '/bin/sh')
        quoted_cmds = []
        for cmnd in commands:
            q = []
            for c in cmnd:
                if ' ' in c:
                    q.append(f'\\"{c}\\"')
                else:
                    q.append(c)
            quoted_cmds.append(' '.join(q))
        cmdstr = ' && '.join(quoted_cmds)
        gen_dict.add_item('shellScript', cmdstr)
        gen_dict.add_item('showEnvVarsInLog', 0)

    def generate_pbx_sources_build_phase(self, objects_dict: PbxDict) -> None:
        for name in self.source_phase:
            phase_dict = PbxDict()
            t = self.build_targets[name]
            objects_dict.add_item(t.buildphasemap[name], phase_dict, 'Sources')
            phase_dict.add_item('isa', 'PBXSourcesBuildPhase')
            phase_dict.add_item('buildActionMask', 2147483647)
            file_arr = PbxArray()
            phase_dict.add_item('files', file_arr)
            for s in self.build_targets[name].sources:
                s = os.path.join(s.subdir, s.fname)
                if not self.environment.is_header(s):
                    file_arr.add_item(self.buildfile_ids[(name, s)], os.path.join(self.environment.get_source_dir(), s))
            generator_id = 0
            for gt in t.generated:
                if isinstance(gt, build.CustomTarget):
                    (srcs, ofilenames, cmd) = self.eval_custom_target_command(gt)
                    for o in ofilenames:
                        file_arr.add_item(self.custom_target_output_buildfile[o],
                                          os.path.join(self.environment.get_build_dir(), o))
                elif isinstance(gt, build.CustomTargetIndex):
                    output_dir = self.get_custom_target_output_dir(gt)
                    for o in gt.get_outputs():
                        o = os.path.join(output_dir, o)
                        file_arr.add_item(self.custom_target_output_buildfile[o],
                                          os.path.join(self.environment.get_build_dir(), o))
                elif isinstance(gt, build.GeneratedList):
                    genfiles = self.generator_buildfile_ids[(name, generator_id)]
                    generator_id += 1
                    for o in genfiles:
                        file_arr.add_item(o)
                else:
                    raise RuntimeError('Unknown input type: ' + str(gt))
            phase_dict.add_item('runOnlyForDeploymentPostprocessing', 0)

    def generate_pbx_target_dependency(self, objects_dict: PbxDict) -> None:
        all_dict = PbxDict()
        objects_dict.add_item(self.build_all_tdep_id, all_dict, 'ALL_BUILD')
        all_dict.add_item('isa', 'PBXTargetDependency')
        all_dict.add_item('target', self.all_id)
        targets = []
        targets.append((self.regen_dependency_id, self.regen_id, 'REGEN', None))
        for t in self.build_targets:
            idval = self.pbx_dep_map[t] # VERIFY: is this correct?
            targets.append((idval, self.native_targets[t], t, self.containerproxy_map[t]))

        for t in self.custom_targets:
            idval = self.pbx_custom_dep_map[t]
            targets.append((idval, self.custom_aggregate_targets[t], t, None)) # self.containerproxy_map[t]))

        # Sort object by ID
        sorted_targets = sorted(targets, key=operator.itemgetter(0))
        for t in sorted_targets:
            t_dict = PbxDict()
            objects_dict.add_item(t[0], t_dict, 'PBXTargetDependency')
            t_dict.add_item('isa', 'PBXTargetDependency')
            t_dict.add_item('target', t[1], t[2])
            if t[3] is not None:
                t_dict.add_item('targetProxy', t[3], 'PBXContainerItemProxy')

    def generate_xc_build_configuration(self, objects_dict: PbxDict) -> None:
        # First the setup for the toplevel project.
        for buildtype in self.buildtypes:
            bt_dict = PbxDict()
            objects_dict.add_item(self.project_configurations[buildtype], bt_dict, buildtype)
            bt_dict.add_item('isa', 'XCBuildConfiguration')
            settings_dict = PbxDict()
            bt_dict.add_item('buildSettings', settings_dict)
            settings_dict.add_item('ARCHS', self.arch)
            settings_dict.add_item('BUILD_DIR', self.environment.get_build_dir())
            settings_dict.add_item('BUILD_ROOT', '$(BUILD_DIR)')
            settings_dict.add_item('ONLY_ACTIVE_ARCH', 'YES')
            settings_dict.add_item('SWIFT_VERSION', '5.0')
            settings_dict.add_item('SDKROOT', 'macosx')
            settings_dict.add_item('OBJROOT', '$(BUILD_DIR)/build')
            bt_dict.add_item('name', buildtype)

        # Then the all target.
        for buildtype in self.buildtypes:
            bt_dict = PbxDict()
            objects_dict.add_item(self.buildall_configurations[buildtype], bt_dict, buildtype)
            bt_dict.add_item('isa', 'XCBuildConfiguration')
            settings_dict = PbxDict()
            bt_dict.add_item('buildSettings', settings_dict)
            warn_array = PbxArray()
            warn_array.add_item('"$(inherited)"')
            settings_dict.add_item('WARNING_CFLAGS', warn_array)

            bt_dict.add_item('name', buildtype)

        # Then the test target.
        for buildtype in self.buildtypes:
            bt_dict = PbxDict()
            objects_dict.add_item(self.test_configurations[buildtype], bt_dict, buildtype)
            bt_dict.add_item('isa', 'XCBuildConfiguration')
            settings_dict = PbxDict()
            bt_dict.add_item('buildSettings', settings_dict)
            warn_array = PbxArray()
            settings_dict.add_item('WARNING_CFLAGS', warn_array)
            warn_array.add_item('"$(inherited)"')
            bt_dict.add_item('name', buildtype)

        # Now finally targets.
        for target_name, target in self.build_targets.items():
            self.generate_single_build_target(objects_dict, target_name, target)

        for target_name, target in self.custom_targets.items():
            bt_dict = PbxDict()
            objects_dict.add_item(self.buildconfmap[target_name][buildtype], bt_dict, buildtype)
            bt_dict.add_item('isa', 'XCBuildConfiguration')
            settings_dict = PbxDict()
            bt_dict.add_item('buildSettings', settings_dict)
            settings_dict.add_item('ARCHS', self.arch)
            settings_dict.add_item('ONLY_ACTIVE_ARCH', 'YES')
            settings_dict.add_item('SDKROOT', 'macosx')
            bt_dict.add_item('name', buildtype)

    def determine_internal_dep_link_args(self, target, buildtype):
        links_dylib = False
        dep_libs = []
        for l in target.link_targets:
            if isinstance(target, build.SharedModule) and isinstance(l, build.Executable):
                continue
            if isinstance(l, build.CustomTargetIndex):
                rel_dir = self.get_custom_target_output_dir(l.target)
                libname = l.get_filename()
            elif isinstance(l, build.CustomTarget):
                rel_dir = self.get_custom_target_output_dir(l)
                libname = l.get_filename()
            else:
                rel_dir = self.get_target_dir(l)
                libname = l.get_filename()
            abs_path = os.path.join(self.environment.get_build_dir(), rel_dir, libname)
            dep_libs.append("'%s'" % abs_path)
            if isinstance(l, build.SharedLibrary):
                links_dylib = True
            if isinstance(l, build.StaticLibrary):
                (sub_libs, sub_links_dylib) = self.determine_internal_dep_link_args(l, buildtype)
                dep_libs += sub_libs
                links_dylib = links_dylib or sub_links_dylib
        return (dep_libs, links_dylib)

    def generate_single_build_target(self, objects_dict, target_name, target) -> None:
        for buildtype in self.buildtypes:
            dep_libs = []
            links_dylib = False
            headerdirs = []
            bridging_header = ""
            is_swift = self.is_swift_target(target)
            langs = set()
            for d in target.include_dirs:
                for sd in d.incdirs:
                    cd = os.path.join(d.curdir, sd)
                    headerdirs.append(os.path.join(self.environment.get_source_dir(), cd))
                    headerdirs.append(os.path.join(self.environment.get_build_dir(), cd))
                for extra in d.extra_build_dirs:
                    headerdirs.append(os.path.join(self.environment.get_build_dir(), extra))
            # Swift can import declarations from C-based code using bridging headers.
            # There can only be one header, and it must be included as a source file.
            for i in target.get_sources():
                if self.environment.is_header(i) and is_swift:
                    relh = i.rel_to_builddir(self.build_to_src)
                    bridging_header = os.path.normpath(os.path.join(self.environment.get_build_dir(), relh))
                    break
            (dep_libs, links_dylib) = self.determine_internal_dep_link_args(target, buildtype)
            if links_dylib:
                dep_libs = ['-Wl,-search_paths_first', '-Wl,-headerpad_max_install_names'] + dep_libs
            dylib_version = None
            if isinstance(target, build.SharedLibrary):
                if isinstance(target, build.SharedModule):
                    ldargs = []
                else:
                    ldargs = ['-dynamiclib']
                ldargs += ['-Wl,-headerpad_max_install_names'] + dep_libs
                install_path = os.path.join(self.environment.get_build_dir(), target.subdir, buildtype)
                dylib_version = target.soversion
            else:
                ldargs = dep_libs
                install_path = ''
            if dylib_version is not None:
                product_name = target.get_basename() + '.' + dylib_version
            else:
                product_name = target.get_basename()
            ldargs += target.link_args
            # Swift is special. Again. You can't mix Swift with other languages
            # in the same target. Thus for Swift we only use
            if is_swift:
                linker, stdlib_args = target.compilers['swift'], []
            else:
                linker, stdlib_args = self.determine_linker_and_stdlib_args(target)
            if not isinstance(target, build.StaticLibrary):
                ldargs += self.build.get_project_link_args(linker, target.subproject, target.for_machine)
                ldargs += self.build.get_global_link_args(linker, target.for_machine)
            cargs = []
            for dep in target.get_external_deps():
                cargs += dep.get_compile_args()
                ldargs += dep.get_link_args()
            for o in target.objects:
                # Add extracted objects to the link line by hand.
                if isinstance(o, build.ExtractedObjects):
                    added_objs = set()
                    for objname_rel in self.determine_ext_objs(o):
                        objname_abs = os.path.join(self.environment.get_build_dir(), o.target.subdir, objname_rel)
                        if objname_abs not in added_objs:
                            added_objs.add(objname_abs)
                            ldargs += [r'\"' + objname_abs + r'\"']
            generator_id = 0
            for o in target.generated:
                if isinstance(o, build.GeneratedList):
                    outputs = self.generator_outputs[target_name, generator_id]
                    generator_id += 1
                    for o_abs in outputs:
                        if o_abs.endswith('.o') or o_abs.endswith('.obj'):
                            ldargs += [r'\"' + o_abs + r'\"']
                else:
                    if isinstance(o, build.CustomTarget):
                        (srcs, ofilenames, cmd) = self.eval_custom_target_command(o)
                        for ofname in ofilenames:
                            if os.path.splitext(ofname)[-1] in LINKABLE_EXTENSIONS:
                                ldargs += [r'\"' + os.path.join(self.environment.get_build_dir(), ofname) + r'\"']
                    elif isinstance(o, build.CustomTargetIndex):
                        for ofname in o.get_outputs():
                            if os.path.splitext(ofname)[-1] in LINKABLE_EXTENSIONS:
                                ldargs += [r'\"' + os.path.join(self.environment.get_build_dir(), ofname) + r'\"']
                    else:
                        raise RuntimeError(o)
            if isinstance(target, build.SharedModule):
                ldargs += linker.get_std_shared_module_link_args(target.get_options())
            elif isinstance(target, build.SharedLibrary):
                ldargs += linker.get_std_shared_lib_link_args()
            ldstr = ' '.join(ldargs)
            valid = self.buildconfmap[target_name][buildtype]
            langargs = {}
            for lang in self.environment.coredata.compilers[target.for_machine]:
                if lang not in LANGNAMEMAP:
                    continue
                compiler = target.compilers.get(lang)
                if compiler is None:
                    continue
                # Start with warning args
                warn_args = compiler.get_warn_args(self.get_target_option(target, 'warning_level'))
                std_args = compiler.get_option_compile_args(target, self.environment, target.subproject)
                std_args += compiler.get_option_std_args(target, self.environment, target.subproject)
                # Add compile args added using add_project_arguments()
                pargs = self.build.projects_args[target.for_machine].get(target.subproject, {}).get(lang, [])
                # Add compile args added using add_global_arguments()
                # These override per-project arguments
                gargs = self.build.global_args[target.for_machine].get(lang, [])
                targs = target.get_extra_args(lang)
                args = warn_args + std_args + pargs + gargs + targs
                if lang == 'swift':
                    # For some reason putting Swift module dirs in HEADER_SEARCH_PATHS does not work,
                    # but adding -I/path to manual args does work.
                    swift_dep_dirs = self.determine_swift_dep_dirs(target)
                    for d in swift_dep_dirs:
                        args += compiler.get_include_args(d, False)
                if args:
                    cti_args = []
                    if compiler and target.implicit_include_directories:
                        # It is unclear what is the cwd when xcode runs. -I. does not seem to
                        # add the root build dir to the search path. So add an absolute path instead.
                        # This may break reproducible builds, in which case patches are welcome.
                        cti_args = self.get_custom_target_dir_include_args(target, compiler, absolute_path=True)
                    # Xcode cannot handle separate compilation flags for C and ObjectiveC. They are both
                    # put in OTHER_CFLAGS. Same with C++ and ObjectiveC++.
                    if lang == 'objc':
                        lang = 'c'
                    elif lang == 'objcpp':
                        lang = 'cpp'
                    langs.add(lang)
                    langname = LANGNAMEMAP[lang]
                    langargs.setdefault(langname, [])
                    langargs[langname] = cargs + cti_args + args
            symroot = os.path.join(self.environment.get_build_dir(), target.subdir).rstrip('/')
            bt_dict = PbxDict()
            objects_dict.add_item(valid, bt_dict, buildtype)
            bt_dict.add_item('isa', 'XCBuildConfiguration')
            settings_dict = PbxDict()
            bt_dict.add_item('buildSettings', settings_dict)
            settings_dict.add_item('COMBINE_HIDPI_IMAGES', 'YES')
            if isinstance(target, build.SharedModule):
                settings_dict.add_item('DYLIB_CURRENT_VERSION', '')
                settings_dict.add_item('DYLIB_COMPATIBILITY_VERSION', '')
            else:
                if dylib_version is not None:
                    settings_dict.add_item('DYLIB_CURRENT_VERSION', str(dylib_version))
            if target.prefix:
                settings_dict.add_item('EXECUTABLE_PREFIX', target.prefix)
            if target.suffix:
                suffix = '.' + target.suffix
                settings_dict.add_item('EXECUTABLE_SUFFIX', suffix)
            settings_dict.add_item('GCC_GENERATE_DEBUGGING_SYMBOLS', BOOL2XCODEBOOL[self.get_target_option(target, 'debug')])
            settings_dict.add_item('GCC_INLINES_ARE_PRIVATE_EXTERN', 'NO')
            opt_flag = OPT2XCODEOPT[self.get_target_option(target, 'optimization')]
            if opt_flag is not None:
                settings_dict.add_item('GCC_OPTIMIZATION_LEVEL', opt_flag)
            if target.has_pch:
                # Xcode uses GCC_PREFIX_HEADER which only allows one file per target/executable. Precompiling various header files and
                # applying a particular pch to each source file will require custom scripts (as a build phase) and build flags per each
                # file. Since Xcode itself already discourages precompiled headers in favor of modules we don't try much harder here.
                pchs = target.get_pch('c') + target.get_pch('cpp') + target.get_pch('objc') + target.get_pch('objcpp')
                # Make sure to use headers (other backends require implementation files like *.c *.cpp, etc; these should not be used here)
                pchs = [pch for pch in pchs if pch.endswith('.h') or pch.endswith('.hh') or pch.endswith('hpp')]
                if pchs:
                    if len(pchs) > 1:
                        mlog.warning(f'Unsupported Xcode configuration: More than 1 precompiled header found "{pchs!s}". Target "{target.name}" might not compile correctly.')
                    relative_pch_path = os.path.join(target.get_subdir(), pchs[0]) # Path relative to target so it can be used with "$(PROJECT_DIR)"
                    settings_dict.add_item('GCC_PRECOMPILE_PREFIX_HEADER', 'YES')
                    settings_dict.add_item('GCC_PREFIX_HEADER', f'$(PROJECT_DIR)/{relative_pch_path}')
            settings_dict.add_item('GCC_PREPROCESSOR_DEFINITIONS', '')
            settings_dict.add_item('GCC_SYMBOLS_PRIVATE_EXTERN', 'NO')
            unquoted_headers = [self.get_target_private_dir_abs(target)]
            if target.implicit_include_directories:
                unquoted_headers.append(os.path.join(self.environment.get_build_dir(), target.get_subdir()))
                unquoted_headers.append(os.path.join(self.environment.get_source_dir(), target.get_subdir()))
            unquoted_headers += headerdirs
            settings_dict.add_item('HEADER_SEARCH_PATHS', self.normalize_header_search_paths(unquoted_headers))
            settings_dict.add_item('INSTALL_PATH', install_path)
            settings_dict.add_item('LIBRARY_SEARCH_PATHS', '')
            if isinstance(target, build.SharedModule):
                settings_dict.add_item('LIBRARY_STYLE', 'BUNDLE')
                settings_dict.add_item('MACH_O_TYPE', 'mh_bundle')
            elif isinstance(target, build.SharedLibrary):
                settings_dict.add_item('LIBRARY_STYLE', 'DYNAMIC')
            self.add_otherargs(settings_dict, langargs)
            settings_dict.add_item('OTHER_LDFLAGS', ldstr)
            settings_dict.add_item('OTHER_REZFLAGS', '')
            settings_dict.add_item('PRODUCT_NAME', product_name)
            settings_dict.add_item('SECTORDER_FLAGS', '')
            if is_swift and bridging_header:
                settings_dict.add_item('SWIFT_OBJC_BRIDGING_HEADER', bridging_header)
                if self.objversion >= 60 and target.uses_swift_cpp_interop():
                    settings_dict.add_item('SWIFT_OBJC_INTEROP_MODE', 'objcxx')
            settings_dict.add_item('BUILD_DIR', symroot)
            settings_dict.add_item('OBJROOT', f'{symroot}/build')
            sysheader_arr = PbxArray()
            # XCode will change every -I flag that points inside these directories
            # to an -isystem. Thus set nothing in it since we control our own
            # include flags.
            settings_dict.add_item('SYSTEM_HEADER_SEARCH_PATHS', sysheader_arr)
            settings_dict.add_item('USE_HEADERMAP', 'NO')
            warn_array = PbxArray()
            settings_dict.add_item('WARNING_CFLAGS', warn_array)
            warn_array.add_item('"$(inherited)"')
            bt_dict.add_item('name', buildtype)

    def normalize_header_search_paths(self, header_dirs) -> PbxArray:
        header_arr = PbxArray()
        for i in header_dirs:
            np = os.path.normpath(i)
            # Make sure Xcode will not split single path into separate entries, escaping space with a slash is not enough
            item = f'"\\\"{np}\\\""' if ' ' in np else f'"{np}"'
            header_arr.add_item(item)
        return header_arr

    def add_otherargs(self, settings_dict, langargs):
        for langname, args in langargs.items():
            if args:
                quoted_args = []
                for a in args:
                    # This works but
                    # a) it's ugly as sin
                    # b) I don't know why it works or why every backslash must be escaped into eight backslashes
                    a = a.replace(chr(92), 8*chr(92)) # chr(92) is backslash, this how we smuggle it in without Python's quoting grabbing it.
                    a = a.replace(r'"', r'\\\"')
                    if ' ' in a or "'" in a:
                        a = r'\"' + a + r'\"'
                    quoted_args.append(a)
                settings_dict.add_item(f'OTHER_{langname}FLAGS', ' '.join(quoted_args))

    def generate_xc_configurationList(self, objects_dict: PbxDict) -> None:
        # FIXME: sort items
        conf_dict = PbxDict()
        objects_dict.add_item(self.project_conflist, conf_dict, f'Build configuration list for PBXProject "{self.build.project_name}"')
        conf_dict.add_item('isa', 'XCConfigurationList')
        confs_arr = PbxArray()
        conf_dict.add_item('buildConfigurations', confs_arr)
        for buildtype in self.buildtypes:
            confs_arr.add_item(self.project_configurations[buildtype], buildtype)
        conf_dict.add_item('defaultConfigurationIsVisible', 0)
        conf_dict.add_item('defaultConfigurationName', self.buildtype)

        # Now the all target
        all_dict = PbxDict()
        objects_dict.add_item(self.all_buildconf_id, all_dict, 'Build configuration list for PBXAggregateTarget "ALL_BUILD"')
        all_dict.add_item('isa', 'XCConfigurationList')
        conf_arr = PbxArray()
        all_dict.add_item('buildConfigurations', conf_arr)
        for buildtype in self.buildtypes:
            conf_arr.add_item(self.buildall_configurations[buildtype], buildtype)
        all_dict.add_item('defaultConfigurationIsVisible', 0)
        all_dict.add_item('defaultConfigurationName', self.buildtype)

        # Test target
        test_dict = PbxDict()
        objects_dict.add_item(self.test_buildconf_id, test_dict, 'Build configuration list for PBXAggregateTarget "RUN_TEST"')
        test_dict.add_item('isa', 'XCConfigurationList')
        conf_arr = PbxArray()
        test_dict.add_item('buildConfigurations', conf_arr)
        for buildtype in self.buildtypes:
            conf_arr.add_item(self.test_configurations[buildtype], buildtype)
        test_dict.add_item('defaultConfigurationIsVisible', 0)
        test_dict.add_item('defaultConfigurationName', self.buildtype)

        # Regen target
        regen_dict = PbxDict()
        objects_dict.add_item(self.regen_buildconf_id, test_dict, 'Build configuration list for PBXAggregateTarget "REGENERATE"')
        regen_dict.add_item('isa', 'XCConfigurationList')
        conf_arr = PbxArray()
        regen_dict.add_item('buildConfigurations', conf_arr)
        for buildtype in self.buildtypes:
            conf_arr.add_item(self.test_configurations[buildtype], buildtype)
        regen_dict.add_item('defaultConfigurationIsVisible', 0)
        regen_dict.add_item('defaultConfigurationName', self.buildtype)

        for target_name in self.build_targets:
            t_dict = PbxDict()
            listid = self.buildconflistmap[target_name]
            objects_dict.add_item(listid, t_dict, f'Build configuration list for PBXNativeTarget "{target_name}"')
            t_dict.add_item('isa', 'XCConfigurationList')
            conf_arr = PbxArray()
            t_dict.add_item('buildConfigurations', conf_arr)
            idval = self.buildconfmap[target_name][self.buildtype]
            conf_arr.add_item(idval, self.buildtype)
            t_dict.add_item('defaultConfigurationIsVisible', 0)
            t_dict.add_item('defaultConfigurationName', self.buildtype)

        for target_name in self.custom_targets:
            t_dict = PbxDict()
            listid = self.buildconflistmap[target_name]
            objects_dict.add_item(listid, t_dict, f'Build configuration list for PBXAggregateTarget "{target_name}"')
            t_dict.add_item('isa', 'XCConfigurationList')
            conf_arr = PbxArray()
            t_dict.add_item('buildConfigurations', conf_arr)
            idval = self.buildconfmap[target_name][self.buildtype]
            conf_arr.add_item(idval, self.buildtype)
            t_dict.add_item('defaultConfigurationIsVisible', 0)
            t_dict.add_item('defaultConfigurationName', self.buildtype)

    def generate_prefix(self, pbxdict: PbxDict) -> PbxDict:
        pbxdict.add_item('archiveVersion', '1')
        pbxdict.add_item('classes', PbxDict())
        pbxdict.add_item('objectVersion', self.objversion)
        objects_dict = PbxDict()
        pbxdict.add_item('objects', objects_dict)

        return objects_dict

    def generate_suffix(self, pbxdict: PbxDict) -> None:
        pbxdict.add_item('rootObject', self.project_uid, 'Project object')
