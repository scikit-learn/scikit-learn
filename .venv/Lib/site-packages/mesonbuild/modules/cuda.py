# SPDX-License-Identifier: Apache-2.0
# Copyright 2017-2025 The Meson development team

from __future__ import annotations

import dataclasses
import re
import typing as T

from ..mesonlib import listify, version_compare
from ..compilers.cuda import CudaCompiler
from ..interpreter.type_checking import NoneType

from . import NewExtensionModule, ModuleInfo

from ..interpreterbase import (
    ContainerTypeInfo, InvalidArguments, KwargInfo, noKwargs, typed_kwargs, typed_pos_args,
)

if T.TYPE_CHECKING:
    from typing_extensions import TypedDict

    from . import ModuleState
    from ..interpreter import Interpreter
    from ..interpreterbase import TYPE_var

    class ArchFlagsKwargs(TypedDict):
        detected: T.Optional[T.List[str]]

    AutoArch = T.Union[str, T.List[str]]


DETECTED_KW: KwargInfo[T.Union[None, T.List[str]]] = KwargInfo('detected', (ContainerTypeInfo(list, str), NoneType), listify=True)


@dataclasses.dataclass
class _CudaVersion:

    meson: str
    windows: str
    linux: str

    def compare(self, version: str, machine: str) -> T.Optional[str]:
        if version_compare(version, f'>={self.meson}'):
            return self.windows if machine == 'windows' else self.linux
        return None


# Copied from: https://docs.nvidia.com/cuda/cuda-toolkit-release-notes/index.html#id7
_DRIVER_TABLE_VERSION: T.List[_CudaVersion] = [
    _CudaVersion('13.0.2', 'unknown', '580.95.05'),
    _CudaVersion('13.0.1', 'unknown', '580.82.07'),
    _CudaVersion('13.0.0', 'unknown', '580.65.06'),
    _CudaVersion('12.9.1', '576.57', '575.57.08'),
    _CudaVersion('12.9.0', '576.02', '575.51.03'),
    _CudaVersion('12.8.1', '572.61', '570.124.06'),
    _CudaVersion('12.8.0', '570.65', '570.26'),
    _CudaVersion('12.6.3', '561.17', '560.35.05'),
    _CudaVersion('12.6.2', '560.94', '560.35.03'),
    _CudaVersion('12.6.1', '560.94', '560.35.03'),
    _CudaVersion('12.6.0', '560.76', '560.28.03'),
    _CudaVersion('12.5.1', '555.85', '555.42.06'),
    _CudaVersion('12.5.0', '555.85', '555.42.02'),
    _CudaVersion('12.4.1', '551.78', '550.54.15'),
    _CudaVersion('12.4.0', '551.61', '550.54.14'),
    _CudaVersion('12.3.1', '546.12', '545.23.08'),
    _CudaVersion('12.3.0', '545.84', '545.23.06'),
    _CudaVersion('12.2.2', '537.13', '535.104.05'),
    _CudaVersion('12.2.1', '536.67', '535.86.09'),
    _CudaVersion('12.2.0', '536.25', '535.54.03'),
    _CudaVersion('12.1.1', '531.14', '530.30.02'),
    _CudaVersion('12.1.0', '531.14', '530.30.02'),
    _CudaVersion('12.0.1', '528.33', '525.85.11'),
    _CudaVersion('12.0.0', '527.41', '525.60.13'),
    _CudaVersion('11.8.0', '522.06', '520.61.05'),
    _CudaVersion('11.7.1', '516.31', '515.48.07'),
    _CudaVersion('11.7.0', '516.01', '515.43.04'),
    _CudaVersion('11.6.1', '511.65', '510.47.03'),  # 11.6.2 is identical
    _CudaVersion('11.6.0', '511.23', '510.39.01'),
    _CudaVersion('11.5.1', '496.13', '495.29.05'),  # 11.5.2 is identical
    _CudaVersion('11.5.0', '496.04', '495.29.05'),
    _CudaVersion('11.4.3', '472.50', '470.82.01'),  # 11.4.4 is identical
    _CudaVersion('11.4.1', '471.41', '470.57.02'),  # 11.4.2 is identical
    _CudaVersion('11.4.0', '471.11', '470.42.01'),
    _CudaVersion('11.3.0', '465.89', '465.19.01'),  # 11.3.1 is identical
    _CudaVersion('11.2.2', '461.33', '460.32.03'),
    _CudaVersion('11.2.1', '461.09', '460.32.03'),
    _CudaVersion('11.2.0', '460.82', '460.27.03'),
    _CudaVersion('11.1.1', '456.81', '455.32'),
    _CudaVersion('11.1.0', '456.38', '455.23'),
    _CudaVersion('11.0.3', '451.82', '450.51.06'),  # 11.0.3.1 is identical
    _CudaVersion('11.0.2', '451.48', '450.51.05'),
    _CudaVersion('11.0.1', '451.22', '450.36.06'),
    _CudaVersion('10.2.89', '441.22', '440.33'),
    _CudaVersion('10.1.105', '418.96', '418.39'),
    _CudaVersion('10.0.130', '411.31', '410.48'),
    _CudaVersion('9.2.148', '398.26', '396.37'),
    _CudaVersion('9.2.88', '397.44', '396.26'),
    _CudaVersion('9.1.85', '391.29', '390.46'),
    _CudaVersion('9.0.76', '385.54', '384.81'),
    _CudaVersion('8.0.61', '376.51', '375.26'),
    _CudaVersion('8.0.44', '369.30', '367.48'),
    _CudaVersion('7.5.16', '353.66', '352.31'),
    _CudaVersion('7.0.28', '347.62', '346.46'),
]

class CudaModule(NewExtensionModule):

    INFO = ModuleInfo('CUDA', '0.50.0', unstable=True)

    def __init__(self, interp: Interpreter):
        super().__init__()
        self.methods.update({
            "min_driver_version": self.min_driver_version,
            "nvcc_arch_flags":    self.nvcc_arch_flags,
            "nvcc_arch_readable": self.nvcc_arch_readable,
        })

    @noKwargs
    def min_driver_version(self, state: 'ModuleState',
                           args: T.List[TYPE_var],
                           kwargs: T.Dict[str, T.Any]) -> str:
        argerror = InvalidArguments('min_driver_version must have exactly one positional argument: ' +
                                    'a CUDA Toolkit version string. Beware that, since CUDA 11.0, ' +
                                    'the CUDA Toolkit\'s components (including NVCC) are versioned ' +
                                    'independently from each other (and the CUDA Toolkit as a whole).')
        if len(args) != 1 or not isinstance(args[0], str):
            raise argerror

        cuda_version = args[0]

        for d in _DRIVER_TABLE_VERSION:
            driver_version = d.compare(cuda_version, state.environment.machines.host.system)
            if driver_version is not None:
                return driver_version
        return 'unknown'

    @typed_pos_args('cuda.nvcc_arch_flags', (str, CudaCompiler), varargs=str)
    @typed_kwargs('cuda.nvcc_arch_flags', DETECTED_KW)
    def nvcc_arch_flags(self, state: 'ModuleState',
                        args: T.Tuple[T.Union[CudaCompiler, str], T.List[str]],
                        kwargs: ArchFlagsKwargs) -> T.List[str]:
        nvcc_arch_args = self._validate_nvcc_arch_args(args, kwargs)
        ret = self._nvcc_arch_flags(*nvcc_arch_args)[0]
        return ret

    @typed_pos_args('cuda.nvcc_arch_readable', (str, CudaCompiler), varargs=str)
    @typed_kwargs('cuda.nvcc_arch_readable', DETECTED_KW)
    def nvcc_arch_readable(self, state: 'ModuleState',
                           args: T.Tuple[T.Union[CudaCompiler, str], T.List[str]],
                           kwargs: ArchFlagsKwargs) -> T.List[str]:
        nvcc_arch_args = self._validate_nvcc_arch_args(args, kwargs)
        ret = self._nvcc_arch_flags(*nvcc_arch_args)[1]
        return ret

    @staticmethod
    def _break_arch_string(s: str) -> T.List[str]:
        s = re.sub('[ \t\r\n,;]+', ';', s)
        return s.strip(';').split(';')

    @staticmethod
    def _detected_cc_from_compiler(c: T.Union[str, CudaCompiler]) -> T.List[str]:
        if isinstance(c, CudaCompiler):
            return [c.detected_cc]
        return []

    def _validate_nvcc_arch_args(self, args: T.Tuple[T.Union[str, CudaCompiler], T.List[str]],
                                 kwargs: ArchFlagsKwargs) -> T.Tuple[str, AutoArch, T.List[str]]:

        compiler = args[0]
        if isinstance(compiler, CudaCompiler):
            cuda_version = compiler.version
        else:
            cuda_version = compiler

        arch_list: AutoArch = args[1]
        arch_list = listify([self._break_arch_string(a) for a in arch_list])
        if len(arch_list) > 1 and not set(arch_list).isdisjoint({'All', 'Common', 'Auto'}):
            raise InvalidArguments('''The special architectures 'All', 'Common' and 'Auto' must appear alone, as a positional argument!''')
        arch_list = arch_list[0] if len(arch_list) == 1 else arch_list

        detected = kwargs['detected'] if kwargs['detected'] is not None else self._detected_cc_from_compiler(compiler)
        detected = [x for a in detected for x in self._break_arch_string(a)]
        if not set(detected).isdisjoint({'All', 'Common', 'Auto'}):
            raise InvalidArguments('''The special architectures 'All', 'Common' and 'Auto' must appear alone, as a positional argument!''')

        return cuda_version, arch_list, detected

    def _filter_cuda_arch_list(self, cuda_arch_list: T.List[str], lo: str, hi: T.Optional[str], saturate: str) -> T.List[str]:
        """
        Filter CUDA arch list (no codenames) for >= low and < hi architecture
        bounds, and deduplicate.
        Architectures >= hi are replaced with saturate.
        """

        filtered_cuda_arch_list = []
        for arch in cuda_arch_list:
            if arch:
                if lo and version_compare(arch, '<' + lo):
                    continue
                if hi and version_compare(arch, '>=' + hi):
                    arch = saturate
                if arch not in filtered_cuda_arch_list:
                    filtered_cuda_arch_list.append(arch)
        return filtered_cuda_arch_list

    def _nvcc_arch_flags(self, cuda_version: str, cuda_arch_list: AutoArch, detected: T.List[str]) -> T.Tuple[T.List[str], T.List[str]]:
        """
        Using the CUDA Toolkit version and the target architectures, compute
        the NVCC architecture flags.
        """

        # Replicates much of the logic of
        #     https://github.com/Kitware/CMake/blob/master/Modules/FindCUDA/select_compute_arch.cmake
        # except that a bug with cuda_arch_list="All" is worked around by
        # tracking both lower and upper limits on GPU architectures.

        cuda_known_gpu_architectures   = []  # noqa: E221
        cuda_common_gpu_architectures  = ['3.0', '3.5', '5.0']           # noqa: E221
        cuda_hi_limit_gpu_architecture = None                            # noqa: E221
        cuda_lo_limit_gpu_architecture = '2.0'                           # noqa: E221
        cuda_all_gpu_architectures     = ['3.0', '3.2', '3.5', '5.0']    # noqa: E221

        # Fermi and Kepler support have been dropped since 12.0
        if version_compare(cuda_version, '<12.0'):
            cuda_known_gpu_architectures.extend(['Fermi', 'Kepler'])

        # Everything older than Turing is dropped by 13.0
        if version_compare(cuda_version, '<13.0'):
            cuda_known_gpu_architectures.append('Maxwell')

            if version_compare(cuda_version, '<7.0'):
                cuda_hi_limit_gpu_architecture = '5.2'

            if version_compare(cuda_version, '>=7.0'):
                cuda_known_gpu_architectures  += ['Kepler+Tegra', 'Kepler+Tesla', 'Maxwell+Tegra']  # noqa: E221
                cuda_common_gpu_architectures += ['5.2']                                            # noqa: E221

                if version_compare(cuda_version, '<8.0'):
                    cuda_common_gpu_architectures += ['5.2+PTX']  # noqa: E221
                    cuda_hi_limit_gpu_architecture = '6.0'        # noqa: E221

            if version_compare(cuda_version, '>=8.0'):
                cuda_known_gpu_architectures  += ['Pascal', 'Pascal+Tegra']  # noqa: E221
                cuda_common_gpu_architectures += ['6.0', '6.1']              # noqa: E221
                cuda_all_gpu_architectures    += ['6.0', '6.1', '6.2']       # noqa: E221

                if version_compare(cuda_version, '<9.0'):
                    cuda_common_gpu_architectures += ['6.1+PTX']  # noqa: E221
                    cuda_hi_limit_gpu_architecture = '7.0'        # noqa: E221

            if version_compare(cuda_version, '>=9.0'):
                cuda_known_gpu_architectures  += ['Volta', 'Xavier'] # noqa: E221
                cuda_common_gpu_architectures += ['7.0']             # noqa: E221
                cuda_all_gpu_architectures    += ['7.0', '7.2']      # noqa: E221
                # https://docs.nvidia.com/cuda/archive/9.0/cuda-toolkit-release-notes/index.html#unsupported-features
                cuda_lo_limit_gpu_architecture = '3.0'               # noqa: E221

                if version_compare(cuda_version, '<10.0'):
                    cuda_common_gpu_architectures += ['7.2+PTX']  # noqa: E221
                    cuda_hi_limit_gpu_architecture = '8.0'        # noqa: E221

        if version_compare(cuda_version, '>=10.0'):
            cuda_known_gpu_architectures  += ['Turing'] # noqa: E221
            cuda_common_gpu_architectures += ['7.5']    # noqa: E221
            cuda_all_gpu_architectures    += ['7.5']    # noqa: E221

            if version_compare(cuda_version, '<11.0'):
                cuda_common_gpu_architectures += ['7.5+PTX']  # noqa: E221
                cuda_hi_limit_gpu_architecture = '8.0'        # noqa: E221

        # need to account for the fact that Ampere is commonly assumed to include
        # SM8.0 and SM8.6 even though CUDA 11.0 doesn't support SM8.6
        cuda_ampere_bin = ['8.0']
        cuda_ampere_ptx = ['8.0']
        if version_compare(cuda_version, '>=11.0'):
            cuda_known_gpu_architectures  += ['Ampere'] # noqa: E221
            cuda_common_gpu_architectures += ['8.0']    # noqa: E221
            cuda_all_gpu_architectures    += ['8.0']    # noqa: E221
            # https://docs.nvidia.com/cuda/archive/11.0/cuda-toolkit-release-notes/index.html#deprecated-features
            cuda_lo_limit_gpu_architecture = '3.5'      # noqa: E221

            if version_compare(cuda_version, '<11.1'):
                cuda_common_gpu_architectures += ['8.0+PTX']  # noqa: E221
                cuda_hi_limit_gpu_architecture = '8.6'        # noqa: E221

        if version_compare(cuda_version, '>=11.1'):
            cuda_ampere_bin += ['8.6'] # noqa: E221
            cuda_ampere_ptx  = ['8.6'] # noqa: E221

            cuda_common_gpu_architectures += ['8.6']             # noqa: E221
            cuda_all_gpu_architectures    += ['8.6']             # noqa: E221

            if version_compare(cuda_version, '<11.8'):
                cuda_common_gpu_architectures += ['8.6+PTX']  # noqa: E221
                cuda_hi_limit_gpu_architecture = '8.7'        # noqa: E221

        if version_compare(cuda_version, '>=11.8'):
            cuda_known_gpu_architectures  += ['Orin', 'Lovelace', 'Hopper']  # noqa: E221
            cuda_common_gpu_architectures += ['8.9', '9.0', '9.0+PTX']       # noqa: E221
            cuda_all_gpu_architectures    += ['8.7', '8.9', '9.0']           # noqa: E221

            if version_compare(cuda_version, '<12'):
                cuda_hi_limit_gpu_architecture = '9.1'        # noqa: E221

        if version_compare(cuda_version, '>=12.0'):
            # https://docs.nvidia.com/cuda/cuda-toolkit-release-notes/index.html#deprecated-features (Current)
            # https://docs.nvidia.com/cuda/archive/12.0/cuda-toolkit-release-notes/index.html#deprecated-features (Eventual?)
            cuda_lo_limit_gpu_architecture = '5.0'            # noqa: E221

            if version_compare(cuda_version, '<13'):
                cuda_hi_limit_gpu_architecture = '10.0'       # noqa: E221

        if version_compare(cuda_version, '>=12.8'):
            cuda_known_gpu_architectures.append('Blackwell')
            cuda_common_gpu_architectures.extend(['10.0', '12.0'])
            cuda_all_gpu_architectures.extend(['10.0', '12.0'])

            if version_compare(cuda_version, '<13'):
                # Yes, 12.8 and 12.9 support 10.1, but 13.0 doesn't
                cuda_common_gpu_architectures.append('10.1')
                cuda_all_gpu_architectures.append('10.1')
                cuda_hi_limit_gpu_architecture = '12.1'

        if version_compare(cuda_version, '>=12.9'):
            cuda_common_gpu_architectures.extend(['10.3', '12.1'])
            cuda_all_gpu_architectures.extend(['10.3', '12.1'])

        if version_compare(cuda_version, '>=13.0'):
            cuda_common_gpu_architectures.append('11.0')
            cuda_all_gpu_architectures.append('11.0')
            cuda_lo_limit_gpu_architecture = '7.5'

            if version_compare(cuda_version, '<14'):
                cuda_hi_limit_gpu_architecture = '12.1'

        if not cuda_arch_list:
            cuda_arch_list = 'Auto'

        if   cuda_arch_list == 'All':     # noqa: E271
            cuda_arch_list = cuda_known_gpu_architectures
        elif cuda_arch_list == 'Common':  # noqa: E271
            cuda_arch_list = cuda_common_gpu_architectures
        elif cuda_arch_list == 'Auto':    # noqa: E271
            if detected:
                if isinstance(detected, list):
                    cuda_arch_list = detected
                else:
                    cuda_arch_list = self._break_arch_string(detected)
                cuda_arch_list = self._filter_cuda_arch_list(cuda_arch_list,
                                                             cuda_lo_limit_gpu_architecture,
                                                             cuda_hi_limit_gpu_architecture,
                                                             cuda_common_gpu_architectures[-1])
            else:
                cuda_arch_list = cuda_common_gpu_architectures
        elif isinstance(cuda_arch_list, str):
            cuda_arch_list = self._break_arch_string(cuda_arch_list)

        cuda_arch_list = sorted(x for x in set(cuda_arch_list) if x)

        cuda_arch_bin: T.List[str] = []
        cuda_arch_ptx: T.List[str] = []
        for arch_name in cuda_arch_list:
            arch_bin: T.Optional[T.List[str]]
            arch_ptx: T.Optional[T.List[str]]
            add_ptx = arch_name.endswith('+PTX')
            if add_ptx:
                arch_name = arch_name[:-len('+PTX')]

            if re.fullmatch('[0-9]+\\.[0-9](\\([0-9]+\\.[0-9]\\))?', arch_name):
                arch_bin, arch_ptx = [arch_name], [arch_name]
            else:
                arch_bin, arch_ptx = {
                    'Fermi':         (['2.0', '2.1(2.0)'], []),
                    'Kepler+Tegra':  (['3.2'],             []),
                    'Kepler+Tesla':  (['3.7'],             []),
                    'Kepler':        (['3.0', '3.5'],      ['3.5']),
                    'Maxwell+Tegra': (['5.3'],             []),
                    'Maxwell':       (['5.0', '5.2'],      ['5.2']),
                    'Pascal':        (['6.0', '6.1'],      ['6.1']),
                    'Pascal+Tegra':  (['6.2'],             []),
                    'Volta':         (['7.0'],             ['7.0']),
                    'Xavier':        (['7.2'],             []),
                    'Turing':        (['7.5'],             ['7.5']),
                    'Ampere':        (cuda_ampere_bin,     cuda_ampere_ptx),
                    'Orin':          (['8.7'],             []),
                    'Lovelace':      (['8.9'],             ['8.9']),
                    'Hopper':        (['9.0'],             ['9.0']),
                    'Blackwell':     (['10.0'],            ['10.0']),
                }.get(arch_name, (None, None))

            if arch_bin is None:
                raise InvalidArguments(f'Unknown CUDA Architecture Name {arch_name}!')

            cuda_arch_bin += arch_bin

            if add_ptx:
                if not arch_ptx:
                    arch_ptx = arch_bin
                cuda_arch_ptx += arch_ptx

        cuda_arch_bin = sorted(set(cuda_arch_bin))
        cuda_arch_ptx = sorted(set(cuda_arch_ptx))

        nvcc_flags = []
        nvcc_archs_readable = []

        for arch in cuda_arch_bin:
            arch, codev = re.fullmatch(
                '([0-9]+\\.[0-9])(?:\\(([0-9]+\\.[0-9])\\))?', arch).groups()

            if version_compare(arch, '<' + cuda_lo_limit_gpu_architecture):
                continue
            if cuda_hi_limit_gpu_architecture and version_compare(arch, '>=' + cuda_hi_limit_gpu_architecture):
                continue

            if codev:
                arch = arch.replace('.', '')
                codev = codev.replace('.', '')
                nvcc_flags += ['-gencode', 'arch=compute_' + codev + ',code=sm_' + arch]
                nvcc_archs_readable += ['sm_' + arch]
            else:
                arch = arch.replace('.', '')
                nvcc_flags += ['-gencode', 'arch=compute_' + arch + ',code=sm_' + arch]
                nvcc_archs_readable += ['sm_' + arch]

        for arch in cuda_arch_ptx:
            arch, codev = re.fullmatch(
                '([0-9]+\\.[0-9])(?:\\(([0-9]+\\.[0-9])\\))?', arch).groups()

            if codev:
                arch = codev

            if version_compare(arch, '<' + cuda_lo_limit_gpu_architecture):
                continue
            if cuda_hi_limit_gpu_architecture and version_compare(arch, '>=' + cuda_hi_limit_gpu_architecture):
                continue

            arch = arch.replace('.', '')
            nvcc_flags += ['-gencode', 'arch=compute_' + arch + ',code=compute_' + arch]
            nvcc_archs_readable += ['compute_' + arch]

        return nvcc_flags, nvcc_archs_readable

def initialize(interp: Interpreter) -> CudaModule:
    return CudaModule(interp)
