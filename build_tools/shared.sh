get_dep() {
    package="$1"
    version="$2"
    if [[ "$version" == "none" ]]; then
        # do not install with none
        echo
    elif [[ "${version%%[^0-9.]*}" ]]; then
        # version number is explicitly passed
        echo "$package==$version"
    elif [[ "$version" == "latest" ]]; then
        # use latest
        echo "$package"
    elif [[ "$version" == "min" ]]; then
        echo "$package==$(python sklearn/_min_dependencies.py $package)"
    fi
}

show_installed_libraries(){
    # Use `conda list` when a usable conda CLI is available (e.g. the GPU build
    # which still relies on miniforge): it shows more info than `pip list`, such
    # as whether OpenBLAS or MKL is installed and their versions. pixi-managed
    # environments do not ship the `conda` executable, so we fall back to
    # `pip list` there.
    if command -v conda > /dev/null 2>&1 && [[ -n "$CONDA_PREFIX" ]]; then
        conda list
    else
        python -m pip list
    fi
}

show_cpu_info() {
    echo "========== CPU information =========="
    if [ -x "$(command -v lscpu)" ] ; then
        lscpu
    elif [ -x "$(command -v system_profiler)" ] ; then
        system_profiler SPHardwareDataType
    elif [ -x "$(command -v powershell)" ] ; then
        powershell -c '$cpu = Get-WmiObject -Class Win32_Processor
            Write-Host "CPU Model: $($cpu.Name)"
            Write-Host "Architecture: $($cpu.Architecture)"
            Write-Host "Physical Cores: $($cpu.NumberOfCores)"
            Write-Host "Logical Processors: $($cpu.NumberOfLogicalProcessors)"
        '
    else
        echo "Could not inspect CPU architecture."
    fi
    echo "====================================="
}

activate_environment() {
    # pixi-managed environments are already activated when the script is run
    # through `pixi run`, so nothing is needed there. Only the apt-based builds
    # (ubuntu_atlas and the debian-32 docker image) rely on a plain virtualenv
    # that we activate here.
    if [[ "$DISTRIB" == "ubuntu" || "$DISTRIB" == "debian-32" || "$DISTRIB" == "scipy-dev" ]]; then
        source $VIRTUALENV/bin/activate
    fi
}
