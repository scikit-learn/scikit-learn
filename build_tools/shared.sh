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
    # use conda list when inside a conda environment. conda list shows more
    # info than pip list, e.g. whether OpenBLAS or MKL is installed as well as
    # the version of OpenBLAS or MKL
    if [[ -n "$CONDA_PREFIX" ]]; then
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
    if [[ "$DISTRIB" =~ ^conda.* ]]; then
        source activate $VIRTUALENV
    elif [[ "$DISTRIB" == "ubuntu" || "$DISTRIB" == "debian-32" ]]; then
        source $VIRTUALENV/bin/activate
    fi
}

create_conda_environment_from_lock_file() {
    ENV_NAME=$1
    LOCK_FILE=$2
    # Because we are using lock-files with the "explicit" format, conda can
    # install them directly, provided the lock-file does not contain pip solved
    # packages. For more details, see
    # https://conda.github.io/conda-lock/output/#explicit-lockfile
    lock_file_has_pip_packages=$(grep -q files.pythonhosted.org $LOCK_FILE && echo "true" || echo "false")
    if [[ "$lock_file_has_pip_packages" == "false" ]]; then
        conda create --quiet --name $ENV_NAME --file $LOCK_FILE
    else
        python -m pip install "$(get_dep conda-lock min)"
        conda-lock install --log-level WARNING --name $ENV_NAME $LOCK_FILE
    fi
}
