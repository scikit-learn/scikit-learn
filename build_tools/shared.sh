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
