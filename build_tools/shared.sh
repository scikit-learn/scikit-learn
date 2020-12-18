get_dep() {
    package="$1"
    version="$2"
    if [[ "$version" == "none" ]]; then
        # do not install with none
        echo
    elif [[ "${version%%[^0-9.]*}" ]]; then
        # version number is explicity passed
        echo "$package==$version"
    elif [[ "$version" == "latest" ]]; then
        # use latest
        echo "$package"
    elif [[ "$version" == "min" ]]; then
        echo "$package==$(python sklearn/_min_dependencies.py $package)"
    fi
}
