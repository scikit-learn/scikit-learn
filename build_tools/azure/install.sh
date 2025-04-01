#!/bin/bash

set -e
set -x

# defines the get_dep and show_installed_libraries functions
source build_tools/shared.sh

UNAMESTR=`uname`
CCACHE_LINKS_DIR="/tmp/ccache"

setup_ccache() {
    CCACHE_BIN=`which ccache || echo ""`
    if [[ "${CCACHE_BIN}" == "" ]]; then
        echo "ccache not found, skipping..."
    elif [[ -d "${CCACHE_LINKS_DIR}" ]]; then
        echo "ccache already configured, skipping..."
    else
        echo "Setting up ccache with CCACHE_DIR=${CCACHE_DIR}"
        mkdir ${CCACHE_LINKS_DIR}
        which ccache
        for name in gcc g++ cc c++ clang clang++ i686-linux-gnu-gcc i686-linux-gnu-c++ x86_64-linux-gnu-gcc x86_64-linux-gnu-c++ x86_64-apple-darwin13.4.0-clang x86_64-apple-darwin13.4.0-clang++; do
        ln -s ${CCACHE_BIN} "${CCACHE_LINKS_DIR}/${name}"
        done
        export PATH="${CCACHE_LINKS_DIR}:${PATH}"
        ccache -M 256M

        # Zeroing statistics so that ccache statistics are shown only for this build
        ccache -z
    fi
}

pre_python_environment_install() {
    if [[ "$DISTRIB" == "ubuntu" ]]; then
        sudo apt-get update
        sudo apt-get install python3-scipy python3-matplotlib \
             libatlas3-base libatlas-base-dev python3-virtualenv ccache

    elif [[ "$DISTRIB" == "debian-32" ]]; then
        apt-get update
        apt-get install -y python3-dev python3-numpy python3-scipy \
                python3-matplotlib libopenblas-dev \
                python3-virtualenv python3-pandas ccache git
    fi
}

check_packages_dev_version() {
    for package in $@; do
        package_version=$(python -c "import $package; print($package.__version__)")
        if [[ $package_version =~ "^[.0-9]+$" ]]; then
            echo "$package is not a development version: $package_version"
            exit 1
        fi
    done
}

python_environment_install_and_activate() {
    if [[ "$DISTRIB" == "conda"* ]]; then
        create_conda_environment_from_lock_file $VIRTUALENV $LOCK_FILE
        activate_environment

    elif [[ "$DISTRIB" == "ubuntu" || "$DISTRIB" == "debian-32" ]]; then
        python3 -m virtualenv --system-site-packages --python=python3 $VIRTUALENV
        activate_environment
        pip install -r "${LOCK_FILE}"

    fi

    # Install additional packages on top of the lock-file in specific cases
    if [[ "$DISTRIB" == "conda-free-threaded" ]]; then
        # TODO: we install scipy with pip. When there is a conda-forge package,
        # we can update build_tools/update_environments_and_lock_files.py and
        # remove the line below
        pip install scipy --only-binary :all:
        # TODO: we install cython 3.1 alpha from pip. When there is a conda-forge package,
        # we can update build_tools/update_environments_and_lock_files.py and
        # remove the line below
        pip install --pre cython --only-binary :all:

    elif [[ "$DISTRIB" == "conda-pip-scipy-dev" ]]; then
        echo "Installing development dependency wheels"
        dev_anaconda_url=https://pypi.anaconda.org/scientific-python-nightly-wheels/simple
        dev_packages="numpy scipy pandas Cython"
        pip install --pre --upgrade --timeout=60 --extra-index $dev_anaconda_url $dev_packages --only-binary :all:

        check_packages_dev_version $dev_packages

        echo "Installing joblib from latest sources"
        pip install https://github.com/joblib/joblib/archive/master.zip
        echo "Installing pillow from latest sources"
        pip install https://github.com/python-pillow/Pillow/archive/main.zip
    fi
}

scikit_learn_install() {
    setup_ccache
    show_installed_libraries

    if [[ "$UNAMESTR" == "Darwin" && "$SKLEARN_TEST_NO_OPENMP" == "true" ]]; then
        # Without openmp, we use the system clang. Here we use /usr/bin/ar
        # instead because llvm-ar errors
        export AR=/usr/bin/ar
        # Make sure omp.h is not present in the conda environment, so that
        # using an unprotected "cimport openmp" will make this build fail. At
        # the time of writing (2023-01-13), on OSX, blas (mkl or openblas)
        # brings in openmp so that you end up having the omp.h include inside
        # the conda environment.
        find $CONDA_PREFIX -name omp.h -delete -print
        # meson >= 1.5 detects OpenMP installed with brew and OpenMP may be installed
        # with brew in CI runner. OpenMP was installed with brew in macOS-12 CI
        # runners which doesn't seem to be the case in macOS-13 runners anymore,
        # but we keep the next line just to be safe ...
        brew uninstall --ignore-dependencies --force libomp
    fi

    if [[ "$UNAMESTR" == "Linux" ]]; then
        # FIXME: temporary fix to link against system libraries on linux
        # https://github.com/scikit-learn/scikit-learn/issues/20640
        export LDFLAGS="$LDFLAGS -Wl,--sysroot=/"
    fi

    if [[ "$PIP_BUILD_ISOLATION" == "true" ]]; then
        # Check that pip can automatically build scikit-learn with the build
        # dependencies specified in pyproject.toml using an isolated build
        # environment:
        pip install --verbose .
    else
        if [[ "$UNAMESTR" == "MINGW64"* ]]; then
           # Needed on Windows CI to compile with Visual Studio compiler
           # otherwise Meson detects a MINGW64 platform and use MINGW64
           # toolchain
           ADDITIONAL_PIP_OPTIONS='-Csetup-args=--vsenv'
        fi
        # Use the pre-installed build dependencies and build directly in the
        # current environment.
        pip install --verbose --no-build-isolation --editable . $ADDITIONAL_PIP_OPTIONS
    fi

    ccache -s || echo "ccache not installed, skipping ccache statistics"
}

main() {
    pre_python_environment_install
    python_environment_install_and_activate
    scikit_learn_install
}

main
