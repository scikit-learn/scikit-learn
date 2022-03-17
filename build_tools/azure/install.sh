#!/bin/bash

set -e
set -x

# defines the get_dep and show_installed_libraries functions
source build_tools/shared.sh

UNAMESTR=`uname`

setup_ccache() {
    echo "Setting up ccache with CCACHE_DIR=${CCACHE_DIR}"
    mkdir /tmp/ccache/
    which ccache
    for name in gcc g++ cc c++ clang clang++ i686-linux-gnu-gcc i686-linux-gnu-c++ x86_64-linux-gnu-gcc x86_64-linux-gnu-c++ x86_64-apple-darwin13.4.0-clang x86_64-apple-darwin13.4.0-clang++; do
      ln -s $(which ccache) "/tmp/ccache/${name}"
    done
    export PATH="/tmp/ccache/:${PATH}"
    ccache -M 256M
}

pre_python_environment_install() {
    if [[ "$DISTRIB" == "ubuntu" ]]; then
        sudo add-apt-repository --remove ppa:ubuntu-toolchain-r/test
        sudo apt-get update
        sudo apt-get install python3-scipy python3-matplotlib \
             libatlas3-base libatlas-base-dev python3-virtualenv ccache

    elif [[ "$DISTRIB" == "debian-32" ]]; then
        apt-get update
        apt-get install -y python3-dev python3-numpy python3-scipy \
                python3-matplotlib libatlas3-base libatlas-base-dev \
                python3-virtualenv python3-pandas ccache

    elif [[ "$DISTRIB" == "conda-pypy3" ]]; then
        # need compilers
        apt-get -yq update
        apt-get -yq install build-essential

    elif [[ "$BUILD_WITH_ICC" == "true" ]]; then
        wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
        sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
        rm GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
        sudo add-apt-repository "deb https://apt.repos.intel.com/oneapi all main"
        sudo apt-get update
        sudo apt-get install intel-oneapi-compiler-dpcpp-cpp-and-cpp-classic
        source /opt/intel/oneapi/setvars.sh
    fi
}

python_environment_install_and_activate() {
    if [[ "$DISTRIB" == "conda"* ]]; then
        conda update -n base conda -y
        # pin conda-lock to latest released version (needs manual update from time to time)
        conda install -c conda-forge conda-lock==1.0.3 -y
        conda-lock install --name $VIRTUALENV $LOCK_FILE
        source activate $VIRTUALENV

    elif [[ "$DISTRIB" == "ubuntu" || "$DISTRIB" == "debian-32" ]]; then
        python3 -m virtualenv --system-site-packages --python=python3 $VIRTUALENV
        source $VIRTUALENV/bin/activate
        pip install -r "${LOCK_FILE}"
    fi

    if [[ "$DISTRIB" == "conda-pip-scipy-dev" ]]; then
        echo "Installing development dependency wheels"
        dev_anaconda_url=https://pypi.anaconda.org/scipy-wheels-nightly/simple
        pip install --pre --upgrade --timeout=60 --extra-index $dev_anaconda_url numpy pandas scipy
        echo "Installing joblib master"
        pip install https://github.com/joblib/joblib/archive/master.zip
        echo "Installing pillow master"
        pip install https://github.com/python-pillow/Pillow/archive/main.zip
    fi
}

scikit_learn_install() {
    setup_ccache
    show_installed_libraries

    # Set parallelism to 3 to overlap IO bound tasks with CPU bound tasks on CI
    # workers with 2 cores when building the compiled extensions of scikit-learn.
    export SKLEARN_BUILD_PARALLEL=3

    if [[ "$UNAMESTR" == "Darwin" && "$SKLEARN_TEST_NO_OPENMP" == "true" ]]; then
        # Without openmp, we use the system clang. Here we use /usr/bin/ar
        # instead because llvm-ar errors
        export AR=/usr/bin/ar
    fi

    if [[ "$UNAMESTR" == "Linux" ]]; then
        # FIXME: temporary fix to link against system libraries on linux
        # https://github.com/scikit-learn/scikit-learn/issues/20640
        export LDFLAGS="$LDFLAGS -Wl,--sysroot=/"
    fi

    if [[ "$BUILD_WITH_ICC" == "true" ]]; then
        # The "build_clib" command is implicitly used to build "libsvm-skl".
        # To compile with a different compiler, we also need to specify the
        # compiler for this command
        python setup.py build_ext --compiler=intelem -i build_clib --compiler=intelem
    fi

    # TODO use a specific variable for this rather than using a particular build ...
    if [[ "$DISTRIB" == "conda-pip-latest" ]]; then
        # Check that pip can automatically build scikit-learn with the build
        # dependencies specified in pyproject.toml using an isolated build
        # environment:
        pip install --verbose --editable .
    else
        # Use the pre-installed build dependencies and build directly in the
        # current environment.
        python setup.py develop
    fi

    ccache -s
}

main() {
    pre_python_environment_install
    python_environment_install_and_activate
    scikit_learn_install
}

main
