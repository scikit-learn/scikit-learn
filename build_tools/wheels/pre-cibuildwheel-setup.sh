set -xe

# OpenMP is not present on macOS by default
if [[ $(uname) == "Darwin" ]]; then
    # Make sure to use a libomp version binary compatible with the oldest
    # supported version of the macos SDK as libomp will be vendored into the
    # scikit-learn wheels for macos.

    if [[ $(uname -m) == "arm64" ]]; then
        # SciPy requires 12.0 on arm to prevent kernel panics
        # https://github.com/scipy/scipy/issues/14688
        # We use the same deployment target to match SciPy.
        echo MACOSX_DEPLOYMENT_TARGET=12.0 >> "$GITHUB_ENV"
        OPENMP_URL="https://anaconda.org/conda-forge/llvm-openmp/11.1.0/download/osx-arm64/llvm-openmp-11.1.0-hf3c4609_1.tar.bz2"
    else
        echo MACOSX_DEPLOYMENT_TARGET=10.9 >> "$GITHUB_ENV"
        OPENMP_URL="https://anaconda.org/conda-forge/llvm-openmp/11.1.0/download/osx-64/llvm-openmp-11.1.0-hda6cdc1_1.tar.bz2"
    fi

    conda create -n build $OPENMP_URL
    PREFIX="$HOME/miniconda3/envs/build"

    echo CC=/usr/bin/clang >> "$GITHUB_ENV"
    echo CXX=/usr/bin/clang++ >> "$GITHUB_ENV"
    echo "CPPFLAGS=$CPPFLAGS -Xpreprocessor -fopenmp" >> "$GITHUB_ENV"
    echo "CFLAGS=$CFLAGS -I$PREFIX/include" >> "$GITHUB_ENV"
    echo "CXXFLAGS=$CXXFLAGS -I$PREFIX/include" >> "$GITHUB_ENV"
    echo "LDFLAGS=$LDFLAGS -Wl,-rpath,$PREFIX/lib -L$PREFIX/lib -lomp" >> "$GITHUB_ENV"
fi
