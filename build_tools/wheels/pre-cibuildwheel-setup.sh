set -xe

# OpenMP is not present on macOS by default
if [[ $(uname) == "Darwin" ]]; then
    # Make sure to use a libomp version binary compatible with the oldest
    # supported version of the macos SDK as libomp will be vendored into the
    # scikit-learn wheels for macos.
    # MACOSX_DEPLOYMENT_TARGET is defined in pyproject.toml
    if [[ $(uname -m) == "arm64" ]]; then
        OPENMP_URL="https://anaconda.org/conda-forge/llvm-openmp/11.1.0/download/osx-arm64/llvm-openmp-11.1.0-hf3c4609_1.tar.bz2"
    else
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
