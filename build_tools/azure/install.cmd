@rem https://github.com/numba/numba/blob/master/buildscripts/incremental/setup_conda_environment.cmd
@rem The cmd /C hack circumvents a regression where conda installs a conda.bat
@rem script in non-root environments.
set CONDA_INSTALL=cmd /C conda install -q -y
set PIP_INSTALL=pip install -q

@echo on

IF "%PYTHON_ARCH%"=="64" (
    @rem Deactivate any environment
    call deactivate
    @rem Clean up any left-over from a previous build
    conda remove --all -q -y -n %VIRTUALENV%
    conda create -n %VIRTUALENV% -q -y python=%PYTHON_VERSION% numpy scipy cython matplotlib wheel pillow joblib

    call activate %VIRTUALENV%

    pip install threadpoolctl

    IF "%PYTEST_VERSION%"=="*" (
        pip install pytest
    ) else (
        pip install pytest==%PYTEST_VERSION%
    )
) else (
    pip install numpy scipy cython pytest wheel pillow joblib threadpoolctl
)

IF "%PYTEST_XDIST%" == "true" (
    pip install pytest-xdist
)

if "%COVERAGE%" == "true" (
    pip install coverage codecov
    REM Use master pytest-cov until 2.10.1 is released that fixes a bug with
    REM  pytest-xdist 2.0.0
    pip install git+https://github.com/pytest-dev/pytest-cov
)
python --version
pip --version

@rem Install the build and runtime dependencies of the project.
python setup.py bdist_wheel bdist_wininst -b doc\logos\scikit-learn-logo.bmp

@rem Install the generated wheel package to test it
pip install --pre --no-index --find-links dist\ scikit-learn

if %errorlevel% neq 0 exit /b %errorlevel%
