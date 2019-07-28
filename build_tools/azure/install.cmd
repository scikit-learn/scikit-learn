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
    conda create -n %VIRTUALENV% -q -y python=%PYTHON_VERSION% numpy scipy cython matplotlib pytest=%PYTEST_VERSION% wheel pillow joblib

    call activate %VIRTUALENV%
    pip install pytest-xdist
) else (
    pip install numpy scipy cython pytest wheel pillow joblib
)
if "%COVERAGE%" == "true" (
    pip install coverage codecov pytest-cov
)
python --version
pip --version

@rem Install the build and runtime dependencies of the project.
python setup.py bdist_wheel bdist_wininst -b doc\logos\scikit-learn-logo.bmp

@rem Install the generated wheel package to test it
pip install --pre --no-index --find-links dist\ scikit-learn

if %errorlevel% neq 0 exit /b %errorlevel%
