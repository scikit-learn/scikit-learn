@rem https://github.com/numba/numba/blob/master/buildscripts/incremental/setup_conda_environment.cmd
@rem The cmd /C hack circumvents a regression where conda installs a conda.bat
@rem script in non-root environments.
set CONDA_INSTALL=cmd /C conda install -q -y
set PIP_INSTALL=pip install -q

@echo on

@rem Deactivate any environment
call deactivate
@rem Display root environment (for debugging)
conda list
@rem Clean up any left-over from a previous build
conda remove --all -q -y -n %VIRTUALENV%
conda env create -n %VIRTUALENV% python=%CONDA_PY% -y

call activate %VIRTUALENV%
conda list

@rem https://github.com/numba/numba/blob/master/buildscripts/incremental/build.cmd

@rem Build numba extensions without silencing compile errors
python setup.py build_ext -q --inplace

@rem Install pandas locally
python -m pip install -e .

if %errorlevel% neq 0 exit /b %errorlevel%
