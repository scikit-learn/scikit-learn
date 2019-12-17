@echo on

@rem Only 64 bit uses conda and uses a python newer than 3.5
IF "%PYTHON_ARCH%"=="64" (
    call activate %VIRTUALENV%
    set PYTEST_ARGS=%PYTEST_ARGS% -n2
)

cp ntpath.py c:\miniconda\envs\testvenv\lib\ntpath.py
mkdir %TMP_FOLDER%
cd %TMP_FOLDER%

if "%CHECK_WARNINGS%" == "true" (
    set PYTEST_ARGS=%PYTEST_ARGS% -Werror::DeprecationWarning -Werror::FutureWarning
)

if "%COVERAGE%" == "true" (
    set PYTEST_ARGS=%PYTEST_ARGS% --cov sklearn
)

pytest -x -k test_import_is_deprecated --junitxml=%JUNITXML% --showlocals --durations=20 %PYTEST_ARGS% --pyargs sklearn
