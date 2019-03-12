@echo on

call activate %VIRTUALENV%

mkdir %TMP_FOLDER%
cd %TMP_FOLDER%

if "%CHECK_WARNINGS%" == "true" (
    set PYTEST_ARGS=%PYTEST_ARGS% -Werror::DeprecationWarning -Werror::FutureWarning
)

if "%COVERAGE%" == "true" (
    set PYTEST_ARGS=%PYTEST_ARGS% --cov sklearn
)

pytest sklearn\linear_model\tests\test_sgd.py::test_sgd_proba  -s --junitxml=%JUNITXML% --showlocals --durations=20
