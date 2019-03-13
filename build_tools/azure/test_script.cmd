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

pytest -s -k "test_seq_dataset_shuffle or test_sgd_proba" --pyargs sklearn
