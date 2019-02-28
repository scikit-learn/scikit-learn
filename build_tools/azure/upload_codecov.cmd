@echo on

call activate %VIRTUALENV%

copy %TEST_DIR%\.coverage %BUILD_REPOSITORY_LOCALPATH%

codecov --root %BUILD_REPOSITORY_LOCALPATH% -t %CODECOV_TOKEN% && (
  echo codecov was successful
) || (
  echo codecov failed
)
