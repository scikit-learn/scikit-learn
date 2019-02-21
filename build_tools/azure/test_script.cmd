call activate %VIRTUALENV%

pip install --timeout=60 --trusted-host 28daf2247a33ed269873-7b1aad3fab3cc330e1fd9d109892382a.r6.cf2.rackcdn.com -r build_tools\appveyor\requirements.txt

mkdir %TMP_FOLDER%
cd %TMP_FOLDER%

pytest --junitxml=%JUNITXML% --showlocals --durations=20 --pyargs sklearn.preprocessing
