set datetimef=%date:~-4%_%date:~3,2%_%date:~0,2%_%time:~0,2%_%time:~3,2%_%time:~6,2%
set datetimef=%date:~-4%_%date:~3,2%_%date:~3,2%_%time:~3,2%_%time:~3,2%_%time:~6,2%
echo %datetimef%
FOR /F "TOKENS=1 eol=/ DELIMS=/ " %%A IN ('DATE/T') DO SET dd=%%A
FOR /F "TOKENS=1,2 eol=/ DELIMS=/ " %%A IN ('DATE/T') DO SET mm=%%B
FOR /F "TOKENS=1,2,3 eol=/ DELIMS=/ " %%A IN ('DATE/T') DO SET yyyy=%%C
SET todaysdate=%yyyy%%mm%%dd%%TIME%
echo %dd%
echo %mm%
echo %yyyy%
echo "The Date Time File Name is:"%datetimef%
echo "Datetimef: " %datetimef%
type nul > %datetimef%
date /t >> %datetimef%
date /t >> %datetimef%
date /t >> %datetimef%
git add .
git commit -m "chore: updating build tasks, package manager configs, ; no production code change made"
git push -f origin master