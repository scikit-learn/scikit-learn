# Sample script to install Python and pip under Windows
# Authors: Olivier Grisel, Jonathan Helmus, Kyle Kastner, and Alex Willmer
# License: CC0 1.0 Universal: http://creativecommons.org/publicdomain/zero/1.0/

$MINICONDA_URL = "http://repo.continuum.io/miniconda/"
$BASE_URL = "https://www.python.org/ftp/python/"
$GET_PIP_URL = "https://bootstrap.pypa.io/get-pip.py"
$GET_PIP_PATH = "C:\get-pip.py"

$PYTHON_PRERELEASE_REGEX = @"
(?x)
(?<major>\d+)
\.
(?<minor>\d+)
\.
(?<micro>\d+)
(?<prerelease>[a-z]{1,2}\d+)
"@


function Download ($filename, $url) {
    $webclient = New-Object System.Net.WebClient

    $basedir = $pwd.Path + "\"
    $filepath = $basedir + $filename
    if (Test-Path $filename) {
        Write-Host "Reusing" $filepath
        return $filepath
    }

    # Download and retry up to 3 times in case of network transient errors.
    Write-Host "Downloading" $filename "from" $url
    $retry_attempts = 2
    for ($i = 0; $i -lt $retry_attempts; $i++) {
        try {
            $webclient.DownloadFile($url, $filepath)
            break
        }
        Catch [Exception]{
            Start-Sleep 1
        }
    }
    if (Test-Path $filepath) {
        Write-Host "File saved at" $filepath
    } else {
        # Retry once to get the error message if any at the last try
        $webclient.DownloadFile($url, $filepath)
    }
    return $filepath
}


function ParsePythonVersion ($python_version) {
    if ($python_version -match $PYTHON_PRERELEASE_REGEX) {
        return ([int]$matches.major, [int]$matches.minor, [int]$matches.micro,
                $matches.prerelease)
    }
    $version_obj = [version]$python_version
    return ($version_obj.major, $version_obj.minor, $version_obj.build, "")
}


function DownloadPython ($python_version, $platform_suffix) {
    $major, $minor, $micro, $prerelease = ParsePythonVersion $python_version

    if (($major -le 2 -and $micro -eq 0) `
        -or ($major -eq 3 -and $minor -le 2 -and $micro -eq 0) `
        ) {
        $dir = "$major.$minor"
        $python_version = "$major.$minor$prerelease"
    } else {
        $dir = "$major.$minor.$micro"
    }

    if ($prerelease) {
        if (($major -le 2) `
            -or ($major -eq 3 -and $minor -eq 1) `
            -or ($major -eq 3 -and $minor -eq 2) `
            -or ($major -eq 3 -and $minor -eq 3) `
            ) {
            $dir = "$dir/prev"
        }
    }

    if (($major -le 2) -or ($major -le 3 -and $minor -le 4)) {
        $ext = "msi"
        if ($platform_suffix) {
            $platform_suffix = ".$platform_suffix"
        }
    } else {
        $ext = "exe"
        if ($platform_suffix) {
            $platform_suffix = "-$platform_suffix"
        }
    }

    $filename = "python-$python_version$platform_suffix.$ext"
    $url = "$BASE_URL$dir/$filename"
    $filepath = Download $filename $url
    return $filepath
}


function InstallPython ($python_version, $architecture, $python_home) {
    Write-Host "Installing Python" $python_version "for" $architecture "bit architecture to" $python_home
    if (Test-Path $python_home) {
        Write-Host $python_home "already exists, skipping."
        return $false
    }
    if ($architecture -eq "32") {
        $platform_suffix = ""
    } else {
        $platform_suffix = "amd64"
    }
    $installer_path = DownloadPython $python_version $platform_suffix
    $installer_ext = [System.IO.Path]::GetExtension($installer_path)
    Write-Host "Installing $installer_path to $python_home"
    $install_log = $python_home + ".log"
    if ($installer_ext -eq '.msi') {
        InstallPythonMSI $installer_path $python_home $install_log
    } else {
        InstallPythonEXE $installer_path $python_home $install_log
    }
    if (Test-Path $python_home) {
        Write-Host "Python $python_version ($architecture) installation complete"
    } else {
        Write-Host "Failed to install Python in $python_home"
        Get-Content -Path $install_log
        Exit 1
    }
}


function InstallPythonEXE ($exepath, $python_home, $install_log) {
    $install_args = "/quiet InstallAllUsers=1 TargetDir=$python_home"
    RunCommand $exepath $install_args
}


function InstallPythonMSI ($msipath, $python_home, $install_log) {
    $install_args = "/qn /log $install_log /i $msipath TARGETDIR=$python_home"
    $uninstall_args = "/qn /x $msipath"
    RunCommand "msiexec.exe" $install_args
    if (-not(Test-Path $python_home)) {
        Write-Host "Python seems to be installed else-where, reinstalling."
        RunCommand "msiexec.exe" $uninstall_args
        RunCommand "msiexec.exe" $install_args
    }
}

function RunCommand ($command, $command_args) {
    Write-Host $command $command_args
    Start-Process -FilePath $command -ArgumentList $command_args -Wait -Passthru
}


function InstallPip ($python_home) {
    $pip_path = $python_home + "\Scripts\pip.exe"
    $python_path = $python_home + "\python.exe"
    if (-not(Test-Path $pip_path)) {
        Write-Host "Installing pip..."
        $webclient = New-Object System.Net.WebClient
        $webclient.DownloadFile($GET_PIP_URL, $GET_PIP_PATH)
        Write-Host "Executing:" $python_path $GET_PIP_PATH
        & $python_path $GET_PIP_PATH
    } else {
        Write-Host "pip already installed."
    }
}


function DownloadMiniconda ($python_version, $platform_suffix) {
    if ($python_version -eq "3.4") {
        $filename = "Miniconda3-3.5.5-Windows-" + $platform_suffix + ".exe"
    } else {
        $filename = "Miniconda-3.5.5-Windows-" + $platform_suffix + ".exe"
    }
    $url = $MINICONDA_URL + $filename
    $filepath = Download $filename $url
    return $filepath
}


function InstallMiniconda ($python_version, $architecture, $python_home) {
    Write-Host "Installing Python" $python_version "for" $architecture "bit architecture to" $python_home
    if (Test-Path $python_home) {
        Write-Host $python_home "already exists, skipping."
        return $false
    }
    if ($architecture -eq "32") {
        $platform_suffix = "x86"
    } else {
        $platform_suffix = "x86_64"
    }
    $filepath = DownloadMiniconda $python_version $platform_suffix
    Write-Host "Installing" $filepath "to" $python_home
    $install_log = $python_home + ".log"
    $args = "/S /D=$python_home"
    Write-Host $filepath $args
    Start-Process -FilePath $filepath -ArgumentList $args -Wait -Passthru
    if (Test-Path $python_home) {
        Write-Host "Python $python_version ($architecture) installation complete"
    } else {
        Write-Host "Failed to install Python in $python_home"
        Get-Content -Path $install_log
        Exit 1
    }
}


function InstallMinicondaPip ($python_home) {
    $pip_path = $python_home + "\Scripts\pip.exe"
    $conda_path = $python_home + "\Scripts\conda.exe"
    if (-not(Test-Path $pip_path)) {
        Write-Host "Installing pip..."
        $args = "install --yes pip"
        Write-Host $conda_path $args
        Start-Process -FilePath "$conda_path" -ArgumentList $args -Wait -Passthru
    } else {
        Write-Host "pip already installed."
    }
}

function main () {
    InstallPython $env:PYTHON_VERSION $env:PYTHON_ARCH $env:PYTHON
    InstallPip $env:PYTHON
}

main
