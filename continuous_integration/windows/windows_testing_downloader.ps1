# Author: Kyle Kastner <kastnerkyle@gmail.com>
# License: BSD 3 clause

#This script is a helper to download the base python, numpy, and scipy
#packages from their respective websites.
#It will also install a version of VC++ and nosetests

#This is a stopgap solution to make Windows testing easier
#until Windows CI issues are resolved.Some of the files are hosted
#on Kyle Kastner's public Dropbox, due to javascript download issues

#Make sure to run:
#set-executionpolicy unrestricted
#so this script can run on Windows Server 2012!

param (
    [string]$python = "None"
)

function DisableInternetExplorerESC {
    #Shamelessly stolen from
    #http://stackoverflow.com/questions/9368305/disable-ie-security-on-windows-server-via-powershell
    $AdminKey = "HKLM:\SOFTWARE\Microsoft\Active Setup\Installed Components\{A509B1A7-37EF-4b3f-8CFC-4F3A74704073}"
    $UserKey = "HKLM:\SOFTWARE\Microsoft\Active Setup\Installed Components\{A509B1A8-37EF-4b3f-8CFC-4F3A74704073}"
    Set-ItemProperty -Path $AdminKey -Name "IsInstalled" -Value 0
    Set-ItemProperty -Path $UserKey -Name "IsInstalled" -Value 0
    Stop-Process -Name Explorer
    Write-Host "IE Enhanced Security Configuration (ESC) has been disabled." -ForegroundColor Green
}

function Python27URLs {
    $urls = @{
        "python" = "https://www.python.org/ftp/python/2.7.7/python-2.7.7.msi"
        "numpy" = "http://sourceforge.net/projects/numpy/files/NumPy/1.8.1/numpy-1.8.1-win32-superpack-python2.7.exe/download"
        "scipy" = "http://sourceforge.net/projects/scipy/files/scipy/0.14.0/scipy-0.14.0-win32-superpack-python2.7.exe/download"
        "nosetests" = "https://dl.dropboxusercontent.com/u/15378192/nose-1.3.3.win32-py2.7.exe"
        "vc++" = "https://dl.dropboxusercontent.com/u/15378192/vcredist_x86.exe"
        "setuptools" = "https://dl.dropboxusercontent.com/u/15378192/setuptools-3.6.win32-py2.7.exe"
        "get-pip" = "https://raw.github.com/pypa/pip/master/contrib/get-pip.py"
    }
    return $urls    
}

function Python34URLs {
    $urls = @{
        "python" = "https://www.python.org/ftp/python/3.4.1/python-3.4.1.msi"
        "numpy" = "http://sourceforge.net/projects/numpy/files/NumPy/1.8.1/numpy-1.8.1-win32-superpack-python3.4.exe/download"
        "scipy" = "http://sourceforge.net/projects/scipy/files/scipy/0.14.0/scipy-0.14.0-win32-superpack-python3.4.exe/download"
        "nosetests" = "https://dl.dropboxusercontent.com/u/15378192/nose-1.3.3.win32-py3.4.exe"
        "vc++" = "https://dl.dropboxusercontent.com/u/15378192/vcredist_x86.exe"
        "setuptools" = "https://dl.dropboxusercontent.com/u/15378192/setuptools-3.6.win32-py3.4.exe"
    }
    return $urls    
}

function main {
    $versions = @{
        "2.7" = Python27URLs
        "3.4" = Python34URLs
    }
    if ($python -eq "None") {
        Write-Host "Python version MUST be passed with -python"
        return $FALSE
    }

    if (!($versions.ContainsKey($python))) {
        Write-Host "This Python version is not yet supported!"
        return $FALSE
    }

    DisableInternetExplorerESC

    $webclient = New-Object System.Net.WebClient
    $pystring = $python -replace "\."
    ForEach ($key in $versions[$python].Keys) {
        $url = $versions[$python][$key]
        $file = $key + "_py" + $pystring
        if ($url -match "(\.exe)") {
            $file = $file + ".exe"
        } elseif ($url -match "(\.msi)") {
            $file = $file + ".msi"
        } else {
            $file = $file + ".py"
        }
        $basedir = $pwd.Path + "\"
        $filepath = $basedir + $file
        Write-Host "Downloading" $file "from" $url
        $webclient.DownloadFile($url, $filepath)
        Write-Host "File saved at" $filepath
    }
    return $TRUE
}

main