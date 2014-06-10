# Author: Kyle Kastner <kastnerkyle@gmail.com>
# License: BSD 3 clause

# This script is a helper to download the base python, numpy, and scipy
# packages from their respective websites.

# This is a stopgap solution to make Windows testing easier
# until Windows CI issues are resolved.

# Rackspace's default Windows VMs have several security features enabled by default.
# The DisableInternetExplorerESC function disables a feature which
# prevents any webpage from opening without explicit permission.
# This is a default setting of Windows VMs on Rackspace, and makes it annoying to
# download other pacakages to test!

# Powershell scripts are also disabled by default. One must run the command:
# set-executionpolicy unrestricted
# from a Powershell terminal with administrator rights to enable scripts.
# To start an administrator Powershell terminal, right click second icon from the left on Windows Server 2012's bottom taskbar

param (
    [string]$python = "None"
)

function DisableInternetExplorerESC {
    # Disables InternetExplorerESC to enable easier manual downloads of testing packages
    # See
    # http://stackoverflow.com/questions/9368305/disable-ie-security-on-windows-server-via-powershell
    $AdminKey = "HKLM:\SOFTWARE\Microsoft\Active Setup\Installed Components\{A509B1A7-37EF-4b3f-8CFC-4F3A74704073}"
    $UserKey = "HKLM:\SOFTWARE\Microsoft\Active Setup\Installed Components\{A509B1A8-37EF-4b3f-8CFC-4F3A74704073}"
    Set-ItemProperty -Path $AdminKey -Name "IsInstalled" -Value 0
    Set-ItemProperty -Path $UserKey -Name "IsInstalled" -Value 0
    Stop-Process -Name Explorer
    Write-Host "IE Enhanced Security Configuration (ESC) has been disabled." -ForegroundColor Green
}

function Python27URLs {
    # Function returns a dictionary of packages to download for Python 2.7
    $urls = @{
        "python" = "https://www.python.org/ftp/python/2.7.7/python-2.7.7.msi"
        "numpy" = "http://sourceforge.net/projects/numpy/files/NumPy/1.8.1/numpy-1.8.1-win32-superpack-python2.7.exe/download"
        "scipy" = "http://sourceforge.net/projects/scipy/files/scipy/0.14.0/scipy-0.14.0-win32-superpack-python2.7.exe/download"
        "get-pip" = "https://raw.github.com/pypa/pip/master/contrib/get-pip.py"
    }
    return $urls    
}

function Python34URLs {
    # Function returns a dictionary of packages to download for Python 3.4
    $urls = @{
        "python" = "https://www.python.org/ftp/python/3.4.1/python-3.4.1.msi"
        "numpy" = "http://sourceforge.net/projects/numpy/files/NumPy/1.8.1/numpy-1.8.1-win32-superpack-python3.4.exe/download"
        "scipy" = "http://sourceforge.net/projects/scipy/files/scipy/0.14.0/scipy-0.14.0-win32-superpack-python3.4.exe/download"
    }
    return $urls    
}

function DownloadPackages ($package_dict, $append_string) {
    $webclient = New-Object System.Net.WebClient

    ForEach ($key in $package_dict.Keys) {
        $url = $package_dict[$key]
        $file = $key + "_py" + $append_string
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
}

function InstallPython($match_string) {
    $pkg_regex = "python" + $match_string + "*"
    $pkg = Get-ChildItem -Filter $pkg_regex -Name
    Invoke-Expression -Command ".\$pkg /qn"
}

function InstallPip($match_string, $python_version) {
    $pkg_regex = "get-pip" + $match_string + "*"
    $py = $python_version -replace "\."
    $pkg = Get-ChildItem -Filter $pkg_regex -Name
    $python_path = "C:\Python" + $py + "\python.exe"
    Invoke-Expression -Command "$python_path $pkg"
}

function PipInstall($pkg_name, $python_version) {
    $py = $python_version -replace "\."
    $pip = "C:\Python" + $py + "\Scripts\pip.exe"
    Invoke-Expression -Command "$pip install $pkg_name"
}

function InstallNose($match_string, $python_version) {
    PipInstall "nose" $python_version 
}

function InstallNumpy($match_string) {
    $pkg_regex = "numpy" + $match_string + "*"
    $pkg = Get-ChildItem -Filter $pkg_regex -Name
    Invoke-Expression -Command ".\$pkg"

}

function InstallScipy($match_string) {
    $pkg_regex = "python" + $match_string + "*"
    $pkg = Get-ChildItem -Filter $pkg_regex -Name
    Invoke-Expression -Command ".\$pkg"
}

function main {
    $versions = @{
        "2.7" = Python27URLs
        "3.4" = Python34URLs
    }

    if (($python -eq "None") -or !($versions.ContainsKey($python))) {
        Write-Host "Python version not recognized!"
        Write-Host "Pass python version with -python"
        Write-Host "Current versions supported are:"
        ForEach ($key in $versions.Keys) {
            Write-Host $key
        }
        return
    }

    DisableInternetExplorerESC
    $pystring = $python -replace "\."
    $pystring = "_py" + $pystring
    $package_dict = $versions[$python]
    DownloadPackages $package_dict $pystring
    InstallPython $pystring
    if ($package_dict.ContainsKey("get-pip")) {
       InstallPip $pystring $python
    }
    InstallNose $pystring $python
    InstallNumpy $pystring
    InstallScipy $pystring
    return
}

main
