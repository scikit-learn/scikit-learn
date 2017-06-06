# Author: Kyle Kastner <kastnerkyle@gmail.com>
# License: BSD 3 clause

# This script is a helper to download the base python, numpy, and scipy
# packages from their respective websites.
# To quickly execute the script, run the following Powershell command:
# powershell.exe -ExecutionPolicy unrestricted "iex ((new-object net.webclient).DownloadString('https://raw.githubusercontent.com/scikit-learn/scikit-learn/master/continuous_integration/windows/windows_testing_downloader.ps1'))"

# This is a stopgap solution to make Windows testing easier
# until Windows CI issues are resolved.

# Rackspace's default Windows VMs have several security features enabled by default.
# The DisableInternetExplorerESC function disables a feature which
# prevents any webpage from opening without explicit permission.
# This is a default setting of Windows VMs on Rackspace, and makes it annoying to
# download other packages to test!

# Powershell scripts are also disabled by default. One must run the command:
# set-executionpolicy unrestricted
# from a Powershell terminal with administrator rights to enable scripts.
# To start an administrator Powershell terminal, right click second icon from the left on Windows Server 2012's bottom taskbar.

param (
    [string]$python = "None",
    [string]$nogit = "False"
)

function DisableInternetExplorerESC {
    # Disables InternetExplorerESC to enable easier manual downloads of testing packages.
    # http://stackoverflow.com/questions/9368305/disable-ie-security-on-windows-server-via-powershell
    $AdminKey = "HKLM:\SOFTWARE\Microsoft\Active Setup\Installed Components\{A509B1A7-37EF-4b3f-8CFC-4F3A74704073}"
    $UserKey = "HKLM:\SOFTWARE\Microsoft\Active Setup\Installed Components\{A509B1A8-37EF-4b3f-8CFC-4F3A74704073}"
    Set-ItemProperty -Path $AdminKey -Name "IsInstalled" -Value 0
    Set-ItemProperty -Path $UserKey -Name "IsInstalled" -Value 0
    Stop-Process -Name Explorer
    Write-Host "IE Enhanced Security Configuration (ESC) has been disabled." -ForegroundColor Green
}

function DownloadPackages ($package_dict, $append_string) {
    $webclient = New-Object System.Net.WebClient

    ForEach ($key in $package_dict.Keys) {
        $url = $package_dict[$key]
        $file = $key + $append_string
        if ($url -match "(\.*exe$)") {
            $file = $file + ".exe"
        } elseif ($url -match "(\.*msi$)") {
            $file = $file + ".msi"
        } else {
            $file = $file + ".py"
        }
        $basedir = $pwd.Path + "\"
        $filepath = $basedir + $file
        Write-Host "Downloading" $file "from" $url

        # Retry up to 5 times in case of network transient errors.
        $retry_attempts = 5
        for($i=0; $i -lt $retry_attempts; $i++){
             try{
                  $webclient.DownloadFile($url, $filepath)
                  break
             }
             Catch [Exception]{
                  Start-Sleep 1
             }
        }
        Write-Host "File saved at" $filepath
    }
}

function InstallPython($match_string) {
    $pkg_regex = "python" + $match_string + "*"
    $pkg = Get-ChildItem -Filter $pkg_regex -Name
    Invoke-Expression -Command "msiexec /qn /i $pkg"
    
    Write-Host "Installing Python"
    Start-Sleep 25
    Write-Host "Python installation complete"
}

function InstallPip($match_string, $python_version) {
    $pkg_regex = "get-pip" + $match_string + "*"
    $py = $python_version -replace "\."
    $pkg = Get-ChildItem -Filter $pkg_regex -Name
    $python_path = "C:\Python" + $py + "\python.exe"
    Invoke-Expression -Command "$python_path $pkg"
}

function EnsurePip($python_version) {
    $py = $python_version -replace "\."
    $python_path = "C:\Python" + $py + "\python.exe"
    Invoke-Expression -Command "$python_path -m ensurepip"
}

function GetPythonHome($python_version) {
    $py = $python_version -replace "\."
    $pypath = "C:\Python" + $py + "\"
    return $pypath
}

function GetPipPath($python_version) {
    $py = $python_version -replace "\."
    $pypath = GetPythonHome $python_version
    if ($py.StartsWith("3")) {
        $pip = $pypath + "Scripts\pip3.exe"
    } else {
        $pip = $pypath + "Scripts\pip.exe"
    }
    return $pip
}

function PipInstall($pkg_name, $python_version, $extra_args) {
    $pip = GetPipPath $python_version
    Invoke-Expression -Command "$pip install $pkg_name"
}

function InstallNose($python_version) {
    PipInstall "nose" $python_version 
}

function WheelInstall($name, $url, $python_version) {
    $pip = GetPipPath $python_version
    $args = "install --use-wheel --no-index"
    Invoke-Expression -Command "$pip $args $url $name"
}

function InstallWheel($python_version) {
    PipInstall "virtualenv" $python_version
    PipInstall "wheel" $python_version
}

function InstallNumpy($package_dict, $python_version) {
    #Don't pass name so we can use URL directly.
    WheelInstall "" $package_dict["numpy"] $python_version
}

function InstallScipy($package_dict, $python_version) {
    #Don't pass name so we can use URL directly.
    WheelInstall "" $package_dict["scipy"] $python_version
}

function InstallGit {
    $pkg_regex = "git*"
    $pkg = Get-ChildItem -Filter $pkg_regex -Name
    $pkg_cmd = $pwd.ToString() + "\" + $pkg + " /verysilent"
    Invoke-Expression -Command $pkg_cmd
    
    Write-Host "Installing Git"
    Start-Sleep 20 
    # Remove the installer - seems to cause weird issues with Git Bash
    Invoke-Expression -Command "rm git.exe"
    Write-Host "Git installation complete"
}

function ReadAndUpdateFromRegistry {
    # http://stackoverflow.com/questions/14381650/how-to-update-windows-powershell-session-environment-variables-from-registry
    foreach($level in "Machine","User") {
    [Environment]::GetEnvironmentVariables($level).GetEnumerator() | % {
       # For Path variables, append the new values, if they're not already in there
       if($_.Name -match 'Path$') { 
          $_.Value = ($((Get-Content "Env:$($_.Name)") + ";$($_.Value)") -split ';' | Select -unique) -join ';'
       }
       $_
       } | Set-Content -Path { "Env:$($_.Name)" }
    }
}

function UpdatePaths($python_version) {
    #This function makes local path updates required in order to install Python and supplementary packages in a single shell.
    $pypath = GetPythonHome $python_version
    $env:PATH = $env:PATH + ";" + $pypath
    $env:PYTHONPATH = $pypath + "DLLs;" + $pypath + "Lib;" + $pypath + "Lib\site-packages"
    $env:PYTHONHOME = $pypath
    Write-Host "PYTHONHOME temporarily set to" $env:PYTHONHOME
    Write-Host "PYTHONPATH temporarily set to" $env:PYTHONPATH
    Write-Host "PATH temporarily set to" $env:PATH
}

function Python27URLs {
    # Function returns a dictionary of packages to download for Python 2.7.
    $urls = @{
        "python" = "https://www.python.org/ftp/python/2.7.7/python-2.7.7.msi"
        "numpy" = "http://28daf2247a33ed269873-7b1aad3fab3cc330e1fd9d109892382a.r6.cf2.rackcdn.com/numpy-1.8.1-cp27-none-win32.whl"
        "scipy" = "http://28daf2247a33ed269873-7b1aad3fab3cc330e1fd9d109892382a.r6.cf2.rackcdn.com/scipy-0.14.0-cp27-none-win32.whl"
        "get-pip" = "https://bootstrap.pypa.io/get-pip.py"
    }
    return $urls    
}

function Python34URLs {
    # Function returns a dictionary of packages to download for Python 3.4.
    $urls = @{
        "python" = "https://www.python.org/ftp/python/3.4.1/python-3.4.1.msi"
        "numpy" = "http://28daf2247a33ed269873-7b1aad3fab3cc330e1fd9d109892382a.r6.cf2.rackcdn.com/numpy-1.8.1-cp34-none-win32.whl"
        "scipy" = "http://28daf2247a33ed269873-7b1aad3fab3cc330e1fd9d109892382a.r6.cf2.rackcdn.com/scipy-0.14.0-cp34-none-win32.whl"
    }
    return $urls    
}

function GitURLs {
    # Function returns a dictionary of packages to download for Git
    $urls = @{
        "git" = "https://github.com/msysgit/msysgit/releases/download/Git-1.9.4-preview20140611/Git-1.9.4-preview20140611.exe"
    }
    return $urls    
}

function main {
    $versions = @{
        "2.7" = Python27URLs
        "3.4" = Python34URLs
    }

    if ($nogit -eq "False") {
        Write-Host "Downloading and installing Gitbash"
        $urls = GitURLs
        DownloadPackages $urls ""
        InstallGit ".exe"
    }

    if (($python -eq "None")) {
        Write-Host "Installing all supported python versions"
        Write-Host "Current versions supported are:"
        ForEach ($key in $versions.Keys) {
            Write-Host $key
            $all_python += @($key)
        }
    } elseif(!($versions.ContainsKey($python))) {
        Write-Host "Python version not recognized!"
        Write-Host "Pass python version with -python"
        Write-Host "Current versions supported are:"
        ForEach ($key in $versions.Keys) {
            Write-Host $key
        }
        return
    } else {
        $all_python += @($python)
    }
    ForEach ($py in $all_python) {
        Write-Host "Installing Python" $py
        DisableInternetExplorerESC
        $pystring = $py -replace "\."
        $pystring = "_py" + $pystring
        $package_dict = $versions[$py]

        # This will download the whl packages as well which is
        # clunky but makes configuration simpler.
        DownloadPackages $package_dict $pystring
        UpdatePaths $py
        InstallPython $pystring
        ReadAndUpdateFromRegistry
        if ($package_dict.ContainsKey("get-pip")) {
           InstallPip $pystring $py
        } else {
           EnsurePip $py 
        }
        InstallNose $py
        InstallWheel $py

        # The installers below here use wheel packages.
        # Wheels were created from CGohlke's installers with
        # wheel convert <exefile>
        # These are hosted in Rackspace Cloud Files.
        InstallNumpy $package_dict $py
        InstallScipy $package_dict $py
    }
    return
}

main
