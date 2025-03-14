# Copied from GStreamer project
# Author: Seungha Yang <seungha.yang@navercorp.com>
#         Xavier Claessens <xclaesse@gmail.com>

$i=1
$ppid=$PID
do {
  $ppid=(Get-CimInstance Win32_Process -Filter "ProcessId=$ppid").parentprocessid
  $pname=(Get-Process -id $ppid).Name
  if($pname -eq "cmd" -Or $pname -eq "powershell" -Or $pname -eq "pwsh") {
    Write-Host ("{0}.exe" -f $pname)
    Break
  }
  # not found yet, find grand parent
  # 10 times iteration seems to be sufficient
  $i++
} while ($i -lt 10)
