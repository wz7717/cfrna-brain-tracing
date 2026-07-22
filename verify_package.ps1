$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$Manifest = Join-Path $Root "SHA256SUMS.txt"
if (-not (Test-Path -LiteralPath $Manifest)) {
    throw "SHA256SUMS.txt is missing"
}

$Failures = 0
Get-Content -LiteralPath $Manifest -Encoding UTF8 | ForEach-Object {
    if (-not $_.Trim()) { return }
    $Parts = $_ -split "  ", 2
    if ($Parts.Count -ne 2) {
        throw "Malformed checksum line: $_"
    }
    $Expected = $Parts[0].ToUpperInvariant()
    $Relative = $Parts[1].Replace("/", [IO.Path]::DirectorySeparatorChar)
    $Path = Join-Path $Root $Relative
    if (-not (Test-Path -LiteralPath $Path)) {
        Write-Host "MISSING $Relative"
        $Failures += 1
        return
    }
    $Actual = (Get-FileHash -Algorithm SHA256 -LiteralPath $Path).Hash
    if ($Actual -ne $Expected) {
        Write-Host "FAILED  $Relative"
        $Failures += 1
    }
}

if ($Failures -ne 0) {
    throw "$Failures package files failed verification"
}
Write-Host "All package files passed SHA256 verification."
