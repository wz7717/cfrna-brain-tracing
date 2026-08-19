[CmdletBinding()]
param(
    [ValidateSet("full", "portable")]
    [string]$Profile = "full",
    [switch]$VerifyOnly,
    [string]$ExternalDataRoot,
    [string]$InputManifest,
    [string]$OutputDir,
    [string[]]$Step,
    [switch]$ReleaseGate
)

# PowerShell is intentionally a thin wrapper.  The Python runner owns the
# canonical external-data resolution, hash verification, command registry and
# fail-closed audit so both entry points enforce exactly the same contract.
$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location -LiteralPath $Root

$Arguments = @("reproduce_all.py", "--profile", $Profile)
if ($VerifyOnly) {
    $Arguments += "--verify-only"
}
if ($ExternalDataRoot) {
    $Arguments += @("--external-data-root", $ExternalDataRoot)
}
if ($InputManifest) {
    $Arguments += @("--input-manifest", $InputManifest)
}
if ($OutputDir) {
    $Arguments += @("--output-dir", $OutputDir)
}
if ($ReleaseGate) {
    $Arguments += "--release-gate"
}
foreach ($RequestedStep in $Step) {
    $Arguments += @("--step", $RequestedStep)
}

& python @Arguments
exit $LASTEXITCODE
