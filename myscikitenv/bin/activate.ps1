# This file must be dot sourced from PoSh; you cannot run it directly. Do this: . ./activate.ps1

$script:THIS_PATH = $myinvocation.mycommand.path
$script:BASE_DIR = split-path (resolve-path "$THIS_PATH/..") -Parent

function global:deactivate([switch] $NonDestructive)
{
    if (test-path variable:_OLD_VIRTUAL_PATH)
    {
        $env:PATH = $variable:_OLD_VIRTUAL_PATH
        remove-variable "_OLD_VIRTUAL_PATH" -scope global
    }

    if (test-path function:_old_virtual_prompt)
    {
        $function:prompt = $function:_old_virtual_prompt
        remove-item function:\_old_virtual_prompt
    }

    if ($env:VIRTUAL_ENV)
    {
        $old_env = split-path $env:VIRTUAL_ENV -leaf
        remove-item env:VIRTUAL_ENV -erroraction silentlycontinue
    }

    if (!$NonDestructive)
    {
        # Self destruct!
        remove-item function:deactivate
        remove-item function:pydoc
    }
}

function global:pydoc
{
    python -m pydoc $args
}

# unset irrelevant variables
deactivate -nondestructive

$VIRTUAL_ENV = $BASE_DIR
$env:VIRTUAL_ENV = $VIRTUAL_ENV

$global:_OLD_VIRTUAL_PATH = $env:PATH
$env:PATH = "$env:VIRTUAL_ENV/bin:" + $env:PATH
if (!$env:VIRTUAL_ENV_DISABLE_PROMPT)
{
    function global:_old_virtual_prompt
    {
        ""
    }
    $function:_old_virtual_prompt = $function:prompt
    function global:prompt
    {
        # Add a prefix to the current prompt, but don't discard it.
        write-host "($( split-path $env:VIRTUAL_ENV -leaf )) " -nonewline
        & $function:_old_virtual_prompt
    }
}
