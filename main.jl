#!/usr/bin/env julia

const WANNIFEST_DIR = @__DIR__
const STARTUP_FILE_DISABLED = 2
const REEXEC_ENV = "WANNIFEST_NO_STARTUP_REEXEC"

function _reexec_without_startup_file_if_needed()
    script_path = abspath(@__FILE__)
    if abspath(PROGRAM_FILE) != script_path
        return nothing
    end
    if Base.JLOptions().startupfile == STARTUP_FILE_DISABLED
        return nothing
    end
    if get(ENV, REEXEC_ENV, "") == "1"
        return nothing
    end

    env = copy(ENV)
    env[REEXEC_ENV] = "1"
    child_cmd = setenv(
        `$(Base.julia_cmd()) --project=$WANNIFEST_DIR --startup-file=no $script_path $ARGS`,
        env,
    )
    proc = run(pipeline(child_cmd; stdout=stdout, stderr=stderr); wait=false)
    wait(proc)
    exit(proc.exitcode)
end

_reexec_without_startup_file_if_needed()

using Wannifest

if abspath(PROGRAM_FILE) == @__FILE__
    try
        exit(Wannifest.run_wannifest(ARGS))
    catch err
        if err isa ArgumentError
            println(stderr, err.msg)
            Wannifest.print_wannifest_usage(stderr)
            exit(1)
        end
        rethrow()
    end
end
