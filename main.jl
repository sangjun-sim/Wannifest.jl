#!/usr/bin/env julia

const WANNIFEST_DIR = @__DIR__
import Pkg
Pkg.activate(WANNIFEST_DIR; io=devnull)

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
