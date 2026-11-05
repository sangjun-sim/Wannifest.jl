module InteractiveWait

export display_available, wait_for_plot_window

# GR window types that render to a file/buffer instead of an on-screen window.
const OFFSCREEN_GKS_TYPES = ("100", "nul", "png", "svg", "pdf", "ps")

function display_available()::Bool
    gks = lowercase(strip(get(ENV, "GKSwstype", "")))
    isempty(gks) || return !(gks in OFFSCREEN_GKS_TYPES)
    (Sys.iswindows() || Sys.isapple()) && return true
    return !isempty(get(ENV, "DISPLAY", "")) || !isempty(get(ENV, "WAYLAND_DISPLAY", ""))
end

stdin_is_terminal()::Bool = stdin isa Base.TTY

"""
    wait_for_plot_window(message) -> Bool

Block until the user presses Enter, so an interactive plot window stays open.
Returns `false` without waiting when no window can be shown (headless session)
or when stdin is not a terminal, because there would be nothing to close and
no way to receive the keypress.
"""
function wait_for_plot_window(message::AbstractString)::Bool
    if !display_available()
        println("No display available; skipping the interactive plot window. Plot files are still written to the output directory.")
        return false
    end
    if !stdin_is_terminal()
        println("stdin is not a terminal; not waiting for the interactive plot window.")
        return false
    end
    println(message)
    readline()
    return true
end

end
