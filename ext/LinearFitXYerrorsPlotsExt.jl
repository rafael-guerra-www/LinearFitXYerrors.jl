# ext/LinearFitXYerrorsPlotsExt.jl — package extension (loaded only if Plots.jl is)
# Contains everything that genuinely needs the Plots backend:
# the display() convenience and the legacy `plot_linfitxy` entry point.
# All graphics themselves live in the backend-free recipe (src/recipes.jl),
# so the core package stays loadable without Plots, and `plot(st)` gains
# Makie support later via a second extension (separate, incompatible API).

module LinearFitXYerrorsPlotsExt

using LinearFitXYerrors: stfitxy
using Plots
using Printf

"""
    plot_linfitxy(st::stfitxy; ratio=1)

Legacy convenience: plot and display a `linearfitxy` result.
Equivalent to `plot(st)` from the backend-free recipe.

`ratio` sets the plot aspect ratio (1 for equal axes, `:auto` otherwise).
"""
function plot_linfitxy(st::stfitxy; ratio = 1)
    p = plot(st; ratio = ratio)
    display(p)
    return p
end

end  # module
