# src/LinearFitXYerrors.jl — module file
# - NO Plots dependency at module level anymore (`using Plots; gr()` removed)
# - stfitxy is now exported (needed for incoming methods / recipes dispatch)
# - isplot=true goes through a Base.get_extension bridge: it works only when
#   the user has loaded Plots (the extension is then active)

module LinearFitXYerrors

using Statistics
using Distributions
using Printf
using RecipesBase          # backend-free: carries the `plot(st)` recipe

include("linearfitxy.jl")   # fit code (Plots-free)
include("recipes.jl")       # @recipe for stfitxy (Plots renders it via RecipesBase)

export linearfitxy, stfitxy   # NEW: stfitxy exported (dispatch target for recipes)

# Bridge for the legacy `isplot` kwarg of linearfitxy: displays the recipe
# through the Plots extension when available.
function _isplay(st::stfitxy; ratio=1)
    ext = Base.get_extension(@__MODULE__, :LinearFitXYerrorsPlotsExt)
    if ext === nothing
        @warn "`isplot=true` requires Plots.jl. Do `using Plots` (before or after " *
              "`using LinearFitXYerrors`) or call `plot(st)` directly. Returning fit only."
        return nothing
    end
    return ext.plot_linfitxy(st; ratio = ratio)
end

end
