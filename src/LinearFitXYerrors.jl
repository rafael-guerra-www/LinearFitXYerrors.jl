module LinearFitXYerrors
  
using Statistics
using Distributions
using Printf
using Measures
using Plots; gr()



plot_font = "Computer Modern";
default(fontfamily=plot_font,framestyle=:axes,yminorgrid=true, legendtitlefontsize=6,fg_color_legend=nothing,
    legendfontsize=6, guidefont=(7,:black),tickfont=(6,:black),size=(600,400),dpi=300, margin=0mm,
    titlefont = (6, plot_font))

include("./linearfitxy.jl")

export linearfitxy
export plot_linfitxy

end 

# TODO: confidence_intervals