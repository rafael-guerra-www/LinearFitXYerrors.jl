module LinearFitXYerrors
  
using Statistics
using Printf
using Measures, Plots; gr()


plot_font = "Computer Modern";
default(fontfamily=plot_font,framestyle=:axes,yminorgrid=true, legendtitlefontsize=6,fg_color_legend=nothing,
    legendfontsize=6, guidefont=(7,:black),tickfont=(6,:black),size=(600,400),dpi=300, margin=0mm,
    titlefont = (6, plot_font))

include("./linearfit_xy_errors.jl")

export linearfit_xy_errors
export plotlinfitxy

end 

# TODO: confidence_intervals