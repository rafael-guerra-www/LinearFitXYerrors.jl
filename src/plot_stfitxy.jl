#TODO: remove Plots.jl dependency and add a Plots recipe

@recipe function f(dt::stfitxy,
            framestyle=:axes,
            yminorgrid=true,
            titefont  = font(6, "Computer Modern"),
            ticfont  = font(6, "Computer Modern"),
            guiefont = font(7, "Computer Modern"),
            legndfont = font(6, "Computer Modern"),
            marershape = :circle,
            marersize = 2,
            marerstrokewidth = 0,
            fg_color_legend = nothing,
            size = (600,400),
            dpi = 300,
            margin = 0mm)

    # Default values
    titlefont  := titlefont
    tickfont  := tickfont
    guidefont := guidefont
    legendfont := legendfont
    markershape := markershape
    markersize := markersize
    markerstrokewidth := markerstrokewidth 
    fg_color_legend := fg_color_legend
    size := size,
    margin := margin,
    dpi := dpi)


    legend := isnothing(legend) ? (size(x,2) <= 4 ? :outertop : :outerright) : legend

    @series begin
        label := permutedims(names_x)
        color_palette := :tab10
        seriescolor := reshape(collect(1:size(x,2)),1,size(x,2))
        seriestype := :line

        ts,x 
    end
end 