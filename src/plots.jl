@recipe function f(sresults::Simulation)
    DF = sresults.sampleddata.DF
    VAF = DF[:VAF]
    freq = DF[:freq]

    if length(sresults.output.clonefreq) > 0
        xint = sresults.output.clonefreq./2 * sresults.input.cellularity

        @series begin
            seriestype --> :bar
            linecolor --> :darkslategrey
            fillcolor --> :darkslategrey
            markerstrokecolor --> :white
            legend --> false
            grid --> false
            VAF, freq
        end

        @series begin
            yaxis --> "Counts"
            xaxis --> "VAF"
            seriestype --> :vline
            legend --> false
            fillcolor --> :darkred
            linewidth --> 3
            xint
        end
    else
        @series begin
            seriestype --> :bar
            yaxis --> "Counts"
            xaxis --> "VAF"
            linecolor --> :darkslategrey
            fillcolor --> :darkslategrey
            markerstrokecolor --> :white
            legend --> false
            grid --> false
            VAF, freq
        end
    end
end
