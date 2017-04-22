
function show(sresult::InputAndAnalysis)

  @printf("Input parameters: \n")
  @printf("\t Mutation rate: %.1f\n", sresult.input.μ)
  @printf("\t Death rate of host population: %.1f\n", sresult.input.d)
  @printf("\t Number of clonal mutation: %d\n", sresult.input.clonalmuts)

  @printf("\t Number of clones: %d\n\n", sresult.input.numclones)
  if sresult.input.numclones > 0
    for i in 1:length(sresult.input.pctfit)
      @printf("Subclone %d \n", i)
      @printf("\tFrequency: %.2f\n", sresult.input.pctfit[i])
      @printf("\tNumber of mutations in subclone: %d\n", sresult.input.clonemuts[i])
      @printf("\tFitness advantage: %.2f\n", sresult.input.s[i])
      @printf("\tTime subclone emerges: %.2f\n", sresult.input.tevent[i])
      @printf("\tPopulation size when subclone emerges: %d\n", sresult.input.cloneN[i])
    end
  else
    @printf("No clones, tumour growth was neutral\n\n")
  end

end

function vafhistogram(sresult; annotateclones = false)

    DF = sresult.sampleddata.DF

    if (annotateclones == true) & (sresult.input.numclones > 0)

      xint = sresult.ouput.pctfit./2
      push!(xint, sum(xint))

      p1 = plot(DF, x="VAF", y="freq",
      Guide.xlabel("Allelic Frequency f"),
      Guide.ylabel("Number of Mutations"),
      xintercept = xint,
      Geom.vline(color=colorant"red"),
      Geom.bar,
      Guide.xticks(ticks = collect(0.0:0.2:1.0)),
      Theme(major_label_font_size = 12pt,
      major_label_font = "Arial",
      minor_label_font_size = 10pt,
      minor_label_font = "Arial",
      key_label_font = "Arial",
      key_label_font_size = 10pt,
      default_color = colormap("blues")[80],
      bar_spacing = -0.05cm))

    else

      p1 = plot(DF, x="VAF", y="freq", Geom.bar,
      Guide.xlabel("Allelic Frequency f"),
      Guide.ylabel("Number of Mutations"),
      Guide.xticks(ticks = collect(0.0:0.2:1.0)),
      Theme(major_label_font_size = 12pt,
      major_label_font = "Arial",
      minor_label_font_size = 10pt,
      minor_label_font = "Arial",
      key_label_font = "Arial",
      key_label_font_size = 10pt,
      default_color = colormap("blues")[80],
      bar_spacing = -0.05cm))

    end


    return p1

end

function makelims(fmax)
  func(x) = "1/$(round(1/(x+(1/fmax)), 3))"
end


function cumulativeplot(sresult)

    DF1 = sresult.output.DF
    DF = getsummary(sresult)

    muin = round(DF[:mu][1], 1)
    mufit = round(DF[:muout][1], 1)
    rsq = round(DF[:rsq][1], 3)

    minrange = minimum(DF1[:v])
    maxrange = maximum(DF1[:v])

    func = makelims(maxrange)

    p2=plot(DF1,

    layer(x = "invf", y = "prediction",Geom.line,
    Theme(default_color = default_color=colormap("reds")[80],
    line_width = 0.08cm)),

    layer(x = "invf", y = "cumsum", Geom.line,
    Theme(line_width = 0.1cm,
    default_color = colormap("blues")[80])),
    Scale.x_continuous(minvalue = minimum(DF1[:invf]),
    maxvalue = maximum(DF1[:invf]), labels = func),
    Scale.y_continuous(minvalue = 0,
    maxvalue = maximum(DF1[:cumsum])),


    Theme(major_label_font_size = 12pt,
    major_label_font = "Arial",
    minor_label_font_size = 10pt,
    minor_label_font = "Arial",
    key_label_font = "Arial",
    key_label_font_size = 10pt),
    Guide.manual_color_key("",
    ["Simulated Data", "Model Fit"],
    [color(colormap("blues")[80]), color(colormap("reds")[80])]),
    Guide.annotation(compose(context(),text(0, maximum(DF1[:cumsum]), "R<sup>2</sup> = $(rsq)"))),
    Guide.annotation(compose(context(),text(0, 0.9*maximum(DF1[:cumsum]), "μ<sub>in</sub> = $(muin)"))),
    Guide.annotation(compose(context(),text(0, 0.8*maximum(DF1[:cumsum]), "μ<sub>fit</sub> = $(mufit)"))),
    Guide.xlabel("Inverse Allelic Frequency 1/f"),
    Guide.ylabel("Cumulative # of Mutations M(f)"),
    Guide.xticks(ticks = [1/maxrange - 1/maxrange, 1/minrange - 1/maxrange]))

end
