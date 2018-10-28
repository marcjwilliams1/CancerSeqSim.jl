"""
    show(sresult::Simulation)

Print out summary of simulation.
"""
function show(io::IO, sresult::Simulation)

  @printf("Input parameters: \n")
  @printf("\t Mutation rate: %.2f\n", sresult.input.μ)
  @printf("\t Death rate of host population: %.2f\n", sresult.input.d)
  @printf("\t Effective mutation rate (μ/β): %.2f\n", sresult.input.μ / ((sresult.input.b-sresult.input.d)/sresult.input.b))
  @printf("\t Number of clonal mutation: %d\n", sresult.input.clonalmutations)

  @printf("\t Number of subclones: %d\n\n", sresult.input.numclones)
  if sresult.input.numclones > 0
    for i in 1:length(sresult.output.clonefreq)
      @printf("Subclone %d \n", i)
      @printf("\tFrequency: %.2f\n", sresult.output.clonefreq[i])
      @printf("\tNumber of mutations in subclone: %d\n", sresult.output.subclonemutations[i])
      @printf("\tFitness advantage: %.2f\n", sresult.input.selection[i])
      @printf("\tTime subclone emerges (population doublings): %.2f\n", log(sresult.output.cloneN[i])/log(2))
      @printf("\tNumber of divisions: %d\n", sresult.output.Ndivisions[i])
      @printf("\tAverage number of divisions per cell: %.2f\n", sresult.output.aveDivisions[i])
      @printf("\tPopulation size when subclone emerges: %d\n", sresult.output.cloneN[i])
      @printf("\tParent of subclone (0 is host): %d\n", sresult.output.clonetype[i])
    end
  else
    @printf("No clones, tumour growth was neutral\n\n")
  end

end

function selection(λ, f, tend, t1)
    #define the equation for selection as above
    s = (λ .* t1 + log.(f ./ (1 - f))) ./ (λ .* (tend - t1))
    return s
end


function selection2clone(λ, f1, f2, tend, t1, t2)
    #define the equation for selection as above

    s1 = zeros(Float64, length(f1))
    s2 = zeros(Float64, length(f1))

    for i in 1:length(f1)
      if (f2[i] .+ f1[i]) < 1.0
        s1[i] = (λ .* t1[i] .+ log.(f1[i] ./ (1 .- f1[i] .- f2[i]))) ./ (λ .* (tend[i] .- t1[i]))
        s2[i] = (λ .* t2[i] .+ log.(f2[i] ./ (1 .- f1[i] .- f2[i]))) ./ (λ .* (tend[i] .- t2[i]))
      else
        s1[i] = (λ .* t1[i] .+ log.((f1[i] .- f2[i]) ./ (1 .- f1[i]))) ./ (λ .* (tend[i] .- t1[i]))
        s2[i] = (λ .* t2[i] .+ log.(f2[i] ./ (1 - f1[i]))) ./ (λ .* (tend[i] - t2[i]))
      end
    end

    return s1, s2
end
