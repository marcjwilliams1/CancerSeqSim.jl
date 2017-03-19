###############################################################################

#Call required packages
using Distributions
using DataFrames
using GLM
using PyCall
using Stats
using HypothesisTests
using Gadfly
###############################################################################
#call samplesim.jl

include("samplesim.jl")
metricp=readtable("julia/metricprobabilities.csv")

###############################################################################
srand(1)
#type definitions

type Priors

    numclones::Array{Int64,1}
    clonalmuts::Array{Int64,1}
    selection::Array{Float64,1}
    μ::Array{Float64,1}
    tevent::Array{Float64,1}
    d::Array{Float64, 1}

end

type Constants

    Nmax::Int64
    det_limit::Float64
    ploidy::Int64
    read_depth::Float64
    fmin::Float64
    fmax::Float64
    b::Float64
    ρ::Float64

end

type PertScale

    clonalmuts::Float64
    selection::Array{Float64, 1}
    μ::Float64
    tevent::Array{Float64, 1}
    d::Float64

end

type PertScaleNeutral

    clonalmuts::Float64
    μ::Float64
    d::Float64

end


DFABC=readtable("/data/home/mpx155/scripts/early_events/ABC-SMCnew/simulations/data/selection_2clone-2_depth.300.csv");
dataout="/data/BCI-EvoCa/marc/early_events/ABC-SMC/samples/selection_highdepth_lowbeta2"

cst = Constants(10^3, 0.017, 2, 300.0, 0.017, 0.75, log(2.0), 0.0 )
Pr = Priors([0,2], [1,1000], [0.0,10.0], [10.0,100.001],[3.0, 8.0], [0.0, 0.0])
ϵ = [1000.0, 0.1]
N=50
runabcsmc(Pr, cst, N, ϵ, metricp, DFABC, dataout, "test", "test1")

DFABC=readtable("/data/home/mpx155/scripts/early_events/ABC-SMCnew/simulations/data/neutral.300.csv");
dataout="/data/BCI-EvoCa/marc/early_events/ABC-SMC/samples/selection_highdepth_lowbeta2"

cst = Constants(10^3, 0.017, 2, 300.0, 0.017, 0.75, log(2.0), 0.0 )
Pr = Priors([0,2], [1,1000], [0.0,10.0], [10.0,100.001],[3.0, 8.0], [0.0, 0.0])
ϵ = [1000.0, 0.1]
N=50
runabcsmc(Pr, cst, N, ϵ, metricp, DFABC, dataout, "test", "test2")
#DFresults, nattempts, adist = runabcsmc(Pr, cst, N, ϵ, metricp, DFABC, "test-things/testfiles", "1","test");
#runabcsmc1(Pr cst, N, ϵ, metricp, DFABC)
# p, nattempts, kdist = runabcsmc1(Pr, cst, N, ϵ, metricp, DFABC);
###############################################################################

function runabcsmc1(P::Priors, cst::Constants, N, ϵ, metricp, DFABC)

    i = 1

    #initialize array of particles to store particles, where the array corresponds to the model indicator
    particles = createparticlearray(P)

    currmodels = createmodelarray(P)

    kdist = Float64[]
    numsims = 0
    sresult = 0

    while i < N + 1
        numsims = numsims + 1

        nc = rand(P.numclones[1] : P.numclones[2])
        correctnc = false

        while correctnc == false

            tevent = sort(rand(Uniform(P.tevent[1], P.tevent[2]), nc))
            clonalmuts = rand(P.clonalmuts[1] : P.clonalmuts[2])
            s = rand(Uniform(P.selection[1], P.selection[2]), nc)
            μ = rand(Uniform(P.μ[1], P.μ[2]))
            if P.d[1] == P.d[2]
                d = P.d[1]
            else
                d = rand(Uniform(P.d[1], P.d[2]))
            end

            IP = InputParameters(nc,
            cst.Nmax,
            cst.det_limit,
            cst.ploidy,
            cst.read_depth,
            cst.fmin,
            cst.fmax,
            clonalmuts,
            s,
            μ,
            cst.b,
            d,
            tevent,
            cst.ρ)

            sresult = finalresults(metricp, DFABC, 10, IP)
            nctemp = sresult.input.numclones
            if nctemp == nc
                correctnc = true
            end
        end

        if sresult.dist < ϵ
            push!(particles[nc + 1], sresult)
            push!(kdist, sresult.dist)
            i = i + 1
        end

    end

    return particles, numsims, kdist

end

function calculatescale(particles, Pr)

    scaletevent = Float64[]
    for i in 1:length(particles[1].parameters.tevent)
        tevents = map(x -> x.parameters.tevent[i], particles)
        scaleteventtemp = (maximum(tevents) - minimum(tevents)) / 2
        if scaleteventtemp == 0.0
            scaleteventtemp = (Pr.tevent[2] - Pr.tevent[1]) / 2
        end
        push!(scaletevent, scaleteventtemp)
    end

    cms = map(x -> x.parameters.clonalmuts, particles)
    scalecm = (maximum(cms) - minimum(cms)) / 2
    if scalecm == 0.0
        scalecm = (Pr.clonalmuts[2] - Pr.clonalmuts[1]) / 2
    end

    scales = Float64[]
    for i in 1:length(particles[1].parameters.selection)
        ss = vcat(map(x -> x.parameters.selection[i], particles)...)
        scalestemp = (maximum(ss) - minimum(ss)) / 2
        if scalestemp == 0.0
            scalestemp = (Pr.selection[2] - Pr.selection[1]) / 2
        end
        push!(scales, scalestemp)
    end

    mus = map(x -> x.parameters.μ, particles)
    scalemu = (maximum(mus) - minimum(mus)) / 2
    if scalemu == 0.0
        scalemu = (Pr.μ[2] - Pr.μ[1]) / 2
    end

    ds = map(x -> x.parameters.d, particles)
    scaled = (maximum(ds) - minimum(ds)) / 2
    if scaled == 0.0
        scaled = (Pr.d[2] - Pr.d[1]) / 2
    end

    return PertScale(scalecm, scales, scalemu, scaletevent, scaled)
end

function calculatescaleneutral(particles, Pr)

    cms = map(x -> x.parameters.clonalmuts, particles)
    scalecm = (maximum(cms) - minimum(cms)) / 2
    if scalecm == 0.0
        scalecm = (Pr.clonalmuts[2] - Pr.clonalmuts[1]) / 2
    end

    mus = map(x -> x.parameters.μ, particles)
    scalemu = (maximum(mus) - minimum(mus)) / 2
    if scalemu == 0.0
        scalemu = (Pr.μ[2] - Pr.μ[1]) / 2
    end

    ds = map(x -> x.parameters.d, particles)
    scaled = (maximum(ds) - minimum(ds)) / 2
    if scaled == 0.0
        scaled = (Pr.d[2] - Pr.d[1]) / 2
    end

    return PertScaleNeutral(scalecm, scalemu, scaled)
end

function kernel_pert(x0, sc, prior)
    if prior[1] == prior[2]
        return x0
    end
    done = false
    while done == false
        x = rand(Uniform(x0 - sc, x0 + sc))
        if (x >= prior[1]) & (x <= prior[2])
            done = true
        return x
        end
    end
end

function kernel_pert2(x0, sc, prior)
    if prior[1] == prior[2]
        return x0
    end
    done = false
    while done == false
        x = int64(rand(Uniform(x0 - sc, x0 + sc)))
        if (x >= prior[1]) & (x <= prior[2])
            done = true
        return x
        end
    end
end

function newinputs(p, scales, prior)

    tevent = zeros(Float64, length(p.parameters.tevent))
    for i in 1:length(p.parameters.tevent)
        tevent[i] = kernel_pert(p.parameters.tevent[i], scales.tevent[i], prior.tevent)
    end

    tevent = sort(tevent)

    s = zeros(Float64, length(p.parameters.selection))
    for i in 1:length(p.parameters.selection)
        s[i] = kernel_pert(p.parameters.selection[i], scales.selection[i], prior.selection)
    end

    μ = kernel_pert(p.parameters.μ, scales.μ, prior.μ)
    d = kernel_pert(p.parameters.d, scales.d, prior.d)
    clonalmuts = kernel_pert2(p.parameters.clonalmuts, scales.clonalmuts, prior.clonalmuts)

    return tevent, s, μ, clonalmuts, d
end

function compute_weights(N, weightsprev, currpar, prevpar, scalesarray, model, survivingmodels, mprob, mfreq)

    N = length(weightsprev)
    denomparticle = zeros(Float64, length(weightsprev))
    denommodel = zeros(Float64, length(weightsprev))
    denommodelmargins = zeros(Float64, length(weightsprev))

    for i in 1:N
        for j in 1:N
            #println("it $j")
            #x = wsample([0, 1, 2], mfreq)
            #println("prev model")
            #println(x)
            #println("curr model")
            #println(model)
            #if model == x
              denomparticle[i] += weightsprev[i] * kernel_prob(currpar[i, :], prevpar[j, :], scalesarray)
            #end

            #println("denomm")
            #println(denomparticle[i])
        end
        #denomparticle[i] = denomparticle[i]./mprob[model + 1]
    end

    for i in 1:N
        denommodelmargins[i] = mprob[model + 1]
    end

    for i in 1:N
        for j in survivingmodels
            denommodel[i] += mprob[j + 1] * getmodelprob(model, j, survivingmodels)
        end
    end

    return denomparticle, denommodel, denommodelmargins
end

function compute_weightsmodel(N, weightsprev, currpar, prevpar, scales, survivingmodels, mprob, pprob, mfreq)

    weights = Array[]
    km = 1

    for k in 1:length(weightsprev)

        if weightsprev[k] == []
            push!(weights, [])
            continue
        end

        model = survivingmodels[km]


        if (survivingmodels[1] == 0) & (k == 1)
            scalesarray = [scales[k].μ, scales[k].clonalmuts, scales[k].d]
        else
            scalesarray = [scales[k].tevent..., scales[k].selection..., scales[k].μ, scales[k].clonalmuts, scales[k].d]
        end

        denomparticle, denommodel, denommodelmargins = compute_weights(N, weightsprev[k], currpar[k], prevpar[k], scalesarray, model, survivingmodels, mprob, mfreq)

        println("particle prior")
        println(pprob[km])
        println("")

        println("denom model")
        println(denommodel)
        println("")

        println("denom particle")
        println(denomparticle)
        println("")

        println("denom model marginal")
        println(denommodelmargins)
        println("")

        W = pprob[km] ./ (denommodel .* (denomparticle ./ denommodelmargins))

        push!(weights, W)

        km = km + 1
    end

    println("weights")
    println(weights)

    sumw = 0.0
    for i in weights
      if i == []
        continue
      end
      sumw = sumw + sum(i)
    end

    for i in 1:length(weights)
      if weights[i] == []
        continue
      end
      weights[i] = weights[i] ./ sumw
    end

    return weights
end

function modelmarginal(weights)

  mmarginal = Float64[]

  for w in weights
    if w == []
        push!(mmarginal, 0.0)
    else
      push!(mmarginal, sum(w))
    end
  end

  return mmarginal

end

function kernel_probold(currpar, prevpar, scalesarray)

    for i in 1:size(currpar, 2)
        if abs(currpar[i] - prevpar[i]) > scalesarray[i]
            return 0
        end
    end

    return 1.0
end

function kernel_prob(currpar, prevpar, scalesarray)

    prob = 1

    for i in 1:size(currpar, 2)
      if (prevpar[i] - scalesarray[i]) < (prevpar[i] + scalesarray[i])

        prob = prob * pdf(Uniform(prevpar[i] - scalesarray[i], prevpar[i] + scalesarray[i]), currpar[i])

      end
    end

    if prob == 0.0
      prob = 10^-5.0
    end

    return prob
end

function appendpararray(p, currpar, prevpar, tevent, s, μ, d, clonalmuts)

    currpar = vcat(currpar, [tevent..., s..., μ, clonalmuts, d]')
    prevpar = vcat(prevpar,[p.parameters.tevent..., p.parameters.selection..., p.parameters.μ, p.parameters.clonalmuts, p.parameters.d]')

    return currpar, prevpar
end


function appendpararray(p, currpar, prevpar, tevent, s, μ, clonalmuts, d, i)

    idx = 1
    for k in 1:length(p.parameters.tevent)
        prevpar[i, idx] = float64(p.parameters.tevent[k])
        currpar[i, idx] = float64(tevent[k])
        idx = idx + 1
    end

    for k in 1:length(p.parameters.selection)
        prevpar[i, idx] = p.parameters.selection[k]
        currpar[i, idx] = s[k]
        idx = idx + 1
    end

    prevpar[i, idx] = p.parameters.μ
    currpar[i, idx] = μ
    idx = idx + 1

    prevpar[i, idx] = float64(p.parameters.clonalmuts)
    currpar[i, idx] = float64(clonalmuts)
    idx = idx + 1

    prevpar[i, idx] = float64(p.parameters.d)
    currpar[i, idx] = float64(d)

    return currpar, prevpar
end

function convertpart2resultsmodel(particles, weights)

    DF = convertpart2results(particles[1], weights[1])

    if length(particles) > 1
        for i in 2:length(particles)
            DF1 = convertpart2results(particles[i], weights[i])
            append!(DF, DF1)
        end
    end

    return DF
end

function convertpart2results(particles, weights)

    DF, AD = getallmetrics(particles[1])
    DF[:draw] = particles[1].dist
    DF[:weight] = weights[1]

    k = 2
    for p in particles[2:end]

        DF1, AD = getallmetrics(p)
        DF1[:draw] = p.dist
        DF1[:weight] = weights[k]
        append!(DF, DF1)
        k = k + 1

    end

    return DF
end

function averagehistogram(particles, N)

    p = InputAndAnalysis[]
    for i in particles
        append!(p, i)
    end

    M = zeros(Int64, 100, N)
    i = 1

    for j in 1:N

        M[:, i] = freqcounts = array(p[i].sampleddata.DF[:freq])
        i = i + 1

    end

    return collect(median(M,2))
end

function averagehistogram2(particles)

    if particles == []
        return 0
    end

    M = zeros(Int64, 100, length(particles))
    i = 1

    for j in 1:length(particles)

        M[:, i] = freqcounts = array(particles[i].sampleddata.DF[:freq])
        i = i + 1

    end

    return collect(median(M,2))
end

function modelprob(particles, N)

    #calculate probabilities of selecting models based on previous population

    probs = Float64[]
    for i in particles
        push!(probs, length(i) / N)
    end

    return probs
end

function createparticlearray(P)

    particle = Array[]
    for k in 1:(maximum(P.numclones) + 1)
        push!(particle, InputAndAnalysis[])
    end

    return particle
end

function createppprobarray(P)

    pprob = Array[]
    for k in 1:(maximum(P.numclones) + 1)
        push!(pprob, Float64[])
    end

    return pprob
end

function createmodelarray(P)

  models = Array[]
  for k in 1:(maximum(P.numclones) + 1)
      push!(models, Int64[])
  end

  return models
end

function getsurvivingmodels(particles)

    survivingmodels = Int64[]

    for i in particles
        if i != []
            push!(survivingmodels, i[1].input.numclones)
        end
    end

    return sort(survivingmodels)
end



function returnscales(particles, Pr, survivingmodels)

    scales = Any[]

    for i in 1:length(particles)
        if particles[i] == []
            push!(scales,[])
            continue
        end
        if i == 1
            push!(scales, calculatescaleneutral(particles[i], Pr))
        else
            push!(scales, calculatescale(particles[i], Pr))
        end
    end

    return scales
end

function prev_curr_parameters(scales, Pr)

    prevpar = Array[]
    currpar = Array[]

    for i in Pr.numclones[1] : Pr.numclones[2]
        push!(prevpar, zeros(Float64, 0, 2i + 3))
        push!(currpar, zeros(Float64, 0, 2i + 3))
    end

    return prevpar, currpar
end

function createweightsprev(P)

    weights = Array[]
    for k in 1:(maximum(P.numclones) + 1)
        push!(weights, Float64[])
    end

    return weights
end

function printresults(DFresults, scales, Pr)

    println(scales)

    println("")

    numclones = sort(unique(DFresults[:numclones]))
    N = length(DFresults[:s])
    println(numclones)

    for i in numclones
        println("Parameters for model with $i clones")
        DF = DFresults[DFresults[:numclones].==i, :]
        println("\tProbability of having $i clones = $(length(DF[:s]) / N)")
        println("\tProbability of having $i clones (marginal) = $(sum(DF[:weight]))")
        println("\tMedian clonal mutations - $(median(DF[:clonalmuts]))")
        println("\tMedian mutation rate - $(median(float64(DF[:mu])))")
        println("\tMedian death rate - $(median(float64(DF[:d])))")

        timeevent = Float64[]
        sevent = Float64[]
        clonemuts = Float64[]
        pctfit = Float64[]

        for k in 1:i

            push!(timeevent, median(map(x -> float64(split(x,",")[k]), DF[:clonetime])))
            push!(sevent, median(map(x -> float64(split(x,",")[k]), DF[:s])))
            push!(clonemuts, median(map(x -> float64(split(x,",")[k]), DF[:clonemuts])))
            push!(pctfit, median(map(x -> float64(split(x,",")[k]), DF[:pctfit])))

        end

        for k in 1:length(timeevent)
            println("\tMedian time of mutations clone $(k): $(timeevent[k])")
            println("\tMedian number of mutations in clone $(k): $(clonemuts[k])")
            println("\tMedian % fit, clone $(k): $(pctfit[k])")
            println("\tMedian selection of clone $(k): $(sevent[k])")

        end

        println("\tScales:")
        println("\t\t$(scales[i + 1])")


        println("")

    end
end


function printresultslog(DFresults, scales, Pr, eps, sname, numsims, dist)

    f = open("log/$(sname).txt","a+")

    write(f, "\n")
    write(f, "######################################\n")
    write(f, "\n")
    write(f, "Finished population with tolerence eps = $(eps)\n")

    write(f, "\n")
    write(f, "Numbers of simulations = $(numsims)\n")
    write(f, "Median distance of all simulations = $(median(dist))\n")
    write(f, "Next epsillon = $(quantile(dist, 0.3))\n")
    write(f, "\n")
    write(f, "Scales:\n")
    write(f, "\t$(scales)\n")
    write(f, "\n")

    numclones = sort(unique(DFresults[:numclones]))
    N = length(DFresults[:s])

    for i in numclones
        write(f, "Parameters for model with $i clones\n")
        DF = DFresults[DFresults[:numclones].==i, :]
        write(f, "\tProbability of having $i clones = $(length(DF[:s]) / N)\n")
        write(f, "\tProbability of having $i clones (marginal) = $(sum(DF[:weight]))\n")
        write(f, "\tMedian clonal mutations - $(median(DF[:clonalmuts]))\n")
        write(f, "\tMedian mutation rate - $(median(float64(DF[:mu])))\n")
        write(f, "\tMedian death rate - $(median(float64(DF[:d])))\n")

        timeevent = Float64[]
        sevent = Float64[]
        clonemuts = Float64[]
        pctfit = Float64[]
        deathrate = Float64[]

        for k in 1:i

            push!(timeevent, median(map(x -> float64(split(x,",")[k]), DF[:clonetime])))
            push!(sevent, median(map(x -> float64(split(x,",")[k]), DF[:s])))
            push!(clonemuts, median(map(x -> float64(split(x,",")[k]), DF[:clonemuts])))
            push!(pctfit, median(map(x -> float64(split(x,",")[k]), DF[:pctfit])))
            push!(deathrate, median(map(x -> float64(split(x,",")[k + 1]), DF[:dr])))

        end

        for k in 1:length(timeevent)
            write(f, "\tMedian time of mutations clone $(k): $(timeevent[k])\n")
            write(f, "\tMedian number of mutations in clone $(k): $(clonemuts[k])\n")
            write(f, "\tMedian % fit, clone $(k): $(pctfit[k])\n")
            write(f, "\tMedian selection of clone $(k): $(sevent[k])\n")
            write(f, "\tMedian deathrate of clone $(k): $(deathrate[k])\n")

        end

    end


    close(f)
end

function modelkernel(currm, survivingmodels)

    prob = 0.7

    mprob = ones(Float64, length(survivingmodels))

    mprob[findin(survivingmodels, currm)] = prob
    mprob[findin(mprob, 1.0)] = (1 - prob) / (length(survivingmodels) - 1)

    wsample(survivingmodels, mprob)

end

function getmodelprob(currmodel, prevmodel, survivingmodels)

  prob = 0.7

  if currmodel == prevmodel
    return prob
  else
    return (1 - prob) / (length(survivingmodels) - 1)
  end

end

function writesimtodisk(particles, survivingmodels, dataout, date, sname,eps)

    for i in survivingmodels
        distances = map(x -> x.dist, particles[i + 1])
        ind = findmin(distances)[2]

        writetable("$(dataout)/individualhist/$(eps).$(particles[i + 1][ind].dist).$(date).$(sname).model$(i).onesimhist.csv", particles[i + 1][ind].sampleddata
        .DF)

        writedlm("$(dataout)/individualhist/$(eps).$(particles[i + 1][ind].dist).$(date).$(sname).model$(i).VAF.txt", particles[i + 1][ind].sampleddata
        .VAF)

    end

end

function priorprob(Pr, tevent, s, μ, clonalmuts, d)

  pprob = 1
  for i in tevent
    pprob = pprob * pdf(Uniform(Pr.tevent[1], Pr.tevent[2]), i)
  end

  for i in s
    pprob = pprob * pdf(Uniform(Pr.selection[1], Pr.selection[2]), i)
  end

  if Pr.μ[1] != Pr.μ[2]
    pprob = pprob * pdf(Uniform(Pr.μ[1], Pr.μ[2]), μ)
  end

  if Pr.clonalmuts[1] != Pr.clonalmuts[2]
    pprob = pprob * pdf(Uniform(Pr.clonalmuts[1], Pr.clonalmuts[2]), clonalmuts)
  end

  if Pr.d[1] != Pr.d[2]
    pprob = pprob * pdf(Uniform(Pr.d[1], Pr.d[2]), d)
  end

  return pprob

end


function runabcsmc(P::Priors, cst::Constants, N, ϵ, metricp, DFABC, dataout, date, sname)

    # run first pass with
    particles, numsims, kdist = runabcsmc1(P, cst, N, ϵ[1], metricp, DFABC)


    survivingmodels = getsurvivingmodels(particles)
    ptemp = Array[]
    for i in survivingmodels
       push!(ptemp, particles[i+1])
    end

    #keep track of number of simulations
    numberattempts = Float64[]
    push!(numberattempts, numsims)

    #keep array of the average eps per populations
    averagedist = Float64[]
    push!(averagedist, mean(kdist))

    #weight of first population are uniform
    weights = Array[]
    for i in particles
        push!(weights, ones(Float64,length(i)) ./ N)
    end

    DFresults = convertpart2resultsmodel(ptemp, weights)
    #calculate model selection probabilities
    mprob = modelmarginal(weights)
    mfreq = modelprob(particles, N)

    #value of epsilon is adaptive, we take the αth quantile
    α = 0.3

    writetable("$(dataout)/eps/$(date).$(sname).$(round(ϵ[1],4)).txt", DFresults, separator = '\t')

    println("Finished population with tolerence ϵ = $(ϵ[1])")
    scales = returnscales(particles, P, survivingmodels)
    printresults(DFresults, scales, P)

    printresultslog(DFresults, scales, P, ϵ[1], sname, numsims, kdist)

    eps = ϵ[1]
    saveeps = Float64[]


    while eps > ϵ[2]

        eps = quantile(kdist, α)
        push!(saveeps, eps)

        kattempts = 0

        #initalize new array to store next set of particles
        #initialize array of particles to store particles, where the array corresponds to the model indicator
        particlesnew = createparticlearray(P)

        modelscurr = createmodelarray(P)

        pprobarray = createppprobarray(P)

        kdist = Float64[]

        scales = returnscales(particles, P, survivingmodels)

        if eps < ϵ[2]
            eps = ϵ[2]
        end

        i = 1
        prevpar, currpar = prev_curr_parameters(scales, P)

        weightsprev = createweightsprev(P)

        tevent, s, μ, clonalmuts, d = 0.0,0.0,0.0,0.0,0.0

        pprob1 = 0.0

        while i < N + 1

            #keep track of number of simulations run
            kattempts = kattempts + 1

            #choose model based on previous population models
            sresult = 0
            mstar = wsample((P.numclones[1] : P.numclones[2]), mprob)
            m = modelkernel(mstar, survivingmodels) + 1
            ind = wsample(1 : length(particles[m]), weights[m])
            p = particles[m][ind]

            correctnc = false
            while correctnc == false

                tevent, s, μ, clonalmuts, d = newinputs(p, scales[m], P)
                pprob1 = priorprob(Pr, tevent, s, μ, clonalmuts, d)

                IP = InputParameters(m - 1,
                cst.Nmax,
                cst.det_limit,
                cst.ploidy,
                cst.read_depth,
                cst.fmin,
                cst.fmax,
                clonalmuts,
                s,
                μ,
                cst.b,
                d,
                tevent,
                cst.ρ)

                sresult = finalresults(metricp, DFABC, 10, IP)
                nctemp = sresult.parameters.numclones
                if nctemp == m - 1
                    correctnc = true
                end
            end

            if pprob1 == 0.0
              continue
            end

            if sresult.dist <= eps

                #println(m)

                push!(particlesnew[m], sresult)
                push!(kdist, sresult.dist)

                currpar[m], prevpar[m] = appendpararray(p, currpar[m], prevpar[m], tevent, s, μ, d, clonalmuts)

                push!(weightsprev[m], weights[m][ind])

                push!(pprobarray[m], pprob1)

                i = i + 1

            end

        end

        push!(numberattempts, kattempts)
        particles = particlesnew
        push!(averagedist , mean(kdist))

        survivingmodels = getsurvivingmodels(particles)


        println("Finished population with tolerence, ϵ = $(eps), number of simulations = $kattempts")
        println("")

        #calculate weights
        #weights = compute_weightsmodel(N, weightsprev, currpar, prevpar, scales, survivingmodels)

        weights = compute_weightsmodel(N, weightsprev, currpar, prevpar, scales, survivingmodels, mprob, pprobarray, mfreq)
        scales = returnscales(particles, P, survivingmodels)
        mprob = modelmarginal(weights)
        mfreq = modelprob(particles, N)

        ptemp = Array[]
        wtemp = Array[]
        for i in survivingmodels
           push!(ptemp, particlesnew[i+1])
           push!(wtemp, weights[i + 1])
        end

        DFresults = convertpart2resultsmodel(ptemp, wtemp)

        printresults(DFresults, scales, P)
        printresultslog(DFresults, scales, P, eps, sname, kattempts, kdist)


        #write results dataframe to file
        writetable("$(dataout)/eps/$(date).$(sname).$(round(eps,4)).txt", DFresults, separator = '\t')
        cp("$(dataout)/eps/$(date).$(sname).$(round(eps,4)).txt",
        "$(dataout)/final/$(sname).txt")

        #write average histogram to file
        writedlm("$(dataout)/eps/$(date).$(sname).$(round(eps,4)).averagehist.txt", averagehistogram(particles, N))
        cp("$(dataout)/eps/$(date).$(sname).$(round(eps,4)).averagehist.txt",
        "$(dataout)/final/$(sname).averagehist.txt")

        for d in survivingmodels
            writedlm("$(dataout)/avehist/$(date).$(sname).$(round(eps,4)).clone$(d).averagehist.txt", averagehistogram2(particles[d + 1]))
            cp("$(dataout)/avehist/$(date).$(sname).$(round(eps,4)).clone$(d).averagehist.txt",
            "$(dataout)/final/$(sname).clone$(d).averagehist.txt")
        end

        #write best models to disk
        writesimtodisk(particles, survivingmodels, dataout, date, sname, round(eps,4))

        epsnew = quantile(kdist, α)
        if abs(epsnew - eps)/eps < 0.01
          f = open("log/$(sname).txt","a+")
          write(f, "% change in eps < 0.5%, end ABC-SMC \t")
          close(f)
          break
        end

    end

    survivingmodels = getsurvivingmodels(particles)
    ptemp = Array[]
    wtemp = Array[]
    for i in survivingmodels
       push!(ptemp, particles[i + 1])
       push!(wtemp, weights[i + 1])
    end
    DFresults = convertpart2resultsmodel(ptemp, wtemp)

    #write best models to disk
    writesimtodisk(particles, survivingmodels, dataout, date, sname, round(eps,4))

    writedlm("$(dataout)/final/$sname.tolerances.txt", saveeps)
    writedlm("$(dataout)/final/$sname.numbersims.txt", numberattempts)

    return DFresults, numberattempts, averagedist

end
