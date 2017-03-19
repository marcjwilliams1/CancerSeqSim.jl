###############################################################################

#Call required packages
using Distributions
using DataFrames
using GLM
#using PyCall
using Stats
using HypothesisTests
using Gadfly
###############################################################################
#call samplesim.jl

include("samplesim.jl")
metricp=readtable("/data/home/mpx155/scripts/early_events/ABC-SMCnew/julia/metricprobabilities.csv")

###############################################################################
#srand(1)
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

type Particle

  nparams::Int64
  modelcurr::Int64
  modelprev::Int64
  weightcurr::Float64
  weightprev::Float64
  paramscurr::Array{Float64, 1}
  paramsprev::Array{Float64, 1}
  pertscale::Array{Float64, 1}
  dist::Float64
  IPA::InputAndAnalysis
  marginscurr::Array{Float64,1}
  marginsprev::Array{Float64,1}
  priorprob::Float64
  priorArray::Array{Float64, 2}

end
#runabcsmcpop1(P, cst, N, 10000, metricp, DFABC)
###############################################################################

#runabcsmcpop1(Pr, cst, N, 10000.0, metricp, DFABC)
function runabcsmcpop1(P::Priors, cst::Constants, N, ϵ, metricp, DFABC, distweight)

    i = 1

    #initialize array of particles to store particles, where the array corresponds to the model indicator
    particles = Array(Particle, N)

    kdist = Float64[]
    numsims = 0
    sresult = 0

    tevent, s, μ, clonalmuts, d = 0.0,0.0,0.0,0.0,0.0

    priorarray = makepriorarray(Pr)

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

            sresult = finalresults(metricp, DFABC, 10, IP, distweight)
            nctemp = sresult.input.numclones
            if nctemp == nc
                correctnc = true
            end
        end
        #println(sresult.dist)

        if sresult.dist < ϵ

            pprob1 = priorprob(Pr, tevent, s, μ, clonalmuts, d)

            p = Particle(3 + 2*nc, #nparams
            nc, #modelcurr
            0, #model_prev
            1/N, #weightcurr
            1/N, #weightprev
            [clonalmuts, μ, d, [tevent s]...], #paramscurr
            [clonalmuts, μ, d, [tevent s]...], #paramsprev
            [0.0, 0.0], #pertscale
            sresult.dist, #distance
            sresult, #simulation output
            [0.0], #marginals curr
            [0.0], #marginals prev
            pprob1, #prior probability
            priorarray #priorarray
            ) #simulation output

            particles[i] = p
            push!(kdist, sresult.dist)
            i = i + 1
        end

    end

    return particles, numsims, kdist, distweight

end



function getvariability(P::Priors, cst::Constants, N, ϵ, metricp, DFABC)
    srand(5)
    i = 1

    #initialize array of particles to store particles, where the array corresponds to the model indicator
    particles = Array(Particle, N)

    kdist = Float64[]
    numsims = 0
    sresult = 0

    tevent, s, μ, clonalmuts, d = 0.0,0.0,0.0,0.0,0.0

    priorarray = makepriorarray(Pr)
    Nsims = 10
    M = zeros(length(DFABC[:cumsum]), Nsims)

    for i in 1:Nsims
        numsims = numsims + 1

        #x = cumulativedist(sampledhist(50*ones(cst.Nmax), cst.Nmax, det_limit = cst.det_limit, ploidy = 2, read_depth = cst.read_depth).VAF, fmin = cst.fmin, fmax = cst.fmax)

        #M[:, i] = map(Float64, x.DF[:cumsum])
        #M[:, i] = float64(sresult.output.DF[:normalized])

    end

    MAD = median(abs(M .- median(M, 2)), 2)
    #MAD = ones(Float64, length(MAD)) - MAD./sum(MAD)
    MAD = ones(Float64, length(MAD))

    writedlm("/data/home/mpx155/scripts/early_events/ABC-SMCnew/weights/weight.$(cst.fmin).$(cst.fmax).$(cst.read_depth).txt", MAD)

    return MAD
end



function getvariability2(P::Priors, cst::Constants, N, ϵ, metricp, DFABC)

    i = 1

    #initialize array of particles to store particles, where the array corresponds to the model indicator
    particles = Array(Particle, N)

    kdist = Float64[]
    numsims = 0
    sresult = 0

    tevent, s, μ, clonalmuts, d = 0.0,0.0,0.0,0.0,0.0

    priorarray = makepriorarray(Pr)
    Nsims = 5000
    M = zeros(length(DFABC[:cumsum]), Nsims)

    for i in 1:Nsims
        numsims = numsims + 1

        nc = rand(0 : 2)
        #nc = 1
        correctnc = false

        while correctnc == false

            tevent = sort(rand(Uniform(P.tevent[1], P.tevent[2]), nc))
            #clonalmuts = int64(round((P.clonalmuts[2] - P.clonalmuts[1])/2))
            clonalmuts = 100
            #clonalmuts = int64(rand(Uniform(P.clonalmuts[1], P.clonalmuts[2])))
            s = rand(Uniform(P.selection[1], P.selection[2]), nc)
            #μ = int64(round((P.μ[2] - P.μ[1])/2))
            μ =20.0
	    #μ = rand(Uniform(P.μ[1], P.μ[2]))
            if P.d[1] == P.d[2]
                d = P.d[1]
            else
                d = rand(Uniform(P.d[1], P.d[2]))
            end

            IP = InputParameters(nc,
            1000,
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

        M[:, i] = map(Float64, sresult.output.DF[:cumsum])
        #M[:, i] = float64(sresult.output.DF[:normalized])

    end

    MAD = median(abs(M .- median(M, 2)), 2)
    MAD = MAD./sum(MAD)
    MAD = ones(Float64, length(MAD))

    writedlm("/data/home/mpx155/scripts/early_events/ABC-SMCnew/weights/weight.$(cst.fmin).$(cst.fmax).$(cst.read_depth).txt", MAD)

    return MAD
end

function getscales(particles)

  #get matrix that has each model, the number of parameters in each model and the number of particles associated with each model
  models = map(x -> [x.modelcurr, x.nparams], particles)
  models = hcat(models...)'
  modelcount = counts(models[:,1], minimum(models[:,1]):maximum(models[:,1]))
  modelcount = modelcount[modelcount.>0]
  models = unique(sort(models, 1), 1)
  models = hcat(models, modelcount)

  scales = Any[]

  for m in 1:size(models, 1)
    x = zeros(Float64, models[m, 3], models[m, 2])
    k = 1

    for p in particles

      if p.modelcurr == models[m, 1]
        x[k, :] = p.paramscurr
        k = k + 1
      end

    end

    #println(x)

    #calculate scale -> (max-min)/2
    if maximum(x, 1) == minimum(x, 1)
      #don't adapt scales if only one particle exists, find particles with correct model and set "new" scale to old scle
      correctmodel = false
      l = 1
      while correctmodel == false
        if particles[l].modelprev == models[m, 1]
          println("Only one particle with model $(models[m, 1]), setting perturbation scale to previous scale")
          push!(scales, particles[l].pertscale)
          correctmodel = true
        end
        l = l + 1
      end

    else
      push!(scales, collect(((maximum(x, 1) - minimum(x, 1))./4)'))
    end

  end

  #use a dictionary to take car of possible dead models
  Dscales = Dict(models[i,1] => scales[i][:] for i in 1:length(models[:,1]))

  for p in particles
    p.pertscale = collect(Dscales[p.modelcurr])
  end

  return particles

end

function getmodelmarginals(particles, Pr, sname)

  models = sort(unique(map(x -> x.modelcurr, particles)))

  modelmarginal = zeros(Float64, Pr.numclones[2] +  1)

  for i in models
    sumweight = 0.0
    for j in particles
      if j.modelcurr == i
        sumweight += j.weightcurr
      end
    end
    modelmarginal[i + 1] = sumweight
  end

  modelmarginal = modelmarginal ./ sum(modelmarginal)

  for i in 1:length(modelmarginal)
    if isnan(modelmarginal[i]) == true
      println("Have NA")
      println(map(x -> x.weightcurr, particles))
      f = open("log/$(sname).txt","a+")
      write(f, "\n")
      write(f, "Have NA\n")
      write(f, "Current weight\n")
      write(f, "$(map(x -> x.weightcurr, particles))\n")
      write(f, "Current model\n")
      write(f, "$(map(x -> x.modelcurr, particles))\n")
      write(f, "Previous model\n")
      write(f, "$(map(x -> x.modelprev, particles))\n")
      write(f, "Current params\n")
      write(f, "$(map(x -> x.paramscurr, particles))\n")
      write(f, "Previous params\n")
      write(f, "$(map(x -> x.paramsprev, particles))\n")
      write(f, "Margins prev\n")
      write(f, "$(map(x -> x.marginsprev, particles))\n")
      write(f, "Margins curr\n")
      write(f, "$(map(x -> x.marginscurr, particles))\n")
      close(f)
      modelmarginal[i] = 1.0
    end
  end

  for p in particles
    p.marginsprev = p.marginscurr
    p.marginscurr = modelmarginal
  end

  return particles, modelmarginal
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

function particleperturbationkernel(x0, scale, Pr)

  if (x0 - scale) >= (x0 + scale)
    return x0
  end

  return rand(Uniform(x0 - scale, x0 + scale))

end

function newinput(particle)

  newin = zeros(Float64, particle.nparams)
  #println(length(newin))
  #println(length(particle.pertscale))
  #println(length(particle.paramscurr))
  #println(length(particle.priorArray))

  for i in 1:particle.nparams
    newin[i] = particleperturbationkernel(particle.paramscurr[i],
                particle.pertscale[i],
                particle.priorArray[i, :])
  end

  clonalmuts = newin[1]
  μ = newin[2]
  d = newin[3]

  tevent = Float64[]
  s = Float64[]
  k = 4
  for i in 1:(particle.modelcurr)
    push!(tevent, newin[k])
    k = k +1
    push!(s, newin[k])
    k = k + 1
  end

  return clonalmuts, μ, d, tevent, s
end

function makepriorarray(Pr)
  #[clonalmuts, μ, d, [tevent s]'[:]]
  x = zeros(Float64, 3 + 2*Pr.numclones[2], 2)
  x[1,:]=Pr.clonalmuts
  x[2,:]=Pr.μ
  x[3,:]=Pr.d

  k = 4
  for i in 1:Pr.numclones[2]
    x[k, :] = Pr.tevent
    k = k + 1
    x[k, :] = Pr.selection
    k = k + 1
  end

  return x

end

function modelkernel(currm, modelmarginal, Pr)

    prob = 0.5

    mprob = ones(Float64, length(modelmarginal))
    mprob[modelmarginal.==0.0] = 0.0

    nsurvivingmodels = sum(mprob)

    mprob[mprob.> 0.0] = (1 - prob) / (nsurvivingmodels - 1)
    mprob[currm + 1] = prob

    wsample(Pr.numclones[1]:Pr.numclones[2], mprob)

end

function getmodelprob(currmodel, prevmodel, modelmarginal)

  prob = 0.5

  if currmodel == prevmodel
    return prob
  elseif sum(modelmarginal.>0.0) > 1
    return (1 - prob) / (sum(modelmarginal.>0.0) - 1)
  else
    return prob
  end

end

function getmodelprob2(currmodel, prevmodel, modelmarginal)

  prob = 0.5

  if currmodel == prevmodel
    return prob
  else
    return (1 - prob) / (sum(modelmarginal.>0.0) - 1)
  end

end

function denom_mfunc(thismodel, modelmarginalsprev, Pr)

  denom_m = 0.0
  for i in Pr.numclones[1]:Pr.numclones[2]
    denom_m = denom_m + modelmarginalsprev[i + 1] * getmodelprob(thismodel, i, modelmarginalsprev)
  end

  return denom_m

end

function kernel_prob(p1, p2)

    prob = 1

    for i in 1:length(p1.paramscurr)

      if (p2.paramsprev[i] - p2.pertscale[i]) < (p2.paramsprev[i] + p2.pertscale[i])

        prob = prob * pdf(Uniform(p2.paramsprev[i] - p2.pertscale[i], p2.paramsprev[i] + p2.pertscale[i]), p1.paramscurr[i])

      end

    end

    return prob
end

function denom_pfunc(p, particles, Pr)

  denom_p = 0.0

  for i in particles

    if p.modelcurr == i.modelprev

      denom_p = denom_p + ((i.weightprev * kernel_prob(p, i)) / i.marginsprev[p.modelcurr + 1])
    end
  end

  return denom_p
end

function getmodelweights(particles, Pr)

  numer = Float64[]
  denom_m = Float64[]
  denom_p = Float64[]

  for p in particles

    push!(numer, p.priorprob)
    push!(denom_m, denom_mfunc(p.modelcurr, p.marginsprev, Pr))
    push!(denom_p, denom_pfunc(p, particles, Pr))

  end

  println("models")
  println(map(x -> x.modelcurr, particles))
  println("")
  println("numer")
  println(numer)
  println("")
  println("denom_m")
  println(denom_m)
  println("")
  println("denom_p")
  println(denom_p)
  println("")

  weights = numer ./ (denom_m .* denom_p)

  weights = weights / sum(weights)

  println("weights")
  println(weights)
  println("")

  i = 1
  for p in particles
    p.weightcurr = weights[i]
    i = i + 1
  end

  return particles, weights, denom_m, denom_p, numer
end

function convertpart2results(particles, denom_m, denom_p, numer)

    DF, AD = getallmetrics(particles[1].IPA)
    DF[:draw] = particles[1].dist
    DF[:weight] = particles[1].weightcurr
    DF[:modelcurr] = particles[1].modelcurr
    DF[:modelprev] = particles[1].modelprev
    DF[:weightprev] = particles[1].weightprev

    k = 2
    for p in particles[2:end]

        DF1, AD = getallmetrics(p.IPA)
        DF1[:draw] = p.dist
        DF1[:weight] = p.weightcurr
        DF1[:modelcurr] = p.modelcurr
        DF1[:modelprev] = p.modelprev
        DF1[:weightprev] = p.weightprev
        append!(DF, DF1)
        k = k + 1

    end

    DF[:denom_m] = denom_m
    DF[:denom_p] = denom_p
    DF[:numer] = numer

    return DF
end


function printresultslog(DFresults, Pr, eps, sname, numsims, dist)

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

    numclones = sort(unique(DFresults[:numclones]))
    N = length(DFresults[:s])

    for i in numclones
        write(f, "Parameters for model with $i clones\n")
        DF = DFresults[DFresults[:numclones].==i, :]
        write(f, "\tProbability of having $i clones = $(length(DF[:s]) / N)\n")
        write(f, "\tProbability of having $i clones (marginal) = $(sum(DF[:weight]))\n")
        write(f, "\tMedian clonal mutations - $(median(DF[:clonalmuts]))\n")
        write(f, "\t95 % Range clonal mutations - ($(quantile(DF[:clonalmuts], 0.025)), $(quantile(DF[:clonalmuts], 0.975)))\n")
        write(f, "\tMedian mutation rate - $(median(map(Float64, DF[:mu])))\n")
        write(f, "\t95 % Range mutation rate - ($(quantile(map(Float64, DF[:mu]), 0.025)), $(quantile(map(Float64, DF[:mu]), 0.975)))\n")
        write(f, "\tMedian death rate - $(median(map(Float64, DF[:d])))\n")
        write(f, "\t95 % Range death rate - ($(quantile(map(Float64, DF[:d]), 0.025)), $(quantile(map(Float64, DF[:d]), 0.975)))\n")

        timeevent = Float64[]
        sevent = Float64[]
        clonemuts = Float64[]
        pctfit = Float64[]
        deathrate = Float64[]

        timeeventrangel = Float64[]
        seventrangel = Float64[]
        clonemutsrangel = Float64[]
        pctfitrangel = Float64[]
        deathraterangel = Float64[]

        timeeventrangeu = Float64[]
        seventrangeu = Float64[]
        clonemutsrangeu = Float64[]
        pctfitrangeu = Float64[]
        deathraterangeu = Float64[]

        for k in 1:i

            push!(timeevent, median(map(x -> parse(Float64, split(x,",")[k]), DF[:clonetime])))
            push!(sevent, median(map(x -> parse(Float64, split(x,",")[k]), DF[:s])))
            push!(clonemuts, median(map(x -> parse(Float64, split(x,",")[k]), DF[:clonemuts])))
            push!(pctfit, median(map(x -> parse(Float64, split(x,",")[k]), DF[:pctfit])))
            push!(deathrate, median(map(x -> parse(Float64, split(x,",")[k + 1]), DF[:dr])))

            push!(timeeventrangeu, quantile(map(Float64, map(x -> parse(Float64, split(x,",")[k]), DF[:clonetime])), 0.975))
            push!(seventrangeu, quantile(map(Float64, map(x -> parse(Float64, split(x,",")[k]), DF[:s])), 0.975))
            push!(clonemutsrangeu, quantile(map(Float64, map(x -> parse(Float64, split(x,",")[k]), DF[:clonemuts])), 0.975))
            push!(pctfitrangeu, quantile(map(Float64, map(x -> parse(Float64, split(x,",")[k]), DF[:pctfit])), 0.975))
            push!(deathraterangeu, quantile(map(Float64, map(x -> parse(Float64, split(x,",")[k + 1]), DF[:dr])), 0.975))

            push!(timeeventrangel, quantile(map(Float64, map(x -> parse(Float64, split(x,",")[k]), DF[:clonetime])), 0.025))
            push!(seventrangel, quantile(map(Float64, map(x -> parse(Float64, split(x,",")[k]), DF[:s])), 0.025))
            push!(clonemutsrangel, quantile(map(Float64, map(x -> parse(Float64, split(x,",")[k]), DF[:clonemuts])), 0.025))
            push!(pctfitrangel, quantile(map(Float64, map(x -> parse(Float64, split(x,",")[k]), DF[:pctfit])), 0.025))
            push!(deathraterangel, quantile(map(Float64, map(x -> parse(Float64, split(x,",")[k + 1]), DF[:dr])), 0.025))

        end

        for k in 1:length(timeevent)
            write(f, "\tMedian time of mutations clone $(k): $(timeevent[k])\n")
            write(f, "\t95% Range - time of mutations clone $(k): ($(timeeventrangel[k]), $(timeeventrangeu[k]))\n")
            write(f, "\tMedian number of mutations in clone $(k): $(clonemuts[k])\n")
            write(f, "\t95% Range -  number of mutations in clone $(k): ($(clonemutsrangel[k]), $(clonemutsrangeu[k]))\n")
            write(f, "\tMedian % fit, clone $(k): $(pctfit[k])\n")
            write(f, "\t95% Range - % fit, clone $(k): ($(pctfitrangel[k]), $(pctfitrangeu[k]))\n")
            write(f, "\tMedian selection of clone $(k): $(sevent[k])\n")
            write(f, "\t95% Range - selection of clone $(k): ($(seventrangel[k]), $(seventrangeu[k]))\n")
            write(f, "\tMedian deathrate of clone $(k): $(deathrate[k])\n")
            write(f, "\t95% Range - deathrate of clone $(k): ($(deathraterangel[k]), $(deathraterangeu[k]))\n")

        end

    end

    write(f, "\n")
    write(f, "######################################\n")
    write(f, "\n")

    close(f)
end

function averagehistogram(particles, N)

    p = InputAndAnalysis[]
    for i in particles
        push!(p, i.IPA)
    end

    M = zeros(Int64, 100, N)
    i = 1

    for j in 1:N

        M[:, i] = freqcounts = convert(Array, p[i].sampleddata.DF[:freq])
        i = i + 1

    end

    mvalues = Float64[]
    meanvalues = Float64[]
    lquant = Float64[]
    lquart = Float64[]
    uquant = Float64[]
    uquart = Float64[]
    sd = Float64[]

    for i in 1:size(M, 1)
      push!(mvalues, median(vec(collect(M[i, :]'))))
      push!(meanvalues, mean(vec(collect(M[i, :]'))))
      push!(lquant, quantile(vec(collect(M[i, :]')), 0.025))
      push!(uquant, quantile(vec(collect(M[i, :]')), 0.975))
      push!(lquart, quantile(vec(collect(M[i, :]')), 0.25))
      push!(uquart, quantile(vec(collect(M[i, :]')), 0.75))
      push!(sd, std(vec(collect(M[i, :]'))))
    end

    DFr = DataFrame(median = mvalues,
                    mean = meanvalues,
                    lowerq = lquant,
                    upperq = uquant,
                    lowerquartile = lquart,
                    upperquartile = uquart,
                    sd = sd)

    return DFr
end

function averagehistogram2(particles, survivingmodel)

    p = InputAndAnalysis[]
    for i in particles
        if i.modelcurr == survivingmodel
          push!(p, i.IPA)
        end
    end

    M = zeros(Int64, 100, length(p))
    i = 1

    for j in 1:length(p)

        M[:, i] = freqcounts = array(p[i].sampleddata.DF[:freq])
        i = i + 1

    end

        mvalues = Float64[]
        meanvalues = Float64[]
        lquant = Float64[]
        lquart = Float64[]
        uquant = Float64[]
        uquart = Float64[]
        sd = Float64[]

        for i in 1:size(M, 1)
          push!(mvalues, median(vec(collect(M[i, :]'))))
          push!(meanvalues, mean(vec(collect(M[i, :]'))))
          push!(lquant, quantile(vec(collect(M[i, :]')), 0.025))
          push!(uquant, quantile(vec(collect(M[i, :]')), 0.975))
          push!(lquart, quantile(vec(collect(M[i, :]')), 0.25))
          push!(uquart, quantile(vec(collect(M[i, :]')), 0.75))
          push!(sd, std(vec(collect(M[i, :]'))))
        end

        DFr = DataFrame(median = mvalues,
                        mean = meanvalues,
                        lowerq = lquant,
                        upperq = uquant,
                        lowerquartile = lquart,
                        upperquartile = uquart,
                        sd = sd)
    return DFr
end


function printresults(DFresults, Pr, eps, sname, numsims, dist)

    println(" ")
    println("######################################")
    println(" ")
    println("Finished population with tolerence eps = $(eps)")

    println(" ")
    println( "Numbers of simulations = $(numsims)")
    println( "Median distance of all simulations = $(median(dist))")
    println( "Next epsillon = $(quantile(dist, 0.3))")


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
        println("\tMedian mutation rate - $(median(map(Float64, DF[:mu])))")
        println("\tMedian death rate - $(median(map(Float64, DF[:d])))")

        timeevent = Float64[]
        sevent = Float64[]
        clonemuts = Float64[]
        pctfit = Float64[]

        for k in 1:i

            push!(timeevent, median(map(x -> parse(Float64, split(x,",")[k]), DF[:clonetime])))
            push!(sevent, median(map(x -> parse(Float64, split(x,",")[k]), DF[:s])))
            push!(clonemuts, median(map(x -> parse(Float64, split(x,",")[k]), DF[:clonemuts])))
            push!(pctfit, median(map(x -> parse(Float64, split(x,",")[k]), DF[:pctfit])))

        end

        for k in 1:length(timeevent)
            println("\tMedian time of mutations clone $(k): $(timeevent[k])")
            println("\tMedian number of mutations in clone $(k): $(clonemuts[k])")
            println("\tMedian % fit, clone $(k): $(pctfit[k])")
            println("\tMedian selection of clone $(k): $(sevent[k])")

        end


        println("")

    end

    println(" ")
    println("######################################")
    println(" ")

end

function modelprob(models, weights, currmodel)

  weightsnew = zeros(Float64, length(models))

  weightsnew[models.==currmodel] = weights[models.==currmodel]

  return weightsnew

end

function writesimtodisk(particles, dataout, date, sname,eps)

        distances = map(x -> x.dist, particles)
        bestsims = sortperm(distances)[1:20]

        if isdir("$(dataout)/individualhist/$(sname)/") == false
          mkdir("$(dataout)/individualhist/$(sname)/")
        end
        rm("$(dataout)/individualhist/$(sname)/", recursive = true)
        mkdir("$(dataout)/individualhist/$(sname)/")

        for i in bestsims

          writetable("$(dataout)/individualhist/$(sname)/$(date).$(i).$(sname).model$(particles[i].modelcurr).eps$(round(particles[i].dist,3)).onesimhist.csv", particles[i].IPA.sampleddata.DF)

          writedlm("$(dataout)/individualhist/$(sname)/$(date).$(i).$(sname).model$(particles[i].modelcurr).eps$(round(particles[i].dist,3)).VAF.txt", particles[i].IPA.sampleddata.VAF)

      end

end

function getmodelfreq(particles, Pr)

  models = map(x -> x.modelcurr, particles)
  freq = Float64[]

  for i in Pr.numclones[1]:Pr.numclones[2]
    push!(freq, sum(models.==i))
  end

  return freq./sum(freq)
end

function runabcsmc(P::Priors, cst::Constants, N, ϵ, metricp, DFABC, dataout, date, sname)

  if isfile("/data/home/mpx155/scripts/early_events/ABC-SMCnew/weights/weight.$(cst.fmin).$(cst.fmax).$(cst.read_depth).txt")
    distweight = collect(readdlm("/data/home/mpx155/scripts/early_events/ABC-SMCnew/weights/weight.$(cst.fmin).$(cst.fmax).$(cst.read_depth).txt"))
    f = open("log/$(sname).txt","a+")
    write(f, "Reading in weights file\n\n")
    close(f)
  else
    f = open("log/$(sname).txt","a+")
    write(f, "Calculate weights\n")
    distweight = getvariability(P, cst, N, ϵ, metricp, DFABC)
    write(f, "Finished calculating weights\n\n")
    close(f)
  end

  #run first populations
  particles, numsims, kdist, distweight = runabcsmcpop1(P, cst, N, 10^6, metricp, DFABC, distweight)

  f = open("log/$(sname).txt","a+")
  write(f, "Finished first pass from prior\n\n")
  close(f)

  #run again with better epsilon start value
  particles, numsims, kdist, distweight = runabcsmcpop1(P, cst, N, quantile(kdist, 0.05), metricp, DFABC, distweight)

  f = open("log/$(sname).txt","a+")
  write(f, "Finished first second pass, with distance = 0.05th quantile of first pass\n Number of simulation = $(numsims)\n\n")
  close(f)

  #calculate scales for perturbation kernels
  particles = getscales(particles)
  particles, modelmarginal = getmodelmarginals(particles, P, sname)
  println("model prob")
  println(modelmarginal)
  DFresults = convertpart2results(particles, zeros(Float64, length(particles)),
  zeros(Float64, length(particles)),
  zeros(Float64, length(particles)))

  writetable("$(dataout)/eps/$(date).$(sname).$(round(ϵ[1],4)).pop1.txt", DFresults, separator = '\t')

  #save model margina probabilities and frequencies
  modelmarginalsave = modelmarginal'
  modelfrequencies = getmodelfreq(particles, P)'


  #print results to log file
  printresultslog(DFresults, P, quantile(kdist, 0.05), sname, numsims, kdist)
  printresults(DFresults, P, quantile(kdist, 0.05), sname, numsims, kdist)

  #value of epsilon is adaptive, we take the αth quantile
  α = 0.3

  #calculate prior array
  priorarray = makepriorarray(Pr)

  eps = ϵ[1]
  saveeps = Float64[]
  push!(saveeps, eps)
  numberattempts = Int64[]
  push!(numberattempts, numsims)

  #keep array of the average eps per populations
  averagedist = Float64[]
  push!(averagedist, median(kdist))

  populationnumber = 1
  eps = quantile(kdist, α)

  while eps > ϵ[2]

      populationnumber = populationnumber + 1

      eps = quantile(kdist, α)
      if eps < ϵ[2]
          eps = ϵ[2]
      end
      push!(saveeps, eps)
      particlesnew = Array(Particle, N)
      weights = map(x -> x.weightcurr, particles)

      #initialize model parameters
      tevent, s, μ, clonalmuts, d = 0.0,0.0,0.0,0.0,0.0

      i = 1
      kattempts = 0

      models = map(x -> x.modelcurr, particles)

      mstarsave = Int64[]
      mstarsave2 = Int64[]
      modelssave = Int64[]

      while i < N + 1

        #keep track of number of simulations run
        kattempts = kattempts + 1

        #choose model based on previous population models
        sresult = 0
        mstar = wsample((P.numclones[1] : P.numclones[2]), modelmarginal)
        push!(mstarsave, mstar)
        m = modelkernel(mstar, modelmarginal, Pr)
        push!(mstarsave2, m)
        #println(m)
        #m = mstar
        correctmodel = false
        p = 0.0
        pprob1 = 0.0
        while correctmodel == false
          ind = wsample(1 : length(particles), modelprob(models, weights, m))
          p = particles[ind]
          if p.modelcurr == m
            correctmodel = true
          end
        end

        #println("model 1st = $m")
        correctnc = false
        j = 1
        while correctnc == false

            clonalmuts, μ, d, tevent, s = newinput(p)
            pprob1 = priorprob(Pr, tevent, s, μ, clonalmuts, d)

            if pprob1 == 0.0
                #println("model = $m")
                #println("cm =  $clonalmuts")
                #println("mu = $μ")
                #println("t = $tevent")
                #println("s = $s")
                #println()
                #sleep(2.0)
                break
            end



            IP = InputParameters(m,
                  cst.Nmax,
                  cst.det_limit,
                  cst.ploidy,
                  cst.read_depth,
                  cst.fmin,
                  cst.fmax,
                  round(Int64, clonalmuts),
                  s,
                  μ,
                  cst.b,
                  d,
                  tevent,
                  cst.ρ)

            sresult = finalresults(metricp, DFABC, 10, IP, distweight)
            nctemp = sresult.parameters.numclones
            #println(nctemp)
            j = j + 1
            if nctemp == m
                correctnc = true
                push!(modelssave, m)
                #println("model 2nd = $m")
            end
        end

        if pprob1 == 0.0

              continue
        end

        #println("model 3rd = $m")
        #println(j)
        #println()
        #sleep(1)


        if sresult.dist <= eps

          if m > 3

            sresult, idx = reorderclones(sresult)

            pnew = Particle(3 + 2*m, #nparams
            m, #modelcurr
            particles[i].modelcurr, #modelprev
            1/N, #weightcurr
            particles[i].weightcurr, #weightprev
            vcat(clonalmuts, μ, d, [tevent[idx] s[idx]]'[:]), #paramscurr
            #[clonalmuts, μ, d, [tevent[idx]  s[idx]]...], #paramscurr
            particles[i].paramscurr, #paramsprev
            particles[i].pertscale, #pertscale
            sresult.dist, #distance
            sresult, #simulation output
            modelmarginal, #marginals curr
            p.marginscurr, #marginals prev
            pprob1, #prior probability
            priorarray #priorarray
            ) #simulation output
          else
            pnew = Particle(3 + 2*m, #nparams
            m, #modelcurr
            particles[i].modelcurr, #modelprev
            1/N, #weightcurr
            particles[i].weightcurr, #weightprev
            vcat(clonalmuts, μ, d, [tevent s]'[:]), #paramscurr
            #[clonalmuts, μ, d, [tevent s]...], #paramscurr
            particles[i].paramscurr, #paramsprev
            particles[i].pertscale, #pertscale
            sresult.dist, #distance
            sresult, #simulation output
            modelmarginal, #marginals curr
            p.marginscurr, #marginals prev
            pprob1, #prior probability
            priorarray #priorarray
            ) #simulation output
          end


          particlesnew[i] = pnew

          i = i + 1


        end

      end

      #calculate scales for perturbation kernels
      particlesnew, weights, denom_m, denom_p, numer = getmodelweights(particlesnew, Pr)
      weights = map(x -> x.weightcurr, particlesnew)
      #println(weights)
      particlesnew = getscales(particlesnew)
      particlesnew, modelmarginal = getmodelmarginals(particlesnew, P, sname)

      particles = deepcopy(particlesnew)
      kdist = map(x -> x.dist, particles)

      #save number of simulations for current epsillon and store median dist
      push!(numberattempts, kattempts)
      push!(averagedist, median(kdist))

      #save model margina probabilities and frequencies
      modelmarginalsave = vcat(modelmarginalsave, modelmarginal')
      modelfrequencies = vcat(modelfrequencies, getmodelfreq(particlesnew, P)')

      #write simulation output
      DFresults = convertpart2results(particlesnew, denom_m, denom_p, numer)
      DFresults = convertpart2results(particles, denom_m, denom_p, numer)
      writetable("$(dataout)/eps/$(date).$(sname).pop$(populationnumber).txt", DFresults, separator = '\t')
      cp("$(dataout)/eps/$(date).$(sname).pop$(populationnumber).txt",
        "$(dataout)/final/$(sname).txt", remove_destination=true)

      #write average histogram to file
      writetable("$(dataout)/eps/$(date).$(sname).pop$(populationnumber).averagehist.txt", averagehistogram(particles, N), separator = '\t')
      cp("$(dataout)/eps/$(date).$(sname).pop$(populationnumber).averagehist.txt",
        "$(dataout)/final/$(sname).averagehist.txt", remove_destination=true)

      #write average histogram for clones to file
      survivingmodels = array(sort(unique(DFresults[:numclones])))
      println(survivingmodels)
      for d in survivingmodels
          writetable("$(dataout)/avehist/$(date).$(sname).clone$(d).averagehist.pop$(populationnumber).txt", averagehistogram2(particles, d), separator = '\t')
          cp("$(dataout)/avehist/$(date).$(sname).clone$(d).averagehist.pop$(populationnumber).txt",
          "$(dataout)/final/$(sname).clone$(d).averagehist.txt", remove_destination=true)
      end

      #write 10best simulations to disk
      writesimtodisk(particles, dataout, date, sname, eps)

      #print results to log file
      printresultslog(DFresults, P, eps, sname, kattempts, kdist)
      printresults(DFresults, P, eps, sname, kattempts, kdist)

      println("sampled probability and observed :sampled frequency")
      println(particles[1].marginsprev)
      println(counts(mstarsave, minimum(mstarsave):maximum(mstarsave))./length(mstarsave))
      println(counts(mstarsave2, minimum(mstarsave2):maximum(mstarsave2))./length(mstarsave2))

      println()
      println("models:")
      println(counts(modelssave, 0:2)./length(modelssave))
      println("equal?")
      println(modelssave == mstarsave)
      #println(modelssave)
      #println(mstarsave)
      println("length modelsave $(length(modelssave))")
      println("length modelstar $(length(mstarsave))")
      sleep(5)
      println()

      f = open("log/$(sname).txt","a+")
      write(f, "\n")
      write(f, "sample frequencies: \n")
      write(f, "model marginals: $(particles[1].marginsprev) \n")
      write(f, "sampled freq $(counts(mstarsave, minimum(mstarsave):maximum(mstarsave))./length(mstarsave) )\n")
      write(f, "sampled freq2 $(counts(mstarsave2, minimum(mstarsave2):maximum(mstarsave2))./length(mstarsave2) )\n")
      close(f)

      writedlm("$(dataout)/eps/$sname.modelmarginalprob.txt", modelmarginalsave)
      writedlm("$(dataout)/eps/$sname.modelfrequency.txt", modelfrequencies)

      epsnew = quantile(kdist, α)

      f = open("log/$(sname).txt","a+")
      write(f, "\n")
      write(f, "% change in eps = $(abs(epsnew - eps)/epsnew)  \n\n")
      close(f)

      println("% change in eps = $(abs(epsnew - eps)/epsnew)")

#      if abs(epsnew - eps) < 0.0001
      if (abs(epsnew - eps)/epsnew <= 0.05) | (sum(numberattempts) > 10^6)
        f = open("log/$(sname).txt","a+")
        write(f, "\n")
        write(f, "% change in eps < 5% or > 10^6 simulations, end ABC-SMC \n")
        close(f)
        break
      end

    end

    #write tolerances and numbers of simulations to disk
    writedlm("$(dataout)/final/$sname.tolerances.txt", saveeps)
    writedlm("$(dataout)/final/$sname.numbersims.txt", numberattempts)
    writedlm("$(dataout)/final/$sname.modelmarginalprob.txt", modelmarginalsave)
    writedlm("$(dataout)/final/$sname.modelfrequency.txt", modelfrequencies)

    return numberattempts, averagedist

end

function reorderclones(sresult)

  idx = sortperm(sresult.input.pctfit, rev = true)

  sresult.input.tevent = sresult.input.tevent[idx]
  sresult.input.s = sresult.input.s[idx]
  sresult.input.pctfit = sresult.input.pctfit[idx]
  sresult.input.clonetime = sresult.input.clonetime[idx]
  sresult.input.clonemuts = sresult.input.clonemuts[idx]
  sresult.input.birthrates = [sresult.input.birthrates[1], sresult.input.birthrates[2:end][idx]]
  sresult.input.deathrates = [sresult.input.deathrates[1], sresult.input.deathrates[2:end][idx]]

  return sresult, idx
end
