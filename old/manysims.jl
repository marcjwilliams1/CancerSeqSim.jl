###############################################################################

#Call required packages
using Distributions
using DataFrames
using GLM
using Stats
using HypothesisTests
using Gadfly
###############################################################################
#call samplesim.jl

include("/data/home/mpx155/scripts/early_events/ABC-SMCnew/julia/samplesim.jl")
metricp=readtable("/data/home/mpx155/scripts/early_events/ABC-SMCnew/julia/metricprobabilities.csv")


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

###############################################################################


function manysims(cst, P, numsims, sname, DFABC)

  tevent = sort(rand(Uniform(P.tevent[1], P.tevent[2]), 1))
  clonalmuts = rand(P.clonalmuts[1] : P.clonalmuts[2])
  s = rand(Uniform(P.selection[1], P.selection[2]), 1)
  μ = rand(Uniform(P.μ[1], P.μ[2]))
  fmin = rand(Uniform(0.05, 0.12))
  det_limit = fmin
  d = 0.0
  nc = 0

  IP = InputParameters(1,
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

  DF, AD = getallmetrics(sresult)

  simname = "temp"

  for i in 1:numsims

    nc = rand(0 : 3)
    correctnc = false

    println("simulation $i, clone $nc")

    while correctnc == false

      tevent = sort(rand(Uniform(P.tevent[1], P.tevent[2]), nc))
      clonalmuts = rand(P.clonalmuts[1] : P.clonalmuts[2])
      s = rand(Uniform(P.selection[1], P.selection[2]), nc)
      μ = rand(Uniform(P.μ[1], P.μ[2]))
      if P.d[1] == P.d[2]
          d = P.d[1]
      else
          d = [0.0, 0.25, 0.5][rand(1:3)]
      end

      tevent = sort(rand(Uniform(P.tevent[1], P.tevent[2]), nc))
      clonalmuts = rand(P.clonalmuts[1] : P.clonalmuts[2])
      s = rand(Uniform(P.selection[1], P.selection[2]), nc)

      datatype = ["WXS", "WGS"][rand(1:2)]

      if datatype == "WGS"
        μ = rand(50.0 : 300.0)
        clonalmuts = rand(2000 : 8000)
      else
        μ = rand(1.0 : 25.0)
        clonalmuts = rand(20 : 300)
      end

      read_depth = [50, 75, 100, 200][rand(1:4)]
      det_limit = 5/read_depth
      fmin = det_limit + 0.05
      ρ = [0.0001, 0.005, 0.01][rand(1:3)]

      simname = "clone$(nc).rd.$(read_depth).rho.$(ρ).mu.$(Int64(μ)).cm.$(clonalmuts).d.$(d).simnum.$(i)"

      IP = InputParameters(nc,
      cst.Nmax,
      det_limit,
      cst.ploidy,
      read_depth,
      IP.fmin,
      cst.fmax,
      clonalmuts,
      s,
      μ,
      cst.b,
      d*log(2),
      tevent,
      ρ)

      sresult = finalresults(metricp, DFABC, 10, IP)
      nctemp = sresult.input.numclones
      if nctemp == nc
          correctnc = true
      end
    end

    DFvaf = DataFrame(counts = sresult.sampleddata.counts, depth = sresult.sampleddata.depth, VAF = sresult.sampleddata.VAF)

    DF1, AD = getallmetrics(sresult)

    writetable("/data/BCI-EvoCa/marc/early_events/manysims/data2/$(simname).csv", DFvaf)

    append!(DF, DF1)

  end

  return DF

end
