using Distributions
using DataFrames
using GLM
using Stats
using HypothesisTests
using Gadfly

include("samplesim.jl")

metricp=readtable("metricprobabilities.csv")

nc = 1
Nmax = 10^3
read_depth = 300.0
ploidy = 2
det_limit = 5./read_depth
fmin = 0.05
fmax = 0.3
clonalmuts = 100
s = [1.0]
μ = 20.0
b = log(2)
d = 0.0*log(2)
tevent = [1.0]
ρ = 0.0

#srand(1)

sresult = simulationfinalresults(metricp, nclones = 1, 0.1, 0.9, read_depth = 300.0, μ = 20.0)

DFvaf = DataFrame(counts = sresult.sampleddata.counts, depth = sresult.sampleddata.depth, VAF = sresult.sampleddata.VAF)
show1(sresult)
plot(DFvaf, x = :VAF, Geom.histogram)
