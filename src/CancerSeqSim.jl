module CancerSeqSim


using Distributions
using DataFrames
using GLM
using Stats
using HypothesisTests
using Gadfly
using Colors
using Compose

import Base.show

export
  # types
  Simulation,
  AllMetrics,
  MetricObj,
  RsqObj,
  AnalysedData,
  SampledData,
  InputParameters,
  SimResult,
  RawOutput,
  cancercell,
  bdprocess,
  SimulationSamples,

  #functions
  simulate,
  getmetrics,
  getsummary,
  vafhistogram,
  cumulativeplot,
  cumulativedist,
  birthdeathprocess,
  getsamples






### source files
include("readdata.jl")
include("runsims.jl")
include("runsimssample.jl")
include("samplesim.jl")
include("util.jl")
include("multiplesamples.jl")
include("stemcellsim.jl")


end
