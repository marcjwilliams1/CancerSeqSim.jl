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
  InputAndAnalysis,
  AllMetrics,
  MetricObj,
  RsqObj,
  AnalysedData,
  SampledData,
  InputParameters,
  SimResult,
  RawOutput,
  cancercell,

  #functions
  simulationfinalresults,
  getmetrics,
  getsummary,
  vafhistogram,
  cumulativeplot,
  cumulativedist






### source files
include("readdata.jl")
include("runsims.jl")
include("samplesim.jl")
include("util.jl")


end
