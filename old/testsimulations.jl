DFABC=readtable("/data/home/mpx155/scripts/early_events/ABC-SMCnew/simulations/data/neutraltest.csv");

cst = Constants(10^3, 0.05, 2, 100.0, 0.1, 0.25, log(2.0), 0.0 )

function simulationtest(cst, Nsims, DFABC)
  srand(1)

    nc = 0
    tevent = [50.0]
    clonalmuts = rand(100:3000)
    #clonalmuts = 250
    s = [0.0]
    μ = rand(Uniform(5.0, 250.0))
    μ = 1.0
    d = 0.0
    distweight = collect(readdlm("weights/neutraltest.txt"))
    readdepth = rand(Uniform(80.0, 200.0))


    IP = InputParameters(nc,
    cst.Nmax,
    cst.det_limit,
    cst.ploidy,
    readdepth,
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
    DFout, AD = getallmetrics(sresult)

    for i in 1:Nsims - 1
      nc = 0
      tevent = [50.0]
      clonalmuts = rand(100:3000)
      #clonalmuts = 250
      s = [0.0]
      μ = rand(Uniform(5.0, 250.0))
      μ = 1.0
      distweight = collect(readdlm("weights/neutraltest.txt"))
      readdepth = rand(Uniform(80.0, 200.0))

      IP = InputParameters(nc,
      cst.Nmax,
      cst.det_limit,
      cst.ploidy,
      readdepth,
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
      DF, AD = getallmetrics(sresult)
      append!(DFout, DF)
    end


    return DFout
end

@time DF = simulationtest(cst, 200, DFABC);
