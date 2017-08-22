
#type definitions

type SampledData

    DF::DataFrame
    VAF::Array{Float64,1}
    counts::Array{Int64,1}
    depth::Array{Int64,1}

end

type AnalysedData

    DF::DataFrame
    VAF::Array{Float64,1}

end

type RsqObj

    metric::Float64
    mu::Float64
    pval::Float64

end

type MetricObj

    metric::Float64
    pval::Float64

end

type AllMetrics

    rsq::RsqObj
    area::MetricObj
    Dk::MetricObj
    meanD::MetricObj

end

type Simulation
  input::InputParameters
  output::SimResult
  sampleddata::SampledData
end

###############################################################################


function sampledhist(AF, cellnum ; detectionlimit = 0.1, ploidy = 2.0, read_depth = 100.0, cellularity = 1.0)

    AF = AF./ploidy
    #read_depth = read_depth * cellularity
    AF = AF .* cellularity
    #detectionlimit = detectionlimit / cellularity
    filter!(x -> x > detectionlimit * cellnum, AF)
    samp_percent = read_depth/cellnum
    depth = rand(Binomial(cellnum,samp_percent), length(AF))
    samp_alleles = map((n, p) -> rand(Binomial(n, p)), depth, AF/cellnum)
    VAF = samp_alleles./depth

    #data for histogram
    x = 0.005:0.01:1.005
    y = fit(Histogram, VAF, x, closed=:right)
    DFhist = DataFrame(VAF = x[1:end-1], freq = y.weights)

    SampledData(DFhist, VAF, samp_alleles, depth)
end

function betabinom(p, n, ρ)

    μ = p * n
    shape1 = (μ / n) * ((1 / ρ) - 1)
    shape2 = n * shape1/μ - shape1

    rand(Binomial(n, rand(Beta(shape1, shape2))))
end

function sampledhist(AF, cellnum, ρ ; detectionlimit = 0.1, ploidy = 2.0, read_depth = 100.0, cellularity = 1.0)

    AF = AF./ploidy
    AF = AF .* cellularity
    #read_depth = read_depth * cellularity
    #detectionlimit = detectionlimit / cellularity
    filter!(x -> x > detectionlimit * cellnum, AF)
    samp_percent = read_depth/cellnum
    depth = rand(Binomial(cellnum, samp_percent), length(AF))
    samp_alleles = map((x, y) -> betabinom(x, y, ρ), AF/cellnum, depth)
    VAF = samp_alleles./depth

    #data for histogram
    x = 0.005:0.01:1.005
    y = fit(Histogram, VAF, x, closed=:right)
    DFhist = DataFrame(VAF = x[1:end-1], freq = y.weights)

    SampledData(DFhist, VAF, samp_alleles, depth)
end


function cumulativedist(sresult; fmin = 0.1, fmax = 0.3)

    VAF = sresult.sampleddata.VAF

    #calculate cumulative sum
    steps = fmax:-0.001:fmin
    cumsum = Array{Int64}(0)
    v = Array{Float64}(0)

    for i in steps
        push!(cumsum, sum(VAF .>= i))
        push!(v, i)
    end
    cumsum = cumsum - cumsum[1]

    DF = DataFrame(cumsum = map(Float64, cumsum), v = v)
    DF[:invf] = 1 ./ DF[:v] - 1 ./ fmax

    DF[:theory] = Mcdf(DF[:v], fmin, fmax)
    DF[:normalized] = DF[:cumsum] ./ maximum(DF[:cumsum])

    lmfit = fit(LinearModel, @formula(cumsum ~ invf + 0), DF)
    DF[:prediction] = predict(lmfit)

    return AnalysedData(DF, VAF)
end

function Mcdf(f,fmin,fmax)

    (1.0./f - 1.0/fmax) ./ (1.0/fmin - 1.0/fmax)

end

function rsq(AD, fmin, fmax, metricp)

    #fit constrained fit
    lmfit = fit(LinearModel, @formula(cumsum ~ invf + 0), AD.DF)
    #calculate R^2 value
    rsqval = 1 - (sum(residuals(lmfit) .^ 2) / sum((AD.DF[:cumsum] - 0) .^ 2))
    #extract coefficient for mutation rate
    mu = coef(lmfit)[1]
    #get pvalue
    pval = metricp[:pval][searchsortedlast(metricp[:invrsqmetric],1 - rsqval)]

    return RsqObj(rsqval, mu, pval)
end

function kolmogorovD(AD, fmin, fmax, metricp)

    xfiltered = AD.VAF[fmin .< AD.VAF .< fmax]

    n = length(xfiltered)

    cdfs = 1 - map(x->Mcdf(x, fmin, fmax), sort(xfiltered))
    δp = maximum((1:n) / n - cdfs)
    δn = -minimum((0:n-1) / n - cdfs)
    δ = max(δn, δp)
    pval = metricp[:pval][searchsortedlast(metricp[:Dkmetric], δ)]

    MetricObj(δ, pval)

end

function Draw(DF1, DF2)

    maximum(abs(DF1[:cumsum] - DF2[:cumsum]))

end

function meanDraw(DF1, DF2)

    #mean(abs(DF1[:cumsum] - DF2[:cumsum]))
    sqrt(sum((DF1[:cumsum] - DF2[:cumsum]).^2))

end

function kolmogorovDmean(AD, fmin, fmax, metricp)

    D = mean(abs(convert(Array, AD.DF[:theory]) - convert(Array, AD.DF[:normalized])))

    pval = metricp[:pval][searchsortedlast(metricp[:meanDmetric], D)]

    return MetricObj(D, pval)
end

function trapz{T<:Number}(x::Vector{T}, y::Vector{T})
    local n = length(x)
    if (length(y) != n)
        error("Vectors 'x', 'y' must be of same length")
    end
    if n == 1; return 0.0; end
    r = 0.0
    for i in 2:n
        r += (x[i] - x[i-1]) * (y[i] + y[i-1])

    end
    r / 2.0
end

function areametric(AD, fmin, fmax, metricp)

    area = abs(trapz(convert(Array, AD.DF[:invf]), convert(Array, AD.DF[:normalized])) - trapz(convert(Array, AD.DF[:invf]), convert(Array, AD.DF[:theory])))
    area = area / (1 / fmin - 1 / fmax)
    pval = metricp[:pval][searchsortedlast(metricp[:areametric], area)]

    return MetricObj(area, pval)
end

function areametricraw(AD, DFABC; fmin = 0.12, fmax = 0.8)

    area = abs(trapz(convert(Array, AD.DF[:v]), convert(Array, AD.DF[:cumsum])) - trapz(convert(Array, AD.DF[:v]), convert(Array, DFABC[:cumsum])))
    area = area / (1 / fmin - 1 / fmax)

    return area
end

function simulate(; nclones = 1, ploidy = 2, read_depth = 100.0, detectionlimit = 5./read_depth, clonalmutations = 100.0, μ = 10.0, d = 0.0, b = log(2), ρ = 0.0, Nmax = 10^3, s = repeat([1.0], inner = nclones), tevent = collect(1.0:0.5:100.0)[1:nclones], cellularity = 1.0, fixedmu = false, timefunction::Function = exptime)

    nclones == length(s) || error("Number of clones is $(nclones), size of selection coefficient array is $(length(s)), these must be the same size ")

    nclones == length(tevent) || error("Number of clones is $(nclones), size of selection coefficient array is $(length(tevent)), these must be the same size ")

    IP = InputParameters(nclones,
    Nmax,
    detectionlimit,
    ploidy,
    read_depth,
    clonalmutations,
    s,
    μ,
    b,
    d,
    tevent,
    ρ,
    cellularity,
    fixedmu,
    timefunction)

    #get simulation data
    simresult, IP = run1simulation(IP, 0.0, 1.0)

    #get sampled VAFs
    if IP.ρ > 0.0
        sampleddata = sampledhist(simresult.trueVAF, IP.Nmax, IP.ρ, detectionlimit = IP.detectionlimit, ploidy = IP.ploidy, read_depth = IP.read_depth, cellularity = IP.cellularity)
    else
        sampleddata = sampledhist(simresult.trueVAF, IP.Nmax, detectionlimit = IP.detectionlimit, ploidy = IP.ploidy, read_depth = IP.read_depth, cellularity = IP.cellularity)
    end

    return Simulation(IP, simresult, sampleddata)
end


function simulate(minclonesize, maxclonesize; nclones = 1, ploidy = 2, read_depth = 100.0, detectionlimit = 5./read_depth, clonalmutations = 100.0, μ = 10.0, d = 0.0, b = log(2), ρ = 0.0, Nmax = 10^3, cellularity = 1.0, fixedmu = false, tmin = 3.0, tmax = 20.0, smin = 0.0, smax = 25.0, timefunction::Function = exptime)

    correctnc = false

    IP, simresult = 0, 0

    while correctnc == false

      tevent = sort(rand(Uniform(tmin, tmax), nclones))
      s = rand(Uniform(smin, smax), nclones)

      IP = InputParameters(nclones,
      Nmax,
      detectionlimit,
      ploidy,
      read_depth,
      clonalmutations,
      s,
      μ,
      b,
      d,
      tevent,
      ρ,
      cellularity,
      fixedmu,
      timefunction)

      simresult, IP = run1simulation(IP, minclonesize, maxclonesize)

      nctemp = IP.numclones
      if nctemp == nclones
          correctnc = true
      end
    end

    #get sampled VAFs
    if IP.ρ > 0.0
        sampleddata = sampledhist(simresult.trueVAF, IP.Nmax, IP.ρ, detectionlimit = IP.detectionlimit, ploidy = IP.ploidy, read_depth = IP.read_depth, cellularity = IP.cellularity)
    else
        sampleddata = sampledhist(simresult.trueVAF, IP.Nmax, detectionlimit = IP.detectionlimit, ploidy = IP.ploidy, read_depth = IP.read_depth, cellularity = IP.cellularity)
    end

    return Simulation(IP, simresult, sampleddata)
end

function simulate(minclonesize, maxclonesize, independentclones::Bool; nclones = 1, ploidy = 2, read_depth = 100.0, detectionlimit = 5./read_depth, clonalmutations = 100.0, μ = 10.0, d = 0.0, b = log(2), ρ = 0.0, Nmax = 10^3, cellularity = 1.0, fixedmu = false, tmin = 3.0, tmax = 20.0, smin = 0.0, smax = 25.0, timefunction::Function = exptime)

  ct = 1
  x = 0.0

  while (ct >= 1)

    x = simulate(minclonesize, maxclonesize; nclones = nclones, ploidy = ploidy, read_depth = read_depth, detectionlimit = detectionlimit, clonalmutations = clonalmutations, μ = μ, d = d, b = b, ρ = ρ, Nmax = Nmax, cellularity = cellularity, fixedmu = fixedmu, tmin = tmin, tmax = tmax, smin = smin, smax = smax, timefunction = timefunction)

    ct = sum(x.output.clonetype)

  end

  return x
end

function simulate(minclonesize, maxclonesize, mindiff::Float64; nclones = 1, ploidy = 2, read_depth = 100.0, detectionlimit = 5./read_depth, clonalmutations = 100.0, μ = 10.0, d = 0.0, b = log(2), ρ = 0.0, Nmax = 10^3, cellularity = 1.0, fixedmu = false, tmin = 3.0, tmax = 20.0, smin = 0.0, smax = 25.0, timefunction::Function = exptime)

  nclones == 2 || error("nclones must be = 2 for this method as it is a function to simulate until we arrive at 2 clones that are greater than mindiff apart")

  ct = 1
  x = 0.0
  freqdiff = false

  while (freqdiff == false)

    x = simulate(minclonesize, maxclonesize; nclones = nclones, ploidy = ploidy, read_depth = read_depth, detectionlimit = detectionlimit, clonalmutations = clonalmutations, μ = μ, d = d, b = b, ρ = ρ, Nmax = Nmax, cellularity = cellularity, fixedmu = fixedmu, tmin = tmin, tmax = tmax, smin = smin, smax = smax, timefunction = timefunction)

    if abs(x.output.clonefreq[2] - x.output.clonefreq[1]) > mindiff
      freqdiff = true
    end

  end

  return x
end

function getmetrics(AD, metricp; fmin = 0.1, fmax = 0.3)

    rsq1 = rsq(AD, fmin, fmax, metricp)
    Dk = kolmogorovD(AD, fmin, fmax, metricp)
    meanD = kolmogorovDmean(AD, fmin, fmax, metricp)
    area = areametric(AD, fmin, fmax, metricp)

    return(AllMetrics(rsq1, area, Dk, meanD))

end

function getsummary(inandout; sname = "", fmin = 0.1, fmax = 0.3)

  simresult = inandout.output
  IP = inandout.input

  sampleddata = inandout.sampleddata

  #get cumulativedistributions
  AD = cumulativedist(inandout, fmin = fmin, fmax = fmax)

  allmetrics = getmetrics(AD, metricp, fmin = fmin, fmax = fmax)

  fitout = rsq(AD, fmin, fmax, metricp)

  DF = DataFrame(
  sname = sname,
  Nmax = IP.Nmax,
  numclones = IP.numclones,
  Nevent = join(simresult.cloneN, ","),
  teventin = join(IP.tevent, ","),
  s = join(IP.selection, ","),
  clonefreq = join(simresult.clonefreq, ","),
  teventtrue = join(simresult.clonetime, ","),
  subclonemutations = join(simresult.subclonemutations, ","),
  clonetype = join(simresult.clonetype, ","),
  br = join(simresult.birthrates, ","),
  dr = join(simresult.deathrates, ","),
  mu = IP.μ,
  muout = fitout.mu,
  b = IP.b,
  d = IP.d,
  tend = simresult.tend,
  clonalmutations = IP.clonalmutations,
  area = allmetrics.area.metric,
  area_pval = allmetrics.area.pval,
  Dk = allmetrics.Dk.metric,
  Dk_pval = allmetrics.Dk.pval,
  meanD = allmetrics.meanD.metric,
  meanD_pval = allmetrics.meanD.pval,
  rsq = allmetrics.rsq.metric,
  rsq_pval = allmetrics.rsq.pval,
  fmin = fmin,
  fmax = fmax,
  detectionlimit = IP.detectionlimit,
  ploidy = IP.ploidy,
  read_depth = IP.read_depth,
  num_muts = AD.DF[:cumsum][end]
  )

  return DF

end
