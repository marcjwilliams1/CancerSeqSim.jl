type SimulationSamples
  input::InputParameters
  #output::SimResult
  sampleddataDF::DataFrame
end


function run1simulationsample(IP::InputParameters, minclonesize, maxclonesize)

  sresult = tumourgrow_birthdeathsample(IP.b, IP.d, IP.Nmax, IP.μ; numclones = IP.numclones, s = IP.selection, tevent = IP.tevent, clonalmutations = IP.clonalmutations);

  return sresult, IP
end

function getsamples(;nsample = 1, nclones = 1, ploidy = 2, read_depth = 100.0, detectionlimit = 5./read_depth, clonalmutations = 100.0, μ = 10.0, d = 0.0, b = log(2), ρ = 0.0, Nmax = 10^3, samplesize = round(Int64, 0.1 * Nmax), s = repeat([1.0], inner = nclones), tevent = collect(1.0:0.5:100.0)[1:nclones], cellularity = 1.0, fixedmu = false)

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
  fixedmu)

  simresult, IP = run1simulationsample(IP, 0.0, 1.0)

  samples = sample(simresult.cells, samplesize * nsample, replace = false)

  samplesvec = []
  for i in 0:(nsample - 1)
      push!(samplesvec, samples[(i * samplesize + 1):(samplesize * (i+ 1))])
  end

  mutations = unique(CancerSeqSim.cellsconvert(samples)[1]);
  DF = DataFrame(mutations = 1:maximum(mutations))

  for i in 1:nsample
      DF[Symbol("sample$i")] = counts(CancerSeqSim.cellsconvert(samplesvec[i])[1], 1:maximum(mutations))
  end

  DFs = sampledhistmultiple(Array(DF[:sample1]), Array(DF[:mutations]), samplesize, 1)
  for i in 2:nsample
      DFtemp = sampledhistmultiple(Array(DF[Symbol("sample$i")]), Array(DF[:mutations]), samplesize, i, detectionlimit = IP.detectionlimit, ploidy = IP.ploidy, read_depth = IP.read_depth, cellularity = IP.cellularity)
      DFs = join(DFs, DFtemp, on = :mutations)
  end

  return SimulationSamples(IP, DFs), simresult, DF

end



function sampledhistmultiple(AF, mutations, samplesize, nsample; detectionlimit = 0.1, ploidy = 2.0, read_depth = 100.0, cellularity = 1.0)

    AF = AF./ploidy

    AF = AF .* cellularity

    ind = AF.>(detectionlimit * samplesize)
    AF = AF[ind]
    mutations = mutations[ind]

    samp_percent = read_depth/samplesize

    depth = rand(Binomial(samplesize,samp_percent), length(AF))

    samp_alleles = map((n, p) -> rand(Binomial(n, p)), depth, AF/samplesize)

    VAF = samp_alleles./depth

    DF = DataFrame(mutations = mutations)
    DF[Symbol("VAF$(nsample)")] = VAF

    return DF
end
