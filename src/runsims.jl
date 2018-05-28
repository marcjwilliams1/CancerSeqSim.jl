#type definitions
@compat abstract type Stem end

type CancerStemCell <: Stem
    mutations::Array{Int64,1}
    fitness::Int64
end

type CancerNonStemCell <: Stem
    mutations::Array{Int64,1}
    fitness::Int64
    Ndivs::Int64
end

type CancerCells
  stemcells::Array{CancerStemCell, 1}
  nonstemcells::Array{CancerNonStemCell, 1}
  ncells::Array{Int64, 1}
  α::Float64
  d::Float64
end

type RawResults
  cells::CancerCells
  Nvec::Array{Float64, 1}
  divisions::Array{Int64, 1}
end

type cancercell
    mutations::Array{Int64,1}
    fitness::Int64
end

type RawOutput
    Nvec::Array{Int64, 1}
    tvec::Array{Float64, 1}
    muts::Array{Int64, 1}
    cells::Array{cancercell, 1}
    birthrates::Array{Float64, 1}
    deathrates::Array{Float64, 1}
    clonetype::Array{Int64, 1}
    clonetime::Array{Float64, 1}
    subclonemutations::Array{Any, 1}
    cloneN::Array{Int64, 1}
    Ndivisions::Array{Int64, 1}
    aveDivisions::Array{Float64, 1}
end

type bdprocess
  N::Array{Int64, 1}
  t::Array{Float64, 1}
  clonefreq::Array{Float64, 1}
end

type SimResult
    clonefreq::Array{Float64,1}
    clonefreqp::Array{Float64,1}
    clonetime::Array{Float64,1}
    subclonemutations::Array{Int64,1}
    birthrates::Array{Float64,1}
    deathrates::Array{Float64,1}
    tend::Float64
    trueVAF::Array{Float64,1}
    cloneN::Array{Int64, 1}
    clonetype::Array{Int64, 1}
    Ndivisions::Array{Int64, 1}
    cells::Array{cancercell, 1}
    aveDivisions::Array{Float64, 1}
end

type InputParameters
    numclones::Int64
    Nmax::Int64
    detectionlimit::Float64
    ploidy::Int64
    read_depth::Float64
    clonalmutations::Int64
    selection::Array{Float64,1}
    μ::Float64
    b::Float64
    d::Float64
    tevent::Array{Float64,1}
    ρ::Float64
    cellularity::Float64
    fixedmu::Bool
    timefunction::Function
    maxclonefreq::Int64
end

type StemCellSimResult
    trueVAF::Array{Float64,1}
    cells::CancerCells
    N::Array{Float64, 1}
    divisions::Array{Float64, 1}
    stemcellfraction::Array{Float64,1}
end

###############################################################################

function newmutations(cancercell, μ, mutID)
    #function to add new mutations to cells based on μ
    numbermutations= 1
    cancercell.mutations = append!(cancercell.mutations, mutID:mutID+numbermutations-1)
    mutID = mutID + numbermutations
    return cancercell, mutID
end

function newmutationsinit(cancercell, μ, mutID)
    numbermutations = 0
    cancercell.mutations = append!(cancercell.mutations,mutID:mutID+numbermutations-1)
    mutID = mutID + numbermutations
    return cancercell, mutID
end

function initializesim(clonalmutations)

    #initialize time to zero
    t = 0.0
    tvec = Float64[]
    push!(tvec,t)

    #population starts with one cell
    N = 1
    Nvec = Int64[]
    push!(Nvec,N)

    #Initialize array of cell type that stores mutations for each cell and their fitness type
    #fitness type of 1 is the host population, lowest fitness
    cells = cancercell[]
    #sizehint!(cells, 10^5)
    push!(cells,cancercell([],1))

    #need to keep track of mutations, assuming infinite sites, new mutations will be unique,
    #we assign each new muation with a unique integer by simply counting up from one
    mutID = 1
    cells[1],mutID = newmutationsinit(cells[1],clonalmutations,mutID)

    return t,tvec,N,Nvec,cells,mutID
end

function birthdeathprocess(b, d, Nmax; nclones = 0, s = repeat([1.0], inner = nclones), tevent = collect(1.0:0.5:100.0)[1:nclones], timefunction::Function = exptime() =  - log(rand()))

  simout = tumourgrow_birthdeath(b, d, Nmax, 0.0; numclones = nclones, clonalmutations = 0.0, s = s, tevent = tevent, timefunction = timefunction)
  M, fitness = cellsconvert(simout.cells)
  pctfit=Float64[]
  for i in 1:nclones push!(pctfit,sum(fitness.==(i+1))/Nmax) end

  return bdprocess(simout.Nvec, simout.tvec, pctfit)
end

exptime() = - log(rand())
meantime() = 1

function copycell(cancercellold::cancercell)
  newcancercell::cancercell = cancercell(copy(cancercellold.mutations), copy(cancercellold.fitness))
end

function tumourgrow_birthdeath(b, d, Nmax, μ; numclones=1, clonalmutations = μ, s = [0.0], tevent=[0.0], maxclonefreq = 200, timefunction::Function = meantime)

    #set array of birthrates
    birthrates = [b]
    deathrates = [d]
    #times = vcat(tevent, 0.0)
    #time is defined in terms of population doublings
    timesN = round.(Int64, vcat(exp.(log(2).*times[1:end-1]), 0.0))

    #depending on number of clones add birthrates to model, fitness is randomly distributed between death and birth rates
    for i in 1:numclones
        push!(deathrates, rand() * deathrates[1])
        push!(birthrates,(1 + s[i]) * (birthrates[1] - deathrates[1]) + deathrates[i + 1])
    end

    #Rmax starts with b + d and changes once a fitter mutant is introduced, this ensures that
    # b and d have correct units
    Rmax = b + d

    #initialize arrays and parameters
    t,tvec,N,Nvec,cells,mutID = initializesim(clonalmutations)
    muts = Int64[]
    push!(muts,mutID)

    #we only want to introduce mutant once so have variable that keeps track of how many mutants have been introduced, keep track of which type of which clone aquires new clone
    fitmutant = 1
    clonetype = Int64[]
    clonetime = Float64[]
    subclonemutations = Any[]
    cloneN = Int64[]
    Ndivisions = Int64[]
    aveDivisions = Float64[]

    clonefreq = zeros(Int64, numclones + 1)
    clonefreq[1] = 1

    executed = false
    changemutrate = broadcast(!, BitArray(numclones + 1))

    while N < Nmax

        #pick a random cell
        randcell = rand(1:N)
        r = rand(Uniform(0,Rmax))
	     Nt = N

        #birth event if r<birthrate, access correct birthrate from cells array
        if r < birthrates[cells[randcell].fitness]

            #population increases by one
            N = N + 1
            #copy cell and mutations for cell that reproduces
            #push!(cells, deepcopy(cells[randcell]))
            push!(cells, copycell(cells[randcell]))
            #add new mutations to both new cells
            if μ > 0.0
              cells[randcell],mutID = newmutations(cells[randcell],μ,mutID)
              cells[end],mutID = newmutations(cells[end],μ,mutID)
            end
            push!(muts,mutID)
            clonefreq[cells[randcell].fitness] = clonefreq[cells[randcell].fitness] + 1
            push!(Nvec, N)
            Δt =  1/(Rmax * Nt) .* timefunction()
            t = t + Δt
            push!(tvec,t)

            #if population time is tevent, cell is mutated into fitter cell
            if N >= timesN[fitmutant]
                if fitmutant != numclones + 1
                    #one mutant turns into another "type" so decreases in frequency

                    clonefreq[cells[randcell].fitness] = clonefreq[cells[randcell].fitness] - 1
                    #keep track of how many clones
                    fitmutant += 1
                    push!(clonetype, cells[randcell].fitness)
                    #change one mutant to fitter mutant
                    cells[randcell].fitness = fitmutant
                    #new type increases in frequency
                    clonefreq[cells[randcell].fitness] = clonefreq[cells[randcell].fitness] + 1

                    #change Rmax given that we now have a new fitter mutant
                    Rmax = maximum(birthrates[1:fitmutant]) + maximum(deathrates[1:fitmutant])

                    push!(clonetime, t)
                    push!(subclonemutations, deepcopy(cells[randcell].mutations))
                    push!(cloneN, N)
                    push!(Ndivisions, length(cells[randcell].mutations))
                    divs = map(x -> length(x.mutations), cells)
                    push!(aveDivisions, mean(divs))

                end
            end
        end

        if (birthrates[cells[randcell].fitness] + deathrates[cells[randcell].fitness]) <= r
          push!(Nvec, N)
          Δt =  1/(Rmax * Nt) * timefunction()
          t = t + Δt
          push!(tvec,t)
        end

        #death event if b<r<b+d
        if (birthrates[cells[randcell].fitness] <= r < birthrates[cells[randcell].fitness] + deathrates[cells[randcell].fitness])

            #population decreases by 1
            N = N - 1
            #frequency of cell type decreases
            clonefreq[cells[randcell].fitness] = clonefreq[cells[randcell].fitness] - 1
            #remove deleted cell
            deleteat!(cells,randcell)
            push!(Nvec,N)
            Δt =  1/(Rmax * Nt) * timefunction()
            t = t + Δt
            push!(tvec,t)
        end

        #every cell dies reinitialize simulation
        if (N == 0)
            t,tvec,N,Nvec,cells,mutID = initializesim(clonalmutations)
            muts = Int64[]
            push!(muts,mutID)
        end

        if (executed == false) && ((clonefreq.>maxclonefreq) == changemutrate)
            #if population of all clones is sufficiently large no new mutations
            #are acquired, can use this approximation as only mutations above 1%
            # frequency can be reliably detected
            μ = 0
            executed = true
        end

    end
    return RawOutput(Nvec, tvec, muts, cells, birthrates, deathrates, clonetype, clonetime, subclonemutations, cloneN, Ndivisions, aveDivisions)
end

function cellsconvert(cells)
    #convert from array of cell types to one array with mutations and one array with cell fitness

    fitness = zeros(Int64,length(cells))
    mutations = Int64[]
    sizehint!(mutations, length(cells) * 10)

    for i in 1:length(cells)
        append!(mutations,cells[i].mutations)
        fitness[i] = cells[i].fitness
    end

    return mutations, fitness
end

function allelefreq(mutations, cellnum)
    #create dictionary that maps mutation ID to allele frequency
    f = counts(mutations,minimum(mutations):maximum(mutations))
    muts = collect(minimum(mutations):maximum(mutations))
    idx = f .> 0.01
    f = map(Float64, f[idx])
    muts = muts[idx]
    Dict{Int64, Float64}(muts[i]::Int64 => f[i]::Float64 for i in 1:length(f))
end

function getresults(tevent::Array{Float64, 1}, s::Array{Float64, 1}, b, d, μ, Nmax; ploidy = 2, clonalmutations = 100, nc = 0, timefunction = exptime, maxclonefreq = 200)

    #Nvec,tvec,mvec,cells,br,dr,ct,clonetime
    sresult = tumourgrow_birthdeath(b, d, Nmax, μ; numclones = nc, s = s, tevent = tevent, clonalmutations = 0, timefunction = timefunction, maxclonefreq = maxclonefreq);
    M,fitness = cellsconvert(sresult.cells)

    return M, fitness, sresult.tvec[end], sresult.clonetime, sresult.subclonemutations, sresult.birthrates, sresult.deathrates, sresult.cloneN, sresult.clonetype, sresult.Ndivisions, sresult.cells, sresult.aveDivisions

end

function allelefreqexpand(AFDict, μ, subclonemutations; fixedmu = false)

  #expand allele frequency given mutation rate and calculate number of mutations in the subclones
  #subclonemutations = convert(Array{Array{Int64,1},1}, subclonemutations)
  if fixedmu == false
    cmuts = zeros(Int64, length(subclonemutations))
    mutfreqs = collect(values(AFDict))
    mutids = collect(keys(AFDict))
    mutations = rand(Poisson(μ), length(mutfreqs))
    AFnew = zeros(Int64, sum(mutations))

    for i in 1:length(cmuts)
      idx = findin(mutids, subclonemutations[i])
      cmuts[i] = sum(mutations[idx])
    end

    j = 0
    for f in 1:length(mutfreqs)
        AFnew[(j + 1): j + mutations[f]] = fill(mutfreqs[f], mutations[f])
        j = j + mutations[f]
    end

  else
    mutfreqs = collect(values(AFDict))
    mutids = collect(keys(AFDict))
    μint = round(Int64, μ)
    mutations = fill(Int64(μ), length(mutfreqs))
    AFnew = zeros(Int64, sum(mutations))
    cmuts = zeros(Int64, length(subclonemutations))

    for i in 1:length(cmuts)
      idx = findin(mutids, subclonemutations[i])
      cmuts[i] = sum(mutations[idx])
    end

    j = 0
    for f in 1:length(mutfreqs)
      AFnew[(j + 1): j + mutations[f]] = fill(mutfreqs[f], mutations[f])
      j = j + mutations[f]
    end
  end

  return AFnew, cmuts
end

function calculateclonefreq(pctfit, cmuts, clonetype)

  cf = deepcopy(pctfit)
  parent = clonetype
  cloneid = collect(1:length(clonetype))
  prev = 0
  for i in reverse(cloneid)[1:end-1]
    if parent[i] == 0
        continue
    end
    cf[parent[i]] = cf[parent[i]] + cf[i]
    cmuts[i] = cmuts[i] - cmuts[parent[i]]
  end

  return cf, cmuts
end


function run1simulation(IP::InputParameters, minclonesize, maxclonesize)

    M, fitness, tend, clonetime, subclonemutations, br, dr, cloneN, clonetype, Ndivisions, cells, aveDivisions = getresults(IP.tevent, IP.selection, IP.b, IP.d, IP.μ, IP.Nmax; ploidy = IP.ploidy, clonalmutations = IP.clonalmutations, nc = IP.numclones, timefunction = IP.timefunction, maxclonefreq = IP.maxclonefreq)

    if length(clonetime)!= IP.numclones

        IP.numclones = length(clonetime)
        IP.tevent = IP.tevent[1:IP.numclones]
        IP.selection = IP.selection[1:IP.numclones]
        br = br[1:IP.numclones]
        dr = dr[1:IP.numclones]
        clonetype = clonetype[1:IP.numclones]

    end

    AF = allelefreq(M, IP.Nmax)
    AF, cmuts = allelefreqexpand(AF, IP.μ, subclonemutations, fixedmu = IP.fixedmu)
    #prepend!(AF, repeat([Float64(IP.Nmax)], inner = IP.clonalmutations))
    prepend!(AF, fill(Float64(IP.Nmax), IP.clonalmutations))

    pctfit=Float64[]
    for i in 1:IP.numclones push!(pctfit,sum(fitness.==(i+1))/IP.Nmax) end

    #remove clones that have frequency < detectionlimit
#    detectableclones = (pctfit.>(IP.detectionlimit)) & (pctfit.<0.95)
    if length(pctfit) > 1
      clonefreq, cmuts = calculateclonefreq(pctfit, cmuts, clonetype - 1)
      !(sum(clonefreq.>1.0) > 0) || error("There is a clone with frequency greater than 1, this should be impossible ($(clonefreq)), $(clonetype), $(pctfit)")
    else
      clonefreq = pctfit
    end

    if VERSION < v"0.6-"
        detectableclones = (clonefreq.>minclonesize) & (clonefreq.<maxclonesize)  # Deprecated as of v0.6
    else
        detectableclones = (clonefreq.>minclonesize) .& (clonefreq.<maxclonesize)
    end

    clonefreq = clonefreq[detectableclones]

    if sum(detectableclones) != IP.numclones
        IP.numclones = sum(detectableclones)
        IP.tevent = IP.tevent[detectableclones]
        IP.selection = IP.selection[detectableclones]
        clonetype = clonetype[detectableclones]
        unshift!(detectableclones, true)
        detectableclones = detectableclones[1:length(br)]
        br = br[detectableclones]
        dr = dr[detectableclones]
    end

    #return SimResults object
    return SimResult(clonefreq, pctfit, clonetime, cmuts, br, dr, tend, AF, cloneN, clonetype - 1, Ndivisions, cells, aveDivisions), IP
end

###############################################################################
