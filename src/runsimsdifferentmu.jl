#type definitions
@compat abstract type Stem end

type RawResults
  cells::CancerCells
  Nvec::Array{Float64, 1}
  divisions::Array{Int64, 1}
end

type cancercellM
    mutationsp::Array{Int64,1}
    mutationsd::Array{Int64,1}
    mutationsneg::Array{Int64,1}
    b::Float64
    d::Float64
    binitial::Float64
    dinitial::Float64
    fitness::Array{Float64, 1}
    fitnessneg::Array{Float64, 1}
    μp::Float64
    μd::Float64
    μneg::Float64
    timedrivers::Array{Float64, 1}
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

type SimResultM
    trueVAF::Array{Float64,1}
    trueVAFp::Array{Float64,1}
    trueVAFd::Array{Float64,1}
    trueVAFneg::Array{Float64,1}
    Ndrivers::Int64
    Npassengers::Int64
    Nnegative::Int64
    timedrivers::Array{Float64, 1}
    cells::Array{cancercellM, 1}
end

type InputParametersM
    Nmax::Int64
    detectionlimit::Float64
    ploidy::Int64
    read_depth::Float64
    clonalmutations::Int64
    μp::Float64
    μd::Float64
    μneg::Float64
    b::Float64
    d::Float64
    ρ::Float64
    cellularity::Float64
    s
    timefunction::Function
end

###############################################################################

function newmutations(cancercell::cancercellM, mutIDp, mutIDd, mutIDneg, Rmax, t, s)

    #function to add new mutations to cells based on μ
    numbermutationsp = rand(Poisson(cancercell.μp))
    numbermutationsd = rand(Poisson(cancercell.μd))
    numbermutationsneg = rand(Poisson(cancercell.μneg))
    cancercell.mutationsd = append!(cancercell.mutationsd, mutIDd:mutIDd+numbermutationsd-1)
    mutIDd = mutIDd + numbermutationsd

    cancercell.mutationsneg = append!(cancercell.mutationsneg, mutIDneg:mutIDneg+numbermutationsneg-1)
    mutIDneg = mutIDneg + numbermutationsneg

    b = cancercell.binitial
    #increase fitness due to driver mutations
    for i in 1:numbermutationsd
      #push!(cancercell.fitness, rand(Exponential(0.1)))
      stemp = s()
      push!(cancercell.fitness, stemp)
      push!(cancercell.timedrivers, t)
    end
    #println(cancercell.fitness)
    for i in cancercell.fitness
      b = b * (1 + i)
    end

    #decrease fitness due to driver mutations
    for i in 1:numbermutationsneg
      #push!(cancercell.fitnessneg, -rand(Exponential(0.01)))
      push!(cancercell.fitness, -s)
    end
    for i in cancercell.fitnessneg
      b = b * (1 + i)
    end

    cancercell.b = b

    cancercell.mutationsp = append!(cancercell.mutationsp, mutIDp:mutIDp+numbermutationsp-1)
    mutIDp = mutIDp + numbermutationsp

    if cancercell.b + cancercell.d > Rmax
      Rmax = cancercell.b + cancercell.d
    end

    return cancercell, mutIDp, mutIDd, mutIDneg, Rmax
end

function initializesim(mup, mud, muneg, b, d)

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
    cells = cancercellM[]
    push!(cells, cancercellM(Int64[], Int64[], Int64[], b, d, b, d, Float64[], Float64[], mup, mud, muneg, []))

    #need to keep track of mutations, assuming infinite sites, new mutations will be unique,
    #we assign each new muation with a unique integer by simply counting up from one
    mutIDp, mutIDd, mutIDneg = 1, 1, 1
    return t, tvec, N, Nvec, cells, mutIDp, mutIDd, mutIDneg
end



function initializesim(mup, mud, muneg, b, d, Ncells)

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
    cells = cancercellM[]
    push!(cells, cancercellM(Int64[], Int64[], Int64[], b, d, b, d, Float64[], Float64[], mup, mud, muneg, []))

    for i in 1:Ncells - 1
        push!(cells, cancercellM(Int64[], Int64[], Int64[], b, d, b, d, Float64[], Float64[], mup, mud, muneg, []))
    end

    #need to keep track of mutations, assuming infinite sites, new mutations will be unique,
    #we assign each new muation with a unique integer by simply counting up from one
    mutIDp, mutIDd, mutIDneg = 1, 1, 1
    return t, tvec, N, Nvec, cells, mutIDp, mutIDd, mutIDneg
end

function tumourmoran(b, d, Nmax, μp, μd, μneg, maxt; clonalmutations = μp, timefunction::Function = exptime, s = 0.1)

    t, tvec, N, Nvec, cells, mutIDp, mutIDd, mutIDneg = initializesim(μp, μd, μneg, b, d, Nmax)
    Rmax = b
    wts = b.*ones(Nmax)

    while t < maxt
      celldie = rand(1:Nmax)
      cellbirth = sample(1:Nmax, weights(wts))
      cells[celldie] = copycell(cells[cellbirth])
      cells[cellbirth], mutIDp, mutIDd, mutIDneg, Rmax = newmutations(cells[cellbirth], mutIDp, mutIDd, mutIDneg, Rmax, t, s)
      cells[celldie], mutIDp, mutIDd, mutIDneg, Rmax = newmutations(cells[celldie], mutIDp, mutIDd, mutIDneg, Rmax, t, s)
      Δt =  1/(Nmax) * timefunction()
      t = Δt + t
      push!(tvec, t)
      wts[celldie] = cells[celldie].b
      wts[cellbirth] = cells[cellbirth].b
    end

    return cells, tvec, Rmax
end


function copycell(cancercellold::cancercellM)
  newcancercell::cancercellM = cancercellM(copy(cancercellold.mutationsp),
  copy(cancercellold.mutationsd),
  copy(cancercellold.mutationsneg),
  copy(cancercellold.b),
  copy(cancercellold.d),
  copy(cancercellold.binitial),
  copy(cancercellold.dinitial),
  copy(cancercellold.fitness),
  copy(cancercellold.fitnessneg),
  copy(cancercellold.μp),
  copy(cancercellold.μd),
  copy(cancercellold.μneg),
  copy(cancercellold.timedrivers))
end

function tumourgrow_birthdeath(b, d, Nmax, μp, μd, μneg; clonalmutations = μp, timefunction::Function = exptime, s = 0.1)

    #Rmax starts with b + d and changes once a fitter mutant is introduced, this ensures that
    # b and d have correct units
    Rmax = b + d

    #initialize arrays and parameters
    t, tvec, N, Nvec, cells, mutIDp, mutIDd, mutIDneg = initializesim(μp, μd, μneg, b, d)
    muts = Int64[]
    push!(muts, mutIDp)

    while N < Nmax

        #pick a random cell
        randcell = rand(1:N)
        r = rand(Uniform(0, Rmax))
	      Nt = N
        Rmaxt = Rmax

        #birth event if r<birthrate, access correct birthrate from cells array
        if r < cells[randcell].b

            #population increases by one
            N = N + 1
            #copy cell and mutations for cell that reproduces
            push!(cells,copycell(cells[randcell]))
            #add new mutations to both new cells
            cells[randcell], mutIDp, mutIDd, mutIDneg, Rmax = newmutations(cells[randcell], mutIDp, mutIDd, mutIDneg, Rmax, t, s)
            cells[end], mutIDp, mutIDd, mutIDneg, Rmax = newmutations(cells[end], mutIDp, mutIDd, mutIDneg, Rmax, t, s)
            push!(Nvec, N)
            Δt =  1/(Rmaxt * Nt) * timefunction()
            t = t + Δt
            push!(tvec,t)
        end

        #nothing if r > b+d
        if ((cells[randcell].b + cells[randcell].d) <= r )
          push!(Nvec, N)
          Δt =  1/(Rmax * Nt) * timefunction()
          t = t + Δt
          push!(tvec,t)
        end

        #death event if b<r<b+d
        if (cells[randcell].b <= r < (cells[randcell].b + cells[randcell].d))
            #population decreases by 1
            N = N - 1
            #frequency of cell type decreases
            #remove deleted cell
            deleteat!(cells,randcell)
            push!(Nvec,N)
            Δt =  1/(Rmax * Nt) * timefunction()
            t = t + Δt
            push!(tvec,t)
        end

        #every cell dies reinitialize simulation
        if (N == 0)
            t, tvec, N, Nvec, cells, mutIDp, mutIDd, mutIDneg = initializesim(μp, μd, μneg, b, d)
            muts = Int64[]
            push!(muts, mutIDp)
        end
    end

    return cells, tvec, Rmax
    #RawOutput(Nvec, tvec, muts, cells, birthrates, deathrates, clonetype, clonetime, subclonemutations, cloneN, Ndivisions, aveDivisions)
end

function cellsconvert(cells::Array{cancercellM, 1})
    #convert from array of cell types to one array with mutations and one array with cell fitness
    fitnessb = zeros(Float64,length(cells))
    mutationsp = Int64[]
    mutationsd = Int64[]
    mutationsneg = Int64[]
    drivertime = Float64[]
    fitness = Float64[]

    for i in 1:length(cells)
        append!(mutationsp,cells[i].mutationsp)
        append!(mutationsd,cells[i].mutationsd)
        append!(mutationsneg,cells[i].mutationsneg)
        fitnessb[i] = cells[i].b
        append!(drivertime, cells[i].timedrivers)
    end

    return mutationsp, mutationsd, mutationsneg, fitnessb, sort(unique(drivertime))
end

function getresults(b, d, μp, μd, μneg, Nmax; ploidy = 2, clonalmutations = 100, timefunction = exptime, s = 0.1)

    #Nvec,tvec,mvec,cells,br,dr,ct,clonetime
    cells, tvec, Rmax = tumourgrow_birthdeath(b, d, Nmax, μp, μd, μneg; clonalmutations = 0, timefunction = timefunction, s = s);
    Mp, Md, Mneg, fitness, timedrivers = cellsconvert(cells)

    return Mp, Md, Mneg, fitness, tvec, cells, timedrivers
end

function getresults(b, d, μp, μd, μneg, Nmax, tmax; ploidy = 2, clonalmutations = 100, timefunction = exptime, s = 0.1)

    #Nvec,tvec,mvec,cells,br,dr,ct,clonetime
    cells, tvec, Rmax = tumourmoran(b, d, Nmax, μp, μd, μneg, tmax; clonalmutations = 0, timefunction = timefunction, s = s);
    Mp, Md, Mneg, fitness, timedrivers = cellsconvert(cells)

    return Mp, Md, Mneg, fitness, tvec, cells, timedrivers
end

function run1simulation(IP::InputParametersM)

    Mp, Md, Mneg, fitness, tvec, cells, timedrivers = getresults(IP.b, IP.d, IP.μp, IP.μd, IP.μneg, IP.Nmax; ploidy = IP.ploidy, clonalmutations = IP.clonalmutations, timefunction = IP.timefunction, s = IP.s)

    if length(Md) > 0
      AFd = counts(Md)
      filter!(x -> x>0, AFd)
    else
      AFd = []
    end

    if length(Mneg) > 0
      AFneg = counts(Mneg)
      filter!(x -> x>0, AFneg)
    else
      AFneg = []
    end

    if length(Mp) > 0
      AFp = counts(Mp)
      filter!(x -> x>0, AFp)
    else
      AFp = []
    end

    AFallmuts = prepend!([AFd; AFp; AFneg], repeat([Float64(IP.Nmax)], inner = IP.clonalmutations))

    return SimResultM(convert(Array{Float64, 1}, AFallmuts),
    convert(Array{Float64, 1}, AFp),
    convert(Array{Float64, 1}, AFd),
    convert(Array{Float64, 1}, AFneg),
    length(AFd),
    length(AFp),
    length(AFneg),
    timedrivers,
    cells), IP
end


function run1simulation(IP::InputParametersM, tmax)

    Mp, Md, Mneg, fitness, tvec, cells, timedrivers = getresults(IP.b, IP.d, IP.μp, IP.μd, IP.μneg, IP.Nmax, tmax; ploidy = IP.ploidy, clonalmutations = IP.clonalmutations, timefunction = IP.timefunction, s = IP.s)

    if length(Md) > 0
      AFd = counts(Md)
      filter!(x -> x>0, AFd)
    else
      AFd = []
    end

    if length(Mneg) > 0
      AFneg = counts(Mneg)
      filter!(x -> x>0, AFneg)
    else
      AFneg = []
    end

    if length(Mp) > 0
      AFp = counts(Mp)
      filter!(x -> x>0, AFp)
    else
      AFp = []
    end

    AFallmuts = prepend!([AFd; AFp; AFneg], repeat([Float64(IP.Nmax)], inner = IP.clonalmutations))

    return SimResultM(convert(Array{Float64, 1}, AFallmuts),
    convert(Array{Float64, 1}, AFp),
    convert(Array{Float64, 1}, AFd),
    convert(Array{Float64, 1}, AFneg),
    length(AFd),
    length(AFp),
    length(AFneg),
    timedrivers,
    cells), IP
end

###############################################################################
