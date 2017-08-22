function newmutations2(cancercell::Stem, μ, mutID)

    #function to add new mutations to cells based on μ

    if μ == 0.0
      return cancercell,mutID
    end

    numbermutations = rand(Poisson(μ))

    cancercell.mutations = append!(cancercell.mutations, mutID:mutID+numbermutations-1)
    mutID = mutID + numbermutations

    return cancercell, mutID
end

function initializesim(α, d)

    #initialize empty arrays and first cell with clonal mutations

    #initialize time to zero
    t = 0.0
    tvec = Float64[]
    push!(tvec,t)

    N = 1
    Nvec = Int64[]
    push!(Nvec, N)

    #Initialize array of cell type that stores mutations for each cell and their fitness type
    #fitness type of 1 is the host population, lowest fitness
    cells = CancerStemCell[]
    push!(cells, CancerStemCell([], 1))

    #need to keep track of mutations, assuming infinite sites, new mutations will be unique,
    #we assign each new muation with a unique integer by simply counting up from one

    mutID = 1

    #cells[1], mutID = newmutationsinit(cells[1],clonalmutations,mutID)

    return t, tvec, N, Nvec, CancerCells(cells, [], [1, 0], α, d), mutID, [0, 0]
end

function newnonstemcell(cells, rcell)
  newcell = CancerNonStemCell(deepcopy(cells.stemcells[rcell].mutations),
            deepcopy(cells.stemcells[rcell].fitness),
            0)
  push!(cells.nonstemcells, newcell)

  return cells
end

function stemcelldivision(cells, μ, mutid, divisions)

  rcell = rand(1:cells.ncells[1])
  r = rand()
  if r < cells.α
    push!(cells.stemcells, deepcopy(cells.stemcells[rcell]))
    cells.stemcells[rcell], mutid = newmutations(cells.stemcells[rcell], μ, mutid)
    cells.stemcells[end], mutid = newmutations(cells.stemcells[end], μ, mutid)
    cells.ncells[1] += 1
    divisions[1] += 1
  else
    cells = newnonstemcell(cells, rcell)
    cells.stemcells[rcell], mutid = newmutations(cells.stemcells[rcell], μ, mutid)
    cells.nonstemcells[end], mutid = newmutations(cells.nonstemcells[end], μ, mutid)
    cells.ncells[2] += 1
    divisions[2] += 1
  end

  return cells, mutid, divisions
end

function nonstemcelldivision(cells, μ, mutid, maxdivisions)

    rcell = rand(1:cells.ncells[2])
    push!(cells.nonstemcells, deepcopy(cells.nonstemcells[rcell]))
    cells.nonstemcells[rcell], mutid = newmutations(cells.nonstemcells[rcell], μ, mutid)
    cells.nonstemcells[end], mutid = newmutations(cells.nonstemcells[end], μ, mutid)

    cells.nonstemcells[rcell].Ndivs += 1
    cells.ncells[2] += 1

    if cells.nonstemcells[rcell].Ndivs >= maxdivisions
      deleteat!(cells.nonstemcells, rcell)
      cells.ncells[2] -= 1
    end

    return cells, mutid
end

function tumourgrow_stemcell(Nmax; α = 0.1, maxdivisions = 5, d = 0.4, μ = 1)

  t, tvec, N, Nvec, cells, mutID, divisions = initializesim(α, d)
  Rmax = 1

  while sum(cells.ncells) < Nmax
    wt = cells.ncells[1] / (cells.ncells[1] + cells.ncells[2])
    celltype = wsample([1, 2], [wt, 1 - wt])
    #println(wt)
    #println(celltype)
    if celltype == 1
      #stem cell division
      cells, mutID, divisions = stemcelldivision(cells, μ, mutID, divisions)
    else
      #nonstemcell division
      cells, mutID = nonstemcelldivision(cells, μ, mutID, maxdivisions)
    end

    push!(Nvec, sum(cells.ncells))
  end

  return RawResults(cells, Nvec, divisions)
end

function run1simulationstem(Nmax; α = 0.1, maxdivisions = 5, d = 0.4, μ = 1, clonalmutations = 100)

  rawresults = tumourgrow_stemcell(Nmax; α = α, maxdivisions = maxdivisions, d = d, μ = μ)

  M, fitness = cellsconvert([rawresults.cells.nonstemcells; rawresults.cells.stemcells])

  AF = allelefreq(M, Nmax)
  AF, cmuts = allelefreqexpand(AF, μ, [])
  #AF = counts(M, minimum(M):maximum(M))
  prepend!(AF, repeat([Float64(Nmax)], inner = Int64(clonalmutations)))

  stemcellfrac = (rawresults.cells.ncells / sum(rawresults.cells.ncells))

  return StemCellSimResult(AF, rawresults.cells, rawresults.Nvec, rawresults.divisions, stemcellfrac)

end
