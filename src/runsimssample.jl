
###############################################################################

function newmutationssample(cancercell, μ, mutID)

    #function to add new mutations to cells based on μ

    if μ == 0.0
      return cancercell,mutID
    end

    numbermutations = rand(Poisson(μ))

    cancercell.mutations = append!(cancercell.mutations, mutID:mutID+numbermutations-1)
    mutID = mutID + numbermutations

    return cancercell, mutID
end

function newmutationsinitsample(cancercell, clonalmutations, mutID)

    numbermutations = clonalmutations

    cancercell.mutations = append!(cancercell.mutations,mutID:mutID+numbermutations-1)
    mutID = mutID + numbermutations


    return cancercell, mutID
end

function initializesimsample(clonalmutations)

    #initialize empty arrays and first cell with clonal mutations

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
    push!(cells,cancercell([],1))

    #need to keep track of mutations, assuming infinite sites, new mutations will be unique,
    #we assign each new muation with a unique integer by simply counting up from one

    mutID = 1

    cells[1],mutID = newmutationsinitsample(cells[1],clonalmutations,mutID)

    return t,tvec,N,Nvec,cells,mutID
end

function tumourgrow_birthdeathsample(b, d, Nmax, μ; numclones=1, clonalmutations = 0, s = [0.0], tevent=[0.0], maxclonefreq = 100)

    #set array of birthrates
    birthrates = [b]
    deathrates = [d]

    times = vcat(tevent, 0.0)

    #map!(x->log(1+x)/(b-d),s)

    #depending on number of clones add birthrates to model
    for i in 1:numclones
        push!(deathrates, rand() * deathrates[1])
        push!(birthrates,(1 + s[i]) * (birthrates[1] - deathrates[1]) + deathrates[i + 1])
    end

    #Rmax starts with b + d and changes once a fitter mutant is introduced, this ensures that
    # b and d have correct units
    Rmax = b + d

    #initialize arrays and parameters
    t,tvec,N,Nvec,cells,mutID = initializesimsample(clonalmutations)
    muts = Int64[]
    push!(muts,mutID)

    #we only want to introduce mutant once so have variable that keeps track of how many mutants have been introduced, keep track of which type of which clone aquires new clone
    fitmutant = 1
    clonetype = Int64[]
    clonetime = Float64[]
    subclonemutations = Any[]
    cloneN = Int64[]
    Ndivisions = Int64[]

    clonefreq = zeros(Int64, numclones + 1)
    clonefreq[1] = 1

    executed = false
    changemutrate = !BitArray(numclones + 1)

    while N < Nmax

        #pick a random cell
        randcell = rand(1:N)

        r = rand(Uniform(0,Rmax))

        #println([r, Rmax, (birthrates[cells[randcell].fitness] + deathrates[cells[randcell].fitness]), randcell])

	      Nt = N

        #birth event if r<birthrate, access correct birthrate from cells array
        if r < birthrates[cells[randcell].fitness]

            #population increases by one
            N = N + 1

            #copy cell and mutations for cell that reproduces
            push!(cells,deepcopy(cells[randcell]))

            #add new mutations to both new cells
            cells[randcell],mutID = newmutationssample(cells[randcell],μ,mutID)
            cells[end],mutID = newmutationssample(cells[end],μ,mutID)

            push!(muts,mutID)

            clonefreq[cells[randcell].fitness] = clonefreq[cells[randcell].fitness] + 1

            push!(Nvec, N)

            Δt =  - 1/(Rmax * Nt) * log(rand())

            t = t + Δt

            push!(tvec,t)

            #if population time is tevent, cell is mutated into fitter cell
            if t > times[fitmutant]
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

                end
            end

        end

        if (birthrates[cells[randcell].fitness] + deathrates[cells[randcell].fitness]) <= r

          push!(Nvec, N)
          Δt =  - 1/(Rmax * Nt) * log(rand())
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

            Δt =  - 1/(Rmax * Nt) * log(rand())

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
            μ = 0.0
            executed = true
        end

    end

    return RawOutput(Nvec, tvec, muts, cells, birthrates, deathrates, clonetype, clonetime, subclonemutations, cloneN, Ndivisions)
end

###############################################################################
