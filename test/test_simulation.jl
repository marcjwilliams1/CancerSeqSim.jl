simdata = simulate(0.4, 0.6, μ = 10.0, clonalmutations = 200)

nclones = [0, 1, 2]
for n in nclones
    global simdata = simulate(nclones = n)
    @test length(simdata.output.clonefreq) == n
end

popsize = [10^2, 10^3, 10^4]
for p in popsize
    global simdata = simulate(nclones = 0, Nmax = p)
    @test length(simdata.output.cells) == p
end

#check mutation rate inference is correct
using GLM
mu = Float64[]
for i in 1:10^3
    global simdata = simulate(nclones = 0, Nmax = 10^3, μ = 12.0, detectionlimit = 0.02, clonalmutations = 0)
    global DF = CancerSeqSim.cumulativedist(simdata, fmin = 0.1, fmax = 0.5).DF
    global lmfit = fit(LinearModel, @formula(cumsum ~ invf + 0), DF)
    push!(mu, coef(lmfit)[1])
end
@test isapprox(mean(mu), 12.0, rtol = 0.1)

#check subclone is returned in correct frequency range
simdata = simulate(0.4, 0.5, nclones = 1)
@test simdata.output.clonefreq[1] > 0.4
@test simdata.output.clonefreq[1] < 0.5
