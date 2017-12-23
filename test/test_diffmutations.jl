srand(2)

mud = 0.01
mup = 10.0
muneg = 0.1

its = 10^7
nmutsp = zeros(Int64, its)
nmutsd = zeros(Int64, its)
nmutsneg = zeros(Int64, its)
Rmax = log(2)
for i in 1:its
  c = CancerSeqSim.cancercellM(Int64[], Int64[], Int64[], log(2), 0.0,log(2), 0.0, Float64[], Float64[], mup, mud, muneg, Float64[])
  mutIDp, mutIDd, mutIDneg = 1, 1, 1
  c, mutIDp, mutIDd, mutIDneg, Rmax = CancerSeqSim.newmutations(c, mutIDp, mutIDd, mutIDneg, Rmax, 1.0, 0.1)
  nmutsp[i] = length(c.mutationsp)
  nmutsd[i] = length(c.mutationsd)
  nmutsneg[i] = length(c.mutationsneg)
end

@test isapprox(mean(nmutsp), mup, rtol = 1e-2)
@test isapprox(mean(nmutsd), mud, rtol = 1e-2)
@test isapprox(mean(nmutsneg), muneg, rtol = 1e-2)

mud = 0.1
mup = 10.0
muneg = 0.5
its = 10^5
Ndivs = 5
nmutsp = zeros(Int64, its)
nmutsd = zeros(Int64, its)
nmutsneg = zeros(Int64, its)
Rmax = log(2)
for i in 1:its
  c = CancerSeqSim.cancercellM(Int64[], Int64[], Int64[], log(2), 0.0,log(2), 0.0, Float64[], Float64[], mup, mud, muneg, Float64[])
  mutIDp, mutIDd, mutIDneg = 1, 1, 1
  for j in 1:Ndivs
      c, mutIDp, mutIDd, mutIDneg, Rmax = CancerSeqSim.newmutations(c, mutIDp, mutIDd, mutIDneg, Rmax, 1.0, 0.1)
  end
  nmutsp[i] = length(c.mutationsp)
  nmutsd[i] = length(c.mutationsd)
  nmutsneg[i] = length(c.mutationsneg)
end

@test isapprox(mean(nmutsp), mup * Ndivs, rtol = 1e-2)
@test isapprox(mean(nmutsd), mud * Ndivs, rtol = 1e-2)

mud = 1.0
mup = 10.0
muneg = 1.0
Rmax = log(2)
c = CancerSeqSim.cancercellM(Int64[], Int64[], Int64[], log(2), 0.0, log(2), 0.0, Float64[], Float64[], mup, mud, muneg, Float64[])
mutIDp, mutIDd, mutIDneg = 1, 1, 1
c, mutIDp, mutIDd, mutIDneg, Rmax = CancerSeqSim.newmutations(c, mutIDp, mutIDd, mutIDneg, Rmax, 1.0, 0.1)
@test isapprox(prod(c.fitness+1) .* prod(c.fitnessneg + 1) * log(2), c.b)
@test Rmax == c.b + c.d

for i in 1:10
  c, mutIDp, mutIDd, mutIDneg, Rmax = CancerSeqSim.newmutations(c, mutIDp, mutIDd, mutIDneg, Rmax, 1.0, 0.1)
end
@test isapprox(prod(c.fitness+1) .* prod(c.fitnessneg + 1) * log(2), c.b)


@time cells, tvec, rm = CancerSeqSim.tumourgrow_birthdeath(log(2), 0.0, 10^3, 10.0, 0.0, 0.0)
@test rm == log(2)

@time cells, tvec, rm = CancerSeqSim.tumourgrow_birthdeath(log(2), 0.0, 10^3, 10.0, 0.1, 0.0)
@test rm == maximum(map(x -> x.b, cells))


IP = CancerSeqSim.InputParametersM(
10^2,
0.1,
2,
200.0,
200,
10,
0.001,
0.0,
log(2),
0.0,
0.0,
1.0,
0.1,
CancerSeqSim.exptime)
#get simulation data
simresult = CancerSeqSim.run1simulation(IP)

# if s = 0.0, mean VAF of passengers and drivers should be equal
meanp = Float64[]
meand = Float64[]
for i in 1:100
  x = simulatedifferentmutations(Nmax = 10^3, μp = 10.0, μd = 10.0, clonalmutations = 0, s = 0.0);
  push!(meand, mean(x.output.trueVAFd))
  push!(meanp, mean(x.output.trueVAFp))
end
@time isapprox(mean(meand), mean(meanp), rtol = 1e-2)

# if s >0.0, mean VAF of drivers > passengers
meanp = Float64[]
meand = Float64[]
for i in 1:100
  x = simulatedifferentmutations(Nmax = 10^3, μp = 0.1, μd = 0.1, clonalmutations = 0, s = 0.1);
  push!(meand, mean(x.output.trueVAFd))
  push!(meanp, mean(x.output.trueVAFp))
end
@test mean(meand) > mean(meanp)
ratio1 = mean(meand) / mean(meanp)

#if s is large above ratio should be larger
meanp = Float64[]
meand = Float64[]
for i in 1:100
  x = simulatedifferentmutations(Nmax = 10^3, μp = 0.1, μd = 0.1, clonalmutations = 0, s = 0.25);
  push!(meand, mean(x.output.trueVAFd))
  push!(meanp, mean(x.output.trueVAFp))
end
@test mean(meand) > mean(meanp)
ratio2 = mean(meand) / mean(meanp)

#if s is large above ratio should be larger
meanp = Float64[]
meand = Float64[]
for i in 1:100
  x = simulatedifferentmutations(Nmax = 10^3, μp = 0.1, μd = 0.1, clonalmutations = 0, s = 0.5);
  push!(meand, mean(x.output.trueVAFd))
  push!(meanp, mean(x.output.trueVAFp))
end
@test mean(meand) > mean(meanp)
ratio3 = mean(meand) / mean(meanp)

#if s is large above ratio should be larger
meanp = Float64[]
meand = Float64[]
for i in 1:100
  x = simulatedifferentmutations(Nmax = 10^3, μp = 0.1, μd = 0.1, clonalmutations = 0, s = 0.75);
  push!(meand, mean(x.output.trueVAFd))
  push!(meanp, mean(x.output.trueVAFp))
end
@test mean(meand) > mean(meanp)
ratio4 = mean(meand) / mean(meanp)

@test ratio4 > ratio3 > ratio2 > ratio1



###
#Moran model
###
μp = 10.0
μd = 0.01
μneg = 0.0
b = log(2)
d = 0.0
Nmax = 10^3
maxt = 5.0
t, tvec, N, Nvec, cells, mutIDp, mutIDd, mutIDneg = CancerSeqSim.initializesim(μp, μd, μneg, b, d, Nmax)

c, t, rm = CancerSeqSim.tumourmoran(b, d, Nmax, μp, μd, μneg, maxt)
Mp, Md, Mneg, fitness, tvec, cells, timedrivers = CancerSeqSim.getresults(b, d, μp, μd, μneg, Nmax, maxt)

@time x = simulatedifferentmutationsmoran(tmax; Nmax = 10^3, μp = 0.1, μd = 0.1, clonalmutations = 0, s = 0.1)

#if s is large above ratio should be larger
tmax = 20.0
meanp = Float64[]
meand = Float64[]
for i in 1:100
  x = simulatedifferentmutationsmoran(tmax, Nmax = 100, μp = 0.1, μd = 0.1, clonalmutations = 0, s = 0.0);
  push!(meand, mean(x.output.trueVAFd))
  push!(meanp, mean(x.output.trueVAFp))
end
@time isapprox(mean(meand), mean(meanp), rtol = 1e-2)
ratio1 = mean(meand) / mean(meanp)

tmax = 20.0
meanp = Float64[]
meand = Float64[]
for i in 1:100
  x = simulatedifferentmutationsmoran(tmax, Nmax = 100, μp = 0.1, μd = 0.1, clonalmutations = 0, s = 0.1);
  push!(meand, mean(x.output.trueVAFd))
  push!(meanp, mean(x.output.trueVAFp))
end
@test mean(meand) > mean(meanp)
ratio1 = mean(meand) / mean(meanp)

tmax = 20.0
meanp = Float64[]
meand = Float64[]
for i in 1:100
  x = simulatedifferentmutationsmoran(tmax, Nmax = 100, μp = 0.1, μd = 0.1, clonalmutations = 0, s = 0.2);
  push!(meand, mean(x.output.trueVAFd))
  push!(meanp, mean(x.output.trueVAFp))
end
@test mean(meand) > mean(meanp)
ratio2 = mean(meand) / mean(meanp)

@test ratio2 > ratio1
