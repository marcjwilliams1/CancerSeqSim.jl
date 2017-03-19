@time include("julia/manysims.jl")
sname = "test2"
Nsims = 10000

cst = Constants(10^3, 0.05, 2, 100.0, 0.02, 0.6, log(2), 0.001)
Pr = Priors([1, 2], [1, 2], [0, 25.0], [1, 2], [3.0, 14.0], [0.0, 0.5])

@time DFABC = readtable("data/AML-primaryPlatinumWGS.csv");

@time DF = manysims(cst, Pr, Nsims, "test", DFABC)

writetable("results/results.$(sname).txt", DF[2:end], separator = '\t')
