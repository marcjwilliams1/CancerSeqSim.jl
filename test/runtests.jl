using CancerSeqSim
using Distributions
using Base.Test

tests = ["1"]

println("Running tests ...")

for t in tests
    fn = "test_$t.jl"
    println("* $fn ...")
    include(fn)
end
