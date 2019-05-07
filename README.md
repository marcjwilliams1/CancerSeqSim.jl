# CancerSeqSim
[![Build Status](https://travis-ci.org/marcjwilliams1/CancerSeqSim.jl.svg?branch=master)](https://travis-ci.org/marcjwilliams1/CancerSeqSim.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/marcjwilliams1/CancerSeqSim.jl?branch=master&svg=true)](https://ci.appveyor.com/project/marcjwilliams1/cancerseqsim-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/marcjwilliams1/CancerSeqSim.jl/badge.svg?branch=master)](https://coveralls.io/github/marcjwilliams1/CancerSeqSim.jl?branch=master)
[![codecov.io](http://codecov.io/github/marcjwilliams1/CancerSeqSim.jl/coverage.svg?branch=master)](http://codecov.io/github/marcjwilliams1/CancerSeqSim.jl?branch=master)


Simulate tumour VAFs with different clonal structures. Package is written in the [Julia](https://julialang.org/) programming language.

Package hosts code related to the simulation framework used in the following publications:
1. [Quantification of subclonal selection in cancer from bulk sequencing data](https://www.nature.com/articles/s41588-018-0128-6)
2. [Identification of neutral tumor evolution across cancer types](https://www.nature.com/articles/ng.3489)

## Getting Started
To download the package, once you're in a Julia session type the following command:
```julia
using Pkg
Pkg.add("CancerSeqSim")
```

## Example
To simulate a tumour and generate synthetic sequencing data simply invoke the `simulate` command. There are many arguments that can be changed, for example it is possible to change the mutation rate or the number of subclones introduced into the population. To see all the optional arguments you can type `?simulate` in a julia session. The command below will simulate a tumour with a single subclone with frequency between 0.4 and 0.6, mutation rate = 10 and 200 clonalmutations.
```julia
simdata = simulate(0.4, 0.6, μ = 10.0, clonalmutations = 200)
```
A summary of the simulation will automatically be printed out.

A VAF histogram can be generated which uses the [Plots.jl](https://github.com/JuliaPlots/Plots.jl) package, so you'll need to install and load the package first:
```julia
Pkg.add("Plots.jl")
using Plots
```

Then using the `plot` function on a simulation object will generate the histogram, the red line shows the frequency of the subclone.
```julia
plot(simdata)
```
![plot](/example/exampleoneclone.png)
