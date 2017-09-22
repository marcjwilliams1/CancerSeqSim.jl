# CancerSeqSim

Simulate tumour VAFs with different clonal structures. Package is written in the [Julia](https://julialang.org/) programming language

## Getting Started
Package has been tested extensively with [Julia](https://julialang.org/) v0.5.1 but should work with later versions. If there any problems please report an issue.

To download the package, once you're in a Julia session type the following command:
```
Pkg.clone("https://github.com/marcjwilliams1/ApproxBayes.jl")
```

## Example
To simulate a tumour and generate synthetic sequencing data simply invoke the `simulate` command. There are many arguments that can be changed, to see all optional arguments type `?simulate` in a julia session. The command below will simulate a tumour with a single subclone with frequency between 0.4 and 0.6, mutation rate = 10 and 200 clonalmutations.
```
simdata = simulate(0.4, 0.6, μ = 10.0, clonalmutations = 200)
```
A summary of the simulation can be printed out and a VAF histogram plotted with the following commands.
```
#print summary
show(simdata)

#plot histogram
vafhistogram(simdata, annotateclones = true)
```
![plot](/example/exampleoneclone.png)

Also included is a function to simulate tumour growth with a stem cell architecture `simulatestemcells`.
