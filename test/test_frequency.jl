# test clone frequency calculations
dummymuts = [1,1,2]
clonetype = [0,1,1]
clonefreq = [0.903,0.058,0.034]
@test isapprox(CancerSeqSim.calculateclonefreq(clonefreq, dummymuts, clonetype)[1], [0.995, 0.058, 0.034])

clonefreq = [0.036,0.013,0.732,0.002]
clonetype = [0,0,1,3]
dummymuts = [1,1,2,4]
@test isapprox(CancerSeqSim.calculateclonefreq(clonefreq, dummymuts, clonetype)[1], [0.77 ,0.013, 0.734, 0.002])

clonefreq = [0.053,0.02,0.002,0.924]
clonetype = [0,1,1,3]
dummymuts = [1,1,2,4]
@test isapprox(CancerSeqSim.calculateclonefreq(clonefreq, dummymuts, clonetype)[1], [0.999, 0.02, 0.926, 0.924])

clonefreq = [0.036,0.003,0.91,0.008,0.004]
clonetype = [0,1,1,1,3]
dummymuts = [1,1,2,4,5]
@test isapprox(CancerSeqSim.calculateclonefreq(clonefreq, dummymuts, clonetype)[1], [0.961, 0.003, 0.914, 0.008, 0.004])
