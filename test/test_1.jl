# test clone frequency calculations
clonetype = [0, 0, 0]
clonefreq = [0.5, 0.4, 0.2]


# test number of mutations in subclone

srand(1)
scmuts = 200
x = simulate(0.3, 0.6, nclones = 1, extrasubclonemutations = [scmuts]);
@test x.output.subclonemutations[1] == scmuts

srand(1)
scmuts = 0
x = simulate(0.3, 0.6, nclones = 1, extrasubclonemutations = [scmuts]);
@test x.output.subclonemutations[1] < 200
