DFABC=readtable("/data/home/mpx155/scripts/early_events/ABC-SMCnew/newdata/data/21lifebreastnew.csv");
sname = "21lifebreastnew"
dataout="/data/BCI-EvoCa/marc/early_events/ABC-SMC/samplesnew/21lifebreastnew"

read_depth=160.0
det_limit=0.03
#fmin=$detlimit
fmin=0.03
fmax=0.7
numclones1=0
t1=1.0
t2=14.0
mu1=1.0
mu2=200.0

clonalm1=1
clonalm2=3000
eps1=10000.0
eps2=0.001
rho=0.004
N=100
db=0.0
Nmax=2000
ploidy=2
numclones2=2
d1=0.0
d2=0.0
s1=0.0
s2=25.0
date="1stfeb"

cst = Constants(Nmax, det_limit, 2, read_depth, fmin, fmax, log(2.0), 0.0 )
Pr = Priors([numclones1, numclones2], [clonalm1, clonalm2], [s1, s2], [mu1, mu2],[t1, t2], [d1, d2])
ϵ = [100000.0, 0.1]

if isfile("log/$sname.txt")
  rm("log/$sname.txt")
end

srand(20)
runabcsmc(Pr, cst, N, ϵ, metricp, DFABC, dataout, date, sname)
