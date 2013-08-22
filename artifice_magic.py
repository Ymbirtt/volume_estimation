from pylab import *

def get_nums(p):
    f = open(p,'r')
    xs = []
    
    for l in f:
        xs.append(double(l))
    return xs
    

xs = get_nums("./additional_magic.txt")
hist(xs)
#xlim(72,80)
title("Distribution of returned volume estimates of [(-1,-1,-1),(1,1,19)]")
ylabel("Density")
xlabel("Estimate")
print len(xs)
print mean(xs)
print std(xs)
print len([x for x in xs if x<72])
show()