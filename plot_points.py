from pylab import *
from scipy.stats import kstest
from scipy.stats import uniform
from scipy.stats import ks_2samp

def read_points(path):
    f = open(path, 'r')
    xs = []
    f.next()
    
    for line in f:
        xs.append(tuple([double(x) for x in line.split(',')]))
    f.close()
    return xs

def plot_points(xs):
    title(str(int(xs[0][0])) + " steps")
    
    for x in xs:
        plot(x[2],x[3], 'x', color = (x[4],x[5],x[6]))
    

#A 2-sample KS test for random vectors would be nice...    
def test_points(xs,ys):
    print xs[0][0], " steps"
    #print ks_2samp(xs[2:],ys[2:])
    print ks_2samp(xs[2],ys[2])
    print ks_2samp(xs[3],ys[3])
    print ks_2samp(xs[4],ys[4])
    print ks_2samp(xs[5],ys[5])
    print ks_2samp(xs[6],ys[6])
    print "======"
    
xs = read_points("results.csv")
ys = read_points("reference.csv")
for x in range(15):
    test_points(xs[1000*x:1000*(x+1)], ys[1000*x:1000*(x+1)])
    plot_points(xs[1000*x:1000*(x+1)])
    xlim([0,1])
    ylim([0,1])
    savefig("./images/hit_and_run_dist" + str((int(xs[1000*x][0]))) + ".png")
    #show()
    cla()