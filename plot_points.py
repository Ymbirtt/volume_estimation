from pylab import *
from scipy.stats import kstest
from scipy.stats import uniform

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
        
def test_points(xs):
    print xs[0][0], " steps"
    print kstest([x[2] for x in xs], 'uniform', args = (0,1))
    print kstest([x[3] for x in xs], 'uniform', args = (0,1))
    print kstest([x[4] for x in xs], 'uniform', args = (0,1))
    print kstest([x[5] for x in xs], 'uniform', args = (0,1))
    print kstest([x[6] for x in xs], 'uniform', args = (0,1))
    print "======"
    
xs = read_points("results.csv")
for x in range(14):
    test_points(xs[1000*x:1000*(x+1)])
    plot_points(xs[1000*x:1000*(x+1)])
    xlim([0,1])
    ylim([0,1])
    savefig("./images/grid_walk_" + str((int(xs[1000*x][0]))) + ".png")
    cla()