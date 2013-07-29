from pylab import *
import subprocess as sub

def dance(executable):
    xs = []
    s = ""
    proc = sub.Popen(executable, stdout = sub.PIPE)
    ss = proc.communicate()[0]
    
    for s in ss.split("\n")[:-1]:
        print s
        xs.append(tuple(float(x) for x in s[1:-1].split(',')))

    print xs
    return xs
    
xs, ys = zip(*dance("./sampling.exe"))
plot(xs, ys)
show()