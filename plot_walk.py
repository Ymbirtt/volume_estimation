from pylab import *
import subprocess as sub

def dance(executable):
    xs = []
    s = ""
    proc = sub.Popen(executable, stdout = sub.PIPE)
    ss = proc.communicate()[0]
    
    for s in ss.split("\n")[:-1]:
        print s
        xs.append(tuple(float(x) for x in s.split(',')))

    print xs
    return xs
    
xs,ys,_,_,_ = zip(*dance("./sampling.exe 5"))

plot(xs,ys)
xlim([0,1])
ylim([0,1])
show()