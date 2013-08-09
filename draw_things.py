from pylab import *
from math  import *

rad_conv = pi/180

fig = gcf()

for radius in range(1,11):
    fig.gca().add_artist(Circle((0,0), sqrt(radius), facecolor = 'none'))

for theta in range(0,180,30):
    x = sqrt(10)*cos(theta*rad_conv)
    y = sqrt(10)*sin(theta*rad_conv)
    plot([x,-x],[y,-y], color = 'black')


axis('equal')
xlim(-sqrt(10), sqrt(10))
ylim(-sqrt(10), sqrt(10))
show()