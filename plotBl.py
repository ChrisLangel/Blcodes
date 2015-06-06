import pylab
from pylab import loadtxt

import matplotlib.pyplot as plt

file = loadtxt('BL_info.txt') 

file2 = loadtxt('n63coords.txt') 

x,H,x1,y = [],[],[],[]

for line in file:
	x.append( line[0] )
        H.append( line[5] )
for line in file2:
	x1.append( line[0] ) 
	y.append( line[1] ) 

plt.plot(x,H,x1,y)
plt.xlim([0.0,0.8])
plt.ylim([0.0,4.0])
plt.show() 
