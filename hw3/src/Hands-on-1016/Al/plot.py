import matplotlib.pyplot as plt
import os,sys
import numpy as np 



if len(sys.argv) < 2 :
    print ("ERROR ev.dat") 
    exit(1)

f = open(sys.argv[1],'r')
lines = f.readlines()
dat = []
for line in lines :
    line = line.split('\n')[0]
    line = line.split(' ')
    dat.append(np.array(line,dtype=np.float))
dat = np.array(dat)

## fitting:
coef = np.polyfit(dat[:,0],dat[:,1],3)
x = np.arange(dat[0,0],dat[-1,0],0.0001)
print (coef)

## get minimum :
d1_coef = np.polyder(coef,1)
a0_min = (-2*coef[1] + np.sqrt(4*coef[1]**2-12*coef[0]*coef[2]) )/(6*coef[0]) 
print ("a0_min : ",a0_min)

## get Bulk modulus :
V = dat[:,0]**3
E = dat[:,1]
coef_2 = np.polyfit(V,E,3)
print coef_2
V_min = (-2*coef_2[1] + np.sqrt(4*coef_2[1]**2-12*coef_2[0]*coef_2[2]) )/(6*coef_2[0]) 
#print V_min
B = (6*coef_2[0]*V_min+2*coef_2[1])*V_min # unit = eV/A^3
print B
B = B * 1.6E+11 # J/m^3 = Pa
B = B * 10**-9
print ("B : ",B,"GPar")

plt.plot(dat[:,0],dat[:,1],'x',label='simulate')
plt.plot(x,np.polyval(coef,x),'-',label = 'fit')
plt.xlabel('a0')
plt.ylabel('E')
plt.grid()
plt.legend()
plt.show()
