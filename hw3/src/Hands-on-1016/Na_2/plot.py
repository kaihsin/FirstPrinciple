import matplotlib.pyplot as plt
import os,sys
import numpy as np 



if len(sys.argv) < 3 :
    print ("ERROR <ev.dat> <ID>") 
    exit(1)

ID = sys.argv[2]

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
x2 = np.arange(V[0],V[-1],0.0001)
print (coef_2)
V_min = (-2*coef_2[1] + np.sqrt(4*coef_2[1]**2-12*coef_2[0]*coef_2[2]) )/(6*coef_2[0]) 
#print V_min
B = (6*coef_2[0]*V_min+2*coef_2[1])*V_min # unit = eV/A^3
#print (B)
B = B * 1.6E+11 # J/m^3 = Pa
B = B * 10**-9
print ("B : ",B,"GPa")

f = open("%s.res"%(ID),'w')
f.write("a0_min = %4.8lf\n"%(a0_min))
f.write("V_min = %4.8lf\n"%(V_min))
f.write("B = %4.8lf GPa\n"%(B))
f.write("[E-a0] coef = ")
for a in coef:
	f.write("%4.10lf "%(a))
f.write('\n')
f.write("[E-V ] coef = ")
for a in coef_2:
    f.write("%4.10lf "%(a))
f.write('\n')
f.close()



plt.figure(1)
plt.plot(dat[:,0],dat[:,1],'x',label='simulate')
plt.plot(x,np.polyval(coef,x),'-',label = 'fit')
plt.xlabel('a0')
plt.ylabel('E')
plt.grid()
plt.legend()
plt.savefig("%s_E-a0"%(ID),format='png',transparent=True,dpi=100)

plt.figure(2)
plt.plot(V,E,'x',label='simulate')
plt.plot(x2,np.polyval(coef_2,x2),'-',label = 'fit')
plt.xlabel('V')
plt.ylabel('E')
plt.grid()
plt.legend()
plt.savefig("%s_E-V"%(ID),format='png',transparent=True,dpi=100)

plt.show()
