import numpy as np 
import os,sys
import matplotlib.pyplot as plt
from matplotlib import rc 
rc('text', usetex=True)

DAT_PATH = "Au_soc_optics/Au_SOC.epsim"
TITLE    = "Au dielectric const. soc [img]"



dat = []

f = open(DAT_PATH,'r')
lines = f.readlines()
f.close()
for line in lines:
	if '#' in line:	
		continue

	line = np.array(line.strip().split())
	line = np.array(line[line!=''],dtype=np.float)
	#print (line)
	#exit(1)
	dat.append(line)

dat = np.array(dat)


plt.figure(1)
## xx
plt.plot(dat[:,0],dat[:,1],label="xx")
## yy
plt.plot(dat[:,0],dat[:,2],label="yy")
## zz
plt.plot(dat[:,0],dat[:,3],label="zz")
## xy
plt.plot(dat[:,0],dat[:,4],label="xy")
## yz
plt.plot(dat[:,0],dat[:,5],label="yz")
## zx
plt.plot(dat[:,0],dat[:,6],label="zx")

plt.title(TITLE)

plt.xlim([0,20])
plt.ylim([-0.5,30])
plt.xlabel("E (eV)")
plt.ylabel(r"$\epsilon_2$")
plt.legend(loc=1)
plt.savefig('Au_soc_epsim.jpg',format='jpeg',dpi=400)

plt.show()



             


