import numpy as np 
import os,sys
import matplotlib.pyplot as plt
from matplotlib import rc 
rc('text', usetex=True)

DAT_PATH = "Au_soc_optics/Au_SOC.alpha1"
TITLE    = "Au absorbtion spectra soc"



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
#plt.ylim([0,6])
plt.xlabel("E (eV)")
plt.ylabel(r"$\alpha$")
plt.legend()
plt.savefig('Au_soc_alpha.jpg',format='jpeg',dpi=400)

plt.show()



             


