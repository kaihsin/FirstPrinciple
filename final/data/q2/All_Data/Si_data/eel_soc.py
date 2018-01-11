import numpy as np 
import os,sys
import matplotlib.pyplot as plt
from matplotlib import rc 
rc('text', usetex=True)

DAT_PATH = "Si_soc_optics/Si_SOC.eels"
TITLE    = "Si energy loss function soc"



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
plt.plot(dat[:,0],dat[:,1])

plt.title(TITLE)

plt.xlim([0,20])
#plt.ylim([-0.5,30])
plt.xlabel("E (eV)")
plt.ylabel("eels")
plt.grid(1)
#plt.legend(loc=1)
plt.savefig('Si_soc_eels.jpg',format='jpeg',dpi=400)

plt.show()



             


