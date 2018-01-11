import numpy as np 
import matplotlib.pyplot as plt
import os,sys

TD_PATH = "super2/THERMO"

f = open(TD_PATH,'r')
lines = f.readlines()
f.close()

dat = []
for line in lines:
	if '#' in line:
		continue
	if 'Nan' in line:
		continue
	line = line.strip()
	line = list(filter(None,line.split(' ')))
	dat.append(np.array(line,dtype=np.float))

dat = np.array(dat)

TD_PATH = "super3/THERMO"

f = open(TD_PATH,'r')
lines = f.readlines()
f.close()

dat2 = []
for line in lines:
	if '#' in line:
		continue
	if 'Nan' in line:
		continue
	line = line.strip()
	line = list(filter(None,line.split(' ')))
	dat2.append(np.array(line,dtype=np.float))

dat2 = np.array(dat2)



## plot E
plt.figure(1)
plt.title("Internal Energy (E)")
plt.plot(dat[::2,0],dat[::2,1],'.',label='[2x2x2]')
plt.plot(dat2[1::2,0],dat2[1::2,1],'.',label='[3x3x3]')
plt.grid(1)
plt.legend(loc=2)
plt.ylabel("E (eV/cell)")
plt.xlabel("T (K)")
plt.savefig("Si_E.jpg",format='jpeg',dpi=400)

## plot F 
plt.figure(2)
plt.title("Free energy (F)")
plt.plot(dat[::2,0],dat[::2,2],'.',label = 'F [2x2x2]')
plt.plot(dat[::2,0],dat[::2,3],'--',label = 'classical lim [2x2x2]')

plt.plot(dat2[1::2,0],dat2[1::2,2],'.',label = 'F [3x3x3]')
plt.plot(dat2[1::2,0],dat2[1::2,3],'--',label = 'classical lim [3x3x3]')

plt.legend()
plt.grid(1)
plt.ylabel("F (eV/cell)")
plt.xlabel("T (K)")
plt.savefig("Si_F.jpg",format='jpeg',dpi=400)

## plot S
plt.figure(3)
plt.title("Entropy (S)")
plt.plot(dat[::2,0],dat[::2,4],'.',label='[2x2x2]')
plt.plot(dat2[1::2,0],dat2[1::2,4],'.',label='[3x3x3]')
plt.grid(1)
plt.legend()
plt.ylabel("S (kB/cell)")
plt.xlabel("T (K)")
plt.savefig("Si_S.jpg",format='jpeg',dpi=400)

## plot Cv:
plt.figure(4)
plt.title("Specific heat (Cv)")
plt.plot(dat[::2,0],dat[::2,5],'.',label='[2x2x2]')
plt.plot(dat2[1::2,0],dat2[1::2,5],'.',label='[3x3x3]')
plt.grid(1)
plt.legend()
plt.ylabel("Cv (kB/cell)")
plt.xlabel("T (K)")
plt.savefig("Si_Cv.jpg",format='jpeg',dpi=400)
plt.show()

