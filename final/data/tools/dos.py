import numpy as np 
import os,sys
import matplotlib.pyplot as plt

paths = []
tit   = []
paths.append("nioaf2.0.u.dat")
tit.append("NiO-AFII dos")
#paths.append("nioaf2.0.d.dat")
#tit.append("NiO-AFII spin-down")

dat = []
lbl = []
for p in paths :
	f = open(p,'r')
	lines = f.readlines()
	f.close()
	tmp = []
	tmpl=[]
	for line in lines:
		if '#' in line:	
			tmpl = line.split('#')[-1]
			tmpl = [x for x in tmpl.split(' ') if x ]
			tmpl = np.array(tmpl) 
		else: 
			tmp.append(np.array(line.strip().split(' '),dtype=np.float))
	dat.append(np.array(tmp))
	lbl.append(np.array(tmpl))
#dat = np.array(dat)

for j in range(len(paths)):
	plt.figure(j)
	for l in range(1):
		plt.plot(dat[j][:,0],dat[j][:,l+1],label=lbl[j][l+1])
	plt.title(tit[j])
	plt.axvline(0,0,6,label='Ef',color='0.6',ls='--')
	plt.xlim([-8,2.5])
	plt.ylim([0,6])
	plt.legend()
	plt.savefig('NiO-totdos.png',format='png',dpi=100,transparent=True)
"""
plt.figure(5,[10,8])
for l in range(len(dat[0][0])-1):
	plt.plot(np.hstack((dat[0][:,0],dat[1][:,0])),np.hstack((dat[0][:,l+1],-dat[1][:,l+1])),label=lbl[j][l+1])

plt.title("NiO-AFII (Ni dos)")
plt.axvline(0,-4.25,4.4,label='Ef')
plt.xlabel("E (eV)")
plt.xlim([-8,2.5])
plt.ylim([-4.25,4.4])
loc,lab = plt.yticks()
plt.yticks(loc,np.abs(loc))
plt.legend()
plt.text(-6,3,"spin-up")
plt.text(-6,-3,"spin-dwn")
#plt.grid()
plt.savefig('NiO-AFII_Ni.png',format='png',transparent=True,dpi=100)
"""
plt.show()



             


