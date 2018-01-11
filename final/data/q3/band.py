import numpy as np 
import matplotlib.pyplot as plt
import os,sys
from matplotlib import rc

band_path  ='band.dat' 
ID = 'Si'

f = open(band_path,'r')
lines= f.readlines()
f.close()

Fixs = [] # <k, ymin , ymax>
for i in np.arange(3,30,3):
	tmp = lines[i].strip()
	tmp = tmp.split('  ')
	tmp2 = lines[i+1].strip()
	tmp2 = tmp2.split('  ')
	Fixs.append(np.array([tmp[0],tmp[1],tmp2[1]],dtype=np.float))

Fixs = np.array(Fixs)	
Fixs = np.vstack((Fixs,np.array([0.0,Fixs[0,1],Fixs[0,2]])))
bands = []
band= []
for l in np.arange(30,len(lines),1):
	lines[l] = lines[l].strip()
	bands.append(np.array(lines[l].split('  '),dtype=np.float))

bands= np.array(bands)

## get points in each band
tmp =  np.sort(np.argwhere(bands[:,0]<1.0e-10))
Npoint = int( (tmp[1] - tmp[0]+1)/2 )

bands = bands.reshape((int(len(bands)/Npoint),Npoint,2))
Nband = len(bands)
print ("N band : %d"%(Nband))
print ("N point per band : %d"%(Npoint))
for b in range(Nband):
	new_band = bands[b,np.argsort(bands[b,:,0])]
	bands[b] = new_band
	

plt.figure(1)
for b in range(Nband):
	plt.plot(bands[b,:,0],bands[b,:,1],'k--')

for ft in Fixs:
	plt.axvline(ft[0], color='0.2', linestyle='-',linewidth=0.2)


new_tick_loc = np.sort(np.unique(Fixs[:,0]))
print (len(new_tick_loc))
new_tick_label = ["W","L",r"$\Gamma$","X","W","K"]
#print plt.xticks()[0]
plt.xticks(new_tick_loc,new_tick_label)


plt.ylabel("E-Ef (eV/Atom)")
plt.xlim([0,1])
plt.ylim([-13,10])
plt.legend()
plt.title(ID + ' GGA-PBE band structure \nEf=5.9080 (eV)')
plt.savefig(ID+'.jpg',format='jpeg',dpi=400)
plt.show()


