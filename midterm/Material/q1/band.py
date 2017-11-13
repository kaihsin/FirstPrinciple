import numpy as np 
import matplotlib.pyplot as plt
import os,sys
from matplotlib import rc

rc('text', usetex=True)

band_path  = "band.dat"
free_path  = "free.dat" 
ID = 'Al'

f = open(band_path,'r')
lines= f.readlines()
f.close()

Fixs = [] # <k, ymin , ymax>
for i in np.arange(3,48,3):
	tmp = lines[i].strip()
	tmp = tmp.split('  ')
	tmp2 = lines[i+1].strip()
	tmp2 = tmp2.split('  ')
	Fixs.append(np.array([tmp[0],tmp[1],tmp2[1]],dtype=np.float))

Fixs = np.array(Fixs)	
Fixs = np.vstack((Fixs,np.array([0.0,Fixs[0,1],Fixs[0,2]])))
bands = []
band= []
for l in np.arange(48,len(lines),1):
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
	

free_dat = []    
f = open(free_path,'r')
lines = f.readlines()
for line in lines:      
    free_dat.append(np.array(line.split(' '),dtype=np.float))
    
f.close()
free_dat = np.array(free_dat)
## plot bands
plt.close('all')
plt.figure(1)
for b in range(Nband):
    plt.plot(bands[b,:,0],bands[b,:,1],'k-')

plt.plot(free_dat[:,0],free_dat[:,1]-11)


for ft in Fixs:
	plt.axvline(ft[0], color='0.2', linestyle='-',linewidth=0.2)


new_tick_loc = np.sort(np.unique(Fixs[:,0]))
#print (len(new_tick_loc))
#W,L,G,X,W,K
new_tick_label = ["L","K","W",r"$\Gamma$","X","W","L",r"$\Gamma$","K"]
#print plt.xticks()[0]
plt.xticks(new_tick_loc,new_tick_label)
plt.xlim([0,1])
plt.ylim([-20,20])
plt.legend()
plt.title(ID.upper() + ' band structure')
plt.savefig(ID+'.png',format='png',dpi=100)
plt.show()


