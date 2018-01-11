import numpy as np 
import matplotlib.pyplot as plt
import os,sys
from matplotlib import rc

band_path  ='super2/dispr/FREQ' 
ID = 'Si'

f = open(band_path,'r')
lines= f.readlines()
f.close()
dat = []
for line in lines:
	line = np.array(line.strip().split(' '))
		
	line = np.array(line[line!=''],dtype=np.float)

	dat.append(line)

dat= np.array(dat)

	
band_path2  ='super3/dispr/FREQ' 

f = open(band_path2,'r')
lines= f.readlines()
f.close()
dat2 = []
for line in lines:
	line = np.array(line.strip().split(' '))
		
	line = np.array(line[line!=''],dtype=np.float)

	dat2.append(line)

dat2= np.array(dat2)






plt.figure(1)
plt.plot(dat[:,0],dat[:,1],'k-')
plt.plot(dat[:,0],dat[:,2],'k-')
plt.plot(dat[:,0],dat[:,3],'k-')
plt.plot(dat[:,0],dat[:,4],'k-')
plt.plot(dat[:,0],dat[:,5],'k-')
plt.plot(dat[:,0],dat[:,6],'k-',label='2x2x2')

plt.plot(dat2[:,0],dat2[:,1],'b--')
plt.plot(dat2[:,0],dat2[:,2],'b--')
plt.plot(dat2[:,0],dat2[:,3],'b--')
plt.plot(dat2[:,0],dat2[:,4],'b--')
plt.plot(dat2[:,0],dat2[:,5],'b--')
plt.plot(dat2[:,0],dat2[:,6],'b--',label='3x3x3')

new_tick_loc = [0.0,0.8603,1.9145,2.6993,3.1962,3.8987]
for i in new_tick_loc:
	plt.axvline(i, color='0.2', linestyle='-',linewidth=0.2)

print (len(new_tick_loc))

new_tick_label = ["L",r"$\Gamma$","K","X","W","L"]
#print plt.xticks()[0]
plt.xticks(new_tick_loc,new_tick_label)


plt.ylabel("Freq (THz)")
plt.xlim([0,new_tick_loc[-1]])
#plt.ylim([-13,10])
plt.legend()
plt.title(ID + ' phonon dispersions')
plt.savefig(ID+'_dispr.jpg',format='jpeg',dpi=400)
plt.show()


