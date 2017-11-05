import numpy as np 
import matplotlib.pyplot as plt
import os,sys
from matplotlib import rc

rc('text', usetex=True)
if len(sys.argv) < 3:
	print ("exec <band> <gap>")
	exit(1)

band_path  = sys.argv[1] 
gap_path   = sys.argv[2]
ID = os.path.basename(gap_path).split('_')[-1]
ID = ID.split('.')[0]

f = open(band_path,'r')
lines= f.readlines()
f.close()

Fixs = [] # <k, ymin , ymax>
for i in np.arange(3,24,3):
	tmp = lines[i].strip()
	tmp = tmp.split('  ')
	tmp2 = lines[i+1].strip()
	tmp2 = tmp2.split('  ')
	Fixs.append(np.array([tmp[0],tmp[1],tmp2[1]],dtype=np.float))

Fixs = np.array(Fixs)	
Fixs = np.vstack((Fixs,np.array([0.0,Fixs[0,1],Fixs[0,2]])))
bands = []
band= []
for l in np.arange(24,len(lines),1):
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
	

## read the LUMO & HOMO 
f = open(gap_path,'r')
lines = f.readlines()
f.close()
fnxt = 0
Egap = 0
lumo_id , homo_id , lumo_E,homo_E = 0,0,0,0
for line in lines:
	if fnxt == 1:
		Egap = np.float(line.strip())
	if 'LUMO' in line:
		line = line.split(' ')
		lumo_id = int(line[2])
		lumo_E = float(line[5])
	
	elif 'HOMO' in line:
		line = line.split(' ')
		homo_id = int(line[2])
		homo_E = float(line[5])
	elif 'E gap' in line:
		fnxt = 1
print ("LUMO : band : %d E = %f"%(lumo_id,lumo_E))
print ("HOMO : band : %d E = %f"%(homo_id,homo_E))
print ("band gap = %f"%(Egap))

plt.figure(1)
for b in range(Nband):
	if b == lumo_id - 1:
		plt.plot(bands[b,:,0],bands[b,:,1],'r-',label='LUMO')
	elif b == homo_id -1:	
		plt.plot(bands[b,:,0],bands[b,:,1],'b-',label='HOMO')
	else:
		plt.plot(bands[b,:,0],bands[b,:,1],'k-')

for ft in Fixs:
	plt.axvline(ft[0], color='0.2', linestyle='-',linewidth=0.2)


new_tick_loc = np.sort(np.unique(Fixs[:,0]))
print (len(new_tick_loc))
new_tick_label = [r"$\Gamma$","X","W","L",r"$\Gamma$"]
#print plt.xticks()[0]
plt.xticks(new_tick_loc,new_tick_label)
plt.xlim([0,1])
plt.legend()
plt.title(ID.upper() + ' band structure')
plt.savefig(ID+'.png',format='png',dpi=100)
plt.show()


