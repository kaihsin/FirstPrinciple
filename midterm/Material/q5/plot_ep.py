import numpy as np
import matplotlib.pyplot as plt
import os,sys
Path = []
Path.append( "euler2o5E-3.dat"	)
Path.append( "predcorr2o5E-3.dat" )
dat	= []
method = []
dT = None

for i in range(len(Path)):
	tmp = []
	mth  = ''
	f = open(Path[i],'r')
	lines = f.readlines()
	for line in lines:
		if '#' in line:
			if 'dT' in line:
				if dT is None:
					dT = float(line.split('=')[-1])
				else :
					if not abs(dT - float(line.split('=')[-1])) < 1.0e-10 :
						print ("dT not consist.")
						exit(1) 

			if 'method' in line:
				mth = line.split('=')[-1].strip()
			continue
		line = line.split('\n')[0]
		line = line.strip()
		line = line.split(' ')
		line = list(filter(None,line))
		#print (line)
		tmp.append(np.array(line,dtype=np.float))
	method.append(mth)
	dat.append( np.array(tmp) )
	f.close()

print (dat)
## plot setting:
Nelem = 100
fig,ax = plt.subplots(nrows=3,ncols=1,figsize=[5,6])

ax[0].set_title("%s-%s \n dT:%4.8f"%(method[0],method[1],dT),loc='left')
for i in range(len(method)):
	ax[0].plot(dat[i][0:Nelem,0],dat[i][0:Nelem,2],':',label=method[i]);

ax[0].set_ylabel('x(t)')
ax[0].grid(1)

for i in range(len(method)):
	ax[1].plot(dat[i][0:Nelem,0],dat[i][0:Nelem,1],':',label=method[i]);
ax[1].set_ylabel('v(t)')
ax[1].grid(1)

for i in range(len(method)):
	ax[2].plot(dat[i][0:Nelem,0],dat[i][0:Nelem,3],':',label=method[i]);
ax[2].set_xlabel('t')
ax[2].set_ylabel('E(t)')
ax[2].grid(1)

ax[0].legend(loc='upper right',bbox_to_anchor=(1.1,1.5))

fig.savefig('%s-%s_%2.4E.png'%(method[0],method[1],dT),format='png',dpi=100,transparent=True)


"""
plt.figure(2)
plt.plot(dat[:,0],dat[:,1],'x');
plt.xlabel('t')
plt.ylabel('v')
plt.title('v-t (dT:%4.8f)'%(dT))

plt.figure(3)
plt.plot(dat[:,0],dat[:,3],'x');
plt.xlabel('t')
plt.ylabel('E')
plt.title('E-t (dT:%4.8f)'%(dT))
"""
plt.show()
