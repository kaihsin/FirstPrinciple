import numpy as np 
import matplotlib.pyplot as plt
import os,sys


f = open(sys.argv[1])
lines = f.readlines()
f.close()

dat = []
dT = 0
method = ''
for line in lines:
	if '#' in line:
		if 'dT' in line:
			dT = float(line.split('=')[-1])
		if 'method' in line:
			method = line.split('=')[-1].strip()
		continue
	line = line.split('\n')[0]
	line = line.strip()
	line = line.split(' ')
	line = list(filter(None,line))
	#print (line)
	dat.append(np.array(line,dtype=np.float))
#print (dat[1])
dat = np.array(dat)

Nelem = 100

fig,ax = plt.subplots(nrows=3,ncols=1,figsize=[5,6])

ax[0].set_title("%s \n dT:%4.8f"%(method,dT))

ax[0].plot(dat[0:Nelem,0],dat[0:Nelem,2],':');
#ax[0].set_xlabel('t')
ax[0].set_ylabel('x(t)')
#ax[0].set_title('x-t (dT:%4.8f)'%(dT))
ax[0].grid(1)

ax[1].plot(dat[0:Nelem,0],dat[0:Nelem,1],':');
#ax[1].set_xlabel('t')
ax[1].set_ylabel('v(t)')
#ax[0].set_title('v-t (dT:%4.8f)'%(dT))
ax[1].grid(1)

ax[2].plot(dat[0:Nelem,0],dat[0:Nelem,3],':');
ax[2].set_xlabel('t')
ax[2].set_ylabel('E(t)')
#ax[2].set_title('E-t (dT:%4.8f)'%(dT))
ax[2].grid(1)

fig.savefig('%s-%2.4E.png'%(method,dT),format='png',dpi=100,transparent=True)


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



