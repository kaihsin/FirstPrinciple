import numpy as np 
import matplotlib.pyplot as plt
import os,sys


f = open(sys.argv[1])
lines = f.readlines()
f.close()

dat = []
dT = 0
for line in lines:
	if '#' in line:
		if 'dT' in line:
			dT = float(line.split('=')[-1])
		continue
	line = line.split('\n')[0]
	line = line.strip()
	line = line.split(' ')
	line = list(filter(None,line))
	#print (line)
	dat.append(np.array(line,dtype=np.float))
#print (dat[1])
dat = np.array(dat)
#print (np.shape(dat))

plt.figure(1)
plt.plot(dat[:,0],dat[:,2],'x');
plt.xlabel('t')
plt.ylabel('x')
plt.title('x-t (dT:%4.8f)'%(dT))

plt.figure(2)
plt.plot(dat[:,0],dat[:,1],'x');
plt.xlabel('t')
plt.ylabel('v')
plt.title('v-t (dT:%4.8f)'%(dT))

plt.figure(3)
plt.plot(dat[:,0],dat[:,3],'x');
plt.xlabel('E')
plt.ylabel('v')
plt.title('E-t (dT:%4.8f)'%(dT))

plt.show()



