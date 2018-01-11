import numpy as np 
import matplotlib.pyplot as plt

def process_line(line):
	if '#' in line:
		return None
	line = line.strip()
	line = np.array(line.split(" "))
	return line[line!=""].astype(np.float)
	


def Readfile(fpath,process_func):
	f = open(fpath,'r')
	lines = f.readlines()
	f.close()
	Dat = []
	for line in lines:
		line = process_func(line)
		if line is None:
			continue

		Dat.append(line)

	return np.array(Dat)





ID = "AuAl2"
pth_Eps_re = "./%s.epsre"%(ID)
pth_Eps_im = "./%s.epsim"%(ID)


# read :
## xx , yy , zz , yz, xz, xy
Eps_re = Readfile(pth_Eps_re,process_line)
Eps_im = Readfile(pth_Eps_im,process_line)

plt.figure(1)
fig,ax = plt.subplots(2,sharey=True)

ax[0].plot(Eps_re[:,0],Eps_re[:,1],'x')
ax[0].set_title("Re")
ax[1].plot(Eps_im[:,0],Eps_im[:,1],'x')
ax[1].set_title("Im")
plt.show()







