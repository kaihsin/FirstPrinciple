import numpy as np 
import os,sys
import matplotlib.pyplot as plt
import scipy as sp
from scipy import constants
def free_dispersion(vec_k):
	"""
		@description : return E = \hbar^2 k^2 /(2*m_e) 
		@input  : k as vector with dim [N x dim(k)]
		@return : E with dim [N]
	"""
	return np.sum(vec_k**2,axis=1)*sp.constants.hbar**2/(2.*sp.constants.m_e)

def vector(x,y,z):
	return np.array([x,y,z],dtype=np.float)

def Construct_line(Charact_pts,Npoint):
	out_k = []
	r_x = []
	r_charact = []
	r_curr = 0
	for i in range(len(Charact_pts)-1):
		dk = (Charact_pts[i+1]-Charact_pts[i])/Npoint
		r_charact.append(r_curr)
		for t in range(Npoint):
			out_k.append(Charact_pts[i] + t*dk)
			r_x.append(r_curr)
			r_curr += np.sqrt(np.sum(dk**2))

	##last one:
	out_k.append(Charact_pts[-1])
	r_x.append(r_curr)
	r_charact.append(r_curr)

	out_k = np.array(out_k)
	r_x   = np.array(r_x)
	r_charact = np.array(r_charact)
	return r_x,out_k,r_charact

## construct character points:
## fcc Al 
## ref : http://lampx.tugraz.at/~hadley/ss1/bzones/fcc.php
a0 = 4.05 * 1.0e-10  # Angstrong

b1 = 2*np.pi/a0 * vector(1 ,-1, 1)
b2 = 2*np.pi/a0 * vector(1 , 1,-1)
b3 = 2*np.pi/a0 * vector(-1, 1, 1)

K = 0.375*b1 + 0.750*b2 + 0.375*b3
L = 0.500*b1 + 0.500*b2 + 0.500*b3
G = 0.000*b1 + 0.000*b2 + 0.000*b3
W = 0.250*b1 + 0.750*b2 + 0.500*b3
X = 0.000*b1 + 0.500*b2 + 0.500*b3
U = 0.250*b1 + 0.625*b2 + 0.625*b3

r_x , ks, r_c = Construct_line([L,K,W,G,X,W,L,G,K],20)
Es = free_dispersion(ks) ## SI unit

Es /= sp.constants.e ## J -> eV


plt.figure(1)
plt.title('free electron band')
plt.ylabel('E (eV)')
plt.plot(r_x,Es,'.')
#plt.show()


f = open('free.dat','w')
for i in range(len(r_x)):
    f.write("%10.11f %10.11f\n"%(r_x[i]/r_x[-1],Es[i]))
f.close()


for r in r_c:
	plt.axvline(r,np.max(Es),color='0.7')

plt.show()







