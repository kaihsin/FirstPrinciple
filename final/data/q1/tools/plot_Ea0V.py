import matplotlib.pyplot as plt
import os,sys
import numpy as np 


def Fit_getMin(x,y):
	coef_2 = np.polyfit(x,y,3)
	xmin   = (-2*coef_2[1] + np.sqrt(4*coef_2[1]**2-12*coef_2[0]*coef_2[2]) )/(6*coef_2[0])
	return xmin, coef_2


def Bulk_mod(xipt,E,is_x_vol=False,verbose=False):

	"""
		@E:
			energy in unit of eV

		@xpit : 
			if is_x_vol == True , xipt is the volumn in unit of (\AA^3)
		            	== False, xipt is the cubic length in unit of (\AA)
		
		@details : 
			using 3rd polyfit with V-E 
	"""

	V = None
	if is_x_vol :
		V = xipt
	else:
		V = xipt**3

	#fitting with 3rd poly and get V_min
	V_min , coef_2  = Fit_getMin(V,E) 
	
	# calc Bulk modulus
	B_eva = (6*coef_2[0]*V_min+2*coef_2[1])*V_min # unit = eV/A^3
	B_pa  = B_eva * 1.6E+11 # J/m^3 = Pa
	B_gpa = B_pa * 10**-9  # GPa

	if verbose:
		print ("=========================")
		print ("Bulk modulus calculation :")
		print ("fit: ax^3 + bx^2 + cx + d")
		print ("a: %10.10f"%(coef_2[0]))
		print ("b: %10.10f"%(coef_2[1]))
		print ("c: %10.10f"%(coef_2[2]))
		print ("d: %10.10f"%(coef_2[3]))
		print ("Minimum: ")
		print ("min_x : %10.10f"%(V_min))
		print ("Bulk Modulo:")
		print ("[eV/A^3] %10.10f"%(B_eva))
		print ("[Pa    ] %10.10f"%(B_pa ))
		print ("[GPa   ] %10.10f"%(B_gpa))
		print ("=========================")
	return V_min , B_gpa , coef_2



if len(sys.argv) < 3 :
    print ("ERROR <ev.dat> <ID>") 
    exit(1)

ID = sys.argv[2]

f = open(sys.argv[1],'r')
lines = f.readlines()
dat = []
for line in lines :
    line = line.split('\n')[0]
    line = line.split(' ')
    dat.append(np.array(line,dtype=np.float))
dat = np.array(dat)

## fitting and get a0_min:
a0_min , coef = Fit_getMin(dat[:,0],dat[:,1])
x = np.arange(dat[0,0],dat[-1,0],0.0001)
print (coef)
print ("a0_min : ",a0_min)

## get Bulk modulus :
V = dat[:,0]**3
E = dat[:,1]
V_min, B , coef_2 = Bulk_Mod(V,E,is_x_vol=True,verbose=True)
x2 = np.arange(V[0],V[-1],0.0001)

f = open("%s.res"%(ID),'w')
f.write("a0_min = %4.8lf\n"%(a0_min))
f.write("V_min = %4.8lf\n"%(V_min))
f.write("B = %4.8lf GPa\n"%(B))
f.write("[E-a0] coef = ")
for a in coef:
	f.write("%4.10lf "%(a))
f.write('\n')
f.write("[E-V ] coef = ")
for a in coef_2:
    f.write("%4.10lf "%(a))
f.write('\n')
f.close()



plt.figure(1)
plt.plot(dat[:,0],dat[:,1],'x',label='simulate')
plt.plot(x,np.polyval(coef,x),'-',label = 'fit')
plt.xlabel('a0')
plt.ylabel('E')
plt.grid()
plt.legend()
#plt.savefig("%s_E-a0"%(ID),format='png',transparent=True,dpi=100)

plt.figure(2)
plt.plot(V,E,'x',label='simulate')
plt.plot(x2,np.polyval(coef_2,x2),'-',label = 'fit')
plt.xlabel('V')
plt.ylabel('E')
plt.grid()
plt.legend()
#plt.savefig("%s_E-V"%(ID),format='png',transparent=True,dpi=100)

plt.show()
