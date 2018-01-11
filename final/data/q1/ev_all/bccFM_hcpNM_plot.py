import matplotlib.pyplot as plt
import matplotlib
import os,sys
import numpy as np 
matplotlib.rcParams['text.usetex'] = True

def Fit_getMin(x,y):
	coef_2 = np.polyfit(x,y,3)
	xmin   = (-2*coef_2[1] + np.sqrt(4*coef_2[1]**2-12*coef_2[0]*coef_2[2]) )/(6*coef_2[0])
	return xmin, coef_2

def Fit_div(x,y):
	coef_2 = np.polyfit(x,y,3)
	coef_2[0] *= 3
	coef_2[1] *= 2
	coef_2 = np.delete(coef_2,3)
	return coef_2

def Interxct(c1,c2):
	
	a1 = -(c1[1]-c2[1]) 
	a2 = a1

	c = np.sqrt( (c1[1] - c2[1])**2 - 4*(c1[0]-c2[0])*(c1[2]-c2[2]) )
	
	a1 += c 
	a2 -= c
	a1 /= (2*(c1[0]-c2[0]))
	a2 /= (2*(c1[0]-c2[0]))


	return a1 , a2

	

def Bulk_Mod(xipt,E,is_x_vol=False,verbose=False):

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



ID = "bcchcp"
FM_DATPATH = "ev.bcc.FM.dat"
NM_DATPATH = "ev.hcp.NM.dat"

f = open(FM_DATPATH,'r')
lines = f.readlines()
datFM = []
for line in lines :
    line = line.split('\n')[0]
    line = line.split(' ')
    datFM.append(np.array(line,dtype=np.float))
datFM = np.array(datFM)
f.close()

## get Bulk modulus :
V_FM = datFM[:,0]**3/2
E_FM = datFM[:,1]
M_FM = datFM[:,2]
V_FM_min, B_FM , coef_FM = Bulk_Mod(V_FM,E_FM,is_x_vol=True,verbose=True)
x_FM = np.arange(V_FM[0],V_FM[-1],0.0001)
coef_div_FM = -Fit_div(V_FM,E_FM)

#------------------------------------------------

f = open(NM_DATPATH,'r')
lines = f.readlines()
datNM = []
for line in lines :
    line = line.split('\n')[0]
    line = line.split(' ')
    datNM.append(np.array(line,dtype=np.float))
datNM = np.array(datNM)
f.close()

## get Bulk modulus :
ca_NM = 1.5860
V_NM = datNM[:,0]**3*(np.sqrt(3)/2*ca_NM)/2
E_NM = datNM[:,1]/2
V_NM_min, B_NM , coef_NM = Bulk_Mod(V_NM,E_NM,is_x_vol=True,verbose=True)
x_NM = np.arange(V_NM[0],V_NM[-1],0.0001)
coef_div_NM = -Fit_div(V_NM,E_NM)


xpt1,xpt2 = Interxct(coef_div_NM,coef_div_FM)
print ("Pc = %f"%(np.polyval(coef_div_FM,xpt2)))



Fsz = 12
plt.figure(1)
plt.title("%s"%(ID),fontsize=Fsz)

plt.plot(V_FM,E_FM,'x',label='bcc-FM')
plt.plot(x_FM,np.polyval(coef_FM,x_FM),'-',label = 'fit bcc-FM')
plt.plot(V_NM,E_NM,'x',label='hcp-NM')
plt.plot(x_NM,np.polyval(coef_NM,x_NM),'-',label = 'fit hcp-NM')
plt.yticks(fontsize=Fsz)
plt.ylabel('E (eV/Atom)',fontsize=Fsz)
plt.legend()
plt.grid()
plt.xlabel(r'V ($\AA^3$/Atom)',fontsize=Fsz)
plt.xticks(fontsize=Fsz)
plt.savefig("%sFMNM.jpg"%(ID),format='jpeg',dpi=200)

plt.figure(2)
plt.title("%s slope"%(ID),fontsize=Fsz)
xx = np.arange(7,15,0.001)
plt.plot(xx,np.polyval(coef_div_FM,xx),'-',label = ' bcc-FM')
plt.plot(xx,np.polyval(coef_div_NM,xx),'-',label = ' hcp-NM')
plt.plot(xpt2,np.polyval(coef_div_FM,xpt2),'ro')
plt.yticks(fontsize=Fsz)
plt.ylabel(r'P (eV/$\AA^3$)',fontsize=Fsz)
plt.legend()
plt.grid()
plt.xlabel(r'V ($\AA^3$/Atom)',fontsize=Fsz)
plt.xticks(fontsize=Fsz)
plt.savefig("%sP-V.jpg"%(ID),format='jpeg',dpi=200)
#plt.savefig("%sFMNM.jpg"%(ID),format='jpeg',dpi=200)
plt.show()
