import numpy as np 
import os,sys
import matplotlib.pyplot as plt

paths = []
paths.append("")

dat = []
for p in paths :
	f = open(p,'r')
	lines = f.readlines()
	for line in lines:
		if '#' in line:	
			continue
		else:
			line = 



