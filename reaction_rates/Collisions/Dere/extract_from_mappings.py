import numpy as np
from astro import constants as c
import math

f = open('/Users/mauerhof/Documents/mappings/lab/data/ionisation/COLLDAT2.txt','r')

lines = f.readlines()
f.close()

Si_lines = lines[139:153]

# ~ print Si_lines
Eion = np.zeros(14)
Tmin = np.zeros(14)
Si_rates = np.zeros((14,20))


for i in range(14):
	test = Si_lines[i].split()
	Eion[i] = float(test[4].replace('D','e'))
	Tmin[i] = float(test[5].replace('D','e'))
	for j in range(20):
		Si_rates[i,j] = float(test[6+j].replace('D','e'))
		
		
x = np.array([0.04, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.90, 0.93])

T = np.zeros_like(Si_rates)
for i in range(14):
	for j in range(20):
		T[i,j] = Eion[i]*c.ev_to_erg/c.kb_cgs*(math.exp(math.log(2)/(1-x[j])) - 2)
	
print T[1,:]
		

# ~ data = open('/Users/mauerhof/Documents/krome/data/Dere.dat','w')

# ~ data.write('# Data for collisional ionizations, from Dere 2007 \n')





# ~ test = np.genfromtxt('/Users/mauerhof/Documents/mappings/lab/data/ionisation/COLLDAT2.txt', skip_header=139, skip_footer=376, deletechars='D')
# ~ print test[:,0].size
# ~ print test[0,:].size
# ~ print test[0,:]
