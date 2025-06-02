import numpy as np

f = open('/Users/mauerhof/Documents/krome/data/database/collisional_ionization.dat','r')

a = f.readlines()

# ~ dE,P,A,X,K = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
coeff = np.zeros((5,403))

#dE, P, A, X, K
for i in range(403):
	for j in range(2,6):
		coeff[j-2,i] = float(a[8+i*7].split()[j].replace(',','').replace('d','e'))
	coeff[4,i] = float(a[8+i*7].split()[6].replace(')','').replace('d','e'))
	
	
	
	
np.savetxt('clist',np.transpose(coeff))


