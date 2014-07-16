import os
import time
import matplotlib.pyplot as plt
import numpy as np

tetra = list()
scwrl = list()
pairwise = list()
with open('tetra-scwrl-errors.out','r') as f:
	a = f.read().split('\n')
	for thing in a:
		if thing == '': continue
		parts = thing.split()
		tetra.append(float(parts[0]))
		scwrl.append(float(parts[1]))
		pairwise.append(float(parts[2]))

tetra.sort()
scwrl.sort()
pairwise.sort()

hist, bins = np.histogram(tetra[:-800], bins=50)
width = 0.7 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2
plt.title('tetra-algorithm error')
plt.xlabel('error (Angstrom)')
plt.ylabel('frequency')
#plt.axis([0,1,0,3000])
plt.bar(center, hist, align='center', width=width)
plt.show()


hist, bins = np.histogram(scwrl[:-500], bins=50)
width = 0.7 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2
plt.title('scwrl-algorithm errors')
plt.xlabel('error (Angstrom)')
plt.ylabel('frequency')
#plt.axis([0,1,0,1000])
plt.bar(center, hist, align='center', width=width)
plt.show()

hist, bins = np.histogram(pairwise[500:], bins=50)
width = 0.1 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2
plt.title('difference of errors between two algorithms')
plt.xlabel('error (Angstrom)')
plt.ylabel('frequency')
plt.axis([-1,1,0,2500])
plt.bar(center, hist, align='center', width=width)
plt.show()
print sum(pairwise[500:])/len(pairwise[500:])
		
