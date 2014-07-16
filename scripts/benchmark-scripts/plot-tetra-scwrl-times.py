import os
import time
import matplotlib.pyplot as plt
import numpy as np

tetra = list()
scwrl = list()
with open('tetra-scwrl-times.out','r') as f:
	a = f.read().split('\n')
	for thing in a:
		if thing == '': continue
		parts = thing.split()
		tetra.append(float(parts[0]))
		scwrl.append(float(parts[1]))

hist, bins = np.histogram(tetra, bins=50)
width = 0.7 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2
plt.title('tetra-algorithm times')
plt.xlabel('time (seconds)')
plt.ylabel('frequency')
#plt.axis([0,5,0,2000])
plt.bar(center, hist, align='center', width=width)
plt.show()


hist, bins = np.histogram(scwrl, bins=50)
width = 0.7 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2
plt.title('scwrl-algorithm times')
plt.xlabel('time (seconds)')
plt.ylabel('frequency')
#plt.axis([0,5,0,2000])
plt.bar(center, hist, align='center', width=width)
plt.show()
		
