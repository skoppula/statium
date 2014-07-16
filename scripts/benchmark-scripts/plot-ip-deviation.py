import os
import time
import matplotlib.pyplot as plt
import numpy as np

my_ips_totals = list()
only_in_mine = list()
same = list()
only_in_joes = list()
joe_ips_totals = list()
with open('ip_deviation.out','r') as f:
	a = f.read().split('\n')
	for thing in a:
		if thing == '': continue
		parts = thing.split()
		my_ips_totals.append(float(parts[0]))
		only_in_mine.append(float(parts[1]))
		same.append(float(parts[2]))
		only_in_joes.append(float(parts[3]))
		joe_ips_totals.append(float(parts[4]))


props = [thing/float(joe_ips_totals[i]) for i,thing in enumerate(same)]
print sum(props)/len(props)
#0.935166190356

print sum(my_ips_totals)/len(my_ips_totals)
print sum(only_in_mine)/len(only_in_mine)
print sum(same)/len(same)
print sum(only_in_joes)/len(only_in_joes)
print sum(joe_ips_totals)/len(joe_ips_totals)

#1148.74
#177.77
#970.97
#75.44
#1046.41

