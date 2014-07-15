from util import *
import os
import pickle

in_dir = 'data/processed_culled_90'
ip_dir = 'data/ip_90_wGLY'

lib_pickles = [os.path.join(in_dir, piql) for piql in os.listdir(in_dir)]
with open('ip_deviation.out', 'w') as outfile:
	for (i, lib_pickle_path) in enumerate(lib_pickles[0:100]):	
		print i, lib_pickle_path
		with open(lib_pickle_path, 'rb') as infile:
			(lib_pdb,lib_ips,lib_distance_matrix) = pickle.load(infile)
			#No CB, no H's, no tetra
			my_ips = set(lib_ips)
			#print sorted(my_ips)
			ip_path = os.path.join(ip_dir, os.path.split(lib_pickle_path)[1].split('.')[0] + '.ip')
			joe_ips = set((int(pair[0]),int(pair[1])) for pair in filelines2deeplist(ip_path,skipComments=True) if pair != [])
			#print sorted(joe_ips)
			d1 = joe_ips.difference(my_ips)
			d2 = my_ips.difference(joe_ips)
			s = my_ips.intersection(joe_ips)
			outfile.write(str(len(my_ips)) + ' ' + str(len(d2)) + ' ' + str(len(s)) + ' ' + str(len(d1)) + str(len(joe_ips)) + '\n')
			print str(len(my_ips)) + ' ' + str(len(d2)) + ' ' + str(len(s)) + ' ' + str(len(d1)) + ' ' + str(len(joe_ips))
