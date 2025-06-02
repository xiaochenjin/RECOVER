import os
import numpy as np
import random
import generate
import find_dK
from RECOVER import read_xyz
from optparse import OptionParser

def write(output,dK_totlist):
	with open(output,'w') as f1:
		for n_shell in range(K_raw):
			dK_list = dK_totlist[n_shell]
			N_dK = len(dK_list)
			if N_dK ==0:
				continue
			f1.write('{0} NN '.format(n_shell+1))
			for i in range(N_dK):
#				dK = dK_list[i]	
				dK = round(dK_list[i],3)
				if i < N_dK - 1 :
					f1.write(str(dK)+'x')
				if i == N_dK- 1:
					f1.write(str(dK)+'\n')
	f1.close()
				

#GET ARGUMENTS
parser = OptionParser()
parser.add_option('--K_raw', type = int,default = 50,help = 'number of shells for computing fk(r) (default: %default)')
parser.add_option('--lattice_constant', type = float,default = 5.86,help = 'direction (default: %default)') #'r','z','xy'
parser.add_option('--direction', type = str,default = 'r',help = 'direction (default: %default)') #r, z, xy
(options, args) = parser.parse_args()
K_raw = options.K_raw
lc = options.lattice_constant
direction = options.direction

input_name = 'test.xyz'
folder_name = '../benchmark-structures/'
output_name = 'KNN-peak-position-direction-{0}.txt'.format(direction)
output = os.path.join(folder_name,output_name)


#read nearest neighbor locations:
with open ('../benchmark-structures/Peak-position-direction-r.txt') as f1:
	neigh_dis = np.loadtxt(f1, delimiter=' ', usecols=(1),unpack=True)

neigh_dis = neigh_dis.tolist()

N_total,species_list,frac_position_totlist,cell_geometry = read_xyz.read(os.path.join(folder_name,input_name))

frac_all_list = []
for index_species in range(len(frac_position_totlist)):
	frac_position_list = frac_position_totlist[index_species]
	for frac_position in frac_position_list:
		frac_all_list.append(frac_position)

dK_totlist = find_dK.compute(frac_all_list,cell_geometry,K_raw,neigh_dis,lc,direction)
write(output,dK_totlist)




