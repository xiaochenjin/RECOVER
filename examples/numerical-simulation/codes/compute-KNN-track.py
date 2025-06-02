import os
import numpy as np
import random
from RECOVER import cal_KNN_track
from optparse import OptionParser
from multiprocessing import Pool
from functools import partial

def extract(file_name):
	infile = open(file_name,'r')
	data = infile.readlines()
	true_KNN_infor_list = []
	for line in data:
		index = int(line.split()[0])
		true_KNN_index_list = []
		for n_shell in range(K_raw):
			if not line.split()[n_shell+1] == 'nan':
				kNN_index = [int(x) for x in line.split()[n_shell+1].split('x')]
			else:
				kNN_index = []
			true_KNN_index_list.append(kNN_index)
		true_KNN_infor = [index,true_KNN_index_list]
		true_KNN_infor_list.append(true_KNN_infor)
	return true_KNN_infor_list

#GET ARGUMENTS
parser = OptionParser()
parser.add_option('--pair_type', type = str,default = 'Sn-Sn',help = 'pair type for computing fk(r) (default: %default)')
parser.add_option('--K_raw', type = int,default = 3,help = 'number of shells for computing fk(r) (default: %default)')
parser.add_option('--lattice_constant', type = float,default = 5.86,help = 'direction (default: %default)')
parser.add_option('--direction', type = str,default = 'r',help = 'direction (default: %default)') #r, z, xy
parser.add_option('--perturb_type', type = str,default = 'iso',help = 'direction (default: %default)') 
parser.add_option('--mu', type = str,default = '0',help = 'mean (default: %default)')
parser.add_option('--sigma', type = str,default = '0.5',help = 'standard deviation (default: %default)')
parser.add_option('--sigma_xy', type = str,default = '0.5' ,help = 'standard deviation (default: %default)')
parser.add_option('--sigma_z', type = str,default = '0.5' ,help = 'standard deviation (default: %default)')
parser.add_option('--collect_efficiency', type = str,default = '1.0',help = 'collection efficiency (default: %default)')
parser.add_option('--N_config', type = int,default = 5,help = 'number of perturbed structures  (default: %default)')
parser.add_option('--N_bin', type = int,default = 10000,help = 'number of bins to compute   (default: %default)')
parser.add_option('--N_thread', type = int,default = 1,help = 'number of threads   (default: %default)')

(options, args) = parser.parse_args()
pair_type = options.pair_type
K_raw = options.K_raw
lc = options.lattice_constant
direction = options.direction
mu = options.mu
perturb_type = options.perturb_type
collect_eff = options.collect_efficiency

N_config = options.N_config
N_thread = options.N_thread
N_bin = options.N_bin

folder_name = '../benchmark-structures/'
true_KNN_file = os.path.join(folder_name,'KNN-info-{0}.txt'.format(pair_type))
true_KNN_infor_list = extract(true_KNN_file)
#print("original info: ")
#print(true_KNN_infor_list[:5])

#read nearest neighbor locations:
with open (os.path.join(folder_name,'structure-info-file.txt')) as f1:
	neigh_dis = np.loadtxt(f1, delimiter=' ', usecols=(1),unpack=True)

neigh_dis = neigh_dis.tolist()

#analysis = '../simulated-structures/efficiency-{0}'.format(collect_eff)

if perturb_type == 'iso':
	sigma = options.sigma
#	analysis = '../../analysis/recover-test-4/test1-64000-atom/efficiency-{0}'.format(collect_eff)
	analysis = '../simulated-structures/isotropic-perturbation/efficiency-{0}'.format(collect_eff)
	directory_name = "mu-"+str(mu)+"-sgm-"+str(sigma)

if perturb_type == 'aniso':
	sigma_z = options.sigma_z
	sigma_xy = options.sigma_xy
#	analysis = '../../analysis/recover-test-5-anisotropic/test1-64000-atom/efficiency-{0}'.format(collect_eff)
	analysis = '../simulated-structures/anisotropic-perturbation/efficiency-{0}'.format(collect_eff)
	directory_name = "mu-"+str(mu)+"-sgm_z-"+str(sigma_z)+"-sgm_xy-"+str(sigma_xy)

folder_name = os.path.join(analysis,directory_name)

index_list = np.arange(N_config)
if __name__ == '__main__':
	with Pool(N_thread) as pool:
		compute_2 = partial(cal_KNN_track.compute,folder_name=folder_name,true_KNN_infor_list=true_KNN_infor_list,K_raw=K_raw,N_bin=N_bin,neigh_dis=neigh_dis,lc=lc,direction=direction,pair_type=pair_type)
		pool.map(compute_2,index_list)





