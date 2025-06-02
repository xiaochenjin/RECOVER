import os
import numpy as np
from RECOVER import cal_PairCorrelation
from multiprocessing import Pool
from functools import partial
from optparse import OptionParser

parser = OptionParser()
parser.add_option('--what_to_compute', type = str,default = 'original',help = 'what systems to compute (default: %default)')
parser.add_option('--direction', type = str,default = 'r',help = 'direction (default: %default)') #r, z, xy
parser.add_option('--pair_type', type = str,default = 'Sn-Sn',help = 'pair type (default: %default)')
parser.add_option('--whether_relax', type = str,default = 'Y',help = 'wheter to relax structure (default: %default)')
parser.add_option('--perturb_type', type = str,default = 'iso',help = 'direction (default: %default)')
parser.add_option('--mu', type = str,default = '0',help = 'mean (default: %default)')
parser.add_option('--sigma', type = str,default = '0.5',help = 'standard deviation (default: %default)')
parser.add_option('--sigma_xy', type = str,default = '0.5' ,help = 'standard deviation (default: %default)')
parser.add_option('--sigma_z', type = str,default = '0.5' ,help = 'standard deviation (default: %default)')
parser.add_option('--collect_efficiency', type = str,default = '1.0',help = 'collection efficiency (default: %default)')
parser.add_option('--N_config', type = int,default = 100,help = 'number of perturbed structures  (default: %default)')
parser.add_option('--N_bin', type = int,default = 10000,help = 'number of bins to compute   (default: %default)')
parser.add_option('--N_thread', type = int,default = 1,help = 'number of threads   (default: %default)')
(options, args) = parser.parse_args()

what_to_compute = options.what_to_compute
whether_relax = options.whether_relax
direction = options.direction
pair_type = options.pair_type
perturb_type = options.perturb_type
mu = options.mu
collect_eff = options.collect_efficiency
N_config = options.N_config
N_bin = options.N_bin
N_thread = options.N_thread

if what_to_compute == 'original':
	if whether_relax == 'Y':
		input_name = 'GeSn-random-relaxed.xyz'
		output_name = 'Relaxed-pair-correlation-{0}'.format(pair_type)+'-direction-{0}.txt'.format(direction)
	if whether_relax == 'N':
		input_name = 'GeSn-random-before-relax.xyz'
		output_name = 'No-relax-pair-correlation-{0}'.format(pair_type)+'-direction-{0}.txt'.format(direction)

	folder_name = '../benchmark-structures/'
	cal_PairCorrelation.compute_original(input_name,folder_name,output_name,pair_type,N_bin,direction)

if what_to_compute == 'perturbed':
	if perturb_type == 'iso':
		sigma = options.sigma
	#   analysis = '../../analysis/recover-test-4/test1-64000-atom/efficiency-{0}'.format(collect_eff)
		analysis = '../simulated-structures/isotropic-perturbation/efficiency-{0}'.format(collect_eff)
		directory_name = "mu-"+str(mu)+"-sgm-"+str(sigma)

	if perturb_type == 'aniso':
		sigma_z = options.sigma_z
		sigma_xy = options.sigma_xy
	#   analysis = '../../analysis/recover-test-5-anisotropic/test1-64000-atom/efficiency-{0}'.format(collect_eff)
		analysis = '../simulated-structures/anisotropic-perturbation/efficiency-{0}'.format(collect_eff)
		directory_name = "mu-"+str(mu)+"-sgm_z-"+str(sigma_z)+"-sgm_xy-"+str(sigma_xy)

	folder_name = os.path.join(analysis,directory_name)

	index_list = np.arange(N_config)
	if __name__ == '__main__':
		with Pool(N_thread) as pool:
			compute_2 = partial(cal_PairCorrelation.compute_perturb,folder_name=folder_name,pair_type=pair_type,N_bin=N_bin,direction=direction)
			pool.map(compute_2,index_list)





