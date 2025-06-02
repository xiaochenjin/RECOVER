import os
import numpy as np
from RECOVER import read_xyz
from RECOVER import write_xyz
from RECOVER import perturb
from optparse import OptionParser


#perturbed data
#GET ARGUMENTS
parser = OptionParser()
parser.add_option('--perturb_type', type = str,default = 'iso',help = 'direction (default: %default)') #r, z, xy
parser.add_option('--mu', type = str,default = '0',help = 'mean (default: %default)')
parser.add_option('--sigma', type = str,default = '0.5',help = 'standard deviation (default: %default)')
parser.add_option('--sigma_xy', type = str,default = '0.5' ,help = 'standard deviation (default: %default)')
parser.add_option('--sigma_z', type = str,default = '0.5' ,help = 'standard deviation (default: %default)')
parser.add_option('--collect_efficiency', type = float,default = 1.0,help = 'collection efficiency (default: %default)')
parser.add_option('--N_config', type = int,default = 100,help = 'number of perturbed structures  (default: %default)')
(options, args) = parser.parse_args()

perturb_type = options.perturb_type
mu = options.mu
collect_eff = options.collect_efficiency
N_config = options.N_config

#basic structure info
input_name = "../benchmark-structures/GeSn-random-relaxed.xyz"
N_total,species_list,old_frac_position_totlist,cell_geometry = read_xyz.read(input_name)
#print(species_list)
#print(len(old_frac_position_totlist[0]),len(old_frac_position_totlist[1]))

old_infor_list = old_frac_position_totlist

#old_infor_list = []
#for index_species in range(len(old_frac_position_totlist)):
#	frac_position_list = old_frac_position_totlist[index_species]
#	for frac_position in frac_position_list:
#		old_infor_list.append(frac_position)

index_list = np.arange(N_config)

analysis = '../simulated-structures'
if not os.path.exists(analysis):
	os.mkdir(analysis)
os.chdir(analysis)

if perturb_type == 'iso':
	next_folder = 'isotropic-perturbation'
	sigma = options.sigma
	directory_name = "mu-"+str(mu)+"-sgm-"+str(sigma)

if perturb_type == 'aniso':
	next_folder = 'anisotropic-perturbation'
	sigma_z = options.sigma_z
	sigma_xy = options.sigma_xy
	directory_name = "mu-"+str(mu)+"-sgm_z-"+str(sigma_z)+"-sgm_xy-"+str(sigma_xy)

if not os.path.exists(next_folder):
	os.mkdir(next_folder)
os.chdir(next_folder)

if not os.path.exists("efficiency-{0}".format(collect_eff)):
	os.mkdir("efficiency-{0}".format(collect_eff))
os.chdir("efficiency-{0}".format(collect_eff))

if not os.path.exists(directory_name):
	os.mkdir(directory_name)
os.chdir(directory_name)

for index in index_list:
	if perturb_type == 'iso':
		r_list,new_infor_list = perturb.iso_perturb(old_infor_list,cell_geometry,mu, sigma,collect_eff)

		with open ("displacement-list-index-{0}.txt".format(index),"w") as f1:
			for i in range(len(r_list)):
				f1.write(str(r_list[i])+'\n')
		f1.close()

	if perturb_type == 'aniso':
		r_z_list,r_xy_list,new_infor_list = perturb.aniso_perturb(old_infor_list,cell_geometry,mu, sigma_z,sigma_xy,collect_eff)

		with open ("xy-displacement-list-index-{0}.txt".format(index),"w") as f1:
			for i in range(len(r_xy_list)):
				f1.write(str(r_xy_list[i])+'\n')
		f1.close()

		with open ("z-displacement-list-index-{0}.txt".format(index),"w") as f1:
			for i in range(len(r_z_list)):
				f1.write(str(r_z_list[i])+'\n')
		f1.close()

	file_name = "perturb-index-{0}.xyz".format(index)
	write_xyz.write_perturb(new_infor_list,cell_geometry,file_name)


