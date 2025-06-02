import os
import numpy as np
import random
import PairCorrelation_AA
import PairCorrelation_AB
from RECOVER import read_xyz

def obtain_pair_index(pair_type,species_list):
	pair_name_list = pair_type.split('-')
	species_1 = pair_name_list[0]
	species_2 = pair_name_list[1]
	pair_name = species_1 + species_2
	index_1 = species_list.index(species_1)
	index_2 = species_list.index(species_2)
	pair_index = [index_1,index_2]
	return pair_name,pair_index

def write(output,R_list,GR_list):
	with open (output,"w") as f1:
		for i in range(len(R_list)):
			f1.write(str(R_list[i])+' '+str(GR_list[i])+'\n')
	f1.close()

def compute(index,data_folder,pair_list,N_bin,direction):
	z_index = index[0]
	xy_index = index[1]
	print("z index: ",z_index)
	print("xy index: ",xy_index)
	folder_name = os.path.join(data_folder,"z-index-{0}".format(z_index))	
	input_name = os.path.join(folder_name,"cube-index-{0}.xyz".format(xy_index))
	N_total,species_list,frac_position_totlist,cell_geometry = read_xyz.read(input_name)

	for pair_type in pair_list:
		if not pair_type == 'all': pair_name,pair_index = obtain_pair_index(pair_type,species_list)
		if pair_type == 'all':
			pair_name = pair_type
			frac_all_list = []
			for index_species in range(len(frac_position_totlist)):
				frac_position_list = frac_position_totlist[index_species]
				for frac_position in frac_position_list:
					frac_all_list.append(frac_position)
			R_list,GR_list = PairCorrelation_AA.compute(frac_all_list,cell_geometry,N_bin,direction)

		else:
			if pair_index[0] == pair_index[1]:
				R_list,GR_list = PairCorrelation_AA.compute(frac_position_totlist[pair_index[0]],cell_geometry,N_bin,direction)

			else:
				R_list,GR_list = PairCorrelation_AB.compute(frac_position_totlist[pair_index[0]],frac_position_totlist[pair_index[1]],cell_geometry,N_bin,direction)

		output = os.path.join(folder_name,'pair-correlation-{0}'.format(pair_name)+'-direction-{0}'.format(direction)+'-perturb-index-{0}.txt'.format(xy_index))
		write(output,R_list,GR_list)



