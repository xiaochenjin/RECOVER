import os
import numpy as np
import random
import KNN_track_AA
import KNN_track_AB
from RECOVER import read_xyz_fkr

def obtain_pair_index(pair_type,species_list):
	pair_name_list = pair_type.split('-')
	species_1 = pair_name_list[0]
	species_2 = pair_name_list[1]
	pair_name = species_1 + species_2
	index_1 = species_list.index(species_1)
	index_2 = species_list.index(species_2)
	pair_index = [index_1,index_2]
	return pair_name,pair_index

def write(output,R_list,GR_trueKNN_totlist,K_raw):
	with open (output,"w") as f1:
		for i in range(len(R_list)):
			f1.write(str(R_list[i])+' ')
			GR_trueKNN_list = GR_trueKNN_totlist[i]
			for n_shell in range(K_raw):
				GR_trueKNN = GR_trueKNN_list[n_shell]
				if n_shell < K_raw - 1:
					f1.write(str(GR_trueKNN)+' ')
				if n_shell == K_raw - 1:
					f1.write(str(GR_trueKNN)+'\n')
	f1.close()

def compute(index,folder_name,true_KNN_infor_list,K_raw,N_bin,neigh_dis,lc,direction,pair_type):
	input_name = os.path.join(folder_name,"perturb-index-{0}.xyz".format(index))
	N_total,species_list,frac_position_totlist,cell_geometry = read_xyz_fkr.read(input_name)

	if pair_type == 'all':
		pair_name = pair_type
		frac_all_list = []
		for index_species in range(len(frac_position_totlist)):
			frac_position_list = frac_position_totlist[index_species]
			for frac_position in frac_position_list:
				frac_all_list.append(frac_position)

		R_list,GR_trueKNN_totlist = KNN_track_AA.compute(frac_all_list,cell_geometry,true_KNN_infor_list,K_raw,N_bin,neigh_dis,lc,direction)

	else:
		pair_name,pair_index = obtain_pair_index(pair_type,species_list)
		if pair_index[0] == pair_index[1]:
			R_list,GR_trueKNN_totlist = KNN_track_AA.compute(frac_position_totlist[pair_index[0]],cell_geometry,true_KNN_infor_list,K_raw,N_bin,neigh_dis,lc,direction)

		else:	
			R_list,GR_trueKNN_totlist  = KNN_track_AB.compute(frac_position_totlist[pair_index[0]],frac_position_totlist[pair_index[1]],cell_geometry,true_KNN_infor_list,K_raw,N_bin,neigh_dis,lc,direction)

	output = os.path.join(folder_name,'pair-correlation-{0}'.format(pair_name)+'-true-KNN-direction-{0}'.format(direction)+'-perturb-index-{0}.txt'.format(index))
	write(output,R_list,GR_trueKNN_totlist,K_raw)



