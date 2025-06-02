import os
import numpy as np
import random
import KNN_info_AA
import KNN_info_AB
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

def write(output,kNN_infor_list,K_raw):
	with open (output,"w") as f1:
		for i in range(len(kNN_infor_list)):
			kNN_infor = kNN_infor_list[i]
			old_index = kNN_infor[0]
			kNN_index = kNN_infor[1]
			f1.write(str(old_index)+' ')
			for i in range(K_raw):
				if not len(kNN_index[i]) == 0:
					kNN_label = 'x'.join([str(x) for x in kNN_index[i]]) #original indices of atoms in 1st shell to atom i
				else:
					kNN_label = 'nan'
				if i < K_raw - 1:
					f1.write(str(kNN_label)+' ')
				if i == K_raw - 1:
					f1.write(str(kNN_label)+'\n')
		f1.close()

def compute(input_name,folder_name,output_name,pair_type,K_raw,neigh_dis,lc):
	#input_name = 'denser-cube-Ge-raw-APT.xyz'
	N_total,species_list,frac_position_totlist,cell_geometry = read_xyz.read(os.path.join(folder_name,input_name))
#	  print(frac_position_Ge_list)

	if pair_type == 'all':
		frac_all_list = []
		for index_species in range(len(frac_position_totlist)):
			frac_position_list = frac_position_totlist[index_species]
			for frac_position in frac_position_list:
				frac_all_list.append(frac_position)

		kNN_infor_list = KNN_info_AA.compute(frac_all_list,cell_geometry,K_raw,neigh_dis,lc)

	else:
		pair_name,pair_index = obtain_pair_index(pair_type,species_list)
		if pair_index[0] == pair_index[1]:
			kNN_infor_list = KNN_info_AA.compute(frac_position_totlist[pair_index[0]],cell_geometry,K_raw,neigh_dis,lc)
		else:
			kNN_infor_list = KNN_info_AB.compute(frac_position_totlist[pair_index[0]],frac_position_totlist[pair_index[1]],cell_geometry,K_raw,neigh_dis,lc)

	output = os.path.join(folder_name,output_name)
	write(output,kNN_infor_list,K_raw)
