import numpy as np
import os
import pandas as pd

def write(R_list,ave_fkr_totlist):
	ave_fkr_totlist_T = np.transpose(ave_fkr_totlist)
	with open (os.path.join(folder_name,output_name),"w") as f1:
		for i in range(len(R_list)):
			f1.write(str(R_list[i])+' ')
			fkr_list = ave_fkr_totlist_T[i]
			for n_shell in range(K_raw):
				fkr = fkr_list[n_shell]
				if n_shell < K_raw - 1:
					f1.write(str(fkr)+' ')
				if n_shell == K_raw - 1:
					f1.write(str(fkr)+'\n')
	f1.close()

def extract(input_name):
	with open (input_name) as f1:
		df = pd.read_csv(input_name,delimiter=' ', header=None)
	R_list = list(df.iloc[:, 0])
	GR_trueKNN_totlist = []
	for n_shell in range(K_raw):
		GR_trueKNN_list = list(df.iloc[:, n_shell+1])
		GR_trueKNN_totlist.append(GR_trueKNN_list)
	return R_list,GR_trueKNN_totlist

def obtain_fkr(R_list,GR_trueKNN_totlist,direction):
#	lx = np.linalg.norm(a)
#	ly = np.linalg.norm(b)
#	lz = np.linalg.norm(c)
	dr = R_list[1] - R_list[0]
	fkr_totlist = []
	for n_shell in range(K_raw):
		Denominator = 0
		GR_trueKNN_list = GR_trueKNN_totlist[n_shell]
		for i in range(len(R_list)):
			R = R_list[i]
			GR_KNN = GR_trueKNN_list[i]
#			if direction == 'x': dV = 2*ly*lz*dr
#			if direction == 'y': dV = 2*lx*lz*dr
#			if direction == 'z': dV = 2*lx*ly*dr
#			if direction == 'xy':dV = 2*3.14*R*lz*dr
			if direction == 'r': dV = 4*3.14*R*R*dr
			Ni = GR_KNN*dV
			Denominator += Ni

		fkr_list = []
		for i in range(len(R_list)):
			R = R_list[i]
			GR_KNN = GR_trueKNN_list[i]
#			if direction == 'x': A_r = 2*ly*lz
#			if direction == 'y': A_r = 2*lx*lz
#			if direction == 'z': A_r = 2*lx*ly
#			if direction == 'xy':A_r = 2*3.14*R*lz
			if direction == 'r': A_r = 4*3.14*R*R
			fkr = A_r*GR_KNN/Denominator
			fkr_list.append(fkr)
		fkr_totlist.append(fkr_list)
	return fkr_totlist

direction = 'r'
pair_type = input("pair: ")
K_raw = int(input("Number of shells to compute fk(r): "))
N_config = int(input("Number of configurations: "))
collect_eff = input("collect efficiency: ")
index_list = np.arange(N_config)
perturb_type = input("perturb type: ")
mu = 0 #input("mean: ")

output_name = 'fkr-file-{0}.txt'.format(pair_type)

if perturb_type == 'iso':
	sigma = input("perturbation length: ")
#   analysis = '../../analysis/recover-test-4/test1-64000-atom/efficiency-{0}'.format(collect_eff)
	analysis = '../simulated-structures/isotropic-perturbation/efficiency-{0}'.format(collect_eff)
	directory_name = "mu-"+str(mu)+"-sgm-"+str(sigma)

if perturb_type == 'aniso':
	sigma_z = input("perturbation length (z): ")
	sigma_xy = input("perturbation length (xy): ")
#   analysis = '../../analysis/recover-test-5-anisotropic/test1-64000-atom/efficiency-{0}'.format(collect_eff)
	analysis = '../simulated-structures/anisotropic-perturbation/efficiency-{0}'.format(collect_eff)
	directory_name = "mu-"+str(mu)+"-sgm_z-"+str(sigma_z)+"-sgm_xy-"+str(sigma_xy)

folder_name = os.path.join(analysis,directory_name)
GR_trueKNN_totlist_all = []
for index in index_list:
	file_name = 'pair-correlation-{0}'.format(pair_type)+'-true-KNN-direction-{0}'.format(direction)+'-perturb-index-{0}.txt'.format(index)
	input_name = os.path.join(folder_name,file_name)
	R_list,GR_trueKNN_totlist = extract(input_name)
	GR_trueKNN_totlist_all.append(GR_trueKNN_totlist)

ave_GRtrueKNN_totlist = np.mean(GR_trueKNN_totlist_all,axis=0)
ave_fkr_totlist = obtain_fkr(R_list,ave_GRtrueKNN_totlist,direction)
write(R_list,ave_fkr_totlist)

#>>> list1 = [[1,2],[3,4]]
#>>> list2 = [[5,6],[7,8]]
#>>> list3 = [list1,list2]
#>>> np.mean(list3,axis=0)
#array([[3., 4.],
#	   [5., 6.]])



