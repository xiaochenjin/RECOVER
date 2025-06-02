import os
import numpy as np

def find_converge(pSRO_list): #find converged pSRO
	for i in reversed(range(len(pSRO_list))):
		diff1 = (pSRO_list[i]-pSRO_list[i-1])/pSRO_list[i]
		diff2 = (pSRO_list[i-1]-pSRO_list[i-2])/pSRO_list[i-1]
		if diff1 <= 0.1 and diff2 <= 0.1:
			pSRO_converged = pSRO_list[i]
			break

	if not 'pSRO_converged' in locals():
		pSRO_converged = pSRO_list[-1]

	return pSRO_converged

def write(infor_list):
	output_file = "{0}NN".format(KNN)+"-{0}-adjusted-SRO-parameter".format(pair)+'-sigma-{0}'.format(str(sigma))+'.txt'

	with open (os.path.join(output_folder1,output_file),'w') as f1:
		for i in range(len(infor_list)):
			infor = infor_list[i]
			for j in range(len(infor)):
				if j < len(infor) - 1:
					f1.write(str(infor[j])+' ')
				else:
					f1.write(str(infor[j])+'\n')
	f1.close()


pair = input('pair: ')
KNN = input("Kth shell: ")
sigma_list = ['2.5','2.7','3.0','3.3','3.5','3.7']

output_folder1 = '../analysis'

for sigma in sigma_list:
	pSRO_new_list = []
	infor_list = []

	file_name = "../example-APT-data/Nz-Nx-Ny.txt"
	with open (file_name) as f1:
		Nz_Nx_Ny = np.loadtxt(f1, delimiter=' ',usecols=(0),dtype=int,unpack=True)
	Nz, Nx, Ny =  Nz_Nx_Ny[0], Nz_Nx_Ny[1], Nz_Nx_Ny[2]
	N_xy = Nx*Ny

	index_z_list = [str(x) for x in np.arange(Nz)]
	index_xy_list = [str(x) for x in np.arange(N_xy)]

	index_list = []
	for index_z in index_z_list:
		for index_xy in index_xy_list:
			index_xyz = [index_z,index_xy]
			index_list.append(index_xyz)

	for index in index_list:
		index_z = index[0]
		index_xy = index[1]
		output_folder2 = 'z-index-{0}'.format(index_z)
		output_folder = os.path.join(output_folder1,output_folder2)

		file_name_all = os.path.join(output_folder,'index-{0}'.format(index_xy)+'-{0}NN-'.format(KNN)+'{0}-SRO-convergence'.format('all')+'-sigma-{0}'.format(str(sigma))+'.txt')
		file_name_pair = os.path.join(output_folder,'index-{0}'.format(index_xy)+'-{0}NN-'.format(KNN)+'{0}-SRO-convergence'.format(pair)+'-sigma-{0}'.format(str(sigma))+'.txt')

		if os.path.exists(file_name_all) and os.path.exists(file_name_pair):
			
			with open (file_name_all) as f1:
				pSRO_all_list = np.loadtxt(f1, delimiter=' ',usecols=(1),unpack=True)
			pSRO_all = find_converge(pSRO_all_list)
#				pSRO_entire_list.append(pSRO_all)
#				print("pSRO all: ",pSRO_all)

			with open (file_name_pair) as f1:
				pSRO_pair_list = np.loadtxt(f1, delimiter=' ',usecols=(1),unpack=True)
			pSRO_pair_old = find_converge(pSRO_pair_list)
#				pSRO_old_list.append(pSRO_pair_old)
#				print("pSRO pair (old): ",pSRO_pair_old)

			pSRO_pair_new = 1 - (1-pSRO_pair_old)/(1-pSRO_all)
			pSRO_new_list.append(pSRO_pair_new)
#				print("pSRO pair (new): ",pSRO_pair_new)

			infor = [index_z,index_xy,pSRO_all,pSRO_pair_old,pSRO_pair_new]
			infor_list.append(infor)		

	if  len(pSRO_new_list) > 0:	
		print("sigma: ", sigma)
		print("# of SRO parameter: ", len(pSRO_new_list))
		print("average SRO parameter: ", round(np.average(pSRO_new_list),3))
		print("standard error: ", round(np.std(pSRO_new_list)/np.sqrt(len(pSRO_new_list)),3))
		write(infor_list)

	else:
		print("Not there yet~")
		break
		
	print('\n')
