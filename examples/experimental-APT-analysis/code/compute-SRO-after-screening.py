import os
import numpy as np

def grep_GR(data_name):
	#obtain by folder/plot-truekNN-proportion-ave.py
	infile =  open(data_name,'r')
	data = infile.readlines()
	infile.close()
	R_totlist = []
	GR_totlist = []
	for line in data:
		R = float(line.split()[0])
		GR = float(line.split()[1])
		if R > float(R_min):
			R_totlist.append(R)
			GR_totlist.append(GR)
#   print(GR_totlist)
	return R_totlist,GR_totlist

def obtain_alpha(screened_index_list,sigma):
	file_name = "{0}NN".format(KNN)+"-{0}-adjusted-SRO-parameter".format(pair)+'-sigma-{0}'.format(sigma)+'.txt'

	infile = open(os.path.join(analysis_folder,file_name),'r')
	data = infile.readlines()
	infile.close()

	alpha_list = []
	screened_infor_list = []
	for line in data:
		index_z,index_xy = line.split()[:2]
		index = [index_z,index_xy]
		if index in screened_index_list:
#			print(index)
			alpha = float(line.split()[-1])
			infor = [index_z,index_xy,alpha]
			screened_infor_list.append(infor)
			alpha_list.append(alpha)
#		print(infor)
	return screened_infor_list,alpha_list



def screen(index_totlist,max_GR_criteria):
	screen_index_list = []
	for index in index_totlist:
		index_z,index_xy = index
		data_dir = os.path.join(data_folder,'z-index-{0}'.format(index_z))
		
		max_GR_list = []
		direction = 'r'
		for pair_name in pair_totlist:
			data_pair = os.path.join(data_dir,'pair-correlation-{0}'.format(pair)+'-direction-{0}'.format(direction)+'-perturb-index-{0}'.format(index_xy)+'-smooth-window-{0}.txt'.format(smooth_window))	

			R_list_pair,GR_list_pair = grep_GR(data_pair)	
			GR_list_selected = []		
			for i in range(len(GR_list_pair)):
				R = R_list_pair[i]
				GR = GR_list_pair[i]
				if R < 1:
					GR_list_selected.append(GR)
				
			max_GR = max(GR_list_selected)
			max_GR_list.append(max_GR)

		if max(max_GR_list) <= float(max_GR_criteria):
			screen_index_list.append(index)

	return screen_index_list

def write_infor(screened_infor_list):
	write_folder = os.path.join(analysis_folder,'screen-max-gr-{0}'.format(max_GR_criteria))
	if not os.path.exists(write_folder):
		os.mkdir(write_folder)
	file_name = "{0}NN".format(KNN)+"-{0}-adjusted-SRO-parameter".format(pair)+'-sigma-{0}'.format(sigma)+'.txt'

	with open (os.path.join(write_folder,file_name),'w') as f1:
		for i in range(len(screened_infor_list)):
			infor = screened_infor_list[i]
			for j in range(len(infor)):
				if j < len(infor) - 1:
					f1.write(str(infor[j])+' ')
				else:
					f1.write(str(infor[j])+'\n')
	f1.close()
	

#pair_list = ['SnSn','GeSn','GeGe']
pair = input('pair: ')
KNN = input("Kth shell: ")
sigma_list = ['2.5','2.7','3.0','3.3','3.5','3.7']
max_GR_criteria = input("maximum gr (r < 1A): ")
#pair_totlist = ['all','Sn-Sn']
pair_totlist = ['all','Sn-Sn','Ge-Sn','Ge-Ge']
smooth_window = input("smooth window: ")
box_length = 100 #10 nm
N_bin = 10000
smooth_window_length = float(smooth_window)*box_length/2/N_bin
R_min = smooth_window_length/2
data_folder = '../example-APT-data'
analysis_folder = '../analysis'


#obtain total index
index_totlist = []
file_name = os.path.join(data_folder,"Nz-Nx-Ny.txt")
with open (file_name) as f1:
	Nz_Nx_Ny = np.loadtxt(f1, delimiter=' ',usecols=(0),dtype=int,unpack=True)
Nz, Nx, Ny =  Nz_Nx_Ny[0], Nz_Nx_Ny[1], Nz_Nx_Ny[2]
N_xy = Nx*Ny

index_z_list = [str(x) for x in np.arange(Nz)]
index_xy_list = [str(x) for x in np.arange(N_xy)]

index_totlist = []
for index_z in index_z_list:
	for index_xy in index_xy_list:
		index_xyz = [index_z,index_xy]
		index_totlist.append(index_xyz)

screened_index_list = screen(index_totlist,max_GR_criteria)
print("Number of data points: ", len(screened_index_list))
if len(screened_index_list) > 0:
	for sigma in sigma_list:
		screened_infor_list,alpha_list = obtain_alpha(screened_index_list,sigma)	
		print("sigma, pSRO: ", sigma, round(np.average(alpha_list),3))
		write_infor(screened_infor_list)
	
		


