import os
import numpy as np
from optparse import OptionParser

def convert(input_name,output_name):
	with open (input_name) as f1:
		R_list,GR_list_old = np.loadtxt(f1, delimiter=' ',usecols=(0,1),unpack=True)

	GR_list_new = []
	N_data = len(R_list)
	for ind in range(N_data):#[:10]:
#		if ind < smooth_half_window or ind > N_data - smooth_half_window:
#			GR_new = GR_list_old[ind]
		if ind < smooth_half_window:
			GR_new = np.average(GR_list_old[0:ind+smooth_half_window])
		if ind > N_data - smooth_half_window:
			GR_new = np.average(GR_list_old[ind-smooth_half_window:N_data])
		if smooth_half_window <= ind <= N_data - smooth_half_window:
			GR_new = np.average(GR_list_old[ind-smooth_half_window:ind+smooth_half_window+1])
		
		GR_list_new.append(GR_new)

	with open (output_name,"w") as f1:
		for i in range(len(R_list)):
			f1.write(str(R_list[i])+' '+str(GR_list_new[i])+'\n')
	f1.close()

#GET ARGUMENTS
parser = OptionParser()
parser.add_option('--smooth_window', type = str,default = '200',help = 'number of bins (default: %default)')
parser.add_option('--direction', type = str,default = 'r',help = 'direction (default: %default)') #r, z, xy
parser.add_option('--N_bin', type = int,default = 10000,help = 'number of bins to compute   (default: %default)')
(options, args) = parser.parse_args()
smooth_window = options.smooth_window
N_bin_old = options.N_bin
direction = options.direction

data_folder = '../example-APT-data' #APT data
with open (os.path.join(data_folder,"Nz-Nx-Ny.txt")) as f1:
	Nz_Nx_Ny = np.loadtxt(f1, delimiter=' ',usecols=(0),dtype=int,unpack=True)
Nz, Nx, Ny =  Nz_Nx_Ny[0], Nz_Nx_Ny[1], Nz_Nx_Ny[2]
			
N_xy = Nx*Ny

#pair = input("pair: ")
pair_list = ['all','SnSn']

#index_z_list = ['0']
index_z_list = np.arange(Nz)
index_xy_list = np.arange(N_xy)

#N_bin_old = 10000
#smooth_window = int(input("smooth window length: "))
smooth_half_window = int(int(smooth_window)/2)

for pair_name in pair_list:
	for index_z in index_z_list:
		data_dir = os.path.join(data_folder,'z-index-{0}'.format(index_z))
		for index_xy in index_xy_list:
			#input_name = os.path.join(data_dir,'rdf-{0}-cube'.format(pair)+'-index-{0}'.format(index_xy)+'-Nbin-{0}.txt'.format(N_bin_old))
			input_name = os.path.join(data_dir,'pair-correlation-{0}'.format(pair_name)+'-direction-{0}'.format(direction)+'-perturb-index-{0}.txt'.format(index_xy))
			output_name = os.path.join(data_dir,'pair-correlation-{0}'.format(pair_name)+'-direction-{0}'.format(direction)+'-perturb-index-{0}'.format(index_xy)+'-smooth-window-{0}.txt'.format(smooth_window))

			convert(input_name,output_name)
