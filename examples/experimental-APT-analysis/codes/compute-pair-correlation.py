import os
import numpy as np
from RECOVER import cal_PairCorrelation_exp
from multiprocessing import Pool
from functools import partial
from optparse import OptionParser

parser = OptionParser()
parser.add_option('--direction', type = str,default = 'r',help = 'direction (default: %default)') #r, z, xy
parser.add_option('--N_bin', type = int,default = 10000,help = 'number of bins to compute   (default: %default)')
parser.add_option('--N_thread', type = int,default = 1,help = 'number of threads   (default: %default)')
(options, args) = parser.parse_args()

direction = options.direction
N_bin = options.N_bin
N_thread = options.N_thread

pair_list = ['all','Sn-Sn','Ge-Ge','Ge-Sn']
data_folder = '../example-APT-data/'

with open (os.path.join(data_folder,"Nz-Nx-Ny.txt")) as f1:
    Nz_Nx_Ny = np.loadtxt(f1, delimiter=' ',usecols=(0),dtype=int,unpack=True)
Nz, Nx, Ny =  Nz_Nx_Ny[0], Nz_Nx_Ny[1], Nz_Nx_Ny[2]
N_xy = Nx*Ny

index_z_list = np.arange(Nz)
index_xy_list = np.arange(N_xy)
index_list = []
for index_z in index_z_list:
    for index_xy in index_xy_list:
        index = [index_z,index_xy]
        index_list.append(index)

if __name__ == '__main__':
	with Pool(N_thread) as pool:
		compute_2 = partial(cal_PairCorrelation_exp.compute,data_folder=data_folder,pair_list=pair_list,N_bin=N_bin,direction=direction)
		pool.map(compute_2,index_list)





