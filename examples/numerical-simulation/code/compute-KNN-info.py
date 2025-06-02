import os
import numpy as np
import random
from RECOVER import generate
from RECOVER import cal_KNN_info
from RECOVER import read_xyz_fkr
from optparse import OptionParser

#GET ARGUMENTS
parser = OptionParser()
parser.add_option('--pair_type', type = str,default = 'Sn-Sn',help = 'pair type for computing fk(r) (default: %default)')
parser.add_option('--K_raw', type = int,default = 3,help = 'number of shells for computing fk(r) (default: %default)')
parser.add_option('--lattice_constant', type = float,default = 5.86,help = 'direction (default: %default)')
(options, args) = parser.parse_args()
K_raw = options.K_raw
lc = options.lattice_constant
pair_type = options.pair_type

folder_name = '../benchmark-structures/'
structure_infor_file = 'structure-info-file.txt'
input_name = 'GeSn-random-relaxed.xyz'
output_name = 'KNN-info-{0}.txt'.format(pair_type)


#read nearest neighbor locations:
with open (os.path.join(folder_name,'structure-info-file.txt')) as f1:
	neigh_dis = np.loadtxt(f1, delimiter=' ', usecols=(1), unpack=True)

neigh_dis = neigh_dis.tolist()

cal_KNN_info.compute(input_name,folder_name,output_name,pair_type,K_raw,neigh_dis,lc)





