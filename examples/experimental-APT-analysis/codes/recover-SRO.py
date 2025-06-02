import numpy as np
import os
from scipy import optimize
from scipy.signal import savgol_filter
from optparse import OptionParser
from multiprocessing import Pool
from functools import partial
from RECOVER import APT_RECOVER

parser = OptionParser()
parser.add_option('--max_iteration_steps', type = int,default = '10000',help = 'maximum iteration steps (default: %default)')
parser.add_option('--R_lower', type = float,default = 0.5,help = 'lower bound of R (default: %default)')
parser.add_option('--R_upper', type = float,default = 4.0,help = 'upper bound of R (default: %default)')
parser.add_option('--tol', type = str,default = '1e-9',help = 'tolerance (default: %default)')
parser.add_option('--p_min', type = float,default = -3,help = 'lower bound for fitted SRO parameters (default: %default)')
parser.add_option('--guess_index', type = str,default = 1,help = 'initial guess (default: %default)')
parser.add_option('--KNN', type = int,default = 1,help = 'shell to obtain SRO (default: %default)')
#parser.add_option('--sigma', type = str,default = 2.5,help = 'window length for smoothing (default: %default)')
parser.add_option('--smooth_window', type = str,default = '200',help = 'number of bins (default: %default)')
parser.add_option('--K_raw', type = int,default = 50,help = 'number of shells for fkr (default: %default)')
#parser.add_option('--R_min',type=float,default=0,help = 'lower bound of R (default: %default)')
parser.add_option('--N_thread', type = int,default = 1,help = 'number of threads (default: %default)')
(options, args) = parser.parse_args()
max_iteration_steps = options.max_iteration_steps
R_lower = options.R_lower
R_upper = options.R_upper
p_min = options.p_min
guess_index = options.guess_index
tol = options.tol
KNN = options.KNN
smooth_window = options.smooth_window
K_raw = options.K_raw
N_thread = options.N_thread

#Inputs 
mu = 0 #mean of atomic perturbation
direction = 'r'
pair_list = ['all','SnSn'] #pairs to obtain 
sigma_list = ['2.5','2.7','3.0','3.3','3.5','3.7']  #ranges of perturbation lengths
atom_species_list = ['Si','Ge','Sn'] #components
lc_list = [5.43,5.6597,6.4892] #corresponding experimental lattice constants


#specify file names and directories for crystal structures and fkr files
infor_file = 'structure-info-file.txt' #for crystal structures
infor_folder = '../../numerical-simulation/benchmark-structures'#for crystal structures
fkr_folder = '../fkr-files' #example fkr files
data_folder = '../example-APT-data' #APT data
analysis_folder = '../analysis'

if not os.path.exists(analysis_folder):
	os.mkdir(analysis_folder)

K_max_list = np.arange(10,K_raw+1,5)  #maximum KNN considered
K_target = int(KNN - 1) #the K  you want to consider

R_range = [R_lower,R_upper]
#part_list = ['lower-part','upper-part']

with open (os.path.join(data_folder,"Nz-Nx-Ny.txt")) as f1:
	Nz_Nx_Ny = np.loadtxt(f1, delimiter=' ',usecols=(0),dtype=int,unpack=True)
Nz, Nx, Ny =  Nz_Nx_Ny[0], Nz_Nx_Ny[1], Nz_Nx_Ny[2]
N_xy = Nx*Ny

index_z_list = [str(x) for x in np.arange(Nz)]
index_xy_list = [str(x) for x in np.arange(N_xy)]

#index_z_list = [0]
#index_xy_list = [0]

index_totlist = []
for index_z in index_z_list:
	for index_xy in index_xy_list:
		index = [index_z,index_xy]
		index_totlist.append(index)


for sigma in sigma_list:
	fkr_file = 'fkr-sigma-{0}.txt'.format(sigma)	

	with Pool(N_thread) as pool:
		fitting_2 = partial(APT_RECOVER.fitting,
		pair_list = pair_list,
		K_target = K_target,
		K_raw = K_raw,
		K_max_list = K_max_list,
		atom_species_list = atom_species_list,
		lc_list = lc_list,
		guess_index = guess_index,
		p_min = p_min,
		max_iteration_steps = max_iteration_steps,
		R_lower = R_lower,
		R_upper = R_upper,
		tol = tol,
		smooth_window = smooth_window,
		mu = mu,
		sigma = sigma,
		direction = direction,
		data_folder = data_folder,
		fkr_folder = fkr_folder,
		fkr_file = fkr_file,
		infor_folder = infor_folder,
		infor_file = infor_file,
		analysis_folder = analysis_folder)

		fitting = pool.map(fitting_2, index_totlist)	



