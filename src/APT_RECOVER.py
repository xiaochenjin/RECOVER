import numpy as np
import os
from scipy import optimize
from scipy.signal import savgol_filter
#import plot_convergence
from optparse import OptionParser
from multiprocessing import Pool

def smooth(original_list):
	window_length = 5
	smoothed_list = savgol_filter(original_list,window_length,polyorder=1)
	return smoothed_list

#obtained by /home/xcjin/Research/other-small-things/generate-perturbed-cell/bulk-GeSn/64000-atom-cell/Sn-0.25/Random/scripts/recover-test-3
#plot-truekNN-fkr-direct.py

def grep_raw_fkr_index(R_list):
	with open (fkr_folder_file) as f1:
		R_raw_totlist = np.loadtxt(f1, delimiter=' ', usecols=(0), unpack=True)

	raw_fkr_index_list = []
	for idx in range(len(R_list)):
		index = min(range(len(R_raw_totlist)), key=lambda i: abs(R_raw_totlist[i]-R_list[idx]))
		raw_fkr_index_list.append(index)
	return raw_fkr_index_list

#obtained by /home/xcjin/Research/other-small-things/generate-perturbed-cell/bulk-GeSn/64000-atom-cell/Sn-0.25/Random/scripts/recover-test-3
#plot-truekNN-fkr-direct.py
def grep_raw_fkr(R_list,n_shell):
	infile = open(fkr_folder_file,'r')
	data = infile.readlines()
	infile.close()
	R_raw_totlist = []
	raw_fkr_totlist = []
	for line in data:
		R_raw = float(line.split()[0])
		raw_fkr = float(line.split()[n_shell+1])
		R_raw_totlist.append(R_raw)
		raw_fkr_totlist.append(raw_fkr)

	raw_fkr_list = []
	for idx in raw_fkr_index_list_SRO:
		raw_fkr = raw_fkr_totlist[idx]
		raw_fkr_list.append(raw_fkr)

	#smooth fkr:
	#print(len(raw_fkr_list))
	raw_fkr_list = smooth(raw_fkr_list)
	return raw_fkr_list

def grep_GR(data_name):
	#obtain by folder/plot-truekNN-proportion-ave.py
	with open (data_name) as f1:
		R_totlist,GR_totlist = np.loadtxt(f1, delimiter=' ', usecols=(0,1),unpack=True)

	#select g(r) in the R_range:
	#https://stackoverflow.com/questions/9706041/finding-index-of-an-item-closest-to-the-value-in-a-list-thats-not-entirely-sort
	R_list = []
	GR_list = []
	R_lower = R_range[0]; R_upper = R_range[1]
	index_lower = min(range(len(R_totlist)), key=lambda i: abs(R_totlist[i]-R_lower))
	index_upper = min(range(len(R_totlist)), key=lambda i: abs(R_totlist[i]-R_upper))
	for R in R_totlist[index_lower:index_upper]:
		R_list.append(R)
	for GR in GR_totlist[index_lower:index_upper]:
		GR_list.append(GR)
	return R_list,GR_list

def obtain_F_kr(R_list):
	with open (infor_folder_file) as f1:
		Nk_list = np.loadtxt(f1, delimiter=' ',usecols=(2),unpack=True)

	F_kr_list = [] #matrix
	for n_shell in range(K_max):
		Nk = Nk_list[n_shell]*collect_efficiency
		f_kr = grep_raw_fkr(R_list,n_shell)
		F_kr = np.multiply(Nk,f_kr)
		F_kr_list.append(F_kr)
	return F_kr_list

def obtain_b_vector(R_list,GR_list): #rescale g(r) to b
	M_list = []
	for i in range(len(R_list)):
		R = R_list[i]
		GR = GR_list[i]
		M = GR*density_0*4*3.14*R**2
		M_list.append(M)
	return M_list

def compute_diff_SRO(alpha_list):
	xk_list = np.subtract(1,alpha_list)
	estimated_b_vector_SRO = [0]*len(R_list_SRO)
	for n_shell in range(K_max):
		F_kr_SRO = F_kr_list_SRO[n_shell]
		xk = xk_list[n_shell]
		scaled_F_kr_SRO = np.multiply(xk,F_kr_SRO)
		estimated_b_vector_SRO = np.add(estimated_b_vector_SRO,scaled_F_kr_SRO)
	total_diff = 0
	for i in range(len(R_list_SRO)):
		R = R_list_SRO[i]
		constant = density_0*4*3.14*R**2
		smoothed_GR = smoothed_GR_list_SRO[i]
		estimated_GR = (estimated_b_vector_SRO[i])/constant
		diff  = estimated_GR-smoothed_GR
		total_diff += diff**2
	RMSE = np.sqrt(total_diff)/len(R_list_SRO)
	return RMSE

def write(K_max_list,alpha_KNN_list_SRO,file_name):
	with open (file_name,'w') as f1:
		for i in range(len(K_max_list)):
			f1.write(str(K_max_list[i])+' '+str(alpha_KNN_list_SRO[i])+'\n')
	f1.close()

def fitting(
	index_xyz,
	pair_list,
	K_target,
	K_raw,
	K_max_list,
	atom_species_list,
	lc_list,
	guess_index,
	p_min,
	max_iteration_steps,
	R_lower,
	R_upper,
	tol,
	smooth_window,
	mu,
	sigma,
	direction,
	data_folder,
	fkr_folder,
	fkr_file,
	infor_folder,
	infor_file,
	analysis_folder):

	global KNN
	KNN = int(K_target + 1) #the K  you want to consider

	global index_z,index_xy
	index_z = index_xyz[0]
	index_xy = index_xyz[1]

	#output folder
	output_folder = os.path.join(analysis_folder,'z-index-{0}'.format(index_z))
	if not os.path.exists(output_folder):
		os.mkdir(output_folder)

	global R_range
	R_range = [R_lower,R_upper]

	global fkr_folder_file
	fkr_folder_file = os.path.join(fkr_folder,fkr_file)

	global infor_folder_file
	infor_folder_file = os.path.join(infor_folder,infor_file)

	global data_dir
	data_dir = os.path.join(data_folder,'z-index-{0}'.format(index_z))
	file_name = os.path.join(data_dir,'cube-index-{0}.xyz'.format(index_xy))
	infile = open(file_name,'r')
	data = infile.readlines()
	infile.close()
	N = int(data[0].split()[0])

	a = float(data[1].split()[0])
	V = a**3
	global density_0 
	density_0 = N/V


	#lattice constant
	N_species = len(atom_species_list)
	N_atom_list = [0]*N_species
	for line in data[2:]:
		atom_species = line.split()[0]
		for species_index in range(N_species):
			if atom_species == atom_species_list[species_index]: N_atom_list[species_index] += 1
	conc_list = [N_atom/N for N_atom in N_atom_list]
	lc = np.sum(np.multiply(lc_list,conc_list))
	#print("concentrations: ", conc_list)
	#print("lattice constant: ",lc)

	#collection efficiency
	V_onecell = lc**3
	N_onecell = 8
	N_max_atoms = V/V_onecell*N_onecell
	global collect_efficiency
	collect_efficiency = N/N_max_atoms
#	  print(collect_efficiency)


	global pair
	for pair in pair_list:
		if smooth_window == '0':
			data_name = os.path.join(data_dir,'pair-correlation-{0}'.format(pair)+'-direction-{0}'.format(direction)+'-perturb-index-{0}.txt'.format(index_xy))
		else:
			data_name = os.path.join(data_dir,'pair-correlation-{0}'.format(pair)+'-direction-{0}'.format(direction)+'-perturb-index-{0}'.format(index_xy)+'-smooth-window-{0}.txt'.format(smooth_window))

		global R_list_SRO,smoothed_GR_list_SRO
		R_list_SRO,GR_list_SRO = grep_GR(data_name)
		smoothed_GR_list_SRO = smooth(GR_list_SRO)
		global raw_fkr_index_list_SRO
		raw_fkr_index_list_SRO =  grep_raw_fkr_index(R_list_SRO)

		alpha_KNN_list_SRO = []
		global K_max
		for K_max in K_max_list:

			#initial guess for SRO parameter:
			if guess_index == '1':initial_guess = [0]*K_max
			if guess_index == '2':initial_guess = [-0.5]*K_max
			if guess_index == '3':initial_guess = [-1]*K_max

		#	   print("{0}NN: ".format(str(K_max)))
			#construct F (r*k)
			global F_kr_list_SRO
			F_kr_list_SRO = obtain_F_kr(R_list_SRO)

			initial_guess = [0]*K_max
			bnds = optimize.Bounds(lb=p_min,ub=1)
			results_SRO = optimize.minimize(compute_diff_SRO,x0=initial_guess,bounds=bnds,method='SLSQP',tol=float(tol),options={'maxiter':int(max_iteration_steps)})
			if results_SRO.success:
				alpha_xk_SRO = results_SRO.x
				alpha_KNN_list_SRO.append(alpha_xk_SRO[K_target]) #K_target = 1NN

			else:
				raise ValueError(results_SRO.message)

		#print("SRO alpha 1NN: ",alpha_KNN_list_SRO[-1])


		file_name = os.path.join(output_folder,'index-{0}'.format(index_xy)+'-{0}NN-'.format(KNN)+'{0}-SRO-convergence'.format(pair)+'-sigma-{0}'.format(str(sigma))+'.txt')
		write(K_max_list,alpha_KNN_list_SRO,file_name)


