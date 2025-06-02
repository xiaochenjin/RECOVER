import numpy as np
import pandas as pd
import os
from scipy import optimize
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
#from scipy.signal import savgol_filter

def grep_raw_fkr(folder):
	file_name =  os.path.join(folder,fkr_file)
	with open (file_name) as f1:
		df = pd.read_csv(file_name,delimiter=' ', header=None)
	R_list = list(df.iloc[:, 0])
	fkr_totlist = []
	for n_shell in range(K_plot):
		fkr_list = list(df.iloc[:, n_shell+1])
		fkr_totlist.append(fkr_list)
	return R_list,fkr_totlist


#		R_raw_list,raw_fkr_1NN_list,raw_fkr_2NN_list,raw_fkr_3NN_list,raw_fkr_4NN_list,raw_fkr_5NN_list = np.loadtxt(f1, delimiter=' ', usecols=(0,1,2,3,4,5),unpack=True)
#	raw_fkr_list = [raw_fkr_1NN_list,raw_fkr_2NN_list,raw_fkr_3NN_list,raw_fkr_4NN_list,raw_fkr_5NN_list]
	f1.close()
	return R_raw_list,raw_fkr_list

def plot():
	plt.rcParams['font.family'] = 'DeJavu Serif'
	plt.rcParams['font.serif'] = ['Times New Roman']

	fontproperties = {'fontweight' : 'bold', 'fontsize' : 10}
	legend_properties = {'weight':'bold','size':11,'style': 'italic'}

#	fig, axs = plt.subplots(1,1, figsize=(6,5))#,dpi=500)
	fig, axs = plt.subplots(1,1, figsize=(4.8,4))#,dpi=500)

	axs.tick_params(which='minor',direction="in")
	axs.tick_params(which='major',direction="in")

	#axs.xaxis.grid(False, which='minor')
	axs.xaxis.set_minor_locator(AutoMinorLocator(2))
	axs.yaxis.set_minor_locator(AutoMinorLocator(2))


	axs.set_ylabel("probability density",fontsize='13',weight = 'bold')
	axs.set_xlabel("r (Å)",fontsize='13',weight = 'bold')
	axs.set_xticklabels(axs.get_xticks(),fontproperties)
	axs.set_yticklabels(axs.get_yticks(),fontproperties)
	axs.yaxis.set_major_formatter(plt.FormatStrFormatter('%.2f')) #https://stackoverflow.com/questions/29188757/matplotlib-specify-format-of-floats-for-tick-labels
	axs.xaxis.set_major_formatter(plt.FormatStrFormatter('%.0f'))
#	axs.set_ylim(-2,1.5)
	axs.set_xlim(0,8)

	color_list = ['r','b','k','g','gold']
	label_list = ['1NN','2NN','3NN','4NN','5NN']
	for n_shell in range(K_plot):
		axs.scatter(R_list,raw_fkr_totlist[n_shell],color = color_list[n_shell],s=1,label=label_list[n_shell])
	axs.set_xlim(0,8)
	axs.set_ylim(top = 1.0)
	plt.legend(prop=legend_properties, frameon=False,loc='best') #bbox_to_anchor=(0.9,-0.1))#loc='best')
	plt.tight_layout()
	plt.show()

mu = 0
K_plot = 5 #plot fkr for first K_plot shell

pair_type = input("pair type: ")  
collect_eff = input("collection efficiency: ")
perturb_type = input("perturb type: ")

fkr_file = 'fkr-file-{0}.txt'.format(pair_type)

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

folder = os.path.join(analysis,directory_name)
R_list,raw_fkr_totlist = grep_raw_fkr(folder)
plot()



