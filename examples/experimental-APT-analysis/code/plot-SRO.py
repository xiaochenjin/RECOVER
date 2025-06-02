import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

def plot():
	plt.rcParams['font.family'] = 'DeJavu Serif'
	plt.rcParams['font.serif'] = ['Times New Roman']

	fontproperties = {'fontweight' : 'bold', 'fontsize' : 14}
	legend_properties = {'weight':'bold','size':12,'style': 'italic'}

	fig, axs = plt.subplots(1,1, figsize=(6,5))

	axs.tick_params(which='minor',direction="in")
	axs.tick_params(which='major',direction="in")

	axs.xaxis.set_minor_locator(AutoMinorLocator(2))
	axs.yaxis.set_minor_locator(AutoMinorLocator(2))

	axs.set_xlabel("Perturbation length "+r'$\sigma$ (Å)',fontsize='15',weight = 'bold')
	axs.set_ylabel("SRO parameter",fontsize='15',weight = 'bold')
	axs.set_xticklabels(sigma_list)
	axs.set_xticklabels(axs.get_xticks(),fontproperties)
	axs.set_yticklabels(axs.get_yticks(),fontproperties)
	axs.yaxis.set_major_formatter(plt.FormatStrFormatter('%.1f')) #https://stackoverflow.com/questions/29188757/matplotlib-specify-format-of-floats-for-tick-labels
	axs.xaxis.set_major_formatter(plt.FormatStrFormatter('%.1f'))

#   X_axis = np.arange(1,len(sigma_list)+1)
	X_axis = [float(x) for x in sigma_list]
	axs.plot(X_axis,[0]*len(X_axis), 'k--',label='random')

	axs.errorbar(X_axis,ave_alpha_list,yerr= ste_alpha_list, marker='s',markersize=4,color='r',capsize=11,label='retreived SRO parameter')

#	axs.set_ylim(-0.3,0.3)

	plt.legend(prop=legend_properties, frameon=False,loc='best') #bbox_to_anchor=(0.8,-0.1))#loc='best')
	plt.tight_layout()
	plt.show()

def obtain_pSRO(sigma,max_GR_criteria):
	file_name = "{0}NN".format(KNN)+"-{0}-adjusted-SRO-parameter".format(pair)+'-sigma-{0}'.format(sigma)+'.txt'

	with open (os.path.join(analysis_folder,'screen-max-gr-{0}'.format(max_GR_criteria),file_name)) as f1:
		alpha_list = np.loadtxt(f1, delimiter=' ', usecols=(-1), unpack=True)	
	return alpha_list


pair = input('pair: ')
KNN = input("Kth shell: ")
sigma_list = ['2.5','2.7','3.0','3.3','3.5','3.7']
max_GR_criteria = input("max GR: ")
analysis_folder = '../analysis'
ave_alpha_list = [];ste_alpha_list = []
for sigma in sigma_list:
	alpha_list = obtain_pSRO(sigma,max_GR_criteria)
	ave_alpha = np.average(alpha_list)
	ste_alpha = np.std(alpha_list)/np.sqrt(len(alpha_list))
	ave_alpha_list.append(ave_alpha)
	ste_alpha_list.append(ste_alpha)

	print("sigma: ", sigma, round(ave_alpha,3),round(ste_alpha,3))

plot()
