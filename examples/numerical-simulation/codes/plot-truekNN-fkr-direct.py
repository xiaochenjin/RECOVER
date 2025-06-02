from __future__ import print_function
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib.legend_handler import HandlerTuple
from scipy.stats import norm
from scipy.optimize import curve_fit
from optparse import OptionParser
import pandas as pd

def obtain_fkr(R_list,GR_truekNN,Nk):
	density_0 = N/V
	Nk_sum = 0
	for idx in range(len(R_list)):
		R = R_list[idx]
		density_truekNN  = density_0*GR_truekNN[idx]
		Ni = (4*3.14*R**2)*density_truekNN*dr
		Nk_sum += Ni

	fkr_list = []
	for idx in range(len(R_list)):
		R = R_list[idx]
		density_truekNN  = density_0*GR_truekNN[idx]
#		fkr = (4*3.14*R**2)*density_truekNN/Nk
		fkr = (4*3.14*R**2)*density_truekNN/Nk_sum
		fkr_list.append(fkr)
	return fkr_list

#https://www.geeksforgeeks.org/python-gaussian-fit/

def gauss(x, A, ave, std):
	return A * np.exp(-(x - ave) ** 2 / (2 * std ** 2))

def normaldist(x,ave,std):
	return 1/(std * np.sqrt(2 * np.pi)) * np.exp(-(x - ave) ** 2 / (2 * std ** 2))

def folded_gauss(x,A,ave,std):
	f1 = A * np.exp(-(x - ave) ** 2 / (2 * std ** 2))
	f2 = A * np.exp(-(x + ave) ** 2 / (2 * std ** 2))
	f = f1+f2
	return f1+f2

def laplace(x,A,ave,b):
	f = 1/(2*b)*np.exp(-abs(x-ave)/b)
	return f

def folded_laplace(x,A,ave,b):
	f1 = A*np.exp(-abs(x-ave)/b)
	f2 = A*np.exp(-abs(x+ave)/b)
	f = f1+f2
	return f

def guess_param(function):
	neigh_dis=[np.sqrt(3/16),np.sqrt(1/2),np.sqrt(11/16),1,np.sqrt(19/16)]
	initial_guess_list = []
	if function == 'gauss' or function == 'folded_gauss':
		for n_shell in range(len(neigh_dis)):
			ave = lc*neigh_dis[n_shell]
			std = 0.5*np.sqrt(2*float(sigma)**2)
			A = 1/(std * np.sqrt(2 * np.pi))
			initial_guess = [A,ave,std]
			initial_guess_list.append(initial_guess)
	if function == 'laplace' or function == 'folded_laplace':
		for n_shell in range(len(neigh_dis)):
			ave = lc*neigh_dis[n_shell]
			std = 0.5*np.sqrt(2*float(sigma)**2)
			A = 1/(2*std)
			initial_guess = [A,ave,std]
			initial_guess_list.append(initial_guess)
	return initial_guess_list
	
	
def fit_to_function(R_ave_list,GR_kNN_list,function,initial_guess_list):
	#https://www.geeksforgeeks.org/python-gaussian-fit/
	#fit to Gaussian
	#remove nan values
	parameters_kNN_list = []
	for n_shell in range(len(GR_kNN_list)):
		GR_kNN = GR_kNN_list[n_shell]
		R_ave_list_fit = []
		GR_kNN_fit = []
		initial_guess = initial_guess_list[n_shell]
		for idx in range(len(R_ave_list)):
			R = R_ave_list[idx]
			GR = GR_kNN[idx]
#		   if not np.isnan(GR):
			R_ave_list_fit.append(R)
			GR_kNN_fit.append(GR)
#		print(np.isnan(GR_kNN_fit).any())

		#https://stackoverflow.com/questions/22021037/convert-string-into-a-function-call/22021058#22021058
		func = globals()[function]
		parameters_kNN, covariance_kNN = curve_fit(func,R_ave_list_fit,GR_kNN_fit,absolute_sigma='False',p0=initial_guess)
		parameters_kNN_list.append(parameters_kNN)
		print("{0}NN parameters: ".format(str(n_shell+1)), parameters_kNN)
		print("initial guess: ", initial_guess)
#		print("Covariance: ",covariance_kNN)
	

#	with open (os.path.join(directory_name,"fitted-gaussian-parameter-normalized-rdf-trueNN.txt"),'w') as f1:
#		   for j in range(3):
#			   f1.write(str(parameters_1NN[j])+' '+str(parameters_2NN[j])+' '+str(parameters_3NN[j])+'\n')
#	f1.close()

	return parameters_kNN_list


def plot(R_ave_list,GR_kNN_list,style):
	if style == 'actual':
		title = 'RDF of true kNN '+directory_name
		y_label = 'g(r) of  true kNN'
	if style == 'normalized':
		title = 'PDF of true kNN '+directory_name
		y_label = 'probability density'

	plt.rcParams['font.family'] = 'DeJavu Serif'
	plt.rcParams['font.serif'] = ['Times New Roman']

	fontproperties = {'fontweight' : 'bold', 'fontsize' : 13}
	legend_properties = {'weight':'bold','size':13,'style': 'italic'}

	fig, axs = plt.subplots(1,1, figsize=(4.8,4))#,dpi=500)

	axs.tick_params(which='minor',direction="in")
	axs.tick_params(which='major',direction="in")

	#axs.xaxis.grid(False, which='minor')
	axs.xaxis.set_minor_locator(AutoMinorLocator(2))
	axs.yaxis.set_minor_locator(AutoMinorLocator(2))

#	axs.plot(R_list,GR_list,'k-')
#	axs.scatter(R_ave_list,GR_kNN_list[0],color='r',marker='x',s=5,label = '1NN',alpha = 0.5)

	color_list = ['r','b','k','g','gold']
	axs.set_xlabel("r (Å)",fontsize='15',weight = 'bold')
	axs.set_ylabel(y_label,fontsize='15',weight = 'bold')
	axs.set_xticklabels(axs.get_xticks(),fontproperties)
	axs.set_yticklabels(axs.get_yticks(),fontproperties)
	axs.yaxis.set_major_formatter(plt.FormatStrFormatter('%.1f')) #https://stackoverflow.com/questions/29188757/matplotlib-specify-format-of-floats-for-tick-labels
	axs.xaxis.set_major_formatter(plt.FormatStrFormatter('%.1f'))

#	axs.plot(R_list,GR_list,'k-')
#	axs.scatter(R_ave_list,GR_kNN_list[0],color='r',marker='x',s=5,label = '1NN',alpha = 0.5)

	color_list = ['r','b','k','g','gold']
	label_list = ['1NN','2NN','3NN','4NN','5NN']
#   function = 'laplace'
#	function = 'gauss' #with amplitute A
#	function = 'normaldist' #without amplitute A
#	function = 'folded_gauss'
	initial_guess_list = guess_param(function)
	parameters_kNN_list = fit_to_function(R_ave_list,GR_kNN_list,function,initial_guess_list)
	for n_shell in range(len(label_list)-2):
#		axs.plot(R_ave_list,GR_kNN_list[n_shell],color=color_list[n_shell],linestyle='--',marker='o',markersize=3,label=label_list[n_shell])
		axs.scatter(R_ave_list,GR_kNN_list[n_shell],color=color_list[n_shell],marker='o',s=1,label=label_list[n_shell]) #,alpha = 0.5)	
		parameters_kNN = parameters_kNN_list[n_shell]
#		parameters_kNN = initial_guess_list[n_shell]
		func = globals()[function]
		if function == 'gauss' or function == 'folded_gauss' or function == 'laplace' or function == 'folded_laplace':
			fitted_GR_kNN = func(R_ave_list,parameters_kNN[0],parameters_kNN[1],parameters_kNN[2])
		if function == 'normaldist':
			fitted_GR_kNN = func(R_ave_list,parameters_kNN[0],parameters_kNN[1])
						
#		axs.plot(R_ave_list,fitted_GR_kNN,linestyle='-',color=color_list[n_shell]) #,label=label_list[n_shell])

	axs.set_xlim(0,8)
	axs.set_ylim(top = 1)
#	plt.legend(prop=legend_properties, frameon=False,loc='best') #bbox_to_anchor=(0.9,-0.1))#loc='best')
#	plt.title(title,fontsize='13',weight = 'bold')
	plt.tight_layout()
	plt.show()

def replace_nan(old_list):
	new_list = []
	for proportion in old_list:
		if np.isnan(proportion):
			proportion = 0
		new_list.append(proportion)
	return new_list

def obtain_rdf_trueNN(sgm):
	if smooth_window == '0':
		file_name = os.path.join(directory_name,'fkr-SnSn-direct.txt')
	else:
		file_name = os.path.join(directory_name,'fkr-SnSn-direct-smooth-window-{0}.txt'.format(smooth_window))
	with open (file_name) as f1:
		df = pd.read_csv(file_name,delimiter=' ', header=None)
		R_list = list(df.iloc[:, 0])
		normalized_GR_truekNN_list = []	
		for n_shell in range(N_plot):
			normalized_GR_truekNN = list(df.iloc[:, n_shell+1])
			normalized_GR_truekNN_list.append(normalized_GR_truekNN)

	plot(R_list,normalized_GR_truekNN_list,'normalized')

#file_name1 = "rdf-SnSn-noperturb"
parser = OptionParser()
parser.add_option('--fit_function', type = str,default = 'folded_laplace',help = 'function for fitting (default: %default)')
(options, args) = parser.parse_args()
function =  options.fit_function
#	function = 'laplace'
#	function = 'folded_laplace'
#	function = 'gauss' #with amplitute A
#	function = 'normaldist' #without amplitute A
#	function = 'folded_gauss'

N = 64000*float(collect_eff)
c_Sn = 0.25
N_Sn = N*c_Sn
a = 117.2
N_bin = 10000
dr = 0.5*a/N_bin
V = a**3
lc = a/20

N_plot = 5
conc_sn = 0.25
collect_eff = input("collect efficiency: ")
mu = 0
sigma = input("sigma: ")
smooth_window = input("smooth window: ")
index_list = np.arange(100)
test_folder = '../../analysis/recover-test-4/test1-64000-atom'
os.chdir(test_folder)
os.chdir("efficiency-{0}".format(collect_eff))

infor_dir = '/home/xcjin/Research/other-small-things/generate-perturbed-cell/bulk-GeSn/64000-atom-cell/Sn-0.25/optimization-test/obtain-kNN-infor'
with open (os.path.join(infor_dir,"kNN-peak-position.txt")) as f1:
	Nk_list = np.loadtxt(f1, delimiter=' ',usecols=(2),unpack=True)

Nk_Sn_list = np.multiply(conc_sn,Nk_list[:N_plot])

with open (os.path.join(infor_dir,"pSRO-SnSn-noperturb.txt")) as f1:
	pSRO_list = np.loadtxt(f1, delimiter=' ',usecols=(2),unpack=True)

Nk_Sn_list_SRO = np.multiply(pSRO_list[:N_plot],Nk_Sn_list)



directory_name = "mu-"+str(mu)+"-sgm-"+str(sigma)
obtain_rdf_trueNN(sigma)



