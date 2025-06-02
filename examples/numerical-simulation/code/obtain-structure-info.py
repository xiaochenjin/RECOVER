import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib.legend_handler import HandlerTuple

def read(file_name):
	infile = open(os.path.join(folder,file_name),'r')
	data = infile.readlines()
	N_total = int(data[0].split()[0])
	a = [float(i) for i in data[1].split()[:3]]
	b = [float(i) for i in data[1].split()[3:6]]
	c = [float(i) for i in data[1].split()[6:]]
	lx = np.linalg.norm(a)
	ly = np.linalg.norm(b)
	lz = np.linalg.norm(c)
	V = np.dot(a,np.cross(b,c))
	return N_total,V,lx,ly,lz

def obtain_structure_info_perfect(file_name):
	with open (os.path.join(folder,file_name)) as f1:
		R_list,GR_list = np.loadtxt(f1, delimiter=' ', usecols=(0,1),unpack=True)

	R_KNN_list = []
	GR_KNN_list = []
	for i in range(len(R_list)):
		R = R_list[i]
		GR = GR_list[i]
		if GR > 0:
			R_KNN_list.append(R)
			GR_KNN_list.append(GR)

	dr = R_list[1] - R_list[0]
	density_0 = N_atom/V
	M_list = []
	for i in range(len(GR_KNN_list)):
		R = R_KNN_list[i]
		GR = GR_KNN_list[i]
		density_R = density_0*GR
		if direction == 'x': dV = 2*ly*lz*dr
		if direction == 'y': dV = 2*lx*lz*dr
		if direction == 'z': dV = 2*lx*ly*dr
		if direction == 'xy':
			dV = 3.14*((R+dr)**2-R**2)*lz
		if direction == 'r':
			dV = 4*3.14*R**2*dr
		M = density_R*dV
		M_list.append(M)

	return R_list,GR_list,R_KNN_list, M_list

def obtain_structure_info_relax(file_name):
	M_list = M_list_perfect

	with open (os.path.join(folder,file_name)) as f1:
		R_list,GR_list = np.loadtxt(f1, delimiter=' ', usecols=(0,1),unpack=True)

	R_KNN_list = []
	for n_shell in range(len(R_KNN_list_perfect)):
		if 0 < n_shell < len(R_KNN_list_perfect)-1:
			length_lower = 0.5*(R_KNN_list_perfect[n_shell]-R_KNN_list_perfect[n_shell-1])
			length_upper = 0.5*(R_KNN_list_perfect[n_shell+1]-R_KNN_list_perfect[n_shell])
			length = min(length_lower,length_upper)
		if n_shell == 0:
			length  = 0.5*(R_KNN_list_perfect[n_shell+1]-R_KNN_list_perfect[n_shell])
		if n_shell == len(R_KNN_list_perfect)-1:
			length = 0.5*(R_KNN_list_perfect[n_shell]-R_KNN_list_perfect[n_shell-1])
		R_upper = R_KNN_list_perfect[n_shell]+length
		R_lower = R_KNN_list_perfect[n_shell]-length
		R_range = [R_lower,R_upper]
		#print(R_range)
		#https://stackoverflow.com/questions/9706041/finding-index-of-an-item-closest-to-the-value-in-a-list-thats-not-entirely-sort
		index_lower = min(range(len(R_list)), key=lambda i: abs(R_list[i]-R_lower))
		index_higher = min(range(len(R_list)), key=lambda i: abs(R_list[i]-R_upper))
		index_KNN_range = [index_lower,index_higher]
		#print(index_KNN_range)
		GR_KNN_range = GR_list[index_lower:index_higher]
		index_KNN = list(GR_list).index(max(GR_KNN_range)) #index of the Sn-Sn peak position
		R_KNN = R_list[index_KNN] #choose peak position
		#print(R_KNN_SnSn)
		R_KNN_list.append(R_KNN)

	return R_list,GR_list,R_KNN_list,M_list

def write(R_KNN_list,M_list,output_name):
	with open ((os.path.join(folder,output_name)),"w") as f1:	
		for i in range(len(R_KNN_list)):
			frac_R = R_KNN_list[i]/lc
			f1.write('{0}NN'.format(i+1)+' '+str(frac_R)+' '+str(round(M_list[i]))+'\n')
	f1.close()

def plot(R_KNN_list, M_list):
	plt.rcParams['font.family'] = 'DeJavu Serif'
	plt.rcParams['font.serif'] = ['Times New Roman']

	fontproperties = {'fontweight' : 'bold', 'fontsize' : 10}
	legend_properties = {'weight':'bold','size':10,'style': 'italic'}

	fig, axs = plt.subplots(1,1, figsize=(6,5))#,dpi=500)

	axs.tick_params(which='minor',direction="in")
	axs.tick_params(which='major',direction="in")

	#axs.xaxis.grid(False, which='minor')
	axs.xaxis.set_minor_locator(AutoMinorLocator(2))
	axs.yaxis.set_minor_locator(AutoMinorLocator(2))

	axs.set_xlabel("r",fontsize='13',weight = 'bold')
	axs.set_ylabel("g(r)",fontsize='13',weight = 'bold')
	axs.set_xticklabels(axs.get_xticks(),fontproperties)
	axs.set_yticklabels(axs.get_yticks(),fontproperties)
	axs.yaxis.set_major_formatter(plt.FormatStrFormatter('%.2f')) #https://stackoverflow.com/questions/29188757/matplotlib-specify-format-of-floats-for-tick-labels
	axs.xaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))

	axs.plot(R_list,GR_list,'k--')
	axs.scatter(R_KNN_list,GR_KNN_list,color='r',marker='s')
	axs.set_xlim(0,3)
#   axs.set_ylim(0,5)
	#plt.legend(prop=legend_properties, frameon=False,bbox_to_anchor=(0.9,-0.1))#loc='best')
	plt.title(file_name,fontsize='13',weight = 'bold')
	plt.show()

direction = 'r'
folder = '../benchmark-structures/'

pair_type = input("pair type: ") #Sn-Sn
lc = float(input("lattice constant: "))
#whether_relax = input("whether relax: ")

structure_file = 'GeSn-random-before-relax.xyz'
input_name_perfect = 'No-relax-pair-correlation-{0}'.format(pair_type)+'-direction-{0}.txt'.format(direction)
input_name_relax = 'Relaxed-pair-correlation-{0}'.format(pair_type)+'-direction-{0}.txt'.format(direction)
#output_name_perfect = 'structure-info-file-no-relax-{0}.txt'.format(pair_type)
output_name_perfect = 'structure-info-file.txt'.format(pair_type)
output_name_relax = 'structure-info-file-relax-{0}.txt'.format(pair_type)

N_atom, V, lx, ly, lz = read(structure_file)

R_list_perfect,GR_list_perfect,R_KNN_list_perfect,M_list_perfect = obtain_structure_info_perfect(input_name_perfect)
write(R_KNN_list_perfect,M_list_perfect,output_name_perfect)

R_list_relax,GR_list_relax,R_KNN_list_relax,M_list_relax = obtain_structure_info_relax(input_name_relax)
write(R_KNN_list_relax,M_list_relax,output_name_relax)

#plot(R_list_relax,GR_list_relax)
	
