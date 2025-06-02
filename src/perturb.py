import numpy as np
import random
from operator import itemgetter

def iso_perturb(old_infor_totlist,cell_geometry,mu, sigma,collect_eff):	
	r_totlist = []
	new_infor_totlist = []
	for old_infor_list in old_infor_totlist:
		#print(old_infor_list)
		N_total = len(old_infor_list)
		N_sample = int(N_total*collect_eff)
		idx_list = random.sample(range(len(old_infor_list)), N_sample)	
		old_collect_list = [old_infor_list[idx] for idx in idx_list]
		r_list = np.random.normal(float(mu),float(sigma),N_sample) #random displacement generated from a Gaussian distribution
		a = cell_geometry[0];b = cell_geometry[1];c = cell_geometry[2]
		H_T = np.array([a,b,c])
		H = H_T.transpose() #CELL MATRIX

		new_infor_list = []
		for i in range(N_sample):
			r = r_list[i]
			ceta = random.uniform(0,2*np.pi)
			phi = random.uniform(0,np.pi)

			dx = r*np.cos(ceta)*np.sin(phi)
			dy = r*np.sin(ceta)*np.sin(phi)
			dz = r*np.cos(phi)
	#		print(r,np.sqrt(dx**2+dy**2+dz**2))

			#convert fraction coordinate to Cartesian coordinate
			old_infor = old_collect_list[i]
	#		print(old_infor)
			old_index = old_infor[0]
			old_species = old_infor[1]
			old_frac_position = old_infor[2]
	#		print("Old: ", old_frac_position)
			old_abs_position = np.matmul(H,old_frac_position)
			
			
			new_abs_position_x = old_abs_position[0] + dx
			new_abs_position_y = old_abs_position[1] + dy
			new_abs_position_z = old_abs_position[2] + dz

			new_abs_position = [new_abs_position_x,new_abs_position_y,new_abs_position_z]
			new_frac_position = list(np.matmul(np.linalg.inv(H),new_abs_position))

			#consider PBC:
			new_frac_position_x = new_frac_position[0] - int(new_frac_position[0])
			new_frac_position_y = new_frac_position[1] - int(new_frac_position[1])
			new_frac_position_z = new_frac_position[2] - int(new_frac_position[2])
			
			new_frac_position = [new_frac_position_x,new_frac_position_y,new_frac_position_z]
			new_infor = [old_index,old_species,new_frac_position]
	#		print(new_infor)
			new_infor_list.append(new_infor)

#		print(r_list)
#		print(r_totlist)
	
		r_totlist += list(r_list)
		new_infor_totlist += new_infor_list

	return r_totlist,new_infor_totlist

def aniso_perturb(old_infor_totlist,cell_geometry,mu, sigma_z,sigma_xy,collect_eff):
	r_z_totlist = []
	r_xy_totlist = []
	new_infor_totlist = []
	for old_infor_list in old_infor_totlist:
		N_total = len(old_infor_list)
		N_sample = int(N_total*collect_eff)
		idx_list = random.sample(range(len(old_infor_list)), N_sample)
		old_collect_list = [old_infor_list[idx] for idx in idx_list]
		r_z_list = np.random.normal(float(mu),float(sigma_z),N_sample) #random displacement generated from a Gaussian distribution
		r_xy_list = np.random.normal(float(mu),float(sigma_xy),N_sample) #random displacement in xy direction generated from a Gaussian distribution
		a = cell_geometry[0];b = cell_geometry[1];c = cell_geometry[2]
		H_T = np.array([a,b,c])
		H = H_T.transpose() #CELL MATRIX

		new_infor_list = []
		for i in range(N_sample):
	#		r = r_list[i]
			r_z = r_z_list[i]
			r_xy = r_xy_list[i]
			ceta = random.uniform(0,2*np.pi)
	#		phi = random.uniform(0,np.pi)

			dz = r_z
			dx = r_xy*np.sin(ceta)
			dy = r_xy*np.cos(ceta)

	#		dx = r*np.cos(ceta)*np.sin(phi)
	#		dy = r*np.sin(ceta)*np.sin(phi)
	#		dz = r*np.cos(phi)
	#	   print(r,np.sqrt(dx**2+dy**2+dz**2))

			#convert fraction coordinate to Cartesian coordinate
			old_infor = old_collect_list[i]
			old_index = old_infor[0]
			old_species = old_infor[1]
			old_frac_position = old_infor[2]
	#	   print("Old: ", old_frac_position)
			old_abs_position = np.matmul(H,old_frac_position)


			new_abs_position_x = old_abs_position[0] + dx
			new_abs_position_y = old_abs_position[1] + dy
			new_abs_position_z = old_abs_position[2] + dz

			new_abs_position = [new_abs_position_x,new_abs_position_y,new_abs_position_z]
			new_frac_position = list(np.matmul(np.linalg.inv(H),new_abs_position))

			#consider PBC:
			new_frac_position_x = new_frac_position[0] - int(new_frac_position[0])
			new_frac_position_y = new_frac_position[1] - int(new_frac_position[1])
			new_frac_position_z = new_frac_position[2] - int(new_frac_position[2])

			new_frac_position = [new_frac_position_x,new_frac_position_y,new_frac_position_z]
			new_infor = [old_index,old_species,new_frac_position]
	#	   print(new_infor)
			new_infor_list.append(new_infor)

		r_xy_totlist += list(r_xy_totlist)
		r_z_totlist += list(r_z_totlist)
		new_infor_totlist += new_infor_list

	#for infor in new_infor_list:
		#print(infor)

	return r_z_totlist,r_xy_totlist,new_infor_totlist
