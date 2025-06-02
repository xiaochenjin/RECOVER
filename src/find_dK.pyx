cimport numpy as np
import numpy as np
from libc.math cimport sqrt

#cpdef list compute(list atom_infor_list,list cell_geometry,int K_raw,list neigh_dis,float lc,str direction):
def compute(atom_infor_list,cell_geometry,K_raw,neigh_dis,lc,direction):
	cdef int N_total 
	N_total = len(atom_infor_list)
	cdef list a,b,c,atom_i,atom_j
	a = cell_geometry[0];b = cell_geometry[1];c = cell_geometry[2]

	cdef float V,density_0,dr,rij
	cdef np.ndarray H, H_T,s_vector,r_vector

	V = np.dot(a,np.cross(b,c))
	H_T = np.array([a,b,c])
	H = H_T.transpose() #CELL MATRIX
	density_0=N_total/V

	cdef list dK_totlist=[]
	for n_shell in range(K_raw):
		dK_totlist.append([])

#	for i in range(N_total):
#		atom_i = atom_infor_list[i][2]	
	atom_i = atom_infor_list[0][2]

	for j in range(1,N_total):
#		for j in range(i+1,N_total):
		atom_j = atom_infor_list[j][2]
		s_xij = atom_j[0] - atom_i[0]
		s_yij = atom_j[1] - atom_i[1]	
		s_zij = atom_j[2] - atom_i[2]
		s_xij=s_xij-round(s_xij)
		s_yij=s_yij-round(s_yij)
		s_zij=s_zij-round(s_zij)
		s_vector=np.array([s_xij,s_yij,s_zij])
		s_vector=np.array([s_xij,s_yij,s_zij])
		r_vector=np.matmul(H,s_vector) #CONVERT VECTOR IN CELL VECTOR DIRECTION TO CARTESIAN
		rij = sqrt(r_vector[0]*r_vector[0]+r_vector[1]*r_vector[1]+r_vector[2]*r_vector[2])

		if direction == 'x': dK = abs(r_vector[0])
		if direction == 'y': dK = abs(r_vector[1])
		if direction == 'z': dK = abs(r_vector[2])
		if direction == 'xy': dK = sqrt(r_vector[0]*r_vector[0]+r_vector[1]*r_vector[1])
		if direction == 'r':dK = sqrt(r_vector[0]*r_vector[0]+r_vector[1]*r_vector[1]+r_vector[2]*r_vector[2])

		n_shell = min(range(K_raw), key=lambda index: abs(rij-lc*neigh_dis[index]))

#		for n_shell in range(K_raw):
#			if abs(rij-lc*neigh_dis[n_shell])<0.25:
		dK_list = dK_totlist[n_shell]
		diff_list = []
		if len(dK_list) > 0:
			for dK_2 in dK_list:
				diff = abs(dK_2 - dK)
				diff_list.append(diff)
			if min(diff_list) < 0.2:
				continue
			dK_totlist[n_shell].append(dK)
		if len(dK_list) == 0:
			dK_totlist[n_shell].append(dK)
#				break

	return dK_totlist

