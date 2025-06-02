cimport numpy as np
import numpy as np
from libc.math cimport sqrt

cpdef tuple compute(list atom_infor_list,list cell_geometry,list true_kNN_infor_list,int K_raw,int N_bin,list neigh_dis,float lc,str direction):
	cdef int N_total 
	N_total = len(atom_infor_list)
	cdef list a,b,c,GR_KNN_list=[],R_list=[],atom_i,atom_j
	a = cell_geometry[0];b = cell_geometry[1];c = cell_geometry[2]

	cdef float V,density_0,dr,rij
	cdef np.ndarray H, H_T,s_vector,r_vector

	V = np.dot(a,np.cross(b,c))
	H_T = np.array([a,b,c])
	H = H_T.transpose() #CELL MATRIX
	density_0=N_total/V

	cdef float lx,ly,lz
	lx = np.linalg.norm(a)
	ly = np.linalg.norm(b)
	lz = np.linalg.norm(c)
	if direction == 'x': dr = (0.5*lx)/N_bin
	if direction == 'y': dr = (0.5*ly)/N_bin
	if direction == 'z': dr = (0.5*lz)/N_bin
	if direction == 'xy': dr = (0.5*min(lx,ly))/N_bin
	if direction == 'r':dr = 0.5*np.linalg.norm(a)/N_bin

	for n in range(N_bin):
		R = dr*(n+1)
		R_list.append(R)

	cdef list old_index_list = []
	cdef list true_kNN_index_list = [] #true kNN index for atom labeled by old index
	for i in range(len(true_kNN_infor_list)):
		true_kNN_infor = true_kNN_infor_list[i]
		old_index = true_kNN_infor[0]
		true_kNN_index = true_kNN_infor[1]
		old_index_list.append(old_index)
		true_kNN_index_list.append(true_kNN_index)
#	print(old_index_list)

	M=[0]*N_bin #list of m
	M_truekNN_totlist = []
	for n_shell in range(K_raw):
		M_truekNN_list = [0]*N_bin #track where true kNN moves
		M_truekNN_totlist.append(M_truekNN_list)

	for i in range(N_total):
		atom_i = atom_infor_list[i][2]	
		index_i = atom_infor_list[i][0]
#		print(index_i)
		true_kNN_index = true_kNN_index_list[old_index_list.index(index_i)]

		for j in range(i+1,N_total):
			atom_j = atom_infor_list[j][2]
			index_j = atom_infor_list[j][0]
			s_xij = atom_j[0] - atom_i[0]
			s_yij = atom_j[1] - atom_i[1]	
			s_zij = atom_j[2] - atom_i[2]
			s_xij=s_xij-round(s_xij)
			s_yij=s_yij-round(s_yij)
			s_zij=s_zij-round(s_zij)
			s_vector=np.array([s_xij,s_yij,s_zij])
			s_vector=np.array([s_xij,s_yij,s_zij])
			r_vector=np.matmul(H,s_vector) #CONVERT VECTOR IN CELL VECTOR DIRECTION TO CARTESIAN
			if direction == 'x': rij = abs(r_vector[0])
			if direction == 'y': rij = abs(r_vector[1])
			if direction == 'z': rij = abs(r_vector[2])
			if direction == 'xy': rij = sqrt(r_vector[0]*r_vector[0]+r_vector[1]*r_vector[1])
			if direction == 'r':rij = sqrt(r_vector[0]*r_vector[0]+r_vector[1]*r_vector[1]+r_vector[2]*r_vector[2])
			k=int((rij-dr)/dr)
			if (k<N_bin):
				for n_shell in range(K_raw):
					if index_j in true_kNN_index[n_shell]: M_truekNN_totlist[n_shell][k] += 2

	cdef list GR_truekNN_list,GR_truekNN_totlist = []
	for l in range(N_bin):
		if direction == 'x': dV = 2*ly*lz*dr
		if direction == 'y': dV = 2*lx*lz*dr
		if direction == 'z': dV = 2*lx*ly*dr
		if direction == 'xy':
			dV = 2*3.14*R_list[l]*lz*dr #same as 3.14*((R_list[l]+dr)**2-R_list[l]**2)*lz if dr is small
		if direction == 'r':
			dV = 4*3.14*R_list[l]*R_list[l]*dr

		GR_truekNN_list = []
		for n_shell in range(K_raw):
			M_truekNN_list = M_truekNN_totlist[n_shell]
			density_R=(M_truekNN_list[l]/N_total)/dV
			GR_truekNN = density_R/density_0
			GR_truekNN_list.append(GR_truekNN)

		GR_truekNN_totlist.append(GR_truekNN_list)
	return R_list,GR_truekNN_totlist

