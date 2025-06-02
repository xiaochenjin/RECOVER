cimport numpy as np
import numpy as np
from libc.math cimport sqrt

def compute(atom_infor_list1,atom_infor_list2,cell_geometry,K_raw,neigh_dis,lc):
	cdef int N1_total,N2_total
	N1_total = len(atom_infor_list1)
	N2_total = len(atom_infor_list2)

	cdef list a,b,c,atom_i,atom_j
	a = cell_geometry[0];b = cell_geometry[1];c = cell_geometry[2]

	cdef float V,dr,rij
	cdef np.ndarray H, H_T,s_vector,r_vector

	V = np.dot(a,np.cross(b,c))
	H_T = np.array([a,b,c])
	H = H_T.transpose() #CELL MATRIX

	cdef list kNN_index_list,kNN_index_totlist=[],NNK_index_totlist

	for i in range(N1_total):
		atom_i = atom_infor_list1[i][2]	
		index_i = atom_infor_list1[i][0]
		NNK_index_totlist = []
		for n_shell in range(K_raw):
			NNK_index_totlist.append([])

		for j in range(N2_total):
			atom_j = atom_infor_list2[j][2]
			index_j = atom_infor_list2[j][0]
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

#			n_shell = min(range(K_raw), key=lambda index: abs(rij-lc*neigh_dis[index]))
#			NNK_index_totlist[n_shell].append(index_j)

			for n_shell in range(K_raw):
				if abs(rij-lc*neigh_dis[n_shell])<0.25:
					NNK_index_totlist[n_shell].append(index_j)
					break

		#indices of kNN atoms for atom i
		kNN_index_list = [index_i,NNK_index_totlist]
		#print(kNN_index_list)
		kNN_index_totlist.append(kNN_index_list)

	return kNN_index_totlist

