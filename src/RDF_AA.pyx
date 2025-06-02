cimport numpy as np
import numpy as np
from libc.math cimport sqrt

cpdef tuple rdf(list atom_infor_list,list cell_geometry,int N_bin):
	cdef int N_total 
	N_total = len(atom_infor_list)
	cdef list a,b,c,GR_list=[],R_list=[],M,atom_i,atom_j
	a = cell_geometry[0];b = cell_geometry[1];c = cell_geometry[2]

	cdef float V,density_0,dr,rij
	cdef np.ndarray H, H_T,s_vector,r_vector

	V = np.dot(a,np.cross(b,c))
	H_T = np.array([a,b,c])
	H = H_T.transpose() #CELL MATRIX
	density_0=N_total/V
	dr = 0.5*np.linalg.norm(a)/N_bin

	for n in range(N_bin):
		R = dr*(n+1)
		R_list.append(R)
	
	M=[0]*N_bin #list of m		
	for i in range(N_total):
		atom_i = atom_infor_list[i][2]	
		for j in range(i+1,N_total):
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
			k=int((rij-dr)/dr)
			if (k<N_bin):
				M[k]=M[k]+2

	for l in range(N_bin):
		density_R=(M[l]/N_total)/(4*3.14*R_list[l]*R_list[l]*dr)
		gr=density_R/density_0
		GR_list.append(gr)

	return R_list,GR_list

