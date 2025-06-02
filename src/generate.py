import random
import os
import numpy as np

def generate(cell_matrix,atom_position_onecell,num_cell,lc,species_list,conc_list):
	a,b,c = cell_matrix[0],cell_matrix[1],cell_matrix[2]
	N_onecell = len(atom_position_onecell)
	atom_position_onecell_T = np.transpose(atom_position_onecell)
	u,v,w = atom_position_onecell_T[0],atom_position_onecell_T[1],atom_position_onecell_T[2]

	#CONSTRUCT FRACTIONAL COORDINATES ALONG EACH DIRECTION OF BASIS VECTOR
	atom_positionlist=[]
	for i in range(num_cell[0]):
		for j in range(num_cell[1]):
			for k in range(num_cell[2]):
				for n in range(N_onecell):
					atom_position=[]
					atom_position.append((u[n]+i)/num_cell[0]) #NORMALIZE THE COORDINATES BY NUMBER OF CELL IN EACH DIRECTION OF BASIS VECTOR
					atom_position.append((v[n]+j)/num_cell[1])
					atom_position.append((w[n]+k)/num_cell[2])
					atom_positionlist.append(atom_position)

	random.shuffle(atom_positionlist)
	
	a_vector = [a[0]*lc*num_cell[0],a[1]*lc*num_cell[0],a[2]*lc*num_cell[0]]
	b_vector = [b[0]*lc*num_cell[1],b[1]*lc*num_cell[1],b[2]*lc*num_cell[1]]
	c_vector = [c[0]*lc*num_cell[2],c[1]*lc*num_cell[2],c[2]*lc*num_cell[2]]
	cell_geometry = [a_vector,b_vector,c_vector]

	#convert to absolute coordinate
	H_T = np.array([a_vector,b_vector,c_vector])
	H = H_T.transpose() #CELL MATRIX

	N_total = len(atom_positionlist)
	N_species_list = [int(N_total*conc) for conc in conc_list]
	atom_species_list = []
	for index_species in range(len(species_list)):
		species = species_list[index_species]
		N_species = N_species_list[index_species]
		atom_species_list += [species for atom in range(N_species)]

	return atom_positionlist,atom_species_list,cell_geometry
