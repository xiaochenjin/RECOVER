import numpy as np

def write_perturb(infor_list,cell_geometry,file_name):
	a = cell_geometry[0]
	b = cell_geometry[1]
	c = cell_geometry[2]
	N_total = len(infor_list)

	H_T = np.array([a,b,c])	
	H = H_T.transpose() #CELL MATRIX

	with open (file_name,"w") as f1:
		f1.write(str(N_total)+'\n')
		f1.write(str(a[0])+' '+str(a[1])+' '+str(a[2])+' '+str(b[0])+' '+str(b[1])+' '+str(b[2])+' '+str(c[0])+' '+str(c[1])+' '+str(c[2])+'\n')
		for i in range(N_total):
			infor = infor_list[i]
#			print(infor)
			old_index = infor[0]
			species = infor[1]
			frac_position = infor[2]
			atom_position = np.matmul(H,frac_position)
			f1.write(str(species)+' '+str(atom_position[0])+' '+str(atom_position[1])+' '+str(atom_position[2])+' '+'#'+' '+str(old_index)+'\n')
	f1.close()

def write_benchmark(atom_position_list,atom_species_list,cell_geometry,file_name):
	a = cell_geometry[0]
	b = cell_geometry[1]
	c = cell_geometry[2]
	N_total = len(atom_position_list)

	H_T = np.array([a,b,c])
	H = H_T.transpose() #CELL MATRIX

	with open (file_name,"w") as f1:
		f1.write(str(N_total)+'\n')
		f1.write(str(a[0])+' '+str(a[1])+' '+str(a[2])+' '+str(b[0])+' '+str(b[1])+' '+str(b[2])+' '+str(c[0])+' '+str(c[1])+' '+str(c[2])+'\n')
		for i in range(N_total):
			frac_position = atom_position_list[i]
			species = atom_species_list[i]
			atom_position = np.matmul(H,frac_position)
			f1.write(str(species)+' '+str(atom_position[0])+' '+str(atom_position[1])+' '+str(atom_position[2])+'\n')
	f1.close()
	
