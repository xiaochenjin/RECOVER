import numpy as np
#from collections import Counter

def read(file_name):
	#obtain number of type of species
#	with open (file_name) as f1:
#		original_species_list = np.loadtxt(f1, delimiter=' ',usecols=(0),dtype='str',unpack=True,skiprows=2)
#	species_list = list(set(original_species_list))
#	print(species_list)
#	N_species = len(species_list)

	species_list = []
	infile = open (file_name,"r")
	data = infile.readlines()
	infile.close()
	
#	infor_species_list = Counter(original_species_list)
#	print(infor_species_list)

	N_total = int(data[0].split()[0])
#	print(N_total)
#	a = [float(i) for i in data[1].split('"')[1].split()][:3]
#	b = [float(i) for i in data[1].split('"')[1].split()][3:6]
#	c = [float(i) for i in data[1].split('"')[1].split()][6:]
	a = [float(i) for i in data[1].split()[:3]]
	b = [float(i) for i in data[1].split()[3:6]]
	c = [float(i) for i in data[1].split()[6:]]

	cell_geometry = [a,b,c]
	H_T = np.array([a,b,c])
	H = H_T.transpose() #CELL MATRIX

	for line in data[2:]:
		atom_species = line.split()[0]
		if not atom_species in species_list:
			species_list.append(atom_species)
#	print(species_list)
	N_species = len(species_list)

#	print(cell_geometry)
#	print(H)
	frac_position_totlist = [[0 for x in range(1)] for y in range(N_species)]
	atom_index = 0
	for line in data[2:]:
		atom_species = line.split()[0]
		abs_position = [float(i) for i in line.split()[1:4]]
#		print(atom_species,abs_position)
		frac_position = list(np.matmul(np.linalg.inv(H),abs_position))
#		atom_infor = [atom_index,atom_species,frac_position]

		#use old index
		old_index = int(line.split()[-1])
		atom_infor = [old_index,atom_species,frac_position]

#		print(frac_position)
		for index_species in range(N_species):
			if atom_species == species_list[index_species]: 
#				frac_position_totlist[index_species].append(frac_position)
#				frac_position_totlist[index_species].append(abs_position)
				frac_position_totlist[index_species].append(atom_infor)
#		atom_index += 1

	for index_species in range(len(frac_position_totlist)):
		frac_position_list = frac_position_totlist[index_species]
		frac_position_list.pop(0)

	return N_total,species_list,frac_position_totlist,cell_geometry


#file_name = 'test.xyz'
#N_total,species_list,frac_position_totlist,cell_geometry=read(file_name)
#for index_species in range(len(frac_position_totlist)):
#	frac_position_list = frac_position_totlist[index_species]
#	for frac_position in frac_position_list:
#		print(frac_position)
