
import os
import numpy as np
from operator import itemgetter

def write_data(folder):
	os.chdir(folder)
	print(folder)
	infile =  open("all-raw-data.txt",'r')
	data = infile.readlines()
	infile.close()

	print("Total number of atoms: ", len(data))

	atom_position_totlist = []
	for species in species_list:
		atom_position_list = []

		for line in data:
			atom_position_x = float(line.split()[0])
			atom_position_y = float(line.split()[1])
			atom_position_z = float(line.split()[2])
			atom_position = [atom_position_x,atom_position_y,atom_position_z]

			if "{0}:1".format(species) in line:
				atom_position_list.append(atom_position)

		atom_position_list = sorted(atom_position_list, key=itemgetter(2))
		atom_position_totlist.append(atom_position_list)				

	with open (output_name,'w') as f1:
		for index_species in range(len(species_list)):
			species = species_list[index_species]
			atom_position_list = atom_position_totlist[index_species]

			for i in range(len(atom_position_list)):
				atom_position = atom_position_list[i]
				f1.write('{0} '.format(species)+str(atom_position[0])+' '+str(atom_position[1])+' '+str(atom_position[2])+'\n')
	f1.close()
	
species_list = ['Ge','Sn']
output_name = "".join(species_list)+'-raw-data.txt'
write_data('../example-APT-data/')
