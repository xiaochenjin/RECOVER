from RECOVER import generate
import os
import numpy as np
from RECOVER import write_xyz

def read_unitcell(file_name):
	with open (file_name) as f1:
		atom_position_x_list,atom_position_y_list,atom_position_z_list = np.loadtxt(f1,usecols=(0,1,2),unpack=True,skiprows=8)
	atom_position_list = list(zip(atom_position_x_list,atom_position_y_list,atom_position_z_list))
	return atom_position_list

#define basis of cell shape
a = [1,0,0]
b = [0,1,0]
c = [0,0,1]
cell_matrix = [a,b,c]

#atom positions in one cell (fractional coordinates)
file_name = '../benchmark-structures/POSCAR-diamond-cubic-unit-cell.poscar'
atom_position_onecell = read_unitcell(file_name)

#specify lattice constant, species, composition, and cell size
lc = 5.86 #lattice constant
species_list = ['Si','Ge','Sn']
conc_list = [0.125,0.625,0.25]
num_cell = [20,20,20] #number of replications of conventional cell along each axis

atom_positionlist,species_list,cell_geometry = generate(cell_matrix,atom_position_onecell,num_cell,lc,species_list,conc_list)

output_name = '../benchmark-structures/test.xyz'
write_xyz.write_benchmark(atom_positionlist,species_list,cell_geometry,output_name)
