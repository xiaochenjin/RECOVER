# Numerical simulation 
Simulate APT GeSn data with perturbed atomic positions and limited collection efficiency. Compute corresponding fkr functions. 

## Step 1: Generate a benchmark structure 

### Code:
   * benchmark-generate.py

### Instructions:
   * Information about the lattice is in _benchmark-structures/POSCAR-diamond-cubic-unit-cell.poscar_
   * Need to specify lattice constant, species, composition, and cell size in the code

## Step 2: Obtain structural information

Suppose you have generated the benchmark structure of GeSn alloy (GeSn-random-before-relax.xyz) and relax the structure using a potential (GeSn-random-relaxed.xyz)

### Code:
	* compute-pair-correlation.py
	* obtain-structure-info.py 

### Examples:
	python3 compute-pair-correlation.py --what_to_compute=original --whether_relax=Y --direction=r --pair_type=Sn-Sn #compute pair correlation before perturbation

### Instructions:
	* compute-pair-correlation.py can compute RDF (--direction=r), z-SDM (--direction=z) and xy-SDM (--direction=xy)
	* obtain-structure-info.py outputs KNN distance and coordination number for Kth shell. 

## Step 3: Compute KNN information of the benchmark structure

For each atom, determine the corresponding atoms at the Kth shell

### Code:
    * compute-KNN-info.py

### Examples:
	python3 compute-KNN-info.py --pair_type=Sn-Sn --K_raw=5 --lattice_constant=5.86
	python3 compute-KNN-info.py --pair_type=all --K_raw=5 --lattice_constant=5.86 #all atoms

### Instructions:
	The code outputs a KNN-info-{pair}.txt file with the output format: 
	atom index, atom indices of 1st shell, atom indices of 2nd shell, ...atom indices of Kth shell

## Step 4: Generated simulated APT measurement

### Code:
    * simulate-APT.py
	* compute-pair-coirrelation.py (parallelized)

### Examples:
    python3 simulate-APT.py --perturb_type=iso --mu=0 --sigma=1.0 --collect_efficiency=1.0 --N_config=50 #isotropic perturbation
    python3 simulate-APT.py --perturb_type=aniso --mu=0 --sigma_z=1.0 --sigma_xy=1.0 --collect_efficiency=1.0 --N_config=50 #anisotropic perturbation
	python3 compute-pair-correlation.py --what_to_compute=perturbed --mu=0 --sigma=1.0 --collect_efficiency=1.0 --N_config=50 --direction=r --pair_type=Sn-Sn --N_thread=2

### Instructions:
	* simulate-APT.py generates perturbed cells in xyz format (last column index of atoms in the original structure)
	* compute-pair-correlation.py compute cooresponding pair correlations  (_no need for only compute fkr_)

## Step 5: Compute fkr functions

### Code:
	* compute-KNN-track.py (parallelized)
    * compute-fkr.py
    * plot-fkr.py

### Examples:
	python3 compute-KNN-track.py --pair_type=Sn-Sn --K_raw=5 --lattice_constant=5.86 --perturb_type=iso --mu=0 --sigma=1.0 --collect_efficiency=1.0 --N_config=50 --N_thread=2
	
### Instructions:
	* compute-KNN-track.py computes RDF for true Kth-nearest-neighbor (KNN) for each Kth shell.
	* compute-fkr.py and plot-fkr.py computes and plot fkr functions, respectively.





