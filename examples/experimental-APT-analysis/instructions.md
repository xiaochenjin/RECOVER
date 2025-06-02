# Experimental APT data analysis
The purpose of this section is to retrieve SRO parameters from the experimental APT data. 
<br> The example APT data is measured based on a GeSn sample grown by MBE. 
<br> fkr files includes fkr functions of 50 shells which are numerically computed from 1000 configurations.

## Step 1: Post-process raw APT data

### Codes:
	* write_data.py 
	* write-data-specific-element.py
	* select-large-box.py
	* select-multiple-cubes.py

### Instructions:
	* write_data.py converts raw APT data format into .txt file
	* write-data-specific-element.py selects atomic information with specified species
	* select-large-box.py selects a box of atomic positions based on the specified range of box dimension in x, y, z direction
	* select-multiple-cubes.py divides the box into 10x10x10 nm^3 cubes. Run python select-multiple-cubes.py directly divides the box, and also generate "Nz-Nx-Ny.txt" which contains number of cubes in z, x, and y direction, respectively.

## Step 2: Compute RDF of specified pairs

### Codes:
	* compute-pair-correlation.py (parallelized)

### Examples:
	python3 compute-pair-correlation.py --direction=r --N_bin=10000 --N_thread=2

### Instructions:
	* compute-pair-correlation.py computes pair correlation of specified pairs

## Step 3: Smooth RDF using  the "moving-average" method

### Codes:
	* convert-rdf-smoothing.py

### Example:
	python3 convert-rdf-smoothing.py  --direction=r --smooth_window=200

### Instructions:
	smoothing length = 0.5*cell_dimension/N_bin*smooth_window

## Step 4: Retreive SRO parameters from raw APT data

### Codes:
	* recover-SRO.py (parallelized)

### Example:
	python3 recover-SRO.py --p_min=-3 --guess_index=1 --tol=1e-9 --KNN=1 --smooth_window='200' --K_raw=50 --R_lower=0.5 --R_upper=4.0 --N_thread=2
	python3 recover-SRO.py --p_min=-3 --guess_index=1 --tol=1e-9 --KNN=3 --smooth_window='200' --K_raw=50 --R_lower=4.0 --R_upper=6.0 --N_thread=2

### Instruction:
	For each data point, this generate the convergence of SRO parameter by increasing number of shells considered. 

## Step 5: Calibrate SRO parameter

### Codes:
	* calibrate-SRO-parameter.py

## Step 6: Compute SRO parameter after data screening, and plot the result

### Codes:
	* compute-SRO-after-screening.py
	* plot-SRO.py




