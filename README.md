### Instructions for pump-probe UV-X-ray spectra of ICN

This directory contains the programs and scripts to generate and run all the necessary files to obtain
the pump-probe UV-Xray spectra of ICN. The procedure consists of generating the Xray response of the
excited state of ICN in a grid and use this grid in a wave-packet propagation of the excited state
to compute the response. 

“Vibronic Dynamics of Photodissociating ICN from Simulations of Ultrafast X-Ray Absorption Spectroscopy,” U. N. Morzan,* P. E. Videla,* M. B. Soley,* E. T. J. Nibbering, V. S. Batista, Angewandte Chemie, 59 (2020) 20044-20048, https://doi.org/10.1002/anie.202007192

##   Generate Xray response on grid

The directory 'GENERATE_GRID/' contains all programs and scripts to generate the Xray response on a grid 
based on QCHEM calculations.

1) In the main directory, the program 'write-input.f90' generates a set of initial coordinates for ICN, 
   distributed according to a grid in R_{IC},theta,r_{CN} (Jacobi coordinates). 
   The variables nptx/npty/nptz control the number of points in the grid, whereas x_min/y_min/z_min and 
   x_max/y_max/z_max control the limits of the grid.

   To compile the program: gfortran write-input.f90 -o write-input 
   To run the program: ./write-input

   The outputs are:
	> qchem_####.coor: contains coordinates for point #### in the grid.
	> movie.xyz: movie of the configurations of the grid. 
	> position.dat: list of file index,R_{IC},theta,r_{CN}. 


2) In the subdirectory 'GENERATE_GRID/GENERATE_RESPONSE/', there are a set of subdirectories 'PROBE_C/', 'PROCE_N/'
   and 'PROBE_I/'. In each directory the Xray response of the corresponding atom is computed. Although
   the same structure is maintained between directories, the QCHEM input options differ in each case.
   The following instructions will focused on the 'PROCE_C/' subdirectory. 


3) In the subdirectory 'PROBE_C/' the script 'cat_script.sh' generates the corresponding QCHEM input files
   to compute the response on a grid. The script takes the coordinates from the files 'qchem_####.coor' 
   previously generated and the QCHEM options from the file 'input-trimed.in'. The QCHEM calculation
   correspond to a MOM calculation to generate the excited state and a CIS calculation to compute the Xray
   response.

   To run the script: ./cat_script.sh

   The outputs are:
	 > qchem_****.in: contains QCHEM input file for point **** in the grid.


   The subdirectory also contains a script 'generate_slurm_script.sh' that generates a set of 'run_script_#.sh'
   files to run QCHEM on Grace. This script uses 'slurm_template.inp' as template for the SLURM queue and divides   
   'num_files' into 'num_proc' (parameters inside 'generate_slurm_script.sh').

   To run the script: ./generate_slurm_script.sh

   The outputs are:
	> run_script_*.sh: contains SLURM script to run on Grace supercomputer for a set of configurations.


4) In the subdirectory 'ENERGY_READER/' there are a set of Python scripts that reads the outputfile of the QCHEM 
   calculations and prepare the data for the grid generation. The script 'run_energy_scripts.sh' runs all the 
   programs of this folder in order. 

   To run the script: ./run_energy_scripts.sh

   The different scripts are: 
   
      (a) energy_reader.sh: read energies from QCHEM outputfiles.

	  On output:
		> energy.dat: contains (in principle) two energies for each QCHEM outputfile corresponding to 
                               ground-state and MOM calculations.
 

      (b) energy_cleaner.py: separate contribution of first ground-state calculation from second MOM calculation.

	  On output:
		> energy_gs.dat: contains energy of ground-state calculation for each QCHEM file.
		> energy_mom.dat: contains energy of MOM calculation for each QCHEM file.

      (c) energy_plot.py (optional): plots energy cuts along R_{IC} variable.

	  On output:
		> energy_*.ps: contains cuts of GS and MOM energies along the R_{IC} variable. Each files 
				correspond to a different value of angle theta. In each graph the different sets 
				of curves correspond to different r_{CN} values. These script is used to visually
				check the amount of points not converged in the QCHEM calculations and the 
				smoothness (or lack of) in the curves obtained.    

      (d) energy_selector.py: performs selection of energies based on threshold. Due to the incorrect convergence
			      of some of the QCHEM calculations, an energy selection to guarantee that the curves
			      obtained are smooth is performed. To this purpouse, some points are excluded from 
			      the grid.    

	  On output:
		> energy_mom_excl.dat: contains the sets of points that will be excluded from the grid.
		> energy_mom_incl.dat: contains the sets of points that are included in the grid.
		> energy_selection_*.ps: contains cuts of GS and MOM with included and excluded energies along the 
					  R_{IC} variable. Each files correspond to a different value of angle theta.
					  In each graph the different sets of curves correspond to different r_{CN} 
					  values. 

5) In the subdirectory 'GENERATE_MATRIX/' the script 'generate_matrix.py' generates the Xray response matrix. 
   It requires the list of points to exclude 'exclusion_list.dat', which is a link to the 'energy_mom_excl.dat' 
   computed in the previous step. The variable 'nlines' inside the scripts determines the number of Xray transitions 
   to be considered.

   To run the script: python generate_matrix.py

   The outputs are:
	> response_mat.dat: contains the Xray response in a grid (in the format theta,R_{IC},r_{CN},energy, intensity)

6) Finally, in the subdirectory 'INTERPOLATION', the program 'interpolator.f90' interpolates the Xray response matrix
   using an Inverse Distance Weighting method to generate a regular grid in the Jacobi variables. The variable 'pvalue'
   inside the program determines the power parameter for the interpolation scheme.

   To compile the program: gfortran interpolator.f90 -o interpolator
   To run the program: ./interpolator

   The outputs are:
	> response_mat_interp.txt: contains the Xray response in a regular grid in the format 
           theta,R_{IC},r_{CN},energy,intensity


# Wavepacket propagation

The directory 'WP_PROPAGATION/' contains the program to run the wave-packet dynamics using the Xray response on a grid
to compute the pump-probe UV-Xray spectra of ICN. The different subdirectories 'PROBE_*/ correspond to the Xray
response of the different atoms (namely, the Xray grid is different). The program uses the 'response_mat_interp.txt' 
file generated in the previous step. See 'main.f' for additional details of the wave-packet propagation.

   To compile the program: make clean
                           make
   To run the program(*): ./3dpropagation

   (*) It is highly recommended to run the program in a cluster. A template for a SLURM queue is provided as the script
       'run_program.sh'

   The most important outputs are:
	> uv-xray-spectrum.dat: contains the time-dependent UV-Xray spectrum computed using the whole wave-packet 
                                 (in the format energy,intensity,time).
	> expect-uv-xray-spect.dat: contains the time-dependent UV-Xray spectrum computed using the average value
				     of the wave-packet (in the format energy,intensity,time).
	> norms.dat: 
  s the norms of wave-packet in each surface.


