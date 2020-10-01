# convergence-test

A simple tool for plane-wave kinetic energy cutoff (encut) and kmesh convergence tests for density functional theory (DFT) calculations.

The main features include:

1. **A script to generate INCAR and KPOINTS files**

   - Using the input parameters provided, INCAR and KPOINTS files are created.
   (**Note:** The ENCUT for the INCAR file will be used for kmesh convergence tests and the kmesh for the KPOINTS file will be used for cutoff/encut convergence tests.)

2. **A script to create folders for cutoff and kmesh convergence, and update the cutoff and kmesh values in each INCAR and KPOINTS file**

3. **Plotting script to create cutoff and kmesh convergence plots.**

   - VASP calculations are imported using Pymatgen (energies extracted from vasprun.xml files)


The code currently primarily supports VASP calculations.


Usage
-------------

The scripts are intended to be used via the command-line. The built-in help (``-h``) option for each command provides a summary of available options.

Currently, the commands provided are:

-``generate_files``: For generating INCAR and KPOINTS input files.<br>
-``generate_folders``: For creating folders for cutoff and kmesh convergence tests, copying the files, and updating the cutoff and kmesh values. (Note: POSCAR (structure file) and POTCAR (pseudopotential file) need to be present in the folder for this command to work.)<br>
-``analyse_conv``: For extracting energies from VASP calculations, creating csv files, and plotting the convergence data. (Note: The command uses vasprun.xml files to extract energies and so complete DFT calculations performed using VASP is a prerequisite.)<br>

Requirements
---------------

The command-line scripts are written in Python and require the following packages:

<li> Click
<li> NumPy
<li> Pandas
<li> Matplotlib
<li> Pymatgen
