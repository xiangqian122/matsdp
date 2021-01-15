
=================
Basic information
=================

MATSDP--The materials simulation and data processing toolkit.

Requirements
------------

- numpy
- scipy
- matplotlib
- scikit-learn

Functions
---------

- Vienna ab-initio simulation package (VASP) analyzing and postprocessing tools 
 
 * Build model by the following methods: atom substitution, atom selection, exfoliation (2D), make supercell, transformation (rotation + translation), reorientation, adding vacuum layer etc.
 * Read and write VASP inputs/outputs.
 * Visualization of model and results: Plot model based on the POSCAR/CONTCAR file, also support color mapping of atom properties); Plot DOS (PDOS, LDOS, TDOS); Plot band structure (including fat band).
 * Analyzing tools: Calculate the nearest neighbor information, perform simple common neighbor analysis, calculate structural energy, overlap peak analyzer of DOS, get band gap.
 * VASP tools: Check VASP errors/warnings and give solutions; Check the job status of multiple jobs; Check lattice parameters of multiple VASP jobs; Conversion of coordinate systems (Fractional/Cartesian);
 * Job management: Processing of multiple VASP jobs automatically; Write task summary/report of multiple VASP jobs.

- Three-dimensional atom probe tomography (APT) postprocessing tools

 * Read the concentration profile .csv file
 * Plot the concentration profile
 
- DVM tools

 * Read the .input, .incar, .otput files
 * Write the .input, .incar, IND.DAT files
 * Write the interatomic energy (IE) files (including the IEs of the first nearest neighbor atoms)
 * The .incar file can also be prepared by atom selection from the vasp_build function in the vasp module 

- PMS tools

 * Job management: Processing of multiple VASP jobs automatically
 * Write task summary (of VASP jobs)

- Others tools

 * file format conversion
 * fig2pdf (converting multiple images to a single .pdf file)
 
Installation
------------

pip install matsdp

Release note
------------

- version 0.2.2

 * Date: 20210116

======
Usage
======

Graphical User Interface (GUI)
------------------------------
- matsdp_gui.exe

Running with Python environment
-------------------------------

- modules that may be imported before using the vasp package

 * from matsdp.vasp import vasp_read
 * from matsdp.vasp import vasp_build
 * from matsdp.vasp import vasp_plot
 * from matsdp.vasp import vasp_analyze
 * from matsdp.vasp import vasp_write

- modules that may be imported before using the apt package

 * from matsdp.apt import apt_read
 * from matsdp.apt import apt_plot
 
- modules that may be imported before using the dvm package

 * from matsdp.dvm import dvm_read
 * from matsdp.dvm import dvm_build
 * from matsdp.dvm import dvm_analyze
 * from matsdp.dvm import dvm_write
 * from matsdp.dvm import dvm_default
 * from matsdp.dvm import dvm_help

- modules that may be imported before using the pms package

 * from matsdp.pms import project_manager
 * from matsdp.pms import task_manager
 
======
Tests
======

To run the tests, please run the runtests.py in the "tests" directory.
