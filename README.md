# matsdp

The materials simulation and data processing toolkit.

## Basic information

MATSDP--The materials simulation and data processing toolkit.

For more information please contact dianwuwang@163.com.

### Requirements

- numpy
- scipy
- matplotlib
- scikit-learn

### Functions

Vienna ab-initio simulation package (VASP) analyzing and postprocessing tools 

 * Build model by atom substitution or atom selection based on a POSCAR file
 * Read information DOSCAR, OUTCAR, POSCAR, and OSZICAR
 * Plot Plot model in the POSCAR/CONTCAR (also support color mapping of atom properties), Required files: POSCAR/CONTCAR or POSCAR with data of atom properties
 * Plot DOS (PDOS, LDOS, TDOS) information. Required files: DOSCAR, OUTCAR, POSCAR
 * Calculate the nearest neighbor information. Required file: POSCAR
 * Perform simple common neighbor analysis
 * Calculate structural energy (E_struct). Required files: CONTCAR, OUTCAR, POSCAR
 * Write atom force information into the POSCAR

Three-dimensional atom probe tomography (APT) postprocessing tools

 * Read the concentration profile .csv file
 * Plot the concentration profile
 
DVM tools

 * Read the .input, .incar, .otput files
 * Write the .input, .incar, IND.DAT files
 * Write the interatomic energy (IE) files (including the IEs of the first nearest neighbor atoms)
 * The .incar file can also be prepared by atom selection from the vasp_build function in the vasp module 

Others tools

 * file format conversion

### Installation

pip install matsdp

### Release note

- version 0.1.9

 * Date: 20191112

## Usage

### Graphical User Interface (GUI)

- matsdp_gui.exe

### Running with Python environment

modules that may be imported before using the vasp package

 * from matsdp.vasp import vasp_read
 * from matsdp.vasp import vasp_build
 * from matsdp.vasp import vasp_plot
 * from matsdp.vasp import vasp_analyze
 * from matsdp.vasp import vasp_write
 * from matsdp.vasp import vasp_default
 * from matsdp.vasp import vasp_help
 
modules that may be imported before using the apt package

 * from matsdp.apt import apt_read
 * from matsdp.apt import apt_plot
 
modules that may be imported before using the dvm package

 * from matsdp.dvm import dvm_read
 * from matsdp.dvm import dvm_build
 * from matsdp.dvm import dvm_analyze
 * from matsdp.dvm import dvm_write
 * from matsdp.dvm import dvm_default
 * from matsdp.dvm import dvm_help

## Tests

To run the tests, please run the runtests.py.
