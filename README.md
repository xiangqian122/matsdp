
# matsdp

The materials simulation and data processing toolkit.

## Basic information

MATSDP--The materials simulation and data processing toolkit.

### Requirements

- numpy
- scipy
- matplotlib
- scikit-learn

### Functions

Vienna ab-initio simulation package (VASP) analyzing and postprocessing tools 

 * Build model by atom substitution or atom selection based on a POSCAR file
 * Read VASP inputs and outputs
 * Plot model in the POSCAR/CONTCAR (also support color mapping of atom properties).
 * Plot DOS (PDOS, LDOS, TDOS), band structure (including fat band).
 * Calculate the nearest neighbor information.
 * Perform simple common neighbor analysis
 * Calculate structural energy.
 * Write atom force information into the POSCAR file.

Three-dimensional atom probe tomography (APT) postprocessing tools

 * Extract information from the concentration profile.csv file.
 * Visualization of the concentration profile across the interface.
 
DVM tools

 * Read the .input, .incar, .otput files
 * Write the .input, .incar, IND.DAT files
 * Ouput the interatomic energy (IE) files (including the IEs of the first nearest neighbor atoms)
 * The .incar file can also be prepared by atom selection from the vasp_build function in the vasp module 

PMS tools

 * Write task summary (of VASP jobs)

Others tools

 * file format conversion
 * fig2pdf (converting multiple images to a single .pdf file)

### Installation

pip install matsdp

### Release note

- version 0.2.1

 * Date: 20201209

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
 * from matsdp.vasp import vasp_tools
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
 
modules that may be imported before using the pms package

 * from matsdp.pms import project_manager
 * from matsdp.pms import task_manager

## Tests

To run the tests, please run the runtests.py in the "tests" directory.
