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

 * vasp_build: Build model by atom substitution or atom selection based on a POSCAR file
 * vasp_read: Read DOSCAR, OUTCAR, POSCAR, and OSZICAR
 * vasp_plot.plot_poscar: Plot POSCAR/CONTCAR model (also support color mapping of atom properties), Required files: POSCAR/CONTCAR or POSCAR with data of atom properties
 * vasp_plot.plot_dos: Plot DOS (PDOS, LDOS, TDOS) information. Required files: DOSCAR, OUTCAR, POSCAR
 * vasp_analyze.nn_map: Calculate the nearest neighbor information. Required file: POSCAR
 * vasp_analyze.simple_cna: Perform simple common neighbor analysis
 * vasp_analyze.estruct: Calculate structural energy (E_struct). Required files: CONTCAR, OUTCAR, POSCAR
 * vasp_write.write_poscar_with_force: Write atom force information into the POSCAR

Three-dimensional atom probe tomography (APT) postprocessing tools

 * apt_read.read_proxigram_csv: Read the concentration profile .csv file
 * apt_plot.plot_proxigram_csv: Plot the concentration profile

### Installation

pip install matsdp

### Release note

- version 0.1.5

 * Date: 20191106

## Usage

### Graphical User Interface (GUI)

- matsdp_gui.exe

### Running with Python environment

modules to import before using the vasp package

 * from matsdp.vasp import vasp_read
 * from matsdp.vasp import vasp_build
 * from matsdp.vasp import vasp_plot
 * from matsdp.vasp import vasp_analyze
 * from matsdp.vasp import vasp_write

modules to import before using the apt package

 * from matsdp.apt import apt_read
 * from matsdp.apt import apt_plot

## Tests

To run the tests, please run the runtests.py.
