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

Vienna ab-initio simulation package (vasp) analyzing and postprocessing tools 

 * vasp_build.substitution: Autosubstitution of element(s) based on a POSCAR file and a .subst file. Required input: POSCAR, .subst file
 * vasp_read: Read DOSCAR, OUTCAR, POSCAR, and OSZICAR
 * vasp_plot.plot_poscar: Plot POSCAR/CONTCAR model, Required input: POSCAR/CONTCAR
 * vasp_plot.plot_dos: Plot DOS (PDOS, LDOS, TDOS) information. Required input: DOSCAR, OUTCAR, POSCAR
 * vasp_analyze.nn_map: Calclate the nearest ngighbor inofmation. Required input: POSCAR
 * vasp_analyze.estruct: Calculate structural energy (E_struct). Required input: CONTCAR, OUTCAR, POSCAR

3-dimensional atom probe tomography (APT) postprocessing tools

 * apt_read.read_proxigram_csv: Read the concentration profile .csv file
 * apt_plot.plot_proxigram_csv: Plot the concentration profile

### Installation

pip install matsdp

### Release note

- version 0.1.0

 * Date: 20190912

## Usage

### Graphical User Interface (GUI)

- matsdp_gui.exe

### Running with Python environment

modules to import before using the vasp package

 * from matsdp.vasp import vasp_read
 * from matsdp.vasp import vasp_build
 * from matsdp.vasp import vasp_plot
 * from matsdp.vasp import vasp_analyze

modules to import before using the apt package

 * from matsdp.apt import apt_read
 * from matsdp.apt import apt_plot

## tests

To run the tests. Please run the runtest.py.
