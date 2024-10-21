# DIFFRACTION PATTERNS BASED BIASED MD
## Installation
  
  1. Download PLUMED (version >= 2.7)
  
  2. Configure PLUMED
     ```bash
       ./configure --prefix=/path/to/plumed PYTHON_BIN=python3 СС=gcc CXX=g++
     ```
     
  3. Copy `src` folder from the current repository to a `plumed/` directory
  4. Compile PLUMED
     ```bash
       make -j N 
     ```
  5. Install PLUMED
     ```bash
       make -j N install
     ```
  6. Turn on PLUMED package in LAMMPS and compile it
     
## Usage
  The method can work in two scenarios: pulling and pushing.
  The sign of a number specified in ```FORCE_COEFF``` keyword determines which algorithm will be used. If it is positive, pushing scenario will be turned on automatically. Otherwise, pulling algorithm will be used. Coefficient zero is not allowed. If you want to turn off the derivation of the diffraction patterns similarity, add the ```NO_DIFF``` keyword to the input file.
  
  ## Calculations, data analysis and pictures
  `calculations` folder in this repository contains the examples of pulling and pushing scenarios application. 
  
  There are 4(5) input files you will need:
  - `plumed.dat`, standard PLUMED input with our keywords
  - `para.lmp`, standard LAMMPS input
  - `*.data`, a file with a configured force-field
  - `names.txt`, a file with atomic numbers of atoms listed in the para.lmp
  - `pattern.dat`, a file with a reference diffraction pattern (only for pull scenario)
  
  To run the calculation you need to run LAMMPS and use `para.lmp` as an input file:
  
  ```bash
    lmp_mpi -in para.lmp
  ```

  There will be 4 output files: 
  - `p.log`, standard PLUMED output
  - `log.lammps`, standard LAMMPS output
  - `COLVAR`, output file with diffraction patterns differences (with pushing scenario, it will be the difference between `i` and `i+1` steps; with pulling scenario, it will be the difference between 1 and `i` steps)
  - `para.lammpstrj`, MD trajectory
  
  `pictures` folder contains the scripts you will need to generate the pictures from the article. To do that, run
  the `generate_picture.py` file.
  
  Push:
   - `diffractions_push` is a file with diffraction patterns generated for each step of a MD trajectory.
   - `global_diff` is a file with a similarities between `i`th step of MD and the first structure.
   - `mono.dat` and `pushed.dat` are the files with generated powder diffraction patterns of the initial and final structures.
   - `mono.png` and `pushed.png` are the images of the structures made in VESTA program.
   - `ref_peaks.txt` and `end_peaks.txt` are the files with the lists of diffraction peaks for the initial and final structures.
  
  Pull:
   - `COLVAR` is a file with a similarities between `i`th step of MD and the first structure.
   - `start.dat`, `pulled.dat` and `ortho.dat` are the files with generated powder diffraction patterns 
     of the initial, final and reference structures.
   - `start.png`, `pulled.png` and `ortho.png` are the images of the structures made in VESTA program.
   - `start_peaks.txt`, `cur_peaks.txt` and `ref_peaks.txt` are the files with the lists of diffraction peaks 
     for the initial, final and reference structures.
  
  Polymorphs:
   - `mono.dat` and `ortho.dat` are the files with generated powder diffraction patterns 
     of the Form I and Form II of paracetamol.
   - `mono.png` and `ortho.png` are the images of the structures made in VESTA program.
  
  The generation and comparison of diffraction patterns is implemented in Critic2 package. With this program, the data in the files
  `diffractions_push` and `global_diff` can be obtained.
