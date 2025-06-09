# ionisation_fractions

This repository provides a full pipeline for computing the ionisation fractions of metals in RAMSES-RT simulation snapshots, using a custom chemical network built with KROME.

## Overview

The goal is to:
- Build a chemical network adapted to the chosen metal elements and ionisation stages.
- Use KROME to generate a Fortran module from this network.
- Compile and run a Fortran code, calling the KROME Fortran module, to post-process RAMSES-RT snapshots and compute the ionisation state of the gas.

## Repository Structure

```
ion_frac_ramses_rt/
├── chemical_network_builder.py       # Builds chemical network text file from element dict
├── generate_krome_module.py          # Runs KROME and prepares module, including compilation
├── prepare_params.py                 # Writes parameter files for the global fortran code
├── reaction_rates/                   # Contains external reaction rate files (e.g., Verner, Voronov)
├── krome_files_to_change/            # Files to be copied into your KROME/src/
├── src_f90/                          # Fortran code that reads RAMSES snapshots and computes ion fractions
├── stellar_seds/                     # Stellar spectra, necessary to compute the cross-section of photoionisation
├── template_params.dat               # Template for the parameter file
├── example_run.ipynb                 # Demo notebook
├── README.md                         # This file
```

## Getting Started

1. **Download and install KROME**

   Download from http://kromepackage.org and extract it somewhere on your machine. This repo does *not* include KROME.

2. **Apply necessary modifications to KROME**

   You must overwrite a few source files in the `KROME/src/` folder. Do:

   ```bash
   cp krome_files_to_change/* /path/to/krome/src/
   ```

3. **Generate a chemical network and compile**

   Run the following, for example, in Python:

   ```python
   from generate_krome_module import generate_krome_module

   elements = {'C': 4, 'Si': 4}
   krome_path = "/absolute/path/to/your/krome"
   generate_krome_module(elements, krome_path)
   ```

   This will:
   - Generate a text file chemical network
   - Run KROME’s python wrapper
   - Move all generated Fortran files into a dedicated folder
   - Compile the final Fortran executable in `src_f90/`
   - Copy `reactions_verbatim.dat` into `src_f90/` (needed at runtime)
   

4. **Prepare parameter file**

Use `prepare_params.py` to generate the global parameter file required to run the Fortran code. This script takes as input the path to a template paramter file (included in this directory), the desired oath to the output files, the path to the RAMSES output, the snapshot number (`snapnum`), the element dictionary, and the path to Ramses-format stellar library, one of which is included in this directory. 

Example usage:

```python
from prepare_params import prepare_basic_param_file

prepare_basic_param_file(template_params, output_path, ramses_repo, snapnum, elements, sed_dir)
```

This creates a `params_krome.dat` and a global parameter file file in output_path, ready for use with the `ion_fractions` executable.


5. **Run the code**

   From inside `src_f90/`, run the `ion_fractions` executable. You need a parameter file, which can be created using `prepare_params.py`.

## Example

See `example_run.ipynb` for a minimal test case.

## Important Notes

- The main Fortran code is in `src_f90/ion_fractions.f90` and modules.
- The pipeline expects to find RAMSES outputs and a configuration file describing where to find fields like temperature, metallicity, HII, HeII, etc.
- The user must manually provide paths and adapt the parameter file.
- The chemical network always includes dummy H and He species to comply with KROME requirements.

## Dependencies

- Python ≥ 3.8
- gfortran with MPI (mpif90)
- KROME installed and modified
- RAMSES-RT output files (not included)

## Acknowledgements

KROME is developed by Tommaso Grassi, Stefano Bovino, and collaborators. See http://kromepackage.org/ for full credits.

This pipeline was created by Valentin Mauerhofer to post-process cosmological simulations with RAMSES-RT.
