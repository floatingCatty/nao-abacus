# op2c: Python Bindings for nao-abacus

`op2c` provides Python bindings for the `Op2c` class from the `nao-abacus` project, enabling the calculation of two-center integrals (such as overlap and position operators) directly from Python.

## Prerequisites

Before installing, ensure the following dependencies are present in your environment:

*   **Compilers**: A C++ compiler with C++14 support (e.g., `g++`, `clang`) and a Fortran compiler (e.g., `gfortran`).
*   **CMake**: Version 3.15 or higher.
*   **Math Libraries**: BLAS and LAPACK libraries (e.g., OpenBLAS, MKL, or Netlib).
*   **Python**: Version 3.8+.
*   **MPI** (Optional): If running in an MPI environment, an MPI implementation (e.g., OpenMPI, MPICH) is recommended.

### Conda Environment Setup (Recommended)
You can easily set up the dependencies using Conda:

```bash
conda create -n op2c python=3.11
conda activate op2c
conda install -y compilers openblas fftw pkg-config cmake make
```

## Installation

Install the package using `pip`:

```bash
pip install .
```

This will compile the C++ extension and install the `op2c` package into your Python environment.

## Functionality

The `Op2C` class allows you to:
*   Initialize orbital and pseudopotential data used by `nao-abacus`.
*   Compute **Overlap** integrals between orbitals.
*   Compute **Position** operator integrals.
*   Compute **Beta** projector integrals.

It handles data structures transparently, converting C++ vectors/matrices to NumPy arrays for easy manipulation in Python.

## Usage Example

Here is a simple example of how to initialize `Op2C` and calculate overlap integrals.

```python
from op2c import Op2C
import numpy as np

# 1. Initialize Op2C
# params: ntype, nspin, lspinorb, orbital_dir, orbit_files, pseudo_dir, pseudo_files, ...
# Note: For serial execution, `comm` can be omitted or set to None.
op = Op2C(
    ntype=1, 
    nspin=1, 
    lspinorb=False,
    orb_dir="./path/to/orbitals/", 
    orb_name=["C_orb_fname.orb"],
    psd_dir="./path/to/pseudos/", 
    psd_name=["C_pseudo_fname.upf"],
    log_file="op2c.log"
)

# 2. Define Inputs
itype = 0  # Atom type index for atom I
jtype = 0  # Atom type index for atom J
Rij = [2.0, 0.0, 0.0] # Vector form atom I to J (in Bohr)

# 3. Compute Overlap
# Returns a 1D flattened numpy array of the overlap matrix block
S_flat = op.overlap(itype, jtype, Rij, is_transpose=False)

# Convert to 2D matrix (assuming you know the dimensions, e.g., 13x13 for d-orbitals)
# S_matrix = S_flat.reshape((n_orbitals_i, n_orbitals_j))

print("Overlap Norm:", np.linalg.norm(S_flat))

# 4. Compute Overlap & Position
S, Sx, Sy, Sz = op.overlap_position(itype, jtype, [0,0,0], [2.0,0,0], False)
```
