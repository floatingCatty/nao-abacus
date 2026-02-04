import numpy as np
from . import _estate

class Op2C:
    """
    Python wrapper for the Op2c C++ class.
    Provides methods to compute quantum operators using 2-center integrals.
    """
    def __init__(self, ntype, nspin, lspinorb, 
                 orb_dir, orb_name, 
                 psd_dir, psd_name, 
                 log_file, comm=None):
        """
        Initialize Op2C.

        Parameters:
        -----------
        ntype : int
            Number of atom types.
        nspin : int
            Number of spins.
        lspinorb : bool
            Whether to use spin-orbit coupling.
        orb_dir : str
            Directory containing orbital files.
        orb_name : list of str
            List of orbital file names.
        psd_dir : str
            Directory containing pseudopotential files.
        psd_name : list of str
            List of pseudopotential file names.
        log_file : str
            Path to log file.
        comm : MPI.Comm, optional
            MPI communicator. If provided, it will be converted to a C MPI_Comm handle.
            Defaults to None (which maps to MPI_COMM_WORLD or Serial in C++).
        """
        self.mpi_handle = 0
        if comm is not None:
             if hasattr(comm, 'py2f'):
                 self.mpi_handle = comm.py2f()
        
        self.op2c = _estate.Op2c(ntype, nspin, lspinorb, 
                                 orb_dir, orb_name, 
                                 psd_dir, psd_name, 
                                 log_file, self.mpi_handle)

    def overlap(self, itype, jtype, Rij, is_transpose=False):
        """
        Compute overlap integral.

        Parameters:
        -----------
        itype, jtype : int
            Atom types.
        Rij : list or array-like of 3 floats
            Vector connecting the two centers.
        is_transpose : bool
            Whether to transpose the result.

        Returns:
        --------
        numpy.ndarray
            The overlap matrix block.
        """
        return self.op2c.overlap(itype, jtype, Rij, is_transpose)

    def overlap_deriv(self, itype, jtype, Rij, is_transpose=False):
        """
        Compute overlap integral and its derivatives.

        Returns:
        --------
        tuple of numpy.ndarray
            (overlap, d_overlap_dx, d_overlap_dy, d_overlap_dz)
        """
        return self.op2c.overlap_deriv(itype, jtype, Rij, is_transpose)

    def overlap_position(self, itype, jtype, Ri, Rj, is_transpose=False):
        """
        Compute overlap and position operator integrals.

        Parameters:
        -----------
        Ri, Rj : list or array-like of 3 floats
            Positions of the two centers.

        Returns:
        --------
        tuple of numpy.ndarray
            (overlap, position_x, position_y, position_z)
        """
        return self.op2c.overlap_position(itype, jtype, Ri, Rj, is_transpose)

    def orb_r_beta(self, itype, ktype, Ri, Rk, is_transpose=False):
        """
        Compute <orb|r|beta> integrals.

        Parameters:
        -----------
        itype : list of int
            List of atom types for 'orb'.
        ktype : int
            Atom type for 'beta'.
        Ri : list of list of 3 floats
            List of positions for 'orb'.
        Rk : list or array-like of 3 floats
            Position for 'beta'.

        Returns:
        --------
        tuple of list of numpy.ndarray
            (ob, oxb, oyb, ozb)
        """
        return self.op2c.orb_r_beta(itype, ktype, Ri, Rk, is_transpose)

    def ncomm_IKJ(self, itype, idx, ktype, jtype, jdx, npol, is_transpose=False):
        """
        Compute non-local commutators? (ncomm_IKJ)

        Returns:
        --------
        tuple
            (ob, oxb, oyb, ozb, vx, vy, vz)
            First 4 are lists of matrices (numpy arrays), last 3 are complex vectors.
        """
        return self.op2c.ncomm_IKJ(itype, idx, ktype, jtype, jdx, npol, is_transpose)
