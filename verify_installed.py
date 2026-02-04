import sys
import os
import numpy as np
try:
    from mpi4py import MPI
except ImportError:
    print("Warning: mpi4py not found. MPI_Init might not be called.")


# Add build/bindings to path
# sys.path.append("/home/zhanghao/softwares/nao-abacus/build/bindings")
# Add source_estate to path to import wrapper
# sys.path.append("/home/zhanghao/softwares/nao-abacus")

from op2c import Op2C

def main():
    orb_dir = "/home/zhanghao/softwares/nao-abacus/test/pporb/"
    psd_dir = "/home/zhanghao/softwares/nao-abacus/test/pporb/"
    orb_name = ["C_gga_7au_100Ry_2s2p1d.orb"]
    psd_name = ["C_ONCV_PBE-1.0.upf"]
    
    print("Initializing Op2C...")
    # ntype=1, nspin=1, lspinorb=False
    try:
        op = Op2C(ntype=1, nspin=1, lspinorb=False,
                  orb_dir=orb_dir, orb_name=orb_name,
                  psd_dir=psd_dir, psd_name=psd_name,
                  log_file="op2c_test.log", comm=None)
    except Exception as e:
        print("Initialization failed:", e)
        sys.exit(1)
    
    print("Op2C initialized.")
    
    # Test overlap
    Rij = [2.0, 0.0, 0.0] # 2.0 Bohr separation
    itype = 0
    jtype = 0
    
    print(f"Computing overlap for Rij={Rij}...")
    try:
        S = op.overlap(itype, jtype, Rij)
        print("Overlap shape:", S.shape)
        print("Overlap matrix norm:", np.linalg.norm(S))
        print("Overlap sample element (0):", S[0])
        
        if S.size == 0:
            print("Error: Overlap matrix is empty!")
            sys.exit(1)
    except Exception as e:
        print("Overlap computation failed:", e)
        sys.exit(1)

    print("Computing overlap_position...")
    try:
        # overlap_position(itype, jtype, Ri, Rj, is_transpose)
        Ri = [0.0, 0.0, 0.0]
        Rj = Rij
        S, Sx, Sy, Sz = op.overlap_position(itype, jtype, Ri, Rj, is_transpose=False)
        print("Sx shape:", Sx.shape)
        print("Sx norm:", np.linalg.norm(Sx))
    except Exception as e:
        print("Overlap_position computation failed:", e)
        sys.exit(1)
    
    print("Verification passed!")

if __name__ == "__main__":
    main()
