#!/usr/bin/env python3
"""
Test script for the new CASPT2 parsing function.
"""

import sys
import os
sys.path.append('.')

from parse_rasscf_logfiles import extract_caspt2_root_energies

def test_caspt2_parsing():
    """Test the CASPT2 parsing function with the log file."""
    
    # Path to the log file
    log_file_path = "log_files/orca_geometry_optimisation_225_120_BP86_def2-SVP_def2J_scf_cas1_ras12i11_cas8.log"
    
    if not os.path.exists(log_file_path):
        print(f"Log file not found: {log_file_path}")
        return
    
    # Read the log file
    with open(log_file_path, 'r') as f:
        log_content = f.read()
    
    # Extract CASPT2 energies
    root_energies = extract_caspt2_root_energies(log_content)
    
    print("CASPT2 Root Energies:")
    print("=" * 50)
    
    for root_num in sorted(root_energies.keys()):
        energies = root_energies[root_num]
        print(f"Root {root_num}:")
        print(f"  Reference Energy: {energies['reference_energy']:.6f} Hartree")
        print(f"  Total Energy:     {energies['total_energy']:.6f} Hartree")
        print(f"  Correlation:      {energies['total_energy'] - energies['reference_energy']:.6f} Hartree")
        print()
    
    print(f"Total roots found: {len(root_energies)}")

if __name__ == "__main__":
    test_caspt2_parsing()
