#!/usr/bin/env python3
"""
Test script to demonstrate the parameterized geometry comparison functionality.
This script shows how to use the updated functions with different geometry combinations.
"""

from plot_comparison_for_single_geometries import create_combined_comparison
from detailed_comparison_analysis import create_comparison_summary

def test_different_geometries():
    """Test the comparison functions with different geometry combinations"""
    
    print("="*80)
    print("TESTING PARAMETERIZED GEOMETRY COMPARISON")
    print("="*80)
    
    # Test 1: Default geometries (Beta-Strand vs Alpha-Helix)
    print("\n1. Testing with default geometries:")
    print("   Beta-Strand (φ=210°, ψ=135°) vs Alpha-Helix (φ=240°, ψ=30°)")
    create_combined_comparison()
    
    # Test 2: Extended vs Turn geometries
    print("\n2. Testing with Extended vs Turn geometries:")
    print("   Extended (φ=0°, ψ=0°) vs Turn (φ=270°, ψ=270°)")
    create_combined_comparison(
        geometry1=(0, 0, "Extended"),
        geometry2=(270, 270, "Turn")
    )
    
    # Test 3: Different Beta-Strand vs Alpha-Helix combinations
    print("\n3. Testing with alternative Beta-Strand vs Alpha-Helix:")
    print("   Beta-Strand (φ=225°, ψ=135°) vs Alpha-Helix (φ=255°, ψ=30°)")
    create_combined_comparison(
        geometry1=(225, 135, "Beta-Strand-Alt"),
        geometry2=(255, 30, "Alpha-Helix-Alt")
    )

def test_analysis_only():
    """Test only the detailed analysis with custom geometries"""
    
    print("\n" + "="*80)
    print("TESTING DETAILED ANALYSIS WITH CUSTOM GEOMETRIES")
    print("="*80)
    
    # Test with different geometry combinations
    print("\nTesting Extended vs Turn analysis:")
    create_comparison_summary(
        geometry1=(0, 0, "Extended"),
        geometry2=(270, 270, "Turn")
    )

if __name__ == "__main__":
    # Test the plotting functionality
    test_different_geometries()
    
    # Test the analysis functionality
    test_analysis_only()
    
    print("\n" + "="*80)
    print("PARAMETERIZED TESTING COMPLETE")
    print("="*80)
    print("\nThe functions now accept custom geometry parameters:")
    print("- geometry1: (phi, psi, name) for first geometry")
    print("- geometry2: (phi, psi, name) for second geometry")
    print("\nExample usage:")
    print("create_combined_comparison(")
    print("    geometry1=(210, 135, 'Beta-Strand'),")
    print("    geometry2=(240, 30, 'Alpha-Helix')")
    print(")") 