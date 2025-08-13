import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def analyze_geometry_differences(csv_file, phi, psi, geometry_name):
    """Perform detailed analysis of MS-CASPT2 vs CASPT2 differences for a specific geometry"""
    
    # Load and filter data
    df = pd.read_csv(csv_file)
    geometry_data = df[(df['phi'] == phi) & (df['psi'] == psi)].copy()
    
    if geometry_data.empty:
        print(f"No data found for geometry phi={phi}, psi={psi}")
        return None
    
    print(f"\n{'='*80}")
    print(f"DETAILED ANALYSIS: {geometry_name} (φ={phi}°, ψ={psi}°)")
    print(f"{'='*80}")
    
    # 1. Overall Statistics
    print(f"\n1. OVERALL STATISTICS:")
    print(f"{'Transition':<20} {'CASPT2 (eV)':<12} {'MS-CASPT2 (eV)':<12} {'Difference (eV)':<15} {'% Difference':<12}")
    print("-" * 80)
    
    transitions = sorted(geometry_data['transition'].unique(), key=lambda x: ('n_to_pi' in x, x))
    
    for transition in transitions:
        trans_data = geometry_data[geometry_data['transition'] == transition]
        if not trans_data.empty:
            caspt2_e = trans_data.iloc[0]['caspt2_energy_difference_ev']
            ms_caspt2_e = trans_data.iloc[0]['ms_caspt2_energy_difference_ev']
            diff = ms_caspt2_e - caspt2_e
            pct_diff = (diff / caspt2_e) * 100
            
            print(f"{transition:<20} {caspt2_e:<12.3f} {ms_caspt2_e:<12.3f} {diff:+<15.3f} {pct_diff:+<12.2f}%")
    
    # 2. State-by-State Analysis
    print(f"\n2. STATE-BY-STATE ANALYSIS:")
    print(f"{'State':<8} {'CASPT2 (Hartree)':<15} {'MS-CASPT2 (Hartree)':<15} {'Difference (Hartree)':<18}")
    print("-" * 65)
    
    states = sorted(geometry_data['end_state'].unique())
    for state in states:
        state_data = geometry_data[geometry_data['end_state'] == state]
        if not state_data.empty:
            caspt2_abs = state_data.iloc[0]['caspt2_energy']
            ms_caspt2_abs = state_data.iloc[0]['ms_caspt2_energy']
            diff_abs = ms_caspt2_abs - caspt2_abs
            
            print(f"{state:<8} {caspt2_abs:<15.6f} {ms_caspt2_abs:<15.6f} {diff_abs:+<18.6f}")
    
    # 3. Transition Dipole Analysis
    print(f"\n3. TRANSITION DIPOLE MOMENT ANALYSIS:")
    print(f"{'Transition':<20} {'Total (D)':<10} {'dx (D)':<10} {'dy (D)':<10} {'dz (D)':<10}")
    print("-" * 70)
    
    for transition in transitions:
        trans_data = geometry_data[geometry_data['transition'] == transition]
        if not trans_data.empty:
            total_dipole = trans_data.iloc[0]['totald_debyes']
            dx = trans_data.iloc[0]['dx_debyes']
            dy = trans_data.iloc[0]['dy_debyes']
            dz = trans_data.iloc[0]['dz_debyes']
            
            print(f"{transition:<20} {total_dipole:<10.3f} {dx:<10.3f} {dy:<10.3f} {dz:<10.3f}")
    
    # 4. Configuration Analysis
    print(f"\n4. ELECTRONIC CONFIGURATION ANALYSIS:")
    print(f"{'Transition':<20} {'Configuration':<15} {'Coefficient':<12} {'Weight':<10}")
    print("-" * 65)
    
    for transition in transitions:
        trans_data = geometry_data[geometry_data['transition'] == transition]
        if not trans_data.empty:
            config = trans_data.iloc[0]['config']
            coeff = trans_data.iloc[0]['coeff']
            weight = trans_data.iloc[0]['weight']
            
            print(f"{transition:<20} {config:<15} {coeff:<12.3f} {weight:<10.3f}")
    
    # 5. Key Insights
    print(f"\n5. KEY INSIGHTS:")
    
    # Calculate average differences
    energy_diffs = []
    for transition in transitions:
        trans_data = geometry_data[geometry_data['transition'] == transition]
        if not trans_data.empty:
            diff = trans_data.iloc[0]['ms_caspt2_energy_difference_ev'] - trans_data.iloc[0]['caspt2_energy_difference_ev']
            energy_diffs.append(diff)
    
    avg_diff = np.mean(energy_diffs)
    std_diff = np.std(energy_diffs)
    max_diff = max(energy_diffs, key=abs)
    min_diff = min(energy_diffs, key=abs)
    
    print(f"   • Average energy difference (MS-CASPT2 - CASPT2): {avg_diff:+.3f} ± {std_diff:.3f} eV")
    print(f"   • Largest energy difference: {max_diff:+.3f} eV")
    print(f"   • Smallest energy difference: {min_diff:+.3f} eV")
    
    # Identify which transitions are most affected
    max_diff_transition = None
    min_diff_transition = None
    for transition in transitions:
        trans_data = geometry_data[geometry_data['transition'] == transition]
        if not trans_data.empty:
            diff = trans_data.iloc[0]['ms_caspt2_energy_difference_ev'] - trans_data.iloc[0]['caspt2_energy_difference_ev']
            if abs(diff) == abs(max_diff):
                max_diff_transition = transition
            if abs(diff) == abs(min_diff):
                min_diff_transition = transition
    
    print(f"   • Most affected transition: {max_diff_transition} ({max_diff:+.3f} eV)")
    print(f"   • Least affected transition: {min_diff_transition} ({min_diff:+.3f} eV)")
    
    # Analyze dipole moments
    dipoles = []
    for transition in transitions:
        trans_data = geometry_data[geometry_data['transition'] == transition]
        if not trans_data.empty:
            dipole = trans_data.iloc[0]['totald_debyes']
            dipoles.append(dipole)
    
    print(f"   • Average transition dipole moment: {np.mean(dipoles):.3f} D")
    print(f"   • Strongest transition: {transitions[np.argmax(dipoles)]} ({max(dipoles):.3f} D)")
    print(f"   • Weakest transition: {transitions[np.argmin(dipoles)]} ({min(dipoles):.3f} D)")
    
    return geometry_data

def display_absolute_energies_table(csv_file, phi, psi, geometry_name):
    """Display a table of absolute energies (CASPT2 and MS-CASPT2) for ground and excited states for a given geometry, with high precision, all in eV. Includes ground state, transition name, and root numbers."""
    df = pd.read_csv(csv_file)
    geometry_data = df[(df['phi'] == phi) & (df['psi'] == psi)].copy()
    if geometry_data.empty:
        print(f"No data found for geometry phi={phi}, psi={psi}")
        return None
    print(f"\n{'='*130}")
    print(f"ABSOLUTE ENERGIES TABLE: {geometry_name} (φ={phi}°, ψ={psi}°)")
    print(f"{'='*130}")
    print(f"{'Start':<7} {'End':<7} {'Transition':<25} {'CASPT2 Root':<12} {'MS-CASPT2 Root':<15} {'CASPT2 Energy (eV)':<25} {'MS-CASPT2 Energy (eV)':<28} {'Difference (eV)':<22}")
    print("-" * 130)
    # Sort by start_state, end_state for clarity
    geometry_data = geometry_data.sort_values(by=['start_state', 'end_state'])
    # Find all unique states (including ground state)
    all_states = set(geometry_data['end_state']).union(set(geometry_data['start_state']))
    # For each unique state, print its absolute energy (ground state: start_state==end_state==1)
    for state in sorted(all_states):
        # Try to find a row where start_state==end_state==state (ground state)
        ground_row = geometry_data[(geometry_data['start_state'] == state) & (geometry_data['end_state'] == state)]
        if not ground_row.empty:
            row = ground_row.iloc[0]
            start_state = row['start_state']
            end_state = row['end_state']
            transition = row['transition'] if 'transition' in row and pd.notnull(row['transition']) else '(ground state)'
            caspt2_root = row['caspt2_root'] if 'caspt2_root' in row and pd.notnull(row['caspt2_root']) else 'N/A'
            ms_caspt2_root = row['ms_caspt2_root'] if 'ms_caspt2_root' in row and pd.notnull(row['ms_caspt2_root']) else 'N/A'
            caspt2_abs_ev = row['caspt2_energy'] * 27.211407953
            ms_caspt2_abs_ev = row['ms_caspt2_energy'] * 27.211407953
            diff_abs_ev = ms_caspt2_abs_ev - caspt2_abs_ev
            caspt2_abs_h = row['caspt2_energy']
            ms_caspt2_abs_h = row['ms_caspt2_energy']
            diff_abs_h = ms_caspt2_abs_h - caspt2_abs_h
            def sig_str(val):
                if val == 0:
                    return f"{0:.10e}"
                return f"{val:.10e}" if abs(val) < 1e-4 or abs(val) > 1e+6 else f"{val:.10g}"
            print(f"{int(start_state):<7} {int(end_state):<7} {transition:<25} {caspt2_root:<12} {ms_caspt2_root:<15} {sig_str(caspt2_abs_ev):<25} {sig_str(ms_caspt2_abs_ev):<28} {sig_str(diff_abs_ev):<22}")
            print(f"{'':<7} {'':<7} {'':<25} {'':<12} {'':<15} {'(Hartree)':<25} {'(Hartree)':<28} {'(Hartree)':<22}")
            print(f"{'':<7} {'':<7} {'':<25} {'':<12} {'':<15} {sig_str(caspt2_abs_h):<25} {sig_str(ms_caspt2_abs_h):<28} {sig_str(diff_abs_h):<22}")
    # Now print all other transitions (where start_state != end_state)
    for _, row in geometry_data.iterrows():
        if row['start_state'] == row['end_state']:
            continue  # already printed above
        start_state = row['start_state']
        end_state = row['end_state']
        transition = row['transition'] if 'transition' in row and pd.notnull(row['transition']) else ''
        caspt2_root = row['caspt2_root'] if 'caspt2_root' in row and pd.notnull(row['caspt2_root']) else 'N/A'
        ms_caspt2_root = row['ms_caspt2_root'] if 'ms_caspt2_root' in row and pd.notnull(row['ms_caspt2_root']) else 'N/A'
        caspt2_abs_ev = row['caspt2_energy'] * 27.211407953
        ms_caspt2_abs_ev = row['ms_caspt2_energy'] * 27.211407953
        diff_abs_ev = ms_caspt2_abs_ev - caspt2_abs_ev
        caspt2_abs_h = row['caspt2_energy']
        ms_caspt2_abs_h = row['ms_caspt2_energy']
        diff_abs_h = ms_caspt2_abs_h - caspt2_abs_h
        def sig_str(val):
            if val == 0:
                return f"{0:.10e}"
            return f"{val:.10e}" if abs(val) < 1e-4 or abs(val) > 1e+6 else f"{val:.10g}"
        print(f"{int(start_state):<7} {int(end_state):<7} {transition:<25} {caspt2_root:<12} {ms_caspt2_root:<15} {sig_str(caspt2_abs_ev):<25} {sig_str(ms_caspt2_abs_ev):<28} {sig_str(diff_abs_ev):<22}")
        print(f"{'':<7} {'':<7} {'':<25} {'':<12} {'':<15} {'(Hartree)':<25} {'(Hartree)':<28} {'(Hartree)':<22}")
        print(f"{'':<7} {'':<7} {'':<25} {'':<12} {'':<15} {sig_str(caspt2_abs_h):<25} {sig_str(ms_caspt2_abs_h):<28} {sig_str(diff_abs_h):<22}")
    print("\n")

def pretty_print_geometry_table(csv_file, phi, psi, geometry_name):
    """Pretty print a summary table for a given geometry from merged_results.csv."""
    try:
        from tabulate import tabulate
        use_tabulate = True
    except ImportError:
        use_tabulate = False
    df = pd.read_csv(csv_file)
    geometry_data = df[(df['phi'] == phi) & (df['psi'] == psi)].copy()
    if geometry_data.empty:
        print(f"No data found for geometry phi={phi}, psi={psi}")
        return None
    print(f"\n{'='*80}")
    print(f"PRETTY PRINTED TABLE: {geometry_name} (φ={phi}°, ψ={psi}°)")
    print(f"{'='*80}")
    # Prepare table data
    table = []
    for _, row in geometry_data.iterrows():
        start_state = int(row['start_state'])
        end_state = int(row['end_state'])
        transition = row['transition'] if 'transition' in row and pd.notnull(row['transition']) else ''
        caspt2_root = row['caspt2_root'] if 'caspt2_root' in row and pd.notnull(row['caspt2_root']) else 'N/A'
        ms_caspt2_root = row['ms_caspt2_root'] if 'ms_caspt2_root' in row and pd.notnull(row['ms_caspt2_root']) else 'N/A'
        caspt2_ev = row['caspt2_energy'] * 27.211407953
        ms_caspt2_ev = row['ms_caspt2_energy'] * 27.211407953
        diff_ev = ms_caspt2_ev - caspt2_ev
        caspt2_h = row['caspt2_energy']
        ms_caspt2_h = row['ms_caspt2_energy']
        diff_h = ms_caspt2_h - caspt2_h
        config = row['config'] if 'config' in row and pd.notnull(row['config']) else ''
        coeff = f"{row['coeff']:.3f}" if 'coeff' in row and pd.notnull(row['coeff']) else ''
        weight = f"{row['weight']:.3f}" if 'weight' in row and pd.notnull(row['weight']) else ''
        table.append([
            start_state, end_state, transition, caspt2_root, ms_caspt2_root,
            f"{caspt2_ev:.6f}", f"{ms_caspt2_ev:.6f}", f"{diff_ev:+.6f}",
            f"{caspt2_h:.8f}", f"{ms_caspt2_h:.8f}", f"{diff_h:+.8f}",
            config, coeff, weight
        ])
    headers = [
        "Start", "End", "Transition", "CASPT2 Root", "MS-CASPT2 Root",
        "CASPT2 (eV)", "MS-CASPT2 (eV)", "Δ (eV)",
        "CASPT2 (Ha)", "MS-CASPT2 (Ha)", "Δ (Ha)",
        "Configuration", "Coeff", "Weight"
    ]
    if use_tabulate:
        print(tabulate(table, headers=headers, tablefmt="fancy_grid", showindex=False))
    else:
        pd_table = pd.DataFrame(table, columns=headers)
        print(pd_table.to_string(index=False))
    print("\n")

def create_comparison_summary(geometry1=None, geometry2=None):
    """Create a comprehensive comparison summary between the two geometries
    
    Parameters:
    -----------
    geometry1 : tuple, optional
        (phi, psi, name) for first geometry. Default: (210, 135, "Beta-Strand")
    geometry2 : tuple, optional
        (phi, psi, name) for second geometry. Default: (240, 30, "Alpha-Helix")
    """
    
    # Set default geometries if not provided
    if geometry1 is None:
        geometry1 = (210, 135, "Beta-Strand")
    if geometry2 is None:
        geometry2 = (240, 30, "Alpha-Helix")
    
    # Analyze both geometries
    phi1, psi1, name1 = geometry1
    phi2, psi2, name2 = geometry2
    
    print(f"Analyzing {name1} (φ={phi1}°, ψ={psi1}°) and {name2} (φ={phi2}°, ψ={psi2}°)")
    
    geometry1_data = analyze_geometry_differences('merged_results.csv', phi1, psi1, name1)
    geometry2_data = analyze_geometry_differences('merged_results.csv', phi2, psi2, name2)
    
    if geometry1_data is None or geometry2_data is None:
        return
    
    print(f"\n{'='*80}")
    print("CROSS-GEOMETRY COMPARISON SUMMARY")
    print(f"{'='*80}")
    
    # Compare energy differences between geometries
    print(f"\nENERGY DIFFERENCE COMPARISON (MS-CASPT2 - CASPT2):")
    print(f"{'Transition':<20} {'Beta-Strand (eV)':<15} {'Alpha-Helix (eV)':<15} {'Difference':<12}")
    print("-" * 70)
    
    transitions = sorted(geometry1_data['transition'].unique(), key=lambda x: ('n_to_pi' in x, x))
    
    for transition in transitions:
        geometry1_trans = geometry1_data[geometry1_data['transition'] == transition]
        geometry2_trans = geometry2_data[geometry2_data['transition'] == transition]
        
        if not geometry1_trans.empty and not geometry2_trans.empty:
            geometry1_diff = geometry1_trans.iloc[0]['ms_caspt2_energy_difference_ev'] - geometry1_trans.iloc[0]['caspt2_energy_difference_ev']
            geometry2_diff = geometry2_trans.iloc[0]['ms_caspt2_energy_difference_ev'] - geometry2_trans.iloc[0]['caspt2_energy_difference_ev']
            geometry_diff = geometry2_diff - geometry1_diff
            
            print(f"{transition:<20} {geometry1_diff:+<15.3f} {geometry2_diff:+<15.3f} {geometry_diff:+<12.3f}")
    
    # Overall conclusions
    print(f"\nOVERALL CONCLUSIONS:")
    
    # Calculate average differences for each geometry
    geometry1_diffs = []
    geometry2_diffs = []
    
    for transition in transitions:
        geometry1_trans = geometry1_data[geometry1_data['transition'] == transition]
        geometry2_trans = geometry2_data[geometry2_data['transition'] == transition]
        
        if not geometry1_trans.empty and not geometry2_trans.empty:
            geometry1_diff = geometry1_trans.iloc[0]['ms_caspt2_energy_difference_ev'] - geometry1_trans.iloc[0]['caspt2_energy_difference_ev']
            geometry2_diff = geometry2_trans.iloc[0]['ms_caspt2_energy_difference_ev'] - geometry2_trans.iloc[0]['caspt2_energy_difference_ev']
            
            geometry1_diffs.append(geometry1_diff)
            geometry2_diffs.append(geometry2_diff)
    
    geometry1_avg = np.mean(geometry1_diffs)
    geometry2_avg = np.mean(geometry2_diffs)
    
    print(f"   • {name1} average MS-CASPT2 vs CASPT2 difference: {geometry1_avg:+.3f} eV")
    print(f"   • {name2} average MS-CASPT2 vs CASPT2 difference: {geometry2_avg:+.3f} eV")
    print(f"   • Geometry effect on MS-CASPT2 vs CASPT2 difference: {geometry2_avg - geometry1_avg:+.3f} eV")
    
    # Identify which geometry shows larger differences
    if abs(geometry1_avg) > abs(geometry2_avg):
        print(f"   • {name1} geometry shows larger MS-CASPT2 vs CASPT2 differences")
    else:
        print(f"   • {name2} geometry shows larger MS-CASPT2 vs CASPT2 differences")
    
    # Analyze which transitions are most sensitive to geometry
    geometry_sensitivities = []
    for transition in transitions:
        geometry1_trans = geometry1_data[geometry1_data['transition'] == transition]
        geometry2_trans = geometry2_data[geometry2_data['transition'] == transition]
        
        if not geometry1_trans.empty and not geometry2_trans.empty:
            geometry1_e = geometry1_trans.iloc[0]['ms_caspt2_energy_difference_ev']
            geometry2_e = geometry2_trans.iloc[0]['ms_caspt2_energy_difference_ev']
            sensitivity = abs(geometry2_e - geometry1_e)
            geometry_sensitivities.append((transition, sensitivity))
    
    geometry_sensitivities.sort(key=lambda x: x[1], reverse=True)
    
    print(f"\n   • Most geometry-sensitive transition: {geometry_sensitivities[0][0]} ({geometry_sensitivities[0][1]:.3f} eV)")
    print(f"   • Least geometry-sensitive transition: {geometry_sensitivities[-1][0]} ({geometry_sensitivities[-1][1]:.3f} eV)")

if __name__ == "__main__":    
    # Example usage with custom geometries
    geometry1 = (225, 120, "Beta-Strand")
    # geometry2 = (240, 30, "Alpha-Helix")
    geometry1 = (300, 300, "Alpha-Helix")

    # Display absolute energies tables for both geometries
    display_absolute_energies_table('merged_results.csv', geometry1[0], geometry1[1], geometry1[2])
    # display_absolute_energies_table('merged_results.csv', geometry2[0], geometry2[1], geometry2[2])

    pretty_print_geometry_table('merged_results.csv', geometry1[0], geometry1[1], geometry1[2])
    # pretty_print_geometry_table('merged_results.csv', geometry2[0], geometry2[1], geometry2[2])
    
    # create_comparison_summary(
    #     geometry1=geometry1,
    #     geometry2=geometry2
    # )
    