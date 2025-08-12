import os
import re
from typing import Dict, List, Tuple, Optional
import pandas as pd
import numpy as np

# Constants for energy conversion
HARTREE_TO_EV = 27.211386245988  # 1 Hartree = 27.211386245988 eV


def hartree_to_ev(hartree_value: float) -> float:
    """
    Convert energy from Hartree to electron volts.
    
    Args:
        hartree_value (float): Energy in Hartree
        
    Returns:
        float: Energy in electron volts
    """
    return hartree_value * HARTREE_TO_EV


def extract_roots(log_content: str) -> Dict[int, Dict]:
    """
    Extract root information (energy, weights, configurations) from RASSCF log file.
    Only processes the section between wave function printout markers.
    
    Args:
        log_content (str): Content of the RASSCF log file
        
    Returns:
        Dict[int, Dict]: Dictionary containing information for each root
    """
    # Define the section markers
    start_marker = "Wave function printout:"
    end_marker = "Natural orbitals and occupation numbers for root"
    
    # Extract the relevant section
    section_pattern = f"{start_marker}.*?{end_marker}"
    section_match = re.search(section_pattern, log_content, re.DOTALL)
    if not section_match:
        return {}
    
    section_content = section_match.group(0)
    root_pattern = re.compile(r'printout of CI-coefficients larger than\s+0\.05 for root\s+(\d+)')
    energy_pattern = re.compile(r'energy=\s+(-?\d+\.\d+)')
    weight_pattern = re.compile(r'^\s+\d+\s+([0-9ud]+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)')

    roots = {}
    current_root = None
    
    # Split content into lines and process
    lines = section_content.split('\n')
    i = 0
    while i < len(lines):
        line = lines[i]
        
        # Check for root marker
        root_match = root_pattern.search(line)
        if root_match:
            current_root = int(root_match.group(1))
            roots[current_root] = {
                'energy': None,
                'configurations': []  # List to store all configurations
            }
            # Look for energy in next line
            if i + 1 < len(lines):
                energy_match = energy_pattern.search(lines[i + 1])
                if energy_match:
                    roots[current_root]['energy'] = float(energy_match.group(1))
                    i += 1  # Skip the energy line
        else:
            # Check for weights
            weight_match = weight_pattern.search(line)
            if weight_match and current_root is not None:
                config = weight_match.group(1)
                coeff = float(weight_match.group(2))
                weight = float(weight_match.group(3))
                roots[current_root]['configurations'].append({
                    'config': config,
                    'coeff': coeff,
                    'weight': weight
                })
        
        i += 1
    return roots


def extract_caspt2_energies(log_content: str) -> Dict[str, Dict[int, float]]:
    """
    Extract CASPT2 and MS-CASPT2 energies for each root from log file.
    
    Args:
        log_content (str): Content of the log file
        
    Returns:
        Dict[str, Dict[int, float]]: Dictionary containing CASPT2 and MS-CASPT2 energies for each root
    """
    # Pattern for CASPT2 energies (must start with CASPT2, not MS-CASPT2)
    caspt2_pattern = re.compile(r'^::\s+CASPT2 Root\s+(\d+)\s+Total energy:\s+([-+]?\d+\.\d+)', re.MULTILINE)
    
    # Pattern for MS-CASPT2 energies (must start with MS-CASPT2)
    ms_caspt2_pattern = re.compile(r'^::\s+MS-CASPT2 Root\s+(\d+)\s+Total energy:\s+([-+]?\d+\.\d+)', re.MULTILINE)
    
    caspt2_energies = {}
    ms_caspt2_energies = {}
    
    # Extract CASPT2 energies
    for match in caspt2_pattern.finditer(log_content):
        root_num = int(match.group(1))
        energy = float(match.group(2))
        caspt2_energies[root_num] = energy
    
    # Extract MS-CASPT2 energies
    for match in ms_caspt2_pattern.finditer(log_content):
        root_num = int(match.group(1))
        energy = float(match.group(2))
        ms_caspt2_energies[root_num] = energy
    
    return {
        'CASPT2': caspt2_energies,
        'MS-CASPT2': ms_caspt2_energies
    }


def extract_caspt2_root_energies(log_content: str, energy_unit: str = 'hartree') -> Dict[int, Dict[str, float]]:
    """
    Extract CASPT2 energies for each root from log file, including reference energy and total energy.
    
    This function parses the 'FINAL CASPT2 RESULT' sections for each root group.
    
    Args:
        log_content (str): Content of the log file
        energy_unit (str): Energy unit for output ('hartree' or 'ev')
        
    Returns:
        Dict[int, Dict[str, float]]: Dictionary containing energy information for each root
                                     Keys: root number, Values: dict with 'reference_energy' and 'total_energy'
    """
    # Pattern to find CASPT2 computation groups
    group_pattern = re.compile(r'\+\+ CASPT2 computation for group\s+(\d+)')
    
    # Pattern to find FINAL CASPT2 RESULT sections
    final_result_pattern = re.compile(
        r'FINAL CASPT2 RESULT:(.*?)(?=\+\+|\Z)',
        re.DOTALL
    )
    
    # Pattern to extract reference energy
    reference_pattern = re.compile(r'Reference energy:\s+([-+]?\d+\.\d+)')
    
    # Pattern to extract total energy
    total_pattern = re.compile(r'Total energy:\s+([-+]?\d+\.\d+)')
    
    root_energies = {}
    
    # Find all CASPT2 computation groups
    group_matches = list(group_pattern.finditer(log_content))
    
    for i, group_match in enumerate(group_matches):
        group_num = int(group_match.group(1))
        
        # Find the corresponding FINAL CASPT2 RESULT section
        if i < len(group_matches) - 1:
            # Extract content between this group and the next
            start_pos = group_match.end()
            end_pos = group_matches[i + 1].start()
            section_content = log_content[start_pos:end_pos]
        else:
            # For the last group, extract to the end of the file
            section_content = log_content[group_match.end():]
        
        # Look for FINAL CASPT2 RESULT in this section
        final_result_match = final_result_pattern.search(section_content)
        if final_result_match:
            result_content = final_result_match.group(1)
            
            # Extract reference energy
            ref_match = reference_pattern.search(result_content)
            reference_energy = None
            if ref_match:
                reference_energy = float(ref_match.group(1))
            
            # Extract total energy
            total_match = total_pattern.search(result_content)
            total_energy = None
            if total_match:
                total_energy = float(total_match.group(1))
            
            # Store energies if both were found
            if reference_energy is not None and total_energy is not None:
                # Convert to eV if requested
                if energy_unit.lower() == 'ev':
                    reference_energy = hartree_to_ev(reference_energy)
                    total_energy = hartree_to_ev(total_energy)
                
                root_energies[group_num] = {
                    'reference_energy': reference_energy,
                    'total_energy': total_energy
                }
    
    return root_energies


def extract_dipole_moments(log_content: str) -> Dict[int, float]:
    """
    Extract dipole moments for each root from RASSCF log file.
    
    Args:
        log_content (str): Content of the RASSCF log file
        
    Returns:
        Dict[int, float]: Dictionary mapping root numbers to their dipole moments
    """
    debye_pattern = re.compile(r'Expectation values of various properties for root number:\s+(\d+)')
    dipole_moments = {}
    current_root = None
    
    for line in log_content.split('\n'):
        debye_match = debye_pattern.search(line)
        if debye_match:
            current_root = int(debye_match.group(1))
        elif current_root is not None:
            debye_match = re.search(r'Total=\s*([-+]?\d+\.\d+E[-+]\d+)', line)
            if debye_match:
                dipole_moments[current_root] = float(debye_match.group(1))
    
    return dipole_moments


def extract_geometry(filename: str) -> Tuple[Optional[int], Optional[int]]:
    """
    Extract phi and psi angles from filename.
    
    Args:
        filename (str): Name of the file containing geometry information
        
    Returns:
        Tuple[Optional[int], Optional[int]]: Phi and psi angles if found, None otherwise
    """
    match = re.search(r'(\d+)_(\d+)', filename)
    if match:
        return int(match.group(1)), int(match.group(2))
    return None, None


def get_section_between_start_and_end(string: str, start_string: str, end_string: str) -> str:
    """
    Extract text between two markers in a string.
    
    Args:
        string (str): Input string to search in
        start_string (str): Starting marker
        end_string (str): Ending marker
        
    Returns:
        str: Extracted text between markers
    """
    block_match = re.search(f"{start_string}(.*?){end_string}", string, re.DOTALL | re.MULTILINE)
    return block_match.group(1) if block_match else ""


def get_dataframe_from_transition_string(
    block: str,
    block_headings: List[str],
    data_pattern: str
) -> pd.DataFrame:
    """
    Convert transition string data to DataFrame.
    
    Args:
        block (str): String containing transition data
        block_headings (List[str]): Column names for the DataFrame
        data_pattern (str): Regex pattern to match data
        
    Returns:
        pd.DataFrame: DataFrame containing transition data
    """
    matches = re.findall(data_pattern, block, re.MULTILINE)
    results = [dict(zip(block_headings, [float(i) for i in match])) for match in matches]
    return pd.DataFrame(results)


def get_transition_data(
    log_content: str,
    start_marker: str,
    end_marker: str,
    block_headings: List[str],
    data_pattern: str = r"^\s*(\d+)\s+(\d+)\s+([-+]?\d+\.\d+E[-+]?\d+)\s+([-+]?\d+\.\d+E[-+]?\d+)\s+([-+]?\d+\.\d+E[-+]?\d+)\s+([-+]?\d+\.\d+E[-+]?\d+)"
) -> pd.DataFrame:
    """
    Extract transition data from log content.
    
    Args:
        log_content (str): Content of the log file
        start_marker (str): Starting marker for transition data
        end_marker (str): Ending marker for transition data
        block_headings (List[str]): Column names for the DataFrame
        data_pattern (str): Regex pattern to match data
        
    Returns:
        pd.DataFrame: DataFrame containing transition data
    """
    block = get_section_between_start_and_end(log_content, start_marker, end_marker)
    transitions_dataframe = get_dataframe_from_transition_string(block, block_headings, data_pattern)
    return transitions_dataframe


def convert_au_columns_to_debye(df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert atomic units to Debye for dipole moment columns.
    
    Args:
        df (pd.DataFrame): Input DataFrame
        
    Returns:
        pd.DataFrame: DataFrame with converted units
    """
    au_debye_conversion_factor = 0.393456
    for col in df.columns:
        if col.endswith('au'):
            df[col.replace('au', 'Debyes')] = df[col] / au_debye_conversion_factor
            df = df.drop(col, axis=1)
    return df


def process_single_log_file(filepath: str) -> pd.DataFrame:
    """
    Process a single RASSCF log file and return the analysis results.
    
    Args:
        filepath (str): Path to the log file
        
    Returns:
        pd.DataFrame: DataFrame containing analysis results
        
    Raises:
        FileNotFoundError: If the file doesn't exist
        ValueError: If geometry information cannot be extracted
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
        
    filename = os.path.basename(filepath)
    phi, psi = extract_geometry(filename)
    
    if phi is None or psi is None:
        raise ValueError(f"Could not extract geometry information from filename: {filename}")
    
    with open(filepath, 'r') as file:
        log_content = file.read()
    
    # Extract data
    roots = extract_roots(log_content)
    dipole_moments = extract_dipole_moments(log_content)
    caspt2_energies = extract_caspt2_energies(log_content)
    transition_dipole_moment = get_transition_data(
        log_content,
        start_marker,
        end_marker,
        block_headings
    )
    
    # Combine data
    data = []
    for root, root_data in roots.items():
        for config_data in root_data['configurations']:
            data.append({
                'Phi': phi,
                'Psi': psi,
                'Root': root,
                'Energy': root_data['energy'],
                'CASPT2_Energy': caspt2_energies['CASPT2'].get(root, None),
                'MS_CASPT2_Energy': caspt2_energies['MS-CASPT2'].get(root, None),
                'Config': config_data['config'],
                'Coeff': config_data['coeff'],
                'Weight': config_data['weight'],
                'Dipole Moment': dipole_moments.get(root, 0)  # Use 0 if no dipole moment found
            })
    
    df = pd.DataFrame(data)

    df = (
        df
        .rename(columns={'Root': 'end_state'})
        .assign(start_state=1)
    )
    df = (
        df
        .merge(
            transition_dipole_moment,
            how='left'
        )
    )
    
    # Calculate energy differences for CASPT2 and MS-CASPT2
    # Group by phi and psi coordinates
    if 'CASPT2_Energy' in df.columns:
        caspt2_diff_hartree = df.groupby(['Phi', 'Psi'])['CASPT2_Energy'].transform(
            lambda x: x - x.min()
        )
        # Convert from Hartree to eV (1 Hartree = 27.2114 eV)
        df['CASPT2_Energy_Difference_eV'] = caspt2_diff_hartree * 27.2114
    
    if 'MS_CASPT2_Energy' in df.columns:
        ms_caspt2_diff_hartree = df.groupby(['Phi', 'Psi'])['MS_CASPT2_Energy'].transform(
            lambda x: x - x.min()
        )
        # Convert from Hartree to eV (1 Hartree = 27.2114 eV)
        df['MS_CASPT2_Energy_Difference_eV'] = ms_caspt2_diff_hartree * 27.2114
    
    return convert_au_columns_to_debye(df)


# Constants
start_marker = r"\+\+ Dipole transition vectors \(spin-free states\):"
end_marker = r"^\s*--\s*$"
block_headings = [
    'start_state',
    'end_state',
    'Dx_au',
    'Dy_au',
    'Dz_au',
    'TotalD_au'
]

def filter_big_weights(df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter out configurations with weights greater than 0.05.
    """
    return df[df['Weight'] > 0.30]


def main():
    """Main function to process all log files in the directory."""
    # Get the directory where the script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Use relative path to log files directory
    log_dir = os.path.join(script_dir, "log_files")
    
    if not os.path.exists(log_dir):
        raise FileNotFoundError(f"Log directory not found: {log_dir}")
        
    results = []
    
    for filename in os.listdir(log_dir):
        if filename.endswith(".log"):
            filepath = os.path.join(log_dir, filename)
            try:
                df = process_single_log_file(filepath)
                results.append(df)
            except (FileNotFoundError, ValueError) as e:
                print(f"Error processing {filename}: {str(e)}")
    
    if results:
        final_df = pd.concat(results, ignore_index=True)
        output_file = os.path.join(script_dir, "transition_analysis_results.csv")
        final_df = filter_big_weights(final_df)
        final_df = final_df.sort_values(by=['Phi', 'Psi', 'start_state', 'end_state'])
        final_df.to_csv(output_file, index=False)

        print(f"Results saved to: {output_file}")
    else:
        print("No results were generated.")


if __name__ == "__main__":
    main()
