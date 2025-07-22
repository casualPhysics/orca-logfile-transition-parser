# RASSCF Log File Analysis Tool

This Python tool analyzes RASSCF (Restricted Active Space Self-Consistent Field) log files to extract important information about electronic states, configurations, and transition properties. It processes multiple log files and combines the results into a single CSV file for further analysis.

## Features

- Extracts root information (energy, configurations, weights) from RASSCF log files
- Processes dipole moments and transition dipole moments
- Handles multiple log files in a directory
- Combines results into a structured CSV format
- Supports geometry information extraction from filenames
- Merges orbital assignments with transition analysis results
- Includes ground state data for complete energy analysis
- Calculates energy differences in electron volts (eV)

## Requirements

- Python 3.6+
- pandas
- numpy

## Installation

1. Clone the repository:
```bash
git clone https://github.com/casualPhysics/orca-logfile-transition-parser.git
cd rasscf-analysis
```

2. Install required packages:
```bash
pip install -r requirements.txt
```

## Usage

### Step 1: Parse RASSCF Log Files
1. Place your RASSCF log files in the `log_files` directory. The filenames should contain geometry information in the format `phi_psi.log` (e.g., `90_180.log`).

2. Run the log file parsing script:
```bash
python parse_rasscf_logfiles.py
```

3. The results will be saved in `transition_analysis_results.csv` in the project root directory.

### Step 2: Merge with Orbital Assignments and Add Ground State
4. **Merge orbital assignments with transition analysis results** to create `merged_results.csv`:
```bash
python construct_transition_data.py
```

This script will:
- Merge transition analysis results with orbital assignments from an external CSV file
- Calculate energy differences in eV (converted from Hartree)
- Filter for specific transitions (n_to_pi_*_L, n_to_pi_*_R, pi_nb_to_pi_*_R, pi_nb_to_pi_*_L)
- **Add ground state rows** (end_state = 1) for each geometry with transition marked as "GROUND"
- Sort results by phi, psi, config, and end_state

**Note**: You'll need to update the `orbital_assignments_path` variable in the script to point to your orbital assignments CSV file.

### Step 3: Analysis and Visualization
5. Run detailed comparison analysis:
```bash
python detailed_comparison_analysis.py
```

6. Generate energy continuity plots:
```bash
python plot_energy_continuity.py
```

## Output Format

The final `merged_results.csv` file contains the following columns:
- `phi`, `psi`: Dihedral angles from filename
- `end_state`: Root number (1 = ground state)
- `start_state`: Initial state (always 1)
- `energy`: Energy of the state in Hartree
- `config`: Configuration string
- `coeff`: CI coefficient
- `weight`: Weight of the configuration
- `dipole_moment`: Total dipole moment
- `dx_debyes`, `dy_debyes`, `dz_debyes`: Transition dipole moments in Debye units
- `totald_debyes`: Total transition dipole moment in Debye units
- `transition`: Transition type (e.g., "n_to_pi_*_L", "GROUND")
- `energy_difference_eV`: Energy difference from ground state in eV

### Example Output

Here's an example of the output data (truncated for clarity):

| phi | psi | end_state | energy | config | coeff | weight | dipole_moment | dx_debyes | dy_debyes | dz_debyes | totald_debyes | transition | energy_difference_eV |
|-----|-----|-----------|---------|---------|--------|---------|---------------|------------|------------|------------|---------------|------------|-------------------|
| 240 | 90 | 1 | -453.959 | 22222200 | 0.963 | 0.927 | 1.429 | - | - | - | - | GROUND | 0.000 |
| 240 | 90 | 2 | -453.743 | 222u220d | 0.947 | 0.898 | 2.516 | -0.065 | -0.068 | 0.173 | 0.197 | n_to_pi_*_R | 5.959 |
| 240 | 90 | 3 | -453.741 | 22u222d0 | -0.894 | 0.800 | 3.013 | -0.244 | -0.208 | -0.036 | 0.323 | n_to_pi_*_L | 6.286 |
| 240 | 90 | 4 | -453.645 | 22222ud0 | -0.830 | 0.689 | 3.273 | -4.238 | 1.524 | -0.824 | 4.578 | pi_nb_to_pi_*_L | 8.327 |
| ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |

## Script Details

### construct_transition_data.py
This script performs the following operations:

1. **Merges orbital assignments** with transition analysis results
2. **Calculates energy differences** in eV (converted from Hartree using 27.2114 eV/Hartree)
3. **Filters transitions** to include only:
   - `n_to_pi_*_L`
   - `n_to_pi_*_R` 
   - `pi_nb_to_pi_*_R`
   - `pi_nb_to_pi_*_L`
4. **Adds ground state data** for each unique (phi, psi, config) geometry
5. **Marks ground states** with transition = "GROUND"
6. **Sorts results** by geometry and state

**Configuration**: Update the `orbital_assignments_path` variable to point to your orbital assignments CSV file.

## Example

For a log file named `90_180.log` containing:
```
Wave function printout:
printout of CI-coefficients larger than 0.05 for root 1
energy= -1234.5678
    1  1u2u3u  -0.1234   0.0152
    2  2u3u4u   0.5678   0.3224
```

The output will include all configurations with their respective coefficients and weights, along with any dipole moment information found in the file, plus ground state data for complete energy analysis.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details. 