# RASSCF Log File Analysis Tool

This Python tool analyzes RASSCF (Restricted Active Space Self-Consistent Field) log files to extract important information about electronic states, configurations, and transition properties. It processes multiple log files and combines the results into a single CSV file for further analysis.

## Features

- Extracts root information (energy, configurations, weights) from RASSCF log files
- Processes dipole moments and transition dipole moments
- Handles multiple log files in a directory
- Combines results into a structured CSV format
- Supports geometry information extraction from filenames

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

1. Place your RASSCF log files in the `log_files` directory. The filenames should contain geometry information in the format `phi_psi.log` (e.g., `90_180.log`).

2. Run the analysis script:
```bash
python analyse_rasscf_logfile.py
```

3. The results will be saved in `transition_analysis_results.csv` in the project root directory.

## Output Format

The output CSV file contains the following columns:
- `Phi`: Phi angle from filename
- `Psi`: Psi angle from filename
- `end_state`: Root number
- `start_state`: Initial state (always 1)
- `Energy`: Energy of the state
- `Config`: Configuration string
- `Coeff`: CI coefficient
- `Weight`: Weight of the configuration
- `Dipole Moment`: Total dipole moment
- `Dx_Debyes`, `Dy_Debyes`, `Dz_Debyes`: Transition dipole moments in Debye units
- `TotalD_Debyes`: Total transition dipole moment in Debye units

### Example Output

Here's an example of the output data (truncated for clarity):

| Phi | Psi | end_state | Energy | Config | Coeff | Weight | Dipole Moment | Dx_Debyes | Dy_Debyes | Dz_Debyes | TotalD_Debyes |
|-----|-----|-----------|---------|---------|--------|---------|---------------|------------|------------|------------|---------------|
| 240 | 90 | 1 | -453.959 | 22222200 | 0.963 | 0.927 | 1.429 | - | - | - | - |
| 240 | 90 | 2 | -453.743 | 222u220d | 0.947 | 0.898 | 2.516 | -0.065 | -0.068 | 0.173 | 0.197 |
| 240 | 90 | 3 | -453.741 | 22u222d0 | -0.894 | 0.800 | 3.013 | -0.244 | -0.208 | -0.036 | 0.323 |
| 240 | 90 | 4 | -453.645 | 22222ud0 | -0.830 | 0.689 | 3.273 | -4.238 | 1.524 | -0.824 | 4.578 |
| ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |

## Example

For a log file named `90_180.log` containing:
```
Wave function printout:
printout of CI-coefficients larger than 0.05 for root 1
energy= -1234.5678
    1  1u2u3u  -0.1234   0.0152
    2  2u3u4u   0.5678   0.3224
```

The output will include all configurations with their respective coefficients and weights, along with any dipole moment information found in the file.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details. 