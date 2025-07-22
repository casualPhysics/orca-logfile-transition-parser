import pandas as pd

def merge_orbital_assignments_and_analysis(orbital_assignments_path, transition_analysis_path):
    orbital_assignments = pd.read_csv(orbital_assignments_path)
    transition_analysis = pd.read_csv(transition_analysis_path)

    # just change everycolumn to lower case here 
    orbital_assignments.columns = orbital_assignments.columns.str.lower()
    transition_analysis.columns = transition_analysis.columns.str.lower()

    # Merge on phi/Phi, psi/Psi, config/Config
    merged = transition_analysis.merge(
        orbital_assignments,
        left_on=["phi", "psi", "config"],
        right_on=["phi", "psi", "config"],
        how="left"
    )

    return merged


def calculate_energy_difference(merged_df):
    # Group by phi and psi coordinates
    energy_diff_hartree = merged_df.groupby(['phi', 'psi'])['energy'].transform(
        lambda x: x - x.min()
    )
    # Convert from Hartree to eV (1 Hartree = 27.2114 eV)
    merged_df['energy_difference_eV'] = energy_diff_hartree * 27.2114
    return merged_df


def filter_transitions(merged_df, desired_transitions = ['n_to_pi_*_L', 'n_to_pi_*_R', 'pi_nb_to_pi_*_R', 'pi_nb_to_pi_*_L']):
    # filter the transitions that have a energy difference of less than 0.01 eV
    merged_df = merged_df[merged_df['transition'].isin(desired_transitions)]
    return merged_df

def main():
    orbital_assignments_path = "/Users/afiqhatta/orbital_labelling/ras12i8_transitions.csv"
    transition_analysis_path = "transition_analysis_results.csv"
    merged_df = merge_orbital_assignments_and_analysis(
        orbital_assignments_path,
        transition_analysis_path
    )
    
    # Calculate energy differences grouped by phi/psi coordinates
    merged_df = calculate_energy_difference(merged_df)
    merged_df = filter_transitions(merged_df)
    merged_df.to_csv("merged_results.csv", index=False)

if __name__ == "__main__":
    main()
