"""
Energy Continuity Analysis Script
================================

This script creates various plots for analyzing energy continuity across different transitions.
All plots are automatically saved to the output folder without displaying them on screen.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend to prevent display
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import os

def create_output_folder(folder_name='output_plots'):
    """Create output folder for saving plots"""
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
        print(f"Created output folder: {folder_name}")
    return folder_name

def load_and_prepare_data(csv_file):
    """Load the CSV data and prepare it for plotting"""
    df = pd.read_csv(csv_file)
    
    # Define the four value columns we want to analyze
    value_columns = ['totald_debyes', 'caspt2_energy_difference_ev', 'ms_caspt2_energy_difference_ev', 'energy_difference_eV']
    
    # Verify the specified columns exist
    missing_columns = [col for col in value_columns if col not in df.columns]
    if missing_columns:
        available_columns = list(df.columns)
        raise ValueError(f"Columns {missing_columns} not found in the CSV file. Available columns: {available_columns}")
    
    # Get unique transitions
    transitions = df['transition'].unique()
    print(f"Available transitions: {transitions}")
    print(f"Will analyze the following columns: {value_columns}")
    
    return df, transitions, value_columns

def transform_phi_psi_coordinates(df):
    """Transform phi-psi coordinates from [0, 360] to [-180, 180] range"""
    df_transformed = df.copy()
    
    # Transform phi coordinates: [0, 360] -> [-180, 180]
    # Map 0-180 to 0 to 180, and 180-360 to -180 to 0
    df_transformed['phi_transformed'] = df_transformed['phi'].apply(
        lambda x: x if x <= 180 else x - 360
    )
    
    # Transform psi coordinates: [0, 360] -> [-180, 180]
    # Map 0-180 to 0 to 180, and 180-360 to -180 to 0
    df_transformed['psi_transformed'] = df_transformed['psi'].apply(
        lambda x: x if x <= 180 else x - 360
    )
    
    return df_transformed

def convert_energy_to_ev(df, energy_columns):
    """Convert energy columns from Hartree to electron volts (eV)"""
    df_converted = df.copy()
    
    # Conversion factor: 1 Hartree = 27.2114 eV
    hartree_to_ev = 27.2114
    
    for col in energy_columns:
        if col in df_converted.columns:
            df_converted[col] = df_converted[col] * hartree_to_ev
    
    return df_converted

def classify_conformation(phi, psi):
    """
    Classify protein conformation based on phi-psi angles.
    
    Args:
        phi (float): Phi angle in degrees
        psi (float): Psi angle in degrees
        
    Returns:
        str: Conformation type ('beta_sheet', 'helix', 'alpha_helix', 'left_handed_helix', 'polyprolineII', or 'other')
    """
    # Transform phi to [-180, 180] range if needed
    if phi > 180:
        phi = phi - 360
    
    # Transform psi to [-180, 180] range if needed  
    if psi > 180:
        psi = psi - 360
    
    # Beta-sheet regions
    if -150 <= phi <= -120 and 90 <= psi <= 180:
        return 'beta_sheet'
    
    # Helix regions
    if -90 <= phi <= -45 and -90 <= psi <= 0:
        return 'helix'
    
    # Alpha-helix regions
    if -120 <= phi <= -90 and -45 <= psi <= 30:
        return 'alpha_helix'
    
    # Left-handed helix regions
    if 45 <= phi <= 75 and 0 <= psi <= 60:
        return 'left_handed_helix'
    
    # Polyproline II regions
    if -105 <= phi <= -45 and 115 <= psi <= 175:
        return 'polyprolineII'
    
    return 'other'

def get_conformation_colors():
    """Get color scheme for different conformations"""
    return {
        'beta_sheet': 'blue',           # Blue circles for beta-sheet conformations
        'helix': 'red',                 # Red circles for helix conformations
        'alpha_helix': 'darkred',       # Dark red circles for alpha-helix conformations
        'left_handed_helix': 'orange',  # Orange circles for left-handed helix conformations
        'polyprolineII': 'purple',      # Purple circles for polyproline II conformations
        'other': 'gray'                 # Gray circles for other conformations
    }

def create_3d_plot(df, transition_type, value_column='totald_debyes', output_folder='output_plots', 
                   save_plot=True, filter_spatial_outliers=False, outlier_threshold_percentile=95):
    """Create a 3D plot for a specific transition type"""
    
    # Filter data for the specific transition
    transition_data = df[df['transition'] == transition_type].copy()
    
    if len(transition_data) == 0:
        print(f"No data found for transition: {transition_type}")
        return
    
    # Transform coordinates to [-180, 180] range
    transition_data = transform_phi_psi_coordinates(transition_data)
    
    # Filter out spatial outliers if requested
    if filter_spatial_outliers:
        # Get geometries to exclude (where dipole > 1)
        high_dipole_geometries = df[
            (df['transition'].isin(['n_to_pi_*_R', 'n_to_pi_*_L'])) & 
            (df['totald_debyes'] > 1)
        ]['geometry'].unique()
        
        # Filter out these geometries from the current transition data
        if len(high_dipole_geometries) > 0:
            transition_data = transition_data[~transition_data['geometry'].isin(high_dipole_geometries)]
            print(f"Filtered out {len(high_dipole_geometries)} geometries with dipole > 1 for {transition_type}")
        
        # Also filter from the main dataframe for this transition to ensure consistency
        df_filtered = df[df['transition'] == transition_type].copy()
        df_filtered = df_filtered[~df_filtered['geometry'].isin(high_dipole_geometries)]
        
        # Update the dipole_data to use filtered data
        dipole_data = df_filtered['totald_debyes']
    else:
        # For other transitions, use original dipole data
        dipole_data = df[df['transition'] == transition_type]['totald_debyes']
    
        if len(dipole_data) > 0:
            # For non-dipole plots, also apply percentile-based filtering if needed
            if value_column != 'totald_debyes':
                # Define upper and lower bounds for outlier filtering
                upper_threshold = np.percentile(dipole_data, outlier_threshold_percentile)
                lower_threshold = np.percentile(dipole_data, 100 - outlier_threshold_percentile)
                
                # Identify outliers (both high and low values)
                upper_outlier_mask = df_filtered['totald_debyes'] >= upper_threshold
                lower_outlier_mask = df_filtered['totald_debyes'] <= lower_threshold
                outlier_mask = upper_outlier_mask | lower_outlier_mask
                
                outlier_indices = df_filtered[outlier_mask].index
                
                # Remove outlier points from transition_data
                transition_data = transition_data[~transition_data.index.isin(outlier_indices)]
                print(f"Filtered out {len(outlier_indices)} spatial outliers based on dipole moment thresholds:")
                print(f"  Upper threshold (95th percentile): {upper_threshold:.3f} Debyes")
                print(f"  Lower threshold (5th percentile): {lower_threshold:.3f} Debyes")
    
    # Create a pivot table for the 3D plot
    # We'll use the first value for each phi-psi combination
    pivot_data = transition_data.groupby(['phi_transformed', 'psi_transformed'])[value_column].first().reset_index()
    
    # Create meshgrid for 3D plotting
    phi_values = sorted(pivot_data['phi_transformed'].unique())
    psi_values = sorted(pivot_data['psi_transformed'].unique())
    
    # Create meshgrid
    phi_mesh, psi_mesh = np.meshgrid(phi_values, psi_values)
    
    # Create value matrix
    value_matrix = np.full((len(psi_values), len(phi_values)), np.nan)
    
    for i, psi in enumerate(psi_values):
        for j, phi in enumerate(phi_values):
            mask = (pivot_data['phi_transformed'] == phi) & (pivot_data['psi_transformed'] == psi)
            if mask.any():
                value_matrix[i, j] = pivot_data[mask][value_column].iloc[0]
    
    # Create the 3D plot
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot the surface
    surf = ax.plot_surface(phi_mesh, psi_mesh, value_matrix, 
                          cmap='viridis', alpha=0.8, linewidth=0, antialiased=True)
    
    # Classify conformations and color-code the scatter points
    conformation_colors = []
    for _, row in pivot_data.iterrows():
        conf_type = classify_conformation(row['phi_transformed'], row['psi_transformed'])
        conformation_colors.append(get_conformation_colors()[conf_type])
    
    # Add scatter points for actual data points with conformation-based colors
    scatter = ax.scatter(pivot_data['phi_transformed'], pivot_data['psi_transformed'], pivot_data[value_column], 
                       c=conformation_colors, s=20, alpha=0.7)
    
    # Add legend for conformation types
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=8, label='Beta-sheet'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=8, label='Helix'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='darkred', markersize=8, label='Alpha-helix'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='orange', markersize=8, label='Left-handed helix'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='purple', markersize=8, label='Polyproline II'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=8, label='Other')
    ]
    ax.legend(handles=legend_elements, loc='upper right')
    
    # Customize the plot
    ax.set_xlabel('Phi (degrees)', fontsize=12)
    ax.set_ylabel('Psi (degrees)', fontsize=12)
    ax.set_zlabel(f'{value_column}', fontsize=12)
    
    # Format the title based on the value column and outlier filtering
    if value_column == 'totald_debyes':
        title = f'Total Dipole Moment (Debyes) Continuity: {transition_type}'
    elif value_column == 'caspt2_energy_difference_ev':
        title = f'CASPT2 Energy Difference (eV) Continuity: {transition_type}'
    elif value_column == 'ms_caspt2_energy_difference_ev':
        title = f'MS-CASPT2 Energy Difference (eV) Continuity: {transition_type}'
    elif value_column == 'energy_difference_eV':
        title = f'Energy Difference (eV) Continuity: {transition_type}'
    else:
        title = f'{value_column} Continuity: {transition_type}'
    
    # Add outlier filtering info to title if applicable
    if filter_spatial_outliers and value_column != 'totald_debyes':
        title += f' (Spatial Outliers Filtered)'
    
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    # Add colorbar
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5, label=value_column)
    
    # Set axis limits to [-180, 180]
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    
    # Add grid
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_plot:
        # Include outlier filtering info in filename if applicable
        filename_suffix = ""
        if filter_spatial_outliers and value_column != 'totald_debyes':
            filename_suffix = "_outliers_filtered"
        
        filename = os.path.join(output_folder, f"{value_column}_3d_continuity_{transition_type.replace('*', 'star').replace(' ', '_')}{filename_suffix}.png")
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"3D plot saved as: {filename}")
    
    # Close the figure instead of showing it
    plt.close(fig)
    
    return fig, ax


def filter_out_erroneous_geometries(df):
    df = df.copy() 
    # remove high dipole and energy difference from above 
    high_dipole_geometries = df[
        (df['transition'].isin(['n_to_pi_*_R', 'n_to_pi_*_L'])) & 
        ((df['totald_debyes'] > 1) | (df['energy_difference_eV'] > 8))
    ]['geometry'].unique()

    # remove low energy dipoles
    low_energy_geometries = df[
        (df['transition'].isin(['pi_nb_to_pi_*_R'])) & 
        (df['energy_difference_eV'] < 3)
    ]['geometry'].unique()

    low_energy_n_pi_star_R_geometries = df[
        (df['transition'].isin(['n_to_pi_*_R'])) & 
        (df['energy_difference_eV'] < 4.5)
    ]['geometry'].unique()

    excluded_set = list(set(list(high_dipole_geometries) + list(low_energy_geometries) + list(low_energy_n_pi_star_R_geometries)))
    
    # Filter out these geometries from the current transition data
    if len(excluded_set) > 0:
        df = df[~df['geometry'].isin(excluded_set)]
        print(f"Filtered out {len(excluded_set)} geometries with dipole > 1")
        return df 
    
    return 


def categorize_geometries(df):
    df = df.copy()
    conformation_colors = []
    conformation_name = []
    for _, row in df.iterrows():
        conf_type = classify_conformation(row['phi_transformed'], row['psi_transformed'])
        conformation_name.append(conf_type)
    df['conformation'] = conformation_name
    return df 


def create_combined_3d_plots(df, value_column, transitions, output_folder='output_plots', 
                            save_plot=True, filter_spatial_outliers=False, outlier_threshold_percentile=95):
    """Create combined 3D plots with all transitions on one page for a specific value column"""
    
    # Filter out GROUND transition
    transitions = [t for t in transitions if t != 'GROUND']
    
    # Create a 2x2 subplot layout
    fig = plt.figure(figsize=(20, 16))

    df = df.copy()
    df = transform_phi_psi_coordinates(df)
    df['geometry'] = df['phi_transformed'].astype(str) + '_' + df['psi_transformed'].astype(str)
    df = categorize_geometries(df)
    df = df.loc[df['conformation'] != 'other']

    if filter_spatial_outliers:
        df = filter_out_erroneous_geometries(df)

    for i, transition_type in enumerate(transitions):
        # Filter data for the specific transition
        transition_data = df[df['transition'] == transition_type].copy()
        
        if len(transition_data) == 0:
            continue
        
        # Create a pivot table for the 3D plot
        pivot_data = transition_data.groupby(['phi_transformed', 'psi_transformed'])[value_column].first().reset_index()
        
        # Create meshgrid for 3D plotting
        phi_values = sorted(pivot_data['phi_transformed'].unique())
        psi_values = sorted(pivot_data['psi_transformed'].unique())
        
        # Create meshgrid
        phi_mesh, psi_mesh = np.meshgrid(phi_values, psi_values)
        
        # Create value matrix
        value_matrix = np.full((len(psi_values), len(phi_values)), np.nan)
        
        for j, psi in enumerate(psi_values):
            for k, phi in enumerate(phi_values):
                mask = (pivot_data['phi_transformed'] == phi) & (pivot_data['psi_transformed'] == psi)
                if mask.any():
                    value_matrix[j, k] = pivot_data[mask][value_column].iloc[0]
        
        # Create subplot
        ax = fig.add_subplot(2, 2, i+1, projection='3d')
        
        # Plot the surface
        surf = ax.plot_surface(phi_mesh, psi_mesh, value_matrix, 
                              cmap='viridis', alpha=0.8, linewidth=0, antialiased=False)
        
        # Classify conformations and color-code the scatter points
        conformation_colors = []
        for _, row in pivot_data.iterrows():
            conf_type = classify_conformation(row['phi_transformed'], row['psi_transformed'])
            conformation_colors.append(get_conformation_colors()[conf_type])
        
        # Add scatter points for actual data points with conformation-based colors
        scatter = ax.scatter(pivot_data['phi_transformed'], pivot_data['psi_transformed'], pivot_data[value_column], 
                           c=conformation_colors, s=20, alpha=0.7)
        
        # Add legend for conformation types
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=8, label='Beta-sheet'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=8, label='Helix'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='darkred', markersize=8, label='Alpha-helix'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='orange', markersize=8, label='Left-handed helix'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='purple', markersize=8, label='Polyproline II'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=8, label='Other')
        ]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=8)
        
        # Customize the subplot
        ax.set_xlabel('Phi (degrees)', fontsize=10)
        ax.set_ylabel('Psi (degrees)', fontsize=10)
        
        # Format the title based on the value column
        if value_column == 'totald_debyes':
            title = f'Total Dipole Moment (Debyes): {transition_type}'
        elif value_column == 'caspt2_energy_difference_ev':
            title = f'CASPT2 Energy Difference (eV): {transition_type}'
        elif value_column == 'ms_caspt2_energy_difference_ev':
            title = f'MS-CASPT2 Energy Difference (eV): {transition_type}'
        elif value_column == 'energy_difference_eV':
            title = f'Energy Difference (eV): {transition_type}'
        else:
            title = f'{value_column}: {transition_type}'
        
        ax.set_title(title, fontsize=12, fontweight='bold')
        
        # Set axis limits
        ax.set_xlim(-180, 180)
        ax.set_ylim(-180, 180)

        # Add grid
        ax.grid(True, alpha=0.3)
    
    # Add overall title
    if value_column == 'totald_debyes':
        overall_title = 'Total Dipole Moment (Debyes) Continuity - All Transitions'
    elif value_column == 'caspt2_energy_difference_ev':
        overall_title = 'CASPT2 Energy Difference (eV) Continuity - All Transitions'
    elif value_column == 'ms_caspt2_energy_difference_ev':
        overall_title = 'MS-CASPT2 Energy Difference (eV) Continuity - All Transitions'
    elif value_column == 'energy_difference_eV':
        overall_title = 'Energy Difference (eV) Continuity - All Transitions'
    else:
        overall_title = f'{value_column} Continuity - All Transitions'
    
    # Add outlier filtering info to title if applicable
    if filter_spatial_outliers:
        overall_title += ' (Spatial Outliers Filtered)'
    
    fig.suptitle(overall_title, fontsize=16, fontweight='bold')
    
    plt.tight_layout()
    
    if save_plot:
        # Include outlier filtering info in filename if applicable
        filename_suffix = ""
        if filter_spatial_outliers:
            filename_suffix = "_outliers_filtered"
        
        filename = os.path.join(output_folder, f"{value_column}_combined_3d_continuity_all_transitions{filename_suffix}.png")
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Combined 3D plot saved as: {filename}")
    
    # Close the figure
    plt.close(fig)
    
    return fig

def create_heatmap_plot(df, transition_type, value_column='totald_debyes', output_folder='output_plots', save_plot=True, filter_spatial_outliers=False):
    """Create a 2D heatmap for a specific transition type"""
    
    # Filter data for the specific transition
    transition_data = df[df['transition'] == transition_type].copy()
    
    if len(transition_data) == 0:
        print(f"No data found for transition: {transition_type}")
        return
    
    # Transform coordinates to [-90, 90] range
    transition_data = transform_phi_psi_coordinates(transition_data)
    
    # Filter out spatial outliers if requested
    if filter_spatial_outliers:
        # Filter out geometries where totald_debyes > 1 for n->pi* L or n->pi* R transitions
        if transition_type in ['n->pi* L', 'n->pi* R']:
            # Get geometries to exclude (where dipole > 1)
            high_dipole_geometries = df[
                (df['transition'] == transition_type) & 
                (df['totald_debyes'] > 1)
            ]['geometry'].unique()
            
            # Filter out these geometries from the current transition data
            if len(high_dipole_geometries) > 0:
                transition_data = transition_data[~transition_data['geometry'].isin(high_dipole_geometries)]
                print(f"Filtered out {len(high_dipole_geometries)} geometries with dipole > 1 for {transition_type}")
    
    # Create a pivot table for the heatmap
    pivot_data = transition_data.groupby(['phi_transformed', 'psi_transformed'])[value_column].mean().reset_index()
    
    # Create pivot table for heatmap
    heatmap_data = pivot_data.pivot(index='psi_transformed', columns='phi_transformed', values=value_column)
    
    # Create the heatmap
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Format the colorbar label based on the value column
    if value_column == 'totald_debyes':
        cbar_label = 'Total Dipole Moment (Debyes)'
    elif value_column == 'caspt2_energy_difference_ev':
        cbar_label = 'CASPT2 Energy Difference (eV)'
    elif value_column == 'ms_caspt2_energy_difference_ev':
        cbar_label = 'MS-CASPT2 Energy Difference (eV)'
    elif value_column == 'energy_difference_eV':
        cbar_label = 'Energy Difference (eV)'
    else:
        cbar_label = value_column
    
    # Create custom colormap
    sns.heatmap(heatmap_data, cmap='viridis', annot=False, cbar_kws={'label': cbar_label})
    
    # Format the title based on the value column
    if value_column == 'totald_debyes':
        title = f'Total Dipole Moment (Debyes) Heatmap: {transition_type}'
    elif value_column == 'caspt2_energy_difference_ev':
        title = f'CASPT2 Energy Difference (eV) Heatmap: {transition_type}'
    elif value_column == 'ms_caspt2_energy_difference_ev':
        title = f'MS-CASPT2 Energy Difference (eV) Heatmap: {transition_type}'
    elif value_column == 'energy_difference_eV':
        title = f'Energy Difference (eV) Heatmap: {transition_type}'
    else:
        title = f'{value_column} Heatmap: {transition_type}'
    
    # Add outlier filtering info to title if applicable
    if filter_spatial_outliers:
        title += ' (Spatial Outliers Filtered)'
    
    plt.title(title, fontsize=14, fontweight='bold')
    plt.xlabel('Phi (degrees)', fontsize=12)
    plt.ylabel('Psi (degrees)', fontsize=12)
    
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45)
    plt.yticks(rotation=0)
    
    plt.tight_layout()
    
    if save_plot:
        filename = os.path.join(output_folder, f"{value_column}_heatmap_{transition_type.replace('*', 'star').replace(' ', '_')}.png")
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Heatmap saved as: {filename}")
    
    # Close the figure instead of showing it
    plt.close(fig)

def create_distribution_plots(df, transition_type, value_column='totald_debyes', output_folder='output_plots', save_plot=True, filter_spatial_outliers=False):
    """Create distribution plots (histogram and box plot) for a specific transition type"""
    
    # Filter data for the specific transition
    transition_data = df[df['transition'] == transition_type].copy()
    
    if len(transition_data) == 0:
        print(f"No data found for transition: {transition_type}")
        return
    
    # Filter out spatial outliers if requested
    if filter_spatial_outliers:
        # Filter out geometries where totald_debyes > 1 for n->pi* L or n->pi* R transitions
        if transition_type in ['n->pi* L', 'n->pi* R']:
            # Get geometries to exclude (where dipole > 1)
            high_dipole_geometries = df[
                (df['transition'] == transition_type) & 
                (df['totald_debyes'] > 1)
            ]['geometry'].unique()
            
            # Filter out these geometries from the current transition data
            if len(high_dipole_geometries) > 0:
                transition_data = transition_data[~transition_data['geometry'].isin(high_dipole_geometries)]
                print(f"Filtered out {len(high_dipole_geometries)} geometries with dipole > 1 for {transition_type}")
    
    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Histogram
    ax1.hist(transition_data[value_column], bins=30, alpha=0.7, color='skyblue', edgecolor='black')
    # Format the x-axis label based on the value column
    if value_column == 'totald_debyes':
        xlabel = 'Total Dipole Moment (Debyes)'
    elif value_column == 'caspt2_energy_difference_ev':
        xlabel = 'CASPT2 Energy Difference (eV)'
    elif value_column == 'ms_caspt2_energy_difference_ev':
        xlabel = 'MS-CASPT2 Energy Difference (eV)'
    elif value_column == 'energy_difference_eV':
        xlabel = 'Energy Difference (eV)'
    else:
        xlabel = value_column
    
    ax1.set_xlabel(xlabel, fontsize=12)
    ax1.set_ylabel('Frequency', fontsize=12)
    # Format the title based on the value column
    if value_column == 'totald_debyes':
        title = f'Total Dipole Moment (Debyes) Distribution: {transition_type}'
    elif value_column == 'caspt2_energy_difference_ev':
        title = f'CASPT2 Energy Difference (eV) Distribution: {transition_type}'
    elif value_column == 'ms_caspt2_energy_difference_ev':
        title = f'MS-CASPT2 Energy Difference (eV) Distribution: {transition_type}'
    elif value_column == 'energy_difference_eV':
        title = f'Energy Difference (eV) Distribution: {transition_type}'
    else:
        title = f'{value_column} Distribution: {transition_type}'
    
    # Add outlier filtering info to title if applicable
    if filter_spatial_outliers:
        title += ' (Spatial Outliers Filtered)'
    
    ax1.set_title(title, fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    # Add statistics text
    mean_val = transition_data[value_column].mean()
    std_val = transition_data[value_column].std()
    min_val = transition_data[value_column].min()
    max_val = transition_data[value_column].max()
    
    stats_text = f'Mean: {mean_val:.3f}\nStd: {std_val:.3f}\nMin: {min_val:.3f}\nMax: {max_val:.3f}'
    ax1.text(0.02, 0.98, stats_text, transform=ax1.transAxes, fontsize=10,
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Box plot
    ax2.boxplot(transition_data[value_column], patch_artist=True, 
                boxprops=dict(facecolor='lightgreen', alpha=0.7))
    ax2.set_ylabel(xlabel, fontsize=12)
    # Format the box plot title based on the value column
    if value_column == 'totald_debyes':
        box_title = f'Total Dipole Moment (Debyes) Box Plot: {transition_type}'
    elif value_column == 'caspt2_energy_difference_ev':
        box_title = f'CASPT2 Energy Difference (eV) Box Plot: {transition_type}'
    elif value_column == 'ms_caspt2_energy_difference_ev':
        box_title = f'MS-CASPT2 Energy Difference (eV) Box Plot: {transition_type}'
    elif value_column == 'energy_difference_eV':
        box_title = f'Energy Difference (eV) Box Plot: {transition_type}'
    else:
        box_title = f'{value_column} Box Plot: {transition_type}'
    
    # Add outlier filtering info to box plot title if applicable
    if filter_spatial_outliers:
        box_title += ' (Spatial Outliers Filtered)'
    
    ax2.set_title(box_title, fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    # Add outliers info
    Q1 = transition_data[value_column].quantile(0.25)
    Q3 = transition_data[value_column].quantile(0.75)
    IQR = Q3 - Q1
    outliers = transition_data[(transition_data[value_column] < Q1 - 1.5*IQR) | 
                              (transition_data[value_column] > Q3 + 1.5*IQR)]
    
    outlier_text = f'Outliers: {len(outliers)} points\n({len(outliers)/len(transition_data)*100:.1f}%)'
    ax2.text(0.02, 0.98, outlier_text, transform=ax2.transAxes, fontsize=10,
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.8))
    
    plt.tight_layout()
    
    if save_plot:
        filename = os.path.join(output_folder, f"{value_column}_distribution_{transition_type.replace('*', 'star').replace(' ', '_')}.png")
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Distribution plots saved as: {filename}")
    
    # Close the figure instead of showing it
    plt.close(fig)
    
    return fig, (ax1, ax2)

def create_scatter_plot(df, transition_type, value_column='totald_debyes', output_folder='output_plots', save_plot=True, filter_spatial_outliers=False):
    """Create a scatter plot showing phi vs psi with color-coded values"""
    
    # Filter data for the specific transition
    transition_data = df[df['transition'] == transition_type].copy()
    
    if len(transition_data) == 0:
        print(f"No data found for transition: {transition_type}")
        return
    
    # Transform coordinates
    transition_data = transform_phi_psi_coordinates(transition_data)
    
    # Filter out spatial outliers if requested
    if filter_spatial_outliers:
        # Filter out geometries where totald_debyes > 1 for n->pi* L or n->pi* R transitions
        if transition_type in ['n->pi* L', 'n->pi* R']:
            # Get geometries to exclude (where dipole > 1)
            high_dipole_geometries = df[
                (df['transition'] == transition_type) & 
                (df['totald_debyes'] > 1)
            ]['geometry'].unique()
            
            # Filter out these geometries from the current transition data
            if len(high_dipole_geometries) > 0:
                transition_data = transition_data[~transition_data['geometry'].isin(high_dipole_geometries)]
                print(f"Filtered out {len(high_dipole_geometries)} geometries with dipole > 1 for {transition_type}")
    
    # Create the scatter plot
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Classify conformations and color-code the scatter points
    conformation_colors = []
    for _, row in transition_data.iterrows():
        conf_type = classify_conformation(row['phi_transformed'], row['psi_transformed'])
        conformation_colors.append(get_conformation_colors()[conf_type])
    
    scatter = ax.scatter(transition_data['phi_transformed'], transition_data['psi_transformed'], 
                        c=conformation_colors, s=50, alpha=0.7)
    
    # Add legend for conformation types
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=10, label='Beta-sheet'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label='Helix'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='darkred', markersize=10, label='Alpha-helix'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='orange', markersize=10, label='Left-handed helix'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='purple', markersize=10, label='Polyproline II'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=10, label='Other')
    ]
    ax.legend(handles=legend_elements, loc='upper right')
    
    # Customize the plot
    ax.set_xlabel('Phi (degrees)', fontsize=12)
    ax.set_ylabel('Psi (degrees)', fontsize=12)
    # Format the title based on the value column
    if value_column == 'totald_debyes':
        title = f'Total Dipole Moment (Debyes) Scatter Plot: {transition_type}'
    elif value_column == 'caspt2_energy_difference_ev':
        title = f'CASPT2 Energy Difference (eV) Scatter Plot: {transition_type}'
    elif value_column == 'ms_caspt2_energy_difference_ev':
        title = f'MS-CASPT2 Energy Difference (eV) Scatter Plot: {transition_type}'
    elif value_column == 'energy_difference_eV':
        title = f'Energy Difference (eV) Scatter Plot: {transition_type}'
    else:
        title = f'{value_column} Scatter Plot: {transition_type}'
    
    # Add outlier filtering info to title if applicable
    if filter_spatial_outliers:
        title += ' (Spatial Outliers Filtered)'
    
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    # Format the colorbar label based on the value column
    if value_column == 'totald_debyes':
        cbar_label = 'Total Dipole Moment (Debyes)'
    elif value_column == 'caspt2_energy_difference_ev':
        cbar_label = 'CASPT2 Energy Difference (eV)'
    elif value_column == 'ms_caspt2_energy_difference_ev':
        cbar_label = 'MS-CASPT2 Energy Difference (eV)'
    elif value_column == 'energy_difference_eV':
        cbar_label = 'Energy Difference (eV)'
    else:
        cbar_label = value_column
    
    # Note: Colorbar removed since we're now using discrete colors for conformations
    
    # Set axis limits
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    
    # Add grid
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_plot:
        filename = os.path.join(output_folder, f"{value_column}_scatter_{transition_type.replace('*', 'star').replace(' ', '_')}.png")
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Scatter plot saved as: {filename}")
    
    # Close the figure instead of showing it
    plt.close(fig)
    
    return fig, ax

def analyze_continuity(df, transition_type, value_column='totald_debyes'):
    """Analyze the continuity of values for a transition"""
    
    transition_data = df[df['transition'] == transition_type].copy()
    
    if len(transition_data) == 0:
        print(f"No data found for transition: {transition_type}")
        return
    
    # Group by phi-psi coordinates and calculate statistics
    stats = transition_data.groupby(['phi', 'psi'])[value_column].agg(['mean', 'std', 'count']).reset_index()
    
    print(f"\n=== Continuity Analysis for {transition_type} - {value_column} ===")
    print(f"Total data points: {len(transition_data)}")
    print(f"Unique phi-psi combinations: {len(stats)}")
    print(f"{value_column} range: {stats['mean'].min():.3f} - {stats['mean'].max():.3f}")
    print(f"Standard deviation range: {stats['std'].min():.3f} - {stats['std'].max():.3f}")
    
    # Check for gaps in the grid
    phi_values = sorted(transition_data['phi'].unique())
    psi_values = sorted(transition_data['psi'].unique())
    
    print(f"Phi range: {min(phi_values)} - {max(phi_values)} degrees")
    print(f"Psi range: {min(psi_values)} - {max(psi_values)} degrees")
    print(f"Expected grid size: {len(phi_values)} x {len(psi_values)} = {len(phi_values) * len(psi_values)}")
    print(f"Actual data points: {len(stats)}")
    
    if len(stats) < len(phi_values) * len(psi_values):
        print("⚠️  Gaps detected in the phi-psi grid!")
    else:
        print("✓ Complete grid coverage")

def analyze_transition_for_all_columns(df, transition_type, value_columns, output_folder='output_plots'):
    """Analyze a single transition for all value columns"""
    
    print(f"\n{'='*80}")
    print(f"Processing transition: {transition_type}")
    print(f"{'='*80}")
    
    for value_column in value_columns:
        print(f"\n--- Analyzing {value_column} ---")
        
        # Analyze continuity
        analyze_continuity(df, transition_type, value_column=value_column)
        
        # Create 3D plot
        print(f"Creating 3D plot for {value_column}...")
        create_3d_plot(df, transition_type, value_column=value_column, output_folder=output_folder, save_plot=True)
        
        # Create filtered 3D plot (without spatial outliers) for energy difference plots
        if value_column:
            print(f"Creating filtered 3D plot for {value_column} (spatial outliers removed)...")
            create_3d_plot(df, transition_type, value_column=value_column, output_folder=output_folder, 
                          save_plot=True, filter_spatial_outliers=True, outlier_threshold_percentile=95)
        
        # # Create heatmap
        # print(f"Creating heatmap for {value_column}...")
        # create_heatmap_plot(df, transition_type, value_column=value_column, output_folder=output_folder, save_plot=True)
        
        # # Create filtered heatmap for energy difference plots
        # if value_column:
        #     print(f"Creating filtered heatmap for {value_column} (spatial outliers removed)...")
        #     create_heatmap_plot(df, transition_type, value_column=value_column, output_folder=output_folder, 
        #                       save_plot=True, filter_spatial_outliers=True)
        
        # # Create distribution plots
        # print(f"Creating distribution plots for {value_column}...")
        # create_distribution_plots(df, transition_type, value_column=value_column, output_folder=output_folder, save_plot=True)
        
        # # Create filtered distribution plots for energy difference plots
        # if value_column:
        #     print(f"Creating filtered distribution plots for {value_column} (spatial outliers removed)...")
        #     create_distribution_plots(df, transition_type, value_column=value_column, output_folder=output_folder, 
        #                            save_plot=True, filter_spatial_outliers=True)
        
        # # Create scatter plot
        # print(f"Creating scatter plot for {value_column}...")
        # create_scatter_plot(df, transition_type, value_column=value_column, output_folder=output_folder, save_plot=True)
        
        # Create filtered scatter plot for energy difference plots
        if value_column:
            print(f"Creating filtered scatter plot for {value_column} (spatial outliers removed)...")
            create_scatter_plot(df, transition_type, value_column=value_column, output_folder=output_folder, 
                              save_plot=True, filter_spatial_outliers=True)
        
        # Create heatmap scatter plot
        print(f"Creating heatmap scatter plot for {value_column}...")
        create_heatmap_scatter_plot(df, transition_type, value_column=value_column, output_folder=output_folder, save_plot=True)
        
        # Create filtered heatmap scatter plot for energy difference plots
        if value_column:
            print(f"Creating filtered heatmap scatter plot for {value_column} (spatial outliers removed)...")
            create_heatmap_scatter_plot(df, transition_type, value_column=value_column, output_folder=output_folder, 
                                      save_plot=True, filter_spatial_outliers=True)
        
        print(f"Completed analysis for {value_column}")
    
    print(f"\nCompleted all analyses for transition: {transition_type}")

def create_heatmap_scatter_plot(df, transition_type, value_column='totald_debyes', output_folder='output_plots', save_plot=True, filter_spatial_outliers=False):
    """Create a heatmap scatter plot for a specific transition type"""
    
    # Filter data for the specific transition
    transition_data = df[df['transition'] == transition_type].copy()
    
    if len(transition_data) == 0:
        print(f"No data found for transition: {transition_type}")
        return
    
    # Transform coordinates to [-180, 180] range
    transition_data = transform_phi_psi_coordinates(transition_data)
    
    # Filter out spatial outliers if requested
    if filter_spatial_outliers:
        # Filter out geometries where totald_debyes > 1 for n->pi* L or n->pi* R transitions
        if transition_type in ['n->pi* L', 'n->pi* R']:
            # Get geometries to exclude (where dipole > 1)
            high_dipole_geometries = df[
                (df['transition'] == transition_type) & 
                (df['totald_debyes'] > 1)
            ]['geometry'].unique()
            
            # Filter out these geometries from the current transition data
            if len(high_dipole_geometries) > 0:
                transition_data = transition_data[~transition_data['geometry'].isin(high_dipole_geometries)]
                print(f"Filtered out {len(high_dipole_geometries)} geometries with dipole > 1 for {transition_type}")
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Create scatter plot with color mapping
    scatter = ax.scatter(transition_data['phi_transformed'], 
                        transition_data['psi_transformed'], 
                        c=transition_data[value_column], 
                        cmap='viridis', 
                        s=60, 
                        alpha=0.8)
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    
    # Format the colorbar label based on the value column
    if value_column == 'totald_debyes':
        cbar_label = 'Total Dipole Moment (Debyes)'
    elif value_column == 'caspt2_energy_difference_ev':
        cbar_label = 'CASPT2 Energy Difference (eV)'
    elif value_column == 'ms_caspt2_energy_difference_ev':
        cbar_label = 'MS-CASPT2 Energy Difference (eV)'
    elif value_column == 'energy_difference_eV':
        cbar_label = 'Energy Difference (eV)'
    else:
        cbar_label = value_column
    
    cbar.set_label(cbar_label, fontsize=12)
    
    # Set axis limits
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    
    # Add grid
    ax.grid(True, alpha=0.3)
    
    # Add conformation regions
    add_conformation_regions(ax)
    
    # Format the title based on the value column
    if value_column == 'totald_debyes':
        title = f'Total Dipole Moment (Debyes) Heatmap Scatter: {transition_type}'
    elif value_column == 'caspt2_energy_difference_ev':
        title = f'CASPT2 Energy Difference (eV) Heatmap Scatter: {transition_type}'
    elif value_column == 'ms_caspt2_energy_difference_ev':
        title = f'MS-CASPT2 Energy Difference (eV) Heatmap Scatter: {transition_type}'
    elif value_column == 'energy_difference_eV':
        title = f'Energy Difference (eV) Heatmap Scatter: {transition_type}'
    else:
        title = f'{value_column} Heatmap Scatter: {transition_type}'
    
    # Add outlier filtering info to title if applicable
    if filter_spatial_outliers:
        title += ' (Spatial Outliers Filtered)'
    
    plt.title(title, fontsize=14, fontweight='bold')
    plt.xlabel('Phi (degrees)', fontsize=12)
    plt.ylabel('Psi (degrees)', fontsize=12)
    
    plt.tight_layout()
    
    if save_plot:
        # Include outlier filtering info in filename if applicable
        filename_suffix = ""
        if filter_spatial_outliers:
            filename_suffix = "_outliers_filtered"
        
        filename = os.path.join(output_folder, f"{value_column}_heatmap_scatter_{transition_type.replace('*', 'star').replace(' ', '_')}{filename_suffix}.png")
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Heatmap scatter plot saved as: {filename}")
    
    # Close the figure
    plt.close(fig)
    
    return fig

def add_conformation_regions(ax):
    """Add protein conformation regions to the plot"""
    # Beta-sheet regions
    ax.axhspan(90, 180, xmin=0.25, xmax=0.42, alpha=0.2, color='blue', label='Beta-sheet')
    
    # Helix regions
    ax.axhspan(-90, 0, xmin=0.25, xmax=0.42, alpha=0.2, color='red', label='Helix')
    
    # Alpha-helix regions
    ax.axhspan(-45, 30, xmin=0.33, xmax=0.5, alpha=0.2, color='darkred', label='Alpha-helix')
    
    # Left-handed helix regions
    ax.axhspan(0, 60, xmin=0.625, xmax=0.71, alpha=0.2, color='orange', label='Left-handed helix')
    
    # Polyproline II regions
    ax.axhspan(115, 175, xmin=0.25, xmax=0.42, alpha=0.2, color='purple', label='Polyproline II')

def create_combined_heatmap_scatter_plots(df, value_column, transitions, output_folder='output_plots', save_plot=True, filter_spatial_outliers=False):
    """Create combined heatmap scatter plots for all transitions on one page"""
    
    # Transform coordinates
    df_transformed = transform_phi_psi_coordinates(df)
    
    # Create figure with subplots
    n_transitions = len(transitions)
    n_cols = 2
    n_rows = (n_transitions + 1) // 2  # Ceiling division
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 8 * n_rows))
    
    # Flatten axes array for easier indexing
    if n_rows == 1:
        axes = [axes] if n_cols == 1 else axes
    else:
        axes = axes.flatten()
    
    # Create a single colorbar for all subplots
    vmin = float('inf')
    vmax = float('-inf')
    
    # filter out 
    df_transformed['geometry'] = df_transformed['phi_transformed'].astype(str) + '_' + df_transformed['psi_transformed'].astype(str)

    if filter_spatial_outliers:
        df_transformed = filter_out_erroneous_geometries(df_transformed)

    # First pass: find global min/max values
    for transition in transitions:
        transition_data = df_transformed[df_transformed['transition'] == transition].copy()        
        if len(transition_data) > 0:
            vmin = min(vmin, transition_data[value_column].min())
            vmax = max(vmax, transition_data[value_column].max())
    
    # Second pass: create the plots
    for i, transition in enumerate(transitions):
        if i >= len(axes):
            break
            
        ax = axes[i]
        transition_data = df_transformed[df_transformed['transition'] == transition].copy()

        if len(transition_data) == 0:
            ax.text(0.5, 0.5, f'No data for {transition}', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'{transition}', fontsize=12, fontweight='bold')
            continue
        
        # Create scatter plot with color mapping
        scatter = ax.scatter(transition_data['phi_transformed'], 
                           transition_data['psi_transformed'], 
                           c=transition_data[value_column], 
                           cmap='viridis', 
                           s=50, 
                           alpha=0.8,
                           vmin=vmin, 
                           vmax=vmax)
        
        # Set axis limits
        ax.set_xlim(-180, 180)
        ax.set_ylim(-180, 180)
        
        # Add grid
        ax.grid(True, alpha=0.3)
        
        # Set title
        ax.set_title(f'{transition}', fontsize=12, fontweight='bold')
        
        # Set labels
        ax.set_xlabel('Phi (degrees)', fontsize=10)
        ax.set_ylabel('Psi (degrees)', fontsize=10)
    
    # Hide unused subplots
    for i in range(len(transitions), len(axes)):
        axes[i].set_visible(False)
    
    # Add colorbar
    cbar = fig.colorbar(scatter, ax=axes, shrink=0.8, aspect=30)
    
    # Format the colorbar label based on the value column
    if value_column == 'totald_debyes':
        cbar_label = 'Total Dipole Moment (Debyes)'
    elif value_column == 'caspt2_energy_difference_ev':
        cbar_label = 'CASPT2 Energy Difference (eV)'
    elif value_column == 'ms_caspt2_energy_difference_ev':
        cbar_label = 'MS-CASPT2 Energy Difference (eV)'
    elif value_column == 'energy_difference_eV':
        cbar_label = 'Energy Difference (eV)'
    else:
        cbar_label = value_column
    
    cbar.set_label(cbar_label, fontsize=12)
    
    # Add overall title
    if value_column == 'totald_debyes':
        overall_title = 'Total Dipole Moment (Debyes) Heatmap Scatter Plots - All Transitions'
    elif value_column == 'caspt2_energy_difference_ev':
        overall_title = 'CASPT2 Energy Difference (eV) Heatmap Scatter Plots - All Transitions'
    elif value_column == 'ms_caspt2_energy_difference_ev':
        overall_title = 'MS-CASPT2 Energy Difference (eV) Heatmap Scatter Plots - All Transitions'
    elif value_column == 'energy_difference_eV':
        overall_title = 'Energy Difference (eV) Heatmap Scatter Plots - All Transitions'
    else:
        overall_title = f'{value_column} Heatmap Scatter Plots - All Transitions'
    
    # Add outlier filtering info to title if applicable
    if filter_spatial_outliers:
        overall_title += ' (Spatial Outliers Filtered)'
    
    fig.suptitle(overall_title, fontsize=16, fontweight='bold')
        
    if save_plot:
        # Include outlier filtering info in filename if applicable
        filename_suffix = ""
        if filter_spatial_outliers:
            filename_suffix = "_outliers_filtered"
        
        filename = os.path.join(output_folder, f"{value_column}_combined_heatmap_scatter_all_transitions{filename_suffix}.png")
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Combined heatmap scatter plot saved as: {filename}")
    
    # Close the figure
    plt.close(fig)
    
    return fig

def main(output_folder='output_plots'):
    """Main function to create all visualizations for all value columns"""
    
    # Create output folder
    output_folder = create_output_folder(output_folder)
    
    # Create combined plots folder
    combined_plots_folder = create_output_folder(os.path.join(output_folder, 'combined_plots'))
    
    # Create subfolders for different plot types
    scatter_plots_folder = create_output_folder(os.path.join(combined_plots_folder, 'scatter_plots'))
    d3_plots_folder = create_output_folder(os.path.join(combined_plots_folder, '3d_plots'))
    
    # Load data
    df, transitions, value_columns = load_and_prepare_data('merged_results.csv')
    
    # Filter out GROUND transition
    transitions = [t for t in transitions if t != 'GROUND']
    
    print(f"\nFound {len(transitions)} transition types (excluding GROUND):")
    for i, transition in enumerate(transitions, 1):
        print(f"{i}. {transition}")
    
    print(f"\nWill analyze {len(value_columns)} value columns:")
    for i, column in enumerate(value_columns, 1):
        print(f"{i}. {column}")
    
    # Create combined 3D plots for each value column (all transitions on one page)
    print(f"\nCreating combined 3D plots for each value column...")
    for value_column in value_columns:
        print(f"Creating combined 3D plot for {value_column}...")
        # Create unfiltered version
        create_combined_3d_plots(df, value_column, transitions, output_folder=d3_plots_folder, save_plot=True)
        
        # Create filtered version (without spatial outliers) for energy difference plots
        print(f"Creating filtered combined 3D plot for {value_column} (spatial outliers removed)...")
        create_combined_3d_plots(df, value_column, transitions, output_folder=d3_plots_folder, 
                                save_plot=True, filter_spatial_outliers=True, outlier_threshold_percentile=95)
    
    # # Create individual plots for each transition and value column
    # for transition in transitions:
    #     analyze_transition_for_all_columns(df, transition, value_columns, output_folder=output_folder)
    
    # Create combined heatmap scatter plots for each value column
    print(f"\nCreating combined heatmap scatter plots for each value column...")
    for value_column in value_columns:
        print(f"Creating combined heatmap scatter plot for {value_column}...")
        # Create unfiltered version
        create_combined_heatmap_scatter_plots(df, value_column, transitions, output_folder=scatter_plots_folder, save_plot=True)
        
        # Create filtered version (without spatial outliers) for energy difference plots
        print(f"Creating filtered combined heatmap scatter plot for {value_column} (spatial outliers removed)...")
        create_combined_heatmap_scatter_plots(df, value_column, transitions, output_folder=scatter_plots_folder, 
                                            save_plot=True, filter_spatial_outliers=True)
    

if __name__ == "__main__":
    OUTPUT_FOLDER = 'output_plots'  # Change this to your desired output folder name
    main(output_folder=OUTPUT_FOLDER)
