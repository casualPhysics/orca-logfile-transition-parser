import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns

def load_and_prepare_data(csv_file, value_column='totald_debyes'):
    """Load the CSV data and prepare it for plotting"""
    df = pd.read_csv(csv_file)
    
    # Verify the specified column exists
    if value_column not in df.columns:
        available_columns = list(df.columns)
        raise ValueError(f"Column '{value_column}' not found in the CSV file. Available columns: {available_columns}")
    
    # Get unique transitions
    transitions = df['transition'].unique()
    print(f"Available transitions: {transitions}")
    print(f"Using column '{value_column}' for analysis")
    
    return df, transitions

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

def create_3d_plot(df, transition_type, value_column='totald_debyes', save_plot=True):
    """Create a 3D plot for a specific transition type"""
    
    # Filter data for the specific transition
    transition_data = df[df['transition'] == transition_type].copy()
    
    if len(transition_data) == 0:
        print(f"No data found for transition: {transition_type}")
        return
    
    # Transform coordinates to [-90, 90] range
    transition_data = transform_phi_psi_coordinates(transition_data)
    
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
    
    # Add scatter points for actual data points
    ax.scatter(pivot_data['phi_transformed'], pivot_data['psi_transformed'], pivot_data[value_column], 
              c='red', s=20, alpha=0.7, label='Data points')
    
    # Customize the plot
    ax.set_xlabel('Phi (degrees)', fontsize=12)
    ax.set_ylabel('Psi (degrees)', fontsize=12)
    ax.set_zlabel(f'{value_column}', fontsize=12)
    ax.set_title(f'{value_column} Continuity: {transition_type}', fontsize=14, fontweight='bold')
    
    # Add colorbar
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5, label=value_column)
    
    # Set axis limits to [-180, 180]
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    
    # Add grid
    ax.grid(True, alpha=0.3)
    
    # Add legend
    ax.legend()
    
    plt.tight_layout()
    
    if save_plot:
        filename = f"energy_continuity_{transition_type.replace('*', 'star').replace(' ', '_')}.png"
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Plot saved as: {filename}")
    
    plt.show()
    
    return fig, ax

def create_heatmap_plot(df, transition_type, value_column='totald_debyes', save_plot=True):
    """Create a 2D heatmap for a specific transition type"""
    
    # Filter data for the specific transition
    transition_data = df[df['transition'] == transition_type].copy()
    
    if len(transition_data) == 0:
        print(f"No data found for transition: {transition_type}")
        return
    
    # Transform coordinates to [-90, 90] range
    transition_data = transform_phi_psi_coordinates(transition_data)
    
    # Create a pivot table for the heatmap
    pivot_data = transition_data.groupby(['phi_transformed', 'psi_transformed'])[value_column].mean().reset_index()
    
    # Create pivot table for heatmap
    heatmap_data = pivot_data.pivot(index='psi_transformed', columns='phi_transformed', values=value_column)
    
    # Create the heatmap
    
    # Create custom colormap
    sns.heatmap(heatmap_data, cmap='viridis', annot=False, cbar_kws={'label': value_column})
    
    plt.title(f'{value_column} Heatmap: {transition_type}', fontsize=14, fontweight='bold')
    plt.xlabel('Phi (degrees)', fontsize=12)
    plt.ylabel('Psi (degrees)', fontsize=12)
    
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45)
    plt.yticks(rotation=0)
    
    plt.tight_layout()
    
    if save_plot:
        filename = f"energy_heatmap_{transition_type.replace('*', 'star').replace(' ', '_')}.png"
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Heatmap saved as: {filename}")
    
    plt.show()

def analyze_continuity(df, transition_type, value_column='totald_debyes'):
    """Analyze the continuity of values for a transition"""
    
    transition_data = df[df['transition'] == transition_type].copy()
    
    if len(transition_data) == 0:
        print(f"No data found for transition: {transition_type}")
        return
    
    # Group by phi-psi coordinates and calculate statistics
    stats = transition_data.groupby(['phi', 'psi'])[value_column].agg(['mean', 'std', 'count']).reset_index()
    
    print(f"\n=== Continuity Analysis for {transition_type} ===")
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

def main(value_column='totald_debyes'):
    """Main function to create all visualizations"""
    
    # Load data
    df, transitions = load_and_prepare_data('merged_results.csv', value_column=value_column)
    
    print(f"\nFound {len(transitions)} transition types:")
    for i, transition in enumerate(transitions, 1):
        print(f"{i}. {transition}")
    
    # Create plots for each transition
    for transition in transitions:
        print(f"\n{'='*60}")
        print(f"Processing transition: {transition}")
        print(f"{'='*60}")
        
        # Analyze continuity
        analyze_continuity(df, transition, value_column=value_column)
        
        # Create 3D plot
        print(f"\nCreating 3D plot for {transition}...")
        create_3d_plot(df, transition, value_column=value_column, save_plot=True)
        
        # Create heatmap
        print(f"\nCreating heatmap for {transition}...")
        create_heatmap_plot(df, transition, value_column=value_column, save_plot=True)
        
        print(f"\nCompleted analysis for {transition}")
    
    print(f"\n{'='*60}")
    print("All visualizations completed!")
    print(f"{'='*60}")

if __name__ == "__main__":
    # You can change this to any column name you want to analyze
    COLUMN_TO_ANALYZE = 'totald_debyes'  # Change this to your desired column name
    main(value_column=COLUMN_TO_ANALYZE) 