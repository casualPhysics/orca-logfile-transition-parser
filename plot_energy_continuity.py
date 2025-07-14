import pandas as pd
import numpy as np
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

def create_3d_plot(df, transition_type, value_column='totald_debyes', output_folder='output_plots', save_plot=True):
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
        filename = os.path.join(output_folder, f"energy_continuity_{transition_type.replace('*', 'star').replace(' ', '_')}.png")
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"3D plot saved as: {filename}")
    
    plt.close(fig)
    
    return fig, ax

def create_heatmap_plot(df, transition_type, value_column='totald_debyes', output_folder='output_plots', save_plot=True):
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
    fig, ax = plt.subplots(figsize=(10, 8))
    
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
        filename = os.path.join(output_folder, f"energy_heatmap_{transition_type.replace('*', 'star').replace(' ', '_')}.png")
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Heatmap saved as: {filename}")
    
    plt.close(fig)

def create_distribution_plots(df, transition_type, value_column='totald_debyes', output_folder='output_plots', save_plot=True):
    """Create distribution plots (histogram and box plot) for a specific transition type"""
    
    # Filter data for the specific transition
    transition_data = df[df['transition'] == transition_type].copy()
    
    if len(transition_data) == 0:
        print(f"No data found for transition: {transition_type}")
        return
    
    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Histogram
    ax1.hist(transition_data[value_column], bins=30, alpha=0.7, color='skyblue', edgecolor='black')
    ax1.set_xlabel(value_column, fontsize=12)
    ax1.set_ylabel('Frequency', fontsize=12)
    ax1.set_title(f'{value_column} Distribution: {transition_type}', fontsize=14, fontweight='bold')
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
    ax2.set_ylabel(value_column, fontsize=12)
    ax2.set_title(f'{value_column} Box Plot: {transition_type}', fontsize=14, fontweight='bold')
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
        filename = os.path.join(output_folder, f"distribution_{transition_type.replace('*', 'star').replace(' ', '_')}.png")
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Distribution plots saved as: {filename}")
    
    plt.close(fig)
    
    return fig, (ax1, ax2)

def create_scatter_plot(df, transition_type, value_column='totald_debyes', output_folder='output_plots', save_plot=True):
    """Create a scatter plot showing phi vs psi with color-coded values"""
    
    # Filter data for the specific transition
    transition_data = df[df['transition'] == transition_type].copy()
    
    if len(transition_data) == 0:
        print(f"No data found for transition: {transition_type}")
        return
    
    # Transform coordinates
    transition_data = transform_phi_psi_coordinates(transition_data)
    
    # Create the scatter plot
    fig, ax = plt.subplots(figsize=(10, 8))
    
    scatter = ax.scatter(transition_data['phi_transformed'], transition_data['psi_transformed'], 
                        c=transition_data[value_column], cmap='viridis', s=50, alpha=0.7)
    
    # Customize the plot
    ax.set_xlabel('Phi (degrees)', fontsize=12)
    ax.set_ylabel('Psi (degrees)', fontsize=12)
    ax.set_title(f'{value_column} Scatter Plot: {transition_type}', fontsize=14, fontweight='bold')
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax, label=value_column)
    
    # Set axis limits
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    
    # Add grid
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_plot:
        filename = os.path.join(output_folder, f"scatter_{transition_type.replace('*', 'star').replace(' ', '_')}.png")
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Scatter plot saved as: {filename}")
    
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

def main(value_column='totald_debyes', output_folder='output_plots'):
    """Main function to create all visualizations"""
    
    # Create output folder
    output_folder = create_output_folder(output_folder)
    
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
        create_3d_plot(df, transition, value_column=value_column, output_folder=output_folder, save_plot=True)
        
        # Create heatmap
        print(f"\nCreating heatmap for {transition}...")
        create_heatmap_plot(df, transition, value_column=value_column, output_folder=output_folder, save_plot=True)
        
        # Create distribution plots
        print(f"\nCreating distribution plots for {transition}...")
        create_distribution_plots(df, transition, value_column=value_column, output_folder=output_folder, save_plot=True)
        
        # Create scatter plot
        print(f"\nCreating scatter plot for {transition}...")
        create_scatter_plot(df, transition, value_column=value_column, output_folder=output_folder, save_plot=True)
        
        print(f"\nCompleted analysis for {transition}")
    
    print(f"\n{'='*60}")
    print(f"All visualizations completed and saved to: {output_folder}")
    print(f"{'='*60}")

if __name__ == "__main__":
    # You can change this to any column name you want to analyze
    COLUMN_TO_ANALYZE = 'totald_debyes'  # Change this to your desired column name
    OUTPUT_FOLDER = 'output_plots'  # Change this to your desired output folder name
    main(value_column=COLUMN_TO_ANALYZE, output_folder=OUTPUT_FOLDER) 