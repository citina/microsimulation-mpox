import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import seaborn as sns

# Set the style for better visualization
plt.style.use('bmh')  # Using a built-in matplotlib style
sns.set_theme()  # This will set a nice default style from seaborn

# Base paths
old_version_path = Path("/Users/liangshiting/Desktop/USCResearch/Mpox_2024/MonteCarloResults")
new_version_path = Path("/Users/liangshiting/Desktop/USCResearch/Mpox_2024/Code_updated/MonteCarloResults")

# Scenarios to compare
scenarios = [
    "mpox2024_S8",
    "mpox2024_S9",
    "mpox2024_S10",
    "mpox2024_S15",
    "mpox2024_S16",
    "mpox2024_S17",
    "mpox2024_S18",
    "mpox2024_S19"
]

# Metrics to compare
metrics = ['ToAware', 'r_t', 'NewInfections']

def plot_comparison(old_data, new_data, metric, scenario, save_path):
    """Create comparison plot for a single metric."""
    plt.figure(figsize=(12, 6))
    
    # Plot both versions
    plt.plot(old_data['Week'], old_data['Mean'], 'b-', label='Old Version', linewidth=2)
    plt.plot(new_data['Week'], new_data['Mean'], 'r--', label='New Version', linewidth=2)
    
    # Add confidence intervals
    plt.fill_between(old_data['Week'], 
                    old_data['LowerBound'], 
                    old_data['UpperBound'], 
                    color='blue', alpha=0.1)
    plt.fill_between(new_data['Week'], 
                    new_data['LowerBound'], 
                    new_data['UpperBound'], 
                    color='red', alpha=0.1)
    
    # Display mean values
    old_mean = old_data['Mean'].mean()
    new_mean = new_data['Mean'].mean()
    plt.text(0.02, 0.98, 
             f'Old Version Mean: {old_mean:.2f}\nNew Version Mean: {new_mean:.2f}',
             transform=plt.gca().transAxes,
             bbox=dict(facecolor='white', alpha=0.8),
             verticalalignment='top')
    
    plt.title(f'{metric} Comparison - {scenario}')
    plt.xlabel('Week')
    plt.ylabel(metric)
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Save plot
    plt.savefig(save_path / f'{metric}_{scenario}_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()

def main():
    # Create output directory for plots
    output_dir = Path("comparison_plots")
    output_dir.mkdir(exist_ok=True)
    
    # Process each scenario and metric
    for scenario in scenarios:
        print(f"\nProcessing scenario: {scenario}")
        
        # Adjust scenario name for new version
        new_scenario = f"test_{scenario}"
        
        for metric in metrics:
            print(f"  Comparing metric: {metric}")
            
            # Construct file paths
            old_file = old_version_path / scenario / f"{metric}_bs_results.csv"
            new_file = new_version_path / new_scenario / f"{metric}_bs_results.csv"
            
            # Check if files exist
            if not old_file.exists() or not new_file.exists():
                print(f"    Warning: Missing files for {metric} in {scenario}")
                continue
            
            # Read data
            try:
                old_data = pd.read_csv(old_file)
                new_data = pd.read_csv(new_file)
                
                # Create comparison plot
                plot_comparison(old_data, new_data, metric, scenario, output_dir)
                
                # Print mean values
                old_mean = old_data['Mean'].mean()
                new_mean = new_data['Mean'].mean()
                print(f"    Old Version Mean: {old_mean:.2f}")
                print(f"    New Version Mean: {new_mean:.2f}")
                
            except Exception as e:
                print(f"    Error processing {metric} in {scenario}: {str(e)}")

if __name__ == "__main__":
    main()