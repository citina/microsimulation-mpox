# microsimulation-mpox-LAC

## Description
This repository hosts the MATLAB code for the microsimulation model used in the study titled **_"Viral introductions and return to baseline sexual behaviors maintain low-level mpox incidence in Los Angeles County, USA, 2023-2024."_** The model was employed to generate the results and figures for the associated journal article. The code is compatible with MATLAB R2024a.

## Installation

### Prerequisites
Ensure you have MATLAB R2024a installed on your machine to run the simulation code. 

### Setup
1. Clone this repository to your local machine: `git clone https://github.com/citina/microsimulation-mpox.git`
2. Navigate to the cloned directory: `cd microsimulation-mpox`
   
## Structure and Usage
The repository contains the following key components:

### Input Files
**input/**: Directory containing all necessary input files for the simulation:
- `Inputs_mpox2024_set2.xlsx`: Main input file containing model parameters and transition probabilities
- Other supporting input files required for the simulation

### Core Simulation Files
#### Shell Scripts
The MATLAB sources now live inside the `Code_updated/` subdirectory. The main shell script for running simulations is therefore located at `Code_updated/Mpox2024_ShellScript.m`. This script supports multiple scenarios; you can configure:
  * Number of Monte Carlo iterations (`num_iterations`)
  * Vaccine effectiveness waning mode (`waning_ve_mode`: 0-4)
    - `0`: No waning
    - `1`: Wanes to 0%
    - `2`: Wanes to 50%
    - `3`: Wanes to 75%
    - `4`: Wanes to 25%
  * Scenarios to run (`scenarios`: array of scenario numbers 0-26)
    - Detailed descriptions of each scenario and their corresponding numbers can be found in `Mpox_2024.pdf`
  * Sensitivity mode (`enable_sensitivity`):
    - `0`: Off (use baseline import schedule)
    - `1`: Double imports (each scheduled import count multiplied by 2)
    - `2`: Halve imports using probabilistic rounding (keeps weekly counts integers; total ≈ 50% in expectation)
  * Import awareness (`import_diagnosed`):
    - `1` (default): Imported cases are diagnosed immediately (sets `pox_aware = 1`)
    - `0`: Imported cases are undiagnosed (sets `pox_aware = 0`)
  * Imported case symptom status (`import_asymptomatic`):
    - `0` (default): Imported cases enter as symptomatic (sets `pox_status = 2`)
    - `1`: Imported cases enter as asymptomatic (sets `pox_status = 1`)
  * Vaccination uptake policy (`enable_vax_policy`):
    - `0` (default): Use baseline `vac_1_new.csv`
    - `1`: Double uptake using `vac_1_new2_2X.csv`
    - `2`: Quadruple uptake using `vac_1_new2_4X.csv`
    - `3`: Tenfold uptake using `vac_1_new2_10X.csv`

#### Main Simulation Scripts:
- `Code_updated/mpox2024_shellMod.m`: Core simulation module
- `Code_updated/mpox2024_parameters.m`: Model parameters configuration
- `Code_updated/infection9.m`: Infection dynamics implementation
- `Code_updated/calc_infec_prob4.m`: Infection probability calculations
- `Code_updated/transition5.m`: State transition logic

#### Helper Functions:
- `Code_updated/find_indices.m`: Utility for index management
- `Code_updated/find_demog_rows.m` and `Code_updated/find_demog_rows2.m`: Demographic data processing
- `Code_updated/create_demog_groups.m`: Population group management
- `Code_updated/cell2csv.m`: Data export utilities
- `Code_updated/col_idx_to_name.m` and `Code_updated/read_table.m`: Data handling utilities

#### Analysis and Visualization:
- `Code_updated/bootstrap_mpox.m`: Bootstrap analysis implementation
- `Code_updated/gen_metric.m`: Metric generation for analysis
- `compare_results.py`: Compare Monte Carlo outputs across scenario folders
- `plot_ve_waning.py`: Plot vaccine waning trajectories across modes and policies
- `visualize_toaware_s20_25.py`: Visualize awareness scenarios (`comparison_plots/`, `vax_comparison_plots/`, and related PNGs contain recent outputs)

### Results and Output
**MonteCarloResults/**: Directory containing simulation results

## Running Simulations
To run a simulation:
1. Open MATLAB
2. Navigate to the repository directory
3. Open `Code_updated/Mpox2024_ShellScript.m` and configure your desired settings:
   - Set `num_iterations` for the number of Monte Carlo simulations
   - Choose `waning_ve_mode` based on your vaccine effectiveness waning scenario
   - Optionally set `enable_sensitivity` to 1 (double imports) or 2 (half imports with probabilistic rounding)
   - Choose imported case symptom status with `import_asymptomatic` (1 = asymptomatic, 0 = symptomatic)
   - Pick a vaccination uptake policy with `enable_vax_policy` if you need alternative uptake schedules
   - Select `scenarios` array with the scenario numbers you want to run
4. Run the script in MATLAB
5. Results will be saved in the `MonteCarloResults` directory, organized by scenario number

## Contact
For queries or collaboration requests, please contact Citina Liang at [citinal@usc.edu](mailto:citinal@usc.edu). Please reach out before using the model to ensure it is applied appropriately. 
