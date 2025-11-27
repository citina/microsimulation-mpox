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
The repository includes a main shell script for running simulations:
- `Mpox2024_ShellScript.m`: Main simulation script that supports multiple scenarios. You can configure:
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
  * Vaccination policy (`enable_vax_policy`):
    - `0` (default): Baseline first-dose uptake schedule
    - `1`: Double uptake (`vac1_2X` transition matrix)
    - `2`: Quadruple uptake (`vac1_4X`)
    - `3`: Ten-fold uptake (`vac1_10X`)

#### Main Simulation Scripts:
- `mpox2024_shellMod.m`: Core simulation module
- `mpox2024_parameters.m`: Model parameters configuration
- `infection9.m`: Infection dynamics implementation
- `calc_infec_prob4.m`: Infection probability calculations
- `transition5.m`: State transition logic

#### Helper Functions:
- `find_indices.m`: Utility for index management
- `find_demog_rows.m` and `find_demog_rows2.m`: Demographic data processing
- `create_demog_groups.m`: Population group management
- `cell2csv.m`: Data export utilities
- `col_idx_to_name.m` and `read_table.m`: Data handling utilities

#### Analysis and Visualization:
- `bootstrap_mpox.m`: Bootstrap analysis implementation that can also produce cumulative metrics and pairwise scenario comparison CSVs under `MonteCarloResults/comparisons/`
- `gen_metric.m`: Metric generation for analysis

### Results and Output
**MonteCarloResults/**: Directory containing simulation results
- `comparisons/`: Optional folder with bootstrap differences between scenario pairs (added when `bootstrap_mpox.m` runs with `comparison_pairs` configured). Each CSV includes week index, mean difference, and 95% bootstrap interval columns.

## Running Simulations
To run a simulation:
1. Open MATLAB
2. Navigate to the repository directory
3. Open `Mpox2024_ShellScript.m` and configure your desired settings:
   - Set `num_iterations` for the number of Monte Carlo simulations
   - Choose `waning_ve_mode` based on your vaccine effectiveness waning scenario
   - Optionally set `enable_sensitivity` to 1 (double imports) or 2 (half imports with probabilistic rounding)
   - Configure import awareness (`import_diagnosed`) and symptom status (`import_asymptomatic`)
   - Select a vaccination uptake policy via `enable_vax_policy` if you need alternative coverage scenarios
   - Select `scenarios` array with the scenario numbers you want to run
4. Run the script in MATLAB
5. Results will be saved in the `MonteCarloResults` directory, organized by scenario number (and, when comparisons are enabled in `bootstrap_mpox.m`, summarized differences will appear under `MonteCarloResults/comparisons/`)

## Contact
For queries or collaboration requests, please contact Citina Liang at [citinal@usc.edu](mailto:citinal@usc.edu). Please reach out before using the model to ensure it is applied appropriately. 
