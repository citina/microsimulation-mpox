%% Mpox2024 Simulation Script
% This script runs Monte Carlo simulations for Mpox transmission modeling
% Author: Citina Liang
% Date: 2025-04-18
% Version: 1.0

%% Initialization
format long
clear all
clear
close all
clc

% Check system capabilities
numCores = feature('numcores');
fprintf('Number of CPU cores: %d\n', numCores);

%% Simulation Configuration Parameters
% Core simulation settings
NUM_ITERATIONS = 3;           % Number of Monte Carlo iterations
TESTING_WEEK = 12;            % Week when "event" people leave LAC (if applicable)

% Vaccine efficacy settings
WANING_VE_MODE = 2;         % 0: No waning
                            % 1: VE wanes to 0 after 12 months
                            % 2: VE wanes to half after 12 months

% Force of infection (FOI) settings
WEEKLY_FOI = 0;              % 1: Use weekly FOI, 0: Use average FOI
FOI_VALUES = [1.8];          % FOI values to test (piso = 0.5) %S1, S2, S3

% Policy intervention settings
ENABLE_POLICY = 0;           % Flag to enable/disable policy interventions
POLICY_START_WEEK = 1000;    % Week when policy starts (week0: Jun26-Jul02)
ENABLE_SENSITIVITY = 0;      % Flag to enable sensitivity analysis

%% File and Directory Configuration
% Input file configuration
SIM_INPUT_FILE = "../input/Inputs_mpox2024_set2.xlsx";

% Output directory configuration
WORKING_DIR = pwd;
OUTPUT_DIR_HEADER = fileparts(WORKING_DIR) + "/MonteCarloResults";
mkdir(OUTPUT_DIR_HEADER);

%% Main Simulation Loop
fprintf('*** Starting Simulation Process ***\n');
tStart = tic;

% Validate input file exists
if ~exist(SIM_INPUT_FILE, 'file')
    error('Input file not found: %s', SIM_INPUT_FILE);
end

% Iterate through FOI values
for foi = FOI_VALUES
    % Create version-specific directory
    testVersion = sprintf('mpox2024_foi_%.2f_iso0.2_test', foi);
    testVerDir = fullfile(OUTPUT_DIR_HEADER, testVersion);
    mkdir(testVerDir);
    
    fprintf('Processing FOI value: %.2f\n', foi);
    
    % Set up parallel processing if multiple iterations
    if NUM_ITERATIONS > 1
        % Calculate maximum number of workers (leave some cores free)
        maxWorkers = max(1, numCores - 4);
        
        % Create a parallel pool if it doesn't exist
        if isempty(gcp('nocreate'))
            parpool('local', min(maxWorkers, NUM_ITERATIONS));
        end
        
        % Run Monte Carlo iterations in parallel
        parfor iter = 1:NUM_ITERATIONS
            fprintf('Running iteration %d/%d\n', iter, NUM_ITERATIONS);
            
            % Create iteration-specific directory
            sim_dataDir = fullfile(testVerDir, sprintf('iter%d/state_matrices/', iter));
            mkdir(sim_dataDir);
            
            % Set random seed for reproducibility
            rng('shuffle');
            
            % Run simulation
            mpox2024_shellMod;
        end
        
        % Clean up parallel pool
        delete(gcp('nocreate'));
    else
        % Single iteration case
        iter = 1;
        fprintf('Running single iteration\n');
        
        % Create iteration-specific directory
        sim_dataDir = fullfile(testVerDir, sprintf('iter%d/state_matrices/', iter));
        mkdir(sim_dataDir);
        
        % Set random seed for reproducibility
        rng('shuffle');
        
        % Run simulation
        mpox2024_shellMod;
    end
    
    % Generate metrics after all iterations
    gen_metric;
end

% Report execution time
tEnd = toc(tStart);
fprintf('Simulation completed in %d hours and %d minutes\n', ...
    floor(tEnd/3600), floor(rem(tEnd,3600)/60));

% Return to original working directory
cd(WORKING_DIR);