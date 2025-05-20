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

% Default configuration structure
config = struct(...
    'workingDir', pwd, ...
    'dataDirHeader', fullfile(pwd, 'MonteCarloResults'), ...
    'inputFile', '../input/Inputs_mpox2024_set2.xlsx', ...
    'num_iterations', 5, ...
    'waning_ve_mode', 2, ...  % 0: No waning, 1: wanes to 0, 2: wanes to 50%, 3: wanes to 75%, 4: wanes to 25%
    'scenarios', [8], ...     % [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22]
    'enable_sensitivity', 0);

% Create output directory if it doesn't exist
if ~exist(config.dataDirHeader, 'dir')
    mkdir(config.dataDirHeader);
end

% Main simulation loop
for scenario = config.scenarios
    % Create version-specific directory
    testVersion = sprintf('test_mpox2024_S%d', scenario);
    testVersionPath = fullfile(config.dataDirHeader, testVersion);
    mkdir(testVersionPath);
    
    % Update config
    config.scenario = scenario;
    
    fprintf('*** Running Simulation: %s ***\n', testVersion);
    tStart = tic;
    
    % Only use parallel processing if we have multiple iterations
    if config.num_iterations > 1
        if isempty(gcp('nocreate'))
            maxWorkers = max(1, numCores - 2);
            if maxWorkers > 1
                parpool('local', maxWorkers);
            end
        end
        
        % Run iterations in parallel
        parfor iters = 1:config.num_iterations
            runIteration(iters, testVersionPath, config);
        end
        
        % Delete the parallel pool after iterations are complete
        delete(gcp('nocreate'));
    else
        % Run single iteration sequentially
        runIteration(1, testVersionPath, config);
    end
    
    % Generate metrics after all iterations complete
    fprintf('Generating metrics for scenario %d...\n', scenario);
    try
        cd(config.workingDir);
        gen_metric;
    catch ME
        fprintf('Error generating metrics for scenario %d: %s\n', scenario, ME.message);
        continue;
    end
    
    % Report execution time
    tEnd = toc(tStart);
    fprintf('Simulation completed in %d hours and %d minutes\n', ...
        floor(tEnd/3600), floor(rem(tEnd,3600)/60));
    
    cd(config.workingDir);
end

function runIteration(iters, testVersionPath, config)
    % Create unique directory for this iteration
    sim_dataDir = fullfile(testVersionPath, sprintf('iter%d/state_matrices/', iters));
    mkdir(sim_dataDir);
    
    % Set random seed for reproducibility
    rng('shuffle');
    
    % Set the current directory to the working directory
    cd(config.workingDir);
    
    % Set global variables for mpox2024_shellMod
    global SIM_INPUT_FILE WANING_VE_MODE ENABLE_SENSITIVITY NUM_ITERATIONS S
    SIM_INPUT_FILE = config.inputFile;
    WANING_VE_MODE = config.waning_ve_mode;
    ENABLE_SENSITIVITY = config.enable_sensitivity;
    NUM_ITERATIONS = config.num_iterations;
    S = config.scenario;
    
    % Run simulation
    try
        mpox2024_shellMod;
        
        % Create tally file path
        tallyPath = fullfile(sim_dataDir, sprintf('Tally_%s.csv', testVersion));
        
        % Save tally data
        if exist('tally', 'var')
            tally_tbl = array2table(tally);
            tally_tbl.Properties.VariableNames(1:end) = tallyHeaders;
            writetable(tally_tbl, tallyPath);
        end
    catch ME
        fprintf('Error in iteration %d for scenario %d: %s\n', iters, config.scenario, ME.message);
    end
end