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

% Store the original working directory
originalDir = pwd;

% Check system capabilities
numCores = feature('numcores');
fprintf('Number of CPU cores: %d\n', numCores);

% Default configuration structure
config = struct(...
    'workingDir', pwd, ...
    'dataDirHeader', fullfile(pwd, 'MonteCarloResults'), ...
    'inputFile', '../input/Inputs_mpox2024_set2.xlsx', ...
    'num_iterations', 20, ...
    'waning_ve_mode', 2, ...  % 0: No waning, 1: wanes to 0, 2: wanes to 50%, 3: wanes to 75%, 4: wanes to 25%
    'scenarios', [17 18 19], ...     % [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22]
    'enable_sensitivity', 0);

% Create output directory if it doesn't exist
if ~exist(config.dataDirHeader, 'dir')
    mkdir(config.dataDirHeader);
end

% Main simulation loop
for scenario = config.scenarios
    % Create version-specific directory
    testVersion = sprintf('new_mpox2024_S%d', scenario);
    testVersionPath = fullfile(config.dataDirHeader, testVersion);
    mkdir(testVersionPath);
    
    % Update config
    config.scenario = scenario;
    
    fprintf('*** Running Simulation: %s ***\n', testVersion);
    tStart = tic;
    
    % Only use parallel processing if we have multiple iterations
    if config.num_iterations > 1
        % Clean up any existing parallel pool first
        if ~isempty(gcp('nocreate'))
            delete(gcp('nocreate'));
        end
        
        maxWorkers = max(1, numCores - 2);
        if maxWorkers > 1
            % Check if running on HPC (SLURM)
            if ~isempty(getenv('SLURM_JOB_ID'))
                % HPC configuration
                pc = parallel.cluster.Local;
                job_folder = fullfile('/scratch1/',getenv('USER'),getenv('SLURM_JOB_ID'));
                if ~exist(job_folder, 'dir')
                    mkdir(job_folder);
                end
                set(pc,'JobStorageLocation',job_folder);
                fprintf('Running on HPC. Using parallel job storage location: %s\n', job_folder);
                
                % Initialize parallel pool with explicit configuration
                parpool(pc, maxWorkers, 'AttachedFiles', {}, 'AutoAddClientPath', false);
            else
                % Local configuration
                fprintf('Running locally. Using default parallel pool configuration.\n');
                parpool('local', maxWorkers);
            end
        end
        
        % Run iterations in parallel
        parfor iters = 1:config.num_iterations
            runIteration(iters, testVersionPath, config, testVersion);
        end
        
        % Delete the parallel pool after iterations are complete
        delete(gcp('nocreate'));
    else
        % Run single iteration sequentially
        runIteration(1, testVersionPath, config, testVersion);
    end
    
    % Generate metrics after all iterations complete
    fprintf('Generating metrics for scenario %d...\n', scenario);
    try
        % Ensure we're in the correct directory
        cd(originalDir);
        % Set global variables for metrics generation
        global NUM_ITERATIONS S testVerDir
        NUM_ITERATIONS = config.num_iterations;
        S = config.scenario;
        testVerDir = testVersionPath;  % Set the directory for memo storage
        gen_metric;
    catch ME
        fprintf('Error generating metrics for scenario %d: %s\n', scenario, ME.message);
        continue;
    end
    
    % Report execution time
    tEnd = toc(tStart);
    fprintf('Simulation completed in %d hours and %d minutes\n', ...
        floor(tEnd/3600), floor(rem(tEnd,3600)/60));
    
    cd(originalDir);
end

function runIteration(iters, testVersionPath, config, testVersion)
    % Create unique directory for this iteration
    sim_dataDir = fullfile(testVersionPath, sprintf('iter%d/state_matrices/', iters));
    mkdir(sim_dataDir);
    
    % Set random seed for reproducibility
    rng('shuffle');
    
    % Set the current directory to the working directory
    cd(config.workingDir);
    
    % Set global variables for mpox2024_shellMod
    global SIM_INPUT_FILE WANING_VE_MODE ENABLE_SENSITIVITY NUM_ITERATIONS S testVerDir testVersion
    SIM_INPUT_FILE = config.inputFile;
    WANING_VE_MODE = config.waning_ve_mode;
    ENABLE_SENSITIVITY = config.enable_sensitivity;
    NUM_ITERATIONS = config.num_iterations;
    S = config.scenario;
    testVerDir = testVersionPath;  % Set the directory for memo storage
    testVersion = testVersion;  % Make testVersion available globally
    
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