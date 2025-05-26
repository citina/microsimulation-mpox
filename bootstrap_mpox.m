clear all
close all
clear
clc

%{
This file creates boostrapped results for selected metrics
Main idea: 
1. create a matrix with x rows and y columns, where x is the number
of bootstraps, and y is the number of iterations we have in total
2. for each row, we randomly from the iterations we got with replacement
3. calculate the avg. for each row - new column named z
4. calculate the avg., UB, LB for column z, and this is our final result
%}

%% bootstrap settings
bs = 500; % should be at least 20 for the code to work
iterations = 1:20;
numWks = 85;

% List of scenarios to process
scenario_names = {
    "new_mpox2024_S17",
    "new_mpox2024_S18",
    "new_mpox2024_S19"
};
%% paths
% Base path for MonteCarloResults
basePath = pwd;
monteCarloPath = fullfile(basePath, "MonteCarloResults");

% Check if MonteCarloResults directory exists
if ~isfolder(monteCarloPath)
    error('MonteCarloResults directory does not exist: %s', monteCarloPath);
end

% Process each scenario
for scenarioIdx = 1:length(scenario_names)
    Scenario_name = scenario_names{scenarioIdx};
    fprintf('\nProcessing scenario %d of %d: %s\n', scenarioIdx, length(scenario_names), Scenario_name);
    
    % Set paths for current scenario
    InPath = fullfile(monteCarloPath, Scenario_name);
    OutPath = InPath;

    % Check if scenario directory exists
    if ~isfolder(InPath)
        warning('Skipping scenario %s: Directory does not exist', Scenario_name);
        continue;
    end

    % Check if scenario directory has required structure
    if ~isfolder(fullfile(InPath, 'iter1'))
        warning('Skipping scenario %s: Missing required directory structure', Scenario_name);
        continue;
    end

    % a matrix to store indecies of iterations
    % e.g.: if 3 iterations, then randomly select 3 samples with replacement
    rng(666);
    smpls = zeros(bs, length(iterations));
    for i = 1:bs
        smpls(i,:) = randsample(iterations, length(iterations), true);
    end

    % metricShelf is a struct where each field is a metric name, and each
    % metric contains subfields for each iteration.
    % Metrics of interest
    metrics_names = {'ToAware', 'ToAware_hiv',...
        'To_aware_b', 'To_aware_h', 'To_aware_w',...
        'ToVax1', 'ToVax2', 'ToVax', 'ToVax1Plwh',...
        'NewInfections', 'newInfect_hiv', 'r_t'};

    % Initialize the struct
    metricShelf = struct();

    for i = 1:length(metrics_names)
        for j = iterations
            metricShelf.(metrics_names{i}).(['Iteration' num2str(j)]) = []; % Initialize each field
        end
    end

    % read in tally for each iteration
    % save all iteration results for each metric
    for i = iterations
        % tally path
        tally_path = fullfile(InPath, sprintf('iter%d', i), 'state_matrices', sprintf('Tally_%s.csv', Scenario_name));
        if ~isfile(tally_path)
            warning('File not found: %s. Skipping scenario.', tally_path);
            continue;
        end
        tally_table = readtable(tally_path);
        for metric_idx = 1:length(metrics_names) 
            metricName = metrics_names{metric_idx};
            metricVector = tally_table.(metricName);
            
            % Store the vector in the corresponding struct field and iteration
            metricShelf.(metricName).(['Iteration' num2str(i)]) = metricVector;
        end    
    end

    % Validate metric dimensions
    validateMetricDimensions(metricShelf, numWks);

    %% bootstrap
    % iterate through metrics
    for metric_idx = 1:length(metrics_names)
        metricName = metrics_names{metric_idx};
        resultShelf = zeros(bs, numWks+1); % +1 b/c weeks start at 0 and go to 85
        for i = 1:bs
            smpls_row = smpls(i,:); 
            tempResults = zeros(length(smpls_row), numWks+1);  % Temporary storage for each iteration's results within a single bootstrap
            
            for k = 1:length(smpls_row)
                j = smpls_row(k);
                % Retrieve the vector for the current iteration and metric
                tempResults(k, :) = metricShelf.(metricName).(['Iteration' num2str(j)]);     
            end

            % Calculate the mean across all selected iterations for this bootstrap
            resultShelf(i, :) = mean(tempResults, 1);  % Dimension 1 averages across the sampled iterations

        end

        % calculate the mean, LB, and UB of the bootstrapped vector
        % Mean across the bootstrap samples for each week
        bootstrapMeans = mean(resultShelf, 1);  % Mean across rows

        % Lower and upper bounds using percentiles
        bootstrapLB = prctile(resultShelf, 2.5, 1);  % 2.5 percentile across rows
        bootstrapUB = prctile(resultShelf, 97.5, 1);  % 97.5 percentile across rows

        % Convert results to a table
        resultsTable = table((0:numWks)', bootstrapMeans', bootstrapLB', bootstrapUB', ...
                              'VariableNames', {'Week', 'Mean', 'LowerBound', 'UpperBound'});

        % Create a filename based on the metric name
        filename = fullfile(OutPath, sprintf('%s_bs_results.csv', metricName));

        % Save the table to a CSV file
        try
            writetable(resultsTable, filename);
        catch ME
            warning('Failed to write results for metric %s: %s', metricName, ME.message);
            continue;
        end
        fprintf('Processing metric %d of %d: %s\n', metric_idx, length(metrics_names), metricName);
        fprintf('Saved %s\n', filename);
    end
end

function validateMetricDimensions(metricShelf, numWks)
    for metric = fieldnames(metricShelf)'
        for iter = fieldnames(metricShelf.(metric{1}))'
            if length(metricShelf.(metric{1}).(iter{1})) ~= numWks + 1
                error('Invalid dimension for metric %s, iteration %s', metric{1}, iter{1});
            end
        end
    end
end

