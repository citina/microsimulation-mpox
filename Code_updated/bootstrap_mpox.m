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
    % "2Xvax_mpox2024_S6",
    % '2Xvax_mpox2024_S21',
    % '2Xvax_mpox2024_S22',
    % '10Xvax_mpox2024_S6',
    % '10Xvax_mpox2024_S21',
    % '10Xvax_mpox2024_S22'
};
%% paths
% Base path for MonteCarloResults
basePath = pwd;
monteCarloPath = fullfile(fileparts(basePath), "MonteCarloResults");

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

    % Define metric configurations so cumulative variants can be added easily
    metrics_definitions = struct( ...
        'source', {'ToAware', 'ToAware_hiv',...
                   'To_aware_b', 'To_aware_h', 'To_aware_w',...
                   'ToVax1', 'ToVax2', 'ToVax', 'ToVax1Plwh',...
                   'NewInfections', 'newInfect_hiv', 'r_t',...
                   'ToAware'}, ...
        'output', {'ToAware', 'ToAware_hiv',...
                   'To_aware_b', 'To_aware_h', 'To_aware_w',...
                   'ToVax1', 'ToVax2', 'ToVax', 'ToVax1Plwh',...
                   'NewInfections', 'newInfect_hiv', 'r_t',...
                   'ToAwareCum'}, ...
        'isCumulative', {false, false,...
                         false, false, false,...
                         false, false, false, false,...
                         false, false, false,...
                         true});

    metrics_names = {metrics_definitions.output};

    % Initialize the struct
    metricShelf = struct();

    % read in tally for each iteration
    % save all iteration results for each metric
    loadedIterations = [];
    for i = iterations
        tally_path = fullfile(InPath, sprintf('iter%d', i), 'state_matrices', sprintf('Tally_%s.csv', Scenario_name));
        if ~isfile(tally_path)
            warning('File not found: %s. Skipping iteration %d for scenario %s.', tally_path, i, Scenario_name);
            continue;
        end
        tally_table = readtable(tally_path);
        metricDataTemp = struct();
        iterationHasAllMetrics = true;
        for metric_idx = 1:length(metrics_definitions)
            metricDef = metrics_definitions(metric_idx);
            if ~ismember(metricDef.source, tally_table.Properties.VariableNames)
                warning('Metric %s not found in %s. Skipping iteration %d.', metricDef.source, tally_path, i);
                iterationHasAllMetrics = false;
                break;
            end
            metricVector = tally_table.(metricDef.source);
            if metricDef.isCumulative
                metricVector = cumsum(metricVector);
            end
            metricDataTemp.(metricDef.output) = metricVector;
        end
        if ~iterationHasAllMetrics
            continue;
        end
        loadedIterations(end+1) = i; %#ok<SAGROW>
        for metricField = fieldnames(metricDataTemp)'
            metricName = metricField{1};
            if ~isfield(metricShelf, metricName)
                metricShelf.(metricName) = struct();
            end
            metricShelf.(metricName).(['Iteration' num2str(i)]) = metricDataTemp.(metricName);
        end
    end

    if isempty(loadedIterations)
        warning('No iterations loaded for scenario %s. Skipping.', Scenario_name);
        continue;
    end

    % a matrix to store indices of iterations actually available
    rng(666);
    numLoadedIterations = numel(loadedIterations);
    smpls = zeros(bs, numLoadedIterations);
    for i = 1:bs
        smpls(i,:) = randsample(loadedIterations, numLoadedIterations, true);
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

%% ========== PAIRWISE SCENARIO COMPARISON ==========
% Compare metrics between two scenarios via bootstrap
% Define pairs: each row is {scenarioA, scenarioB, metric_source, metric_output, isCumulative}
comparison_pairs = {
    % 'new_mpox2024_S6', '2Xvax_mpox2024_S6', 'ToAware', 'ToAware', false;
    'new_mpox2024_S6', 'new_mpox2024_S26', 'ToAware', 'ToAwareCum', true;
};

if ~isempty(comparison_pairs)
    fprintf('\n===== Pairwise Scenario Comparisons =====\n');
    for pairIdx = 1:size(comparison_pairs, 1)
        scenarioA = comparison_pairs{pairIdx, 1};
        scenarioB = comparison_pairs{pairIdx, 2};
        metricSource = comparison_pairs{pairIdx, 3};
        metricOutput = comparison_pairs{pairIdx, 4};
        isCumulative = comparison_pairs{pairIdx, 5};
        
        fprintf('\nComparing %s vs %s for metric %s\n', scenarioA, scenarioB, metricOutput);
        
        % Load data for both scenarios
        [dataA, iterA] = loadScenarioMetric(monteCarloPath, scenarioA, metricSource, isCumulative, iterations);
        [dataB, iterB] = loadScenarioMetric(monteCarloPath, scenarioB, metricSource, isCumulative, iterations);
        
        numIterA = numel(iterA);
        numIterB = numel(iterB);
        
        % Bootstrap the difference
        rng(666);
        diffShelf = zeros(bs, numWks+1);
        for b = 1:bs
            % Sample iterations with replacement for each scenario
            sampA = randsample(iterA, numIterA, true);
            sampB = randsample(iterB, numIterB, true);
            
            % Compute mean for each scenario in this bootstrap
            tempA = zeros(numIterA, numWks+1);
            tempB = zeros(numIterB, numWks+1);
            for k = 1:numIterA
                tempA(k,:) = dataA.(['Iteration' num2str(sampA(k))]);
            end
            for k = 1:numIterB
                tempB(k,:) = dataB.(['Iteration' num2str(sampB(k))]);
            end
            meanA = mean(tempA, 1);
            meanB = mean(tempB, 1);
            diffShelf(b,:) = meanA - meanB;
        end
        
        % Summarize
        diffMean = mean(diffShelf, 1);
        diffLB = prctile(diffShelf, 2.5, 1);
        diffUB = prctile(diffShelf, 97.5, 1);
        
        % Save results
        resultsTable = table((0:numWks)', diffMean', diffLB', diffUB', ...
            'VariableNames', {'Week', 'MeanDiff', 'LowerBound', 'UpperBound'});
        outDir = fullfile(monteCarloPath, 'comparisons');
        if ~isfolder(outDir), mkdir(outDir); end
        filename = fullfile(outDir, sprintf('%s_vs_%s_%s_diff.csv', scenarioA, scenarioB, metricOutput));
        writetable(resultsTable, filename);
        fprintf('Saved %s\n', filename);
    end
end

function [metricData, loadedIters] = loadScenarioMetric(basePath, scenarioName, metricSource, isCumulative, iterations)
    metricData = struct();
    loadedIters = [];
    scenarioPath = fullfile(basePath, scenarioName);
    for i = iterations
        tallyPath = fullfile(scenarioPath, sprintf('iter%d', i), 'state_matrices', sprintf('Tally_%s.csv', scenarioName));
        if ~isfile(tallyPath), continue; end
        tbl = readtable(tallyPath);
        vec = tbl.(metricSource);
        if isCumulative, vec = cumsum(vec); end
        metricData.(['Iteration' num2str(i)]) = vec;
        loadedIters(end+1) = i; %#ok<AGROW>
    end
end

