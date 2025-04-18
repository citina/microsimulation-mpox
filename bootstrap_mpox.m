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
iterations = 1:10;
numWks = 85;
Scenario_name = "mpox2024_foi_2iso0.2";

% a matrix to store indecies of iterations
% e.g.: if 3 iterations, then randomly select 3 samples with replacement
rng(666);
smpls = zeros(bs, length(iterations));
for i = 1:bs
    smpls(i,:) = randsample(iterations, length(iterations), true);
end

%% paths
% Tally path of each iteration
InPath = fileparts(pwd) + "/MonteCarloResults/" + Scenario_name;
% InPath = fileparts(pwd) + "/create files/" + Scenario_name;

% metric matrix output path
% OutPath = fileparts(pwd) + "/MonteCarloResults/" + Scenario_name;
OutPath = InPath;

% memo for bootstrap
Tabl = table(length(iterations), bs, datetime('now'));
writetable(Tabl, OutPath+'/memo.txt', 'Delimiter','|')

% metricShelf is a struct where each field is a metric name, and each
% metric contains subfields for each iteration.
% Metrics of interest
metrics_names = {'ToAware', 'ToAware_hiv',...
    'To_aware_b', 'To_aware_h', 'To_aware_w',...
    'ToVax1', 'ToVax2', 'ToVax', 'ToVax1Plwh',...
    'NewInfections', 'newInfect_hiv', 'r_t', 'infected_per_infectious'};

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
    tally_path = InPath + '/iter' + num2str(i) + '/state_matrices' + '/Tally_' + Scenario_name + '.csv';
    tally_table = readtable(tally_path);
    for metric_idx = 1:length(metrics_names) 
        metricName = metrics_names{metric_idx};
        metricVector = tally_table.(metricName);
        
        % Store the vector in the corresponding struct field and iteration
        metricShelf.(metricName).(['Iteration' num2str(i)]) = metricVector;
    end    
end

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
    filename = sprintf('/%s_bs_results.csv', metricName);

    % Save the table to a CSV file
    writetable(resultsTable, OutPath + filename);
    fprintf('Saved %s\n', filename);

end

