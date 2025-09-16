% format long
% clear all
% close all
% clc
% testVersion = 'mpox2024_S0';
% iterations = 1:10;
% numWks = 85;

%% set paths (provided by shell script)
global NUM_ITERATIONS testVerDir testVersion
InPath = testVerDir; 
OutPath = testVerDir; 

%% Load tally
% each cell element holds tally data for a specified iteration
talllyShelf = [];  

% import the tally.csv
for i = 1:NUM_ITERATIONS
    % the path to the data must be modified by the iteration number and policy. 
    dataPath = fullfile(InPath, sprintf('iter%d/state_matrices/Tally_%s.csv', i, testVersion));
    if ~exist(dataPath,'file')
        error('gen_metric:MissingTally', 'Missing tally file: %s', dataPath);
    end
    dataStruct = readtable(dataPath, 'PreserveVariableNames',true);
    talllyShelf(:,:,i) = dataStruct.Variables;
end

%% Average metrics
n_wks = size(dataStruct,1);
n_metrics = size(dataStruct,2);
metric_names = dataStruct.Properties.VariableNames;

cd(OutPath)
AvgTally = round(mean(talllyShelf,3), 2);
avgtally_tbl = array2table(AvgTally);
avgtally_tbl.Properties.VariableNames = metric_names;
writetable(avgtally_tbl, ['AvgTally_',char(testVersion),'.csv'])

% cd(OutPath)
% save AvgTally.mat resultShelf
% bootstrap;
% calibration_graph;











