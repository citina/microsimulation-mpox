format long
clear all
clear
close all
clc

% things to check:
% 1. waning vaccine efficacy? waning_ve
% 2. time-varying foi? infectCalib_
% 3. foi values? infectCalib_
% 4. policy on? policy
% 5. new input people? starts from when? leaving? ends till when? infect_id

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
runSimulations = 1;
iterations = 1:10;
waning_ve = 2; % 0-waning ve off, 1-ve wanes to 0 after 12 months, 2-ve wanes to half after 12 months
scenarios = [15, 16, 17, 18, 19]; %[11, 15, 16, 17, 18, 19];

% Week when intervention starts
policy = 0; % turn on or off policy
p = '';
policy_wk = 100; % ploicy starting week, week0 Jun26-Jul02, week14 Oct2-Oct8
sens = 0; % turn on or off sens analysis

% defining the input file name
inputHeader = "input/Inputs_mpox2024_set2";
fileFormat = ".xlsx";
sim_inputFile = strcat(inputHeader, fileFormat);

% create list of output locations
workingDir = pwd;
dataDirHeader = fileparts(workingDir) + "/MonteCarloResults";
mkdir(dataDirHeader)

for scenario = scenarios
    
    % randomly select weeks to assign low FoI 
    rng(123)
    if scenario==17
        selected_wks = randsample(1:T, 22);
    elseif scenario==18
        selected_wks = randsample(1:T, 13);
    elseif scenario==19
        selected_wks = randsample(1:T, 9);
    end
 
    testVersion = strcat('mpox2024_S', num2str(scenario));% S11, S15, S16

    % update the complete dataDir header
    testVerDir = strcat(dataDirHeader,"/",testVersion);
    mkdir(testVerDir)

    if runSimulations == 1
        disp("*** Running Simulation ***")

        % start timer 
        tStart = tic;

        disp(sim_inputFile)

        for iters = iterations
            sim_dataDir = strcat(testVerDir, "/iter", num2str(iters), "/state_matrices/"); 
            mkdir(sim_dataDir)

            rng('shuffle')
            disp(strcat("Iteration Number: ", num2str(find(iterations == iters)),"/", num2str(length(iterations))))
            mpox2024_shellMod4;
        end

        gen_metric;

        % end timer 
        tEnd = toc(tStart);

        fprintf('Simulation time %d hours and %d minutes\n', floor(tEnd/3600), floor(rem(tEnd,3600)/60));            
    end   
    cd(workingDir)
end