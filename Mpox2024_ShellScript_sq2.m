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
testing_wk = 12; % the week added "event" people leave LAC, not used if "added" column not used
waning_ve = 2; % 0-waning ve off, 1-ve wanes to 0 after 12 months, 2-ve wanes to half after 12 months
weekly_foi = 0; % 1-if weekly foi, 0-if average foi
foi = 2.2;
piso_input = 0.8;%[0, 0.2, 0.3, 0.4]; %S7, S6, S5, S4, S14

% Week when intervention starts
policy = 0; % turn on or off policy
p = '';
policy_wk = 100; % ploicy starting week, week0 Jun26-Jul02, week14 Oct2-Oct8
sens = 0; % turn on or off sens analysis

% defining the input file name
inputHeader = "../input/Inputs_mpox2024_set2";
fileFormat = ".xlsx";
sim_inputFile = strcat(inputHeader, fileFormat);

% create list of output locations
workingDir = pwd;
dataDirHeader = fileparts(workingDir) + "/MonteCarloResults";
mkdir(dataDirHeader)

for piso = piso_input
    for tweek = testing_wk %%debug
        
        testVersion = strcat('mpox2024_foi_', num2str(foi), 'iso', num2str(piso));%, string(tweek));  %'mpox_'
    
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
                mpox2024_shellMod2;
            end
    
            gen_metric;
    
            % end timer 
            tEnd = toc(tStart);
    
            fprintf('Simulation time %d hours and %d minutes\n', floor(tEnd/3600), floor(rem(tEnd,3600)/60));            
        end   
    end
    cd(workingDir)
end