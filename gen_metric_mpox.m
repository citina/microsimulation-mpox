% clear all
% close all
% clear
% clc
% 
% %{
% The following code geneartes a file named result.mat
% 
% result.mat is a N by 1 matrix
% * N is the number of metrics generated
% * in each cell among the N cells, there stores a matrix
% * those matricies can be 2D or 3D denpend on if the metric is aggregated or by race or age
% * For each 2D matrix, number of rows equals to number of time cycles, number of columns equals to number of iterations
% * For each 3D matrix, number of rows equals to number of time cycles, number of columns equals to number of breakout buckets, 
% number of the 3rd dimension equals to the number of iterations
% * Each element represents the corresponding measurement as its file name
% %}
% 
% %% set inputs and outputs
% init_wk = 0;
% wksPassIn = (0:85);
% iterations = 1:10;
% 
% %% set paths
% % state matrices path
% dataFolder = "/mpox2024_S10";
% InPath = fileparts(pwd) + "/MonteCarloResults" + dataFolder; 
% 
% % initial population path
% initPop = "init_2024_new3.csv";
% 
% % metric matrix output path
% OutPath = fileparts(pwd) + "/MonteCarloResults" + dataFolder;
% 
% % memo
% iter = iterations(end);
% date = datetime('now');
% Tabl = table(iter, date);
% writetable(Tabl, OutPath+'/results_memo.txt', 'Delimiter', '|')
% 
% %% Load basic inputs state matrices
% % State matrices
% numWks = size(wksPassIn,2);
% wksKey = zeros(numWks, 2);
% 
% for wkIndex = 1:numWks
%     wksKey(wkIndex,:) = [wksPassIn(wkIndex), wkIndex-1];
% end
% 
% % each cell element holds EOW population data for a specified iteration
% % row: weeks, column: iterations, cell: state matrix for that week
% eowShelf = cell(numWks, length(iterations));  
% 
% % fill in state matrices
% for EOWindx = 1:numWks  %iterate over week index
%     for i = 1:length(iterations)
%         EOWnum = EOWindx-1;
%         % the path to the data must be modified by the iteration number and policy. 
%         dataPath = strcat(InPath, "/iter", num2str(i), "/state_matrices");
%         dataStruct = load(fullfile(dataPath, num2str(EOWnum)));
%         eowShelf{EOWindx, i} = dataStruct.state_matrix; %load all state matries
%     end
% end
% 
% % Initial population variable names
% inputfolder = pwd + "/input/";
% initialPopSpring = readtable(fullfile(inputfolder, initPop));
% EOWvarNames = initialPopSpring.Properties.VariableNames; %column names of init pop file
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % a shlef to store all results
% % each cell is the number for below metric of all years
% resultShelf = cell(4, 1);
% resultShelf_metrics_names = {'aware', 'aware_plwh', 'incidence', 'inc_plwh'};
% 
% agg_indx = 1:4;
%     %1: Mpox new diagnosed cases
%     %2: Mpox new diagnosed cases | PLWH 
%     %3: Mpox Incident cases
%     %4: Mpox Incident cases | PLWH
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % an empty matrix for each metric
% [new_diag, new_diag_plwh, incidence, inc_plwh] = deal(zeros(numWks,length(iterations)));
% 
% %% metrics calculation
% for iter = iterations 
%     for EOWindx = 0:numWks-1
% 
%         % read in a state matrix 
%         if EOWindx == 0
%             continue;
%         else
%             EOWmat_prev = eowShelf{EOWindx, iter};
%             EOWmat = eowShelf{EOWindx+1, iter};
%             alive = EOWmat(:,strcmp(EOWvarNames, 'alive'));
%             hiv_awareVec = EOWmat(:,strcmp(EOWvarNames, 'hiv_aware'));
%             hiv_statusVec = EOWmat(:,strcmp(EOWvarNames, 'hiv_status'));
%             pox_statusVec = EOWmat(:,strcmp(EOWvarNames, 'pox_status')); 
%             pox_awareVec = EOWmat(:,strcmp(EOWvarNames, 'pox_aware'));
%             pox_statusVec_prev = EOWmat_prev(:,strcmp(EOWvarNames, 'pox_status')); 
%             pox_awareVec_prev = EOWmat_prev(:,strcmp(EOWvarNames, 'pox_aware'));
%         end
% 
%         %% Mpox incidence
%         if EOWindx == 0 % the initial week before simulation starts
%             % indecies of people whose pox status == 0 the 1st week
%             not_infected_indx = find(pox_statusVec==0 & alive==1);
%             incidence(EOWindx+1, iter) = nan; % assign nan since no inc in the initial pop
%             inc_plwh(EOWindx+1, iter) = nan;
%         else
%             infected_indx = find((pox_statusVec==1 & alive==1) | (pox_statusVec==2 & alive==1));
%             inc_bool = (ismember(infected_indx, not_infected_indx));
%             incidence(EOWindx+1, iter) = sum(inc_bool);    
% 
%             % BY HIV aware status
%             infected_plwh_indx = find((pox_statusVec==1 & alive==1 & hiv_statusVec ~= 0) | (pox_statusVec==2 & alive==1 & hiv_statusVec ~= 0));
%             inc_plwh_bool = (ismember(infected_plwh_indx, not_infected_indx));
%             inc_plwh(EOWindx+1, iter) = sum(inc_plwh_bool);
% 
%             % indecies of people who are pox suspectible this year, will be used as comparison in the next year
%             not_infected_indx = find(pox_statusVec==0);
% 
%             % indecies of people who are mpox infected this year
%             prev_infected_indx = find(pox_statusVec==1 | pox_statusVec==2);
%         end
% 
%         %% New Mpox diagnosis % not yet solved
%         % new diagnosis: 
%             % 1. still infected: pox_status==1 or 2, pox_aware_prev==0, pox_aware_now==1
%             % 2. recovered: pox_status_prev==1 pox_aware_prev==0, pox_aware_now==0
%         if EOWindx == 0
%             % indecies of people who were infected but undiagnosed in the initial week
%             undiagnosed_indx = find(pox_awareVec==0 & (pox_statusVec==1 | pox_statusVec==2) & alive==1);
%             sus_indx = find(pox_statusVec==0 & alive==1);
%             new_diag(EOWindx+1, iter) = nan;
%             new_diag_plwh(EOWindx+1, iter) = nan;
%         else                   
%             % also need to count those people whose aware status been turned off bc recovered
%             diagnosed_indx = find(pox_awareVec==1 & alive==1); 
%             recovered_indx = find(pox_statusVec==3 & alive==1);
%             find((pox_statusVec==1 | pox_statusVec==2) & alive==1)
%             new_diag(EOWindx+1, iter) = sum(ismember(diagnosed_indx, undiagnosed_indx))+...
%                                         sum(ismember(recovered_indx, sus_indx));             
% 
%             % BY HIV aware status
%             diagnosed_plwh_indx = find(pox_awareVec==1 & hiv_awareVec==1);
%             new_diag_plwh_bool = (ismember(diagnosed_plwh_indx, undiagnosed_indx));        
%             new_diag_plwh(EOWindx+1, iter) = sum(new_diag_plwh_bool);
% 
%             undiagnosed_indx = find(pox_awareVec==0);
%             sus_indx = find(pox_statusVec~=0);
%         end      
%     end
% 
%     resultShelf{1} = new_diag;
%     resultShelf{2} = new_diag_plwh;
%     resultShelf{3} = incidence;
%     resultShelf{4} = inc_plwh;
% end
% 
% cd(OutPath)
% save result.mat resultShelf
% cd([fileparts(fileparts(pwd)),'/Code'])
% bootstrap_mpox;
% 
% % disp("you might want to run bootstrap next")