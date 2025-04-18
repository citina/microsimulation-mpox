
% Input parameters
mpox2024_parameters;

if sens
    asym2sym_transition_path = asym2sym_transition_path_9;
end

% the state matrix directory is different. It should also be listed as a char not a string
state_matrices_path = convertStringsToChars(sim_dataDir);

% Read in the initial population from a csv
% Assumption: the csv file has all the people pre-populated with characteristics
[state_matrix, StateMatCols] = read_table(init_pop_file);

% Read in all possible demographic groups
% age min and max MUST be listed first in excel sheet
[demog_table, DemogTblCols] = create_demog_groups(demog_var_def_file);

% Read in demographic groups specified in the mixing matrix
% age min and max MUST be listed first in excel sheet
[mixing_table, MixingTblCols] = create_demog_groups(mixing_mat_def_file);

% Read in the mixing matrix
mixing_matrix = importdata(mixing_mat_file);

% read in number of partners vector
numPartners = importdata(num_partners_file);

% an empty matrix to store tally
finalTransitionTally = [];

% genereate final output from state matrices
workingDir_home = pwd;

% clear current matrices in folder
cd(sim_dataDir);
delete *.mat
cd(workingDir_home)

% save iniital population as a state matrix
matrix_name = strcat(state_matrices_path, int2str(0));
save(matrix_name, 'state_matrix');

tallyHeaders = {'Week', 'New Births', 'New Infections', 'infected_per_infectious',...
    'r_t', 'max(infectious_wks)'...
    'new infect|hiv', 'new infect|vax1', 'new infect|vax2',...
    'new infect b', 'new infect h', 'new infect w',...
    'new infect a1', 'new infect a2', 'new infect a3',...
    'new infect a4', 'new infect a5',...
    'To Vax 1', 'To Vax 2', 'To Vax', 'To Vax1 plwh',...
    'To Aware (asym)', 'To Aware (sym)', 'To Aware', 'To aware|hiv',...
    'To Aware b|hiv', 'To Aware h|hiv', 'To Aware w|hiv',...
    'To Aware ls 40|hiv', 'To Aware ge 40|hiv',...
    'To_aware_b', 'To_aware_h', 'To_aware_w',...
    'To_aware_a1', 'To_aware_a2', 'To_aware_a3', 'To_aware_a4', 'To_aware_a5',...
    'To_vax1_b', 'To_vax1_h', 'To_vax1_w',...
    'To_vax1_a1', 'To_vax1_a2', 'To_vax1_a3', 'To_vax1_a4', 'To_vax1_a5',... 
    'To_vax2_b', 'To_vax2_h', 'To_vax2_w',...
    'To_vax2_a1', 'To_vax2_a2', 'To_vax2_a3', 'To_vax2_a4', 'To_vax2_a5',...
    'To_vax_b', 'To_vax_h', 'To_vax_w',...
    'To_vax_a1', 'To_vax_a2', 'To_vax_a3', 'To_vax_a4', 'To_vax_a5',...
    'Asymp to Symp', 'Symp to Recover', 'Asymp to Recover',...
    'To on Treatment', 'To Death', 'To isolation',...
    'Alive MSM', 'Sus', 'Asymp', 'Symp', 'recovered',...
    'aware', 'aware_asym', 'aware_sym',...
    'vax_1', 'vax_2',...
    'vac1_sus', 'vac1_asy', 'vac1_sym', 'vac1_rec',...
    'vac2_sus', 'vac2_asy', 'vac2_sym', 'vac2_rec',...
    'vac1_b', 'vac1_h', 'vac1_w',...
    'vac2_b', 'vac2_h', 'vac2_w',...
    'vac1_hiv', 'vac2_hiv',...
    'on trt', 'on iso', 'hiv aware', 'pox + hiv aware'};
tally = zeros(T+1, length(tallyHeaders)); % +1 cuz include the initial wk as 0

% tally for the initial state matrix
alive = state_matrix(:, StateMatCols.alive);
pox_status = state_matrix(:, StateMatCols.pox_status);
pox_aware = state_matrix(:, StateMatCols.pox_aware);
vaccinated = state_matrix(:, StateMatCols.vaccinated);
treatment = state_matrix(:, StateMatCols.treatment);
isolation = state_matrix(:, StateMatCols.isolation);
race = state_matrix(:, StateMatCols.race);
hiv_aware = state_matrix(:, StateMatCols.hiv_aware);
hiv_status = state_matrix(:, StateMatCols.hiv_status);
vax_wk = state_matrix(:, StateMatCols.vax_wk);
ve = state_matrix(:, StateMatCols.ve);
infectious_wks = state_matrix(:, StateMatCols.infectious_wks);

tally(1,:) = [zeros(1,68),...
              sum(alive),...
              sum(pox_status==0), sum(pox_status==1),...
              sum(pox_status==2), sum(pox_status==3),...
              sum(pox_aware==1),...
              sum(pox_aware==1 & pox_status==1), sum(pox_aware==1 & pox_status==2),...
              sum(vaccinated==1), sum(vaccinated==2),...
              sum(vaccinated==1 & pox_status==0), sum(vaccinated==1 & pox_status==1), sum(vaccinated==1 & pox_status==2), sum(vaccinated==1 & pox_status==3),...
              sum(vaccinated==2 & pox_status==0), sum(vaccinated==2 & pox_status==1), sum(vaccinated==2 & pox_status==2), sum(vaccinated==2 & pox_status==3),...
              sum(vaccinated==1 & race==0), sum(vaccinated==1 & race==1), sum(vaccinated==1 & race==2),...
              sum(vaccinated==2 & race==0), sum(vaccinated==2 & race==1), sum(vaccinated==2 & race==2),...
              sum(vaccinated==1 & hiv_status~=0), sum(vaccinated==2 & hiv_status~=0),...
              sum(treatment==1), sum(isolation==1), sum(hiv_aware==1), sum(hiv_aware==1 & pox_aware==1)];      

%%%%% 2024 set2 updates %%%%%
% weekly input schedule: the week and its corresponding number of imports
import_wk = [1 3 5 7 8 12 13 15 16 20 21 22 23 24 26 27 28 30 31 32 33 35 37 39 41 45 46 50 52 53 54 57 59 60 61 63 64 65 67 68 69 70];
import_num = [1 1 2 2 1 2 2 1 1 3 1 2 1 2 2 1 3 1 2 4 2 2 1 1 1 1 1 1 1 1 2 3 1 1 2 1 1 1 1 1 1 1];
% import_num = import_num*2;

% For each week
for t = 1:T
    
    [infectionTally, vac1Tally, vac2Tally, awareTally1, awareTally2,...
     asym2symTally, sym2recoverTally, asym2recoverTally, ontrtTally, deathTally,...
     onisoTally] = deal(0);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% update ve  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     # immunocompetent: one dose ve = 0.721, two doses ve = 0.878
%     # immunocompromised: one dose ve = 0.51, two doses ve = 0.702
%     # reach one dose ve after 2 weeks of 1st dose and reach full ve after 2 weeks of 2nd dose
%     # reference paper link: https://www.nejm.org/doi/full/10.1056/NEJMoa1817307
%     # Assumption: ve follows a piecewise function (x is # of weeks 1st dose vaccinated , y is ve):
%
%     # For immunocompetent individuals
%     # week 0-2: y = 0.3605x
%     # week 2-4: y = 0.721
%     # --> if didn't get 2nd dose, ve drops from dose 1 ve to 0 after 52 weeks: y = -0.014x+0.777
%     # --> if got 2nd dose, ve develop to dose 2 level between 4-6: y = 0.079x+0.405				
%     #     --> ve gradually drops to 0 after 52 weeks: y = -0.017x + 0.986
%
%     # For immunocompromised individuals - 39% of PLWH (from 2022 LAC surveillance)
%     # week 0-2: y = 0.255x
%     # week 2-4: y = 0.51
%     # --> if didn't get 2nd dose, ve drops from dose 1 ve to 0 after 52 weeks: y = -0.01x + 0.56
%     # --> if got 2nd dose, ve develop to dose 2 level between 4-6: y = 0.096x + 0.126				
%     #     --> ve gradually drops to 0 after 52 weeks: y = -0.014x + 0.786			
    
    if waning_ve == 0  % no waning ve       
        plwh_vax1_id = hiv_status~=0 & vaccinated==1 & vax_wk>=2 & alive==1; 
        plwh_vax2_id = hiv_status~=0 & vaccinated==2 & vax_wk>=6  & alive==1;
        state_matrix(plwh_vax1_id, StateMatCols.ve) = vac1_plwh;
        state_matrix(plwh_vax2_id, StateMatCols.ve) = vac2_plwh;
        
        normal_vax1_id = hiv_status==0 & vaccinated==1 & vax_wk>=2 & alive==1;
        normal_vax2_id = hiv_status==0 & vaccinated==2 & vax_wk>=6 & alive==1;
        state_matrix(normal_vax1_id, StateMatCols.ve) = vac1_normal;        
        state_matrix(normal_vax2_id, StateMatCols.ve) = vac2_normal;
    end

    if waning_ve ~= 0    
        % for immunocompetent individuals (now using plwh)
        idx1 = hiv_status==0 & vaccinated==1 & vax_wk>=0 & vax_wk<=2 & alive==1;
        idx2 = hiv_status==0 & vaccinated==1 & vax_wk>2 & vax_wk<=4 & alive==1;
        idx3 = hiv_status==0 & vaccinated==1 & vax_wk>4 & alive==1; % if not getting 2nd dose
        idx4 = hiv_status==0 & vaccinated==2 & vax_wk>=4 & vax_wk<=6 & alive==1; % if got 2nd dose

        state_matrix(idx1, StateMatCols.ve) = 0.3605*vax_wk(idx1);
        state_matrix(idx2, StateMatCols.ve) = 0.721;
        state_matrix(idx3, StateMatCols.ve) = max(round(-0.014*vax_wk(idx3)+0.777, 3), 0); 
        state_matrix(idx4, StateMatCols.ve) = max(round(0.079*vax_wk(idx4)+0.405, 3), 0);	

        % for PLWH
        idx1_plwh = hiv_status~=0 & vaccinated==1 & vax_wk>=0 & vax_wk<=2 & alive==1;
        idx2_plwh = hiv_status~=0 & vaccinated==1 & vax_wk>2 & vax_wk<=4 & alive==1;
        idx3_plwh = hiv_status~=0 & vaccinated==1 & vax_wk>4 & alive==1; % if not getting 2nd dose
        idx4_plwh = hiv_status~=0 & vaccinated==2 & vax_wk>=4 & vax_wk<=6 & alive==1; % if got 2nd dose
        
        state_matrix(idx1_plwh, StateMatCols.ve) = 0.255*vax_wk(idx1_plwh);        
        state_matrix(idx2_plwh, StateMatCols.ve) = 0.51;
        state_matrix(idx3_plwh, StateMatCols.ve) = max(round(-0.01*vax_wk(idx3_plwh)+0.56, 3), 0); 
        state_matrix(idx4_plwh, StateMatCols.ve) = max(round(0.096*vax_wk(idx4_plwh)+0.126, 3), 0);	
        
        if waning_ve == 2 % remian half of full ve after 12 months
            idx5 = hiv_status==0 & vaccinated==2 & vax_wk>6 & vax_wk<=58 & alive==1;
            idx6 = hiv_status==0 & vaccinated==2 & vax_wk>58 & alive==1;
            idx7 = vaccinated==1 & vax_wk>56 & hiv_status==0 & alive==1;
            idx5_plwh = hiv_status~=0 & vaccinated==2 & vax_wk>6 & vax_wk<=58 & alive==1;
            idx6_plwh = hiv_status~=0 & vaccinated==2 & vax_wk>58 & alive==1;
            idx7_plwh = vaccinated==1 & vax_wk>56 & hiv_status~=0 & alive==1;
            
            state_matrix(idx5, StateMatCols.ve) = max(round(-0.017*vax_wk(idx5)+0.986, 3),0);            
            state_matrix(idx6, StateMatCols.ve) = 0.439;         
            state_matrix(idx7, StateMatCols.ve) = 0.3605;
            state_matrix(idx5_plwh, StateMatCols.ve) = max(round(-0.014*vax_wk(idx5_plwh) + 0.786, 3),0);
            state_matrix(idx6_plwh, StateMatCols.ve) = 0.351;
            state_matrix(idx7_plwh, StateMatCols.ve) = 0.255;
            
        elseif waning_ve == 1 % ve wanes to 0 after 12 months
            idx5 = vaccinated==2 & vax_wk>6 & hiv_status==0 & alive==1;
            idx6 = vaccinated==2 & vax_wk>58 & hiv_status==0 & alive==1;
            idx7 = vaccinated==1 & vax_wk>56 & hiv_status==0 & alive==1;
            idx5_plwh = vaccinated==2 & vax_wk>6 & hiv_status~=0 & alive==1;
            idx6_plwh = vaccinated==2 & vax_wk>58 & hiv_status~=0 & alive==1;
            idx7_plwh = vaccinated==1 & vax_wk>56 & hiv_status~=0 & alive==1;

            state_matrix(idx5, StateMatCols.ve) = max(round(-0.017*vax_wk(idx5)+0.986, 3),0);      
            state_matrix(idx6, StateMatCols.ve) = 0;   
            state_matrix(idx7, StateMatCols.ve) = 0; 
            state_matrix(idx5_plwh, StateMatCols.ve) = max(round(-0.014*vax_wk(idx5_plwh) + 0.786, 3),0);		         
            state_matrix(idx6_plwh, StateMatCols.ve) = 0;
            state_matrix(idx7_plwh, StateMatCols.ve) = 0;
        end
    end
    
    %%%%%%%%%%%%%% mpox 2024: assign time varying params %%%%%%%%%%%%%%%%%    
    infectCalib_ = 2.2;

    % scenario S11, S15, S16, S17, S18, S19
    if scenario==11 % S11: week 4 – 16 (Apr-Jun) and week 26 – 34 (Sept-Oct), FOI = 0.7 
        if t>=4 && t<=16
            infectCalib_ = 0.7;
        elseif t>=26 && t<=34
            infectCalib_ = 0.7; 
        end
    elseif scenario==15 % S15: For week 4 – 16 (Apr-Jun), FOI = 0.7
        if t>=4 && t<=16
            infectCalib_ = 0.7;
        end
    elseif scenario==16 % S16: For week 26 – 34 (Sept-Oct), FOI = 0.7
        if t>=26 && t<=34
            infectCalib_ = 0.7; 
        end
    elseif scenario==17 % S17: Randomly select 22 weeks with FoI = 0.7
        if ismember(t, selected_wks)
            infectCalib_ = 0.7; 
        end
    elseif scenario==18 % S18: Randomly select 13 weeks with FoI = 0.7
        if ismember(t, selected_wks)
            infectCalib_ = 0.7; 
        end
    elseif scenario==19 % S19: Randomly select 9 weeks with FoI = 0.7
        if ismember(t, selected_wks)
            infectCalib_ = 0.7; 
        end
    end
    
    %%%%%%%%%%%%% 1. Update vax_wk and vax2_wk if vaccinated and alive %%%%%%%%%%%%%%%
    % for people die at the end of the week vax_wk is updated here before they die
    % for new birth this week, their vax_wk is 0
    vax_wk_id = find(alive & vaccinated ~= 0);
    state_matrix(vax_wk_id, StateMatCols.vax_wk) = state_matrix(vax_wk_id, StateMatCols.vax_wk)+1;
    
    vax2_wk_id = find(alive & vaccinated == 2);
    state_matrix(vax2_wk_id, StateMatCols.vax2_wk) = state_matrix(vax2_wk_id, StateMatCols.vax2_wk)+1;
 
    %%%%%%%%%%%%%%%% 2024 update: infection every week %%%%%%%%%%%%%%%%%%%%%
   % % converting X people from sus to symptomatic every week starting from Jul 30 2023
   %  infect_id = find(alive & pox_status==0);
   %  infect_id = randsample(infect_id, 15);
   %  state_matrix(infect_id, StateMatCols.pox_status) = 2;
   % state_matrix(infect_id, StateMatCols.pox_aware) = 1;

    % converting X people from sus to symptomatic according to import_wk and import_num
    for wk = import_wk
        if t == wk
            tmp_i = find(t == import_wk);
            infect_id = find(alive==1 & pox_status==0);
            infect_id = randsample(infect_id, import_num(tmp_i));
            state_matrix(infect_id, StateMatCols.pox_status) = 2;
            state_matrix(infect_id, StateMatCols.pox_aware) = 1;
        end
    end 
    
    alive_beginning = state_matrix(:, StateMatCols.alive);
    pox_status_beginning = state_matrix(:, StateMatCols.pox_status);
    %%%%%%%%%%%%%%%%%%%%%%%%%% 2. New births %%%%%%%%%%%%%%%%%%%%%%%%
    % a bunch of 14 year-old enter the model
    [state_matrix, birthTally] = birth2(state_matrix, StateMatCols, age_def, race_def, race_prop, inflow);
    
    % update columns from the state matrix
    alive = state_matrix(:, StateMatCols.alive);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% 3. INFECTION %%%%%%%%%%%%%%%%%%%%%%%%%%
    % create a state matrix copy to do infections because we want newly infected people 
    % start to infect others in the next time period (not this time period).
    state_matrix_bfInfected = state_matrix;
    
    % For each demographic group
    for demog_group_def = demog_table'
        
        % Get state matrix row indices for people in this demographic group    
        demog_group_idx = find_demog_rows2(state_matrix, StateMatCols, demog_group_def, DemogTblCols);

        % If no one in the state matrix in this demographic group, continue to the next
        if isempty(demog_group_idx)
            continue
        end

        [state_matrix, infectionTallytmp] = infection9(demog_group_def, DemogTblCols,...
                                    mixing_matrix, mixing_table, ...
                                    MixingTblCols, state_matrix, StateMatCols, ...
                                    numPartners, infectCalib_, hivInfect,...
                                    youngInfect, midInfect,...
                                    oldInfect, raceInfect, ...
                                    isolation_adherence);

        % tally for all the demographic groups summed for this week
        infectionTally = infectionTally + infectionTallytmp;
    end
    
    % people newly to infection
    new_infect = unique(find(state_matrix_bfInfected(:,StateMatCols.pox_status)~=state_matrix(:, StateMatCols.pox_status)));
    
    % number of people do not infect any new people


    % number of people in each group turn to infected
    infect_hiv = sum(state_matrix(new_infect, StateMatCols.hiv_status)~=0);
    
    infect_vax1 = sum(state_matrix(new_infect, StateMatCols.vaccinated)==1);
    infect_vax2 = sum(state_matrix(new_infect, StateMatCols.vaccinated)==2);
    
    infect_b = sum(state_matrix(new_infect, StateMatCols.race)==0);
    infect_h = sum(state_matrix(new_infect, StateMatCols.race)==1);
    infect_w = sum(state_matrix(new_infect, StateMatCols.race)==2);
    
    infect_a1 = sum(state_matrix(new_infect, StateMatCols.age_bucket)==1);
    infect_a2 = sum(state_matrix(new_infect, StateMatCols.age_bucket)==2);
    infect_a3 = sum(state_matrix(new_infect, StateMatCols.age_bucket)==3);
    infect_a4 = sum(state_matrix(new_infect, StateMatCols.age_bucket)==4);
    infect_a5 = sum(state_matrix(new_infect, StateMatCols.age_bucket)==5);
    
    % save the state matrix before aware to track aware and vax1 by race
    sm_before_aware = state_matrix;
    sm_before_vax1 = state_matrix;
    
    %%%%%%%%%%%%%%%%%%%% vaccination under policies %%%%%%%%%%%%%%%%%%%%%% 
    % 95% of people aware of their pox status and asymp will be prioritized getting vax1
    vax1_pri = find(state_matrix(:,StateMatCols.alive)==1 &...
                 state_matrix(:,StateMatCols.pox_status)==1 &...
                 state_matrix(:,StateMatCols.vaccinated)==0 &...
                 state_matrix(:,StateMatCols.pox_aware)==1);
                 
    if strcmp(p, 'vax12') % fixed vax size, all for PLWH
        % find those eligible for 1st dose        
        vax1_eligible_Ind = find(state_matrix(:,StateMatCols.alive)==1 &...
                     state_matrix(:,StateMatCols.hiv_aware)==1 &...
                     state_matrix(:,StateMatCols.vaccinated)==0 &...
                     state_matrix(:,StateMatCols.pox_aware)==0);
        vax1_other_eligible_Ind = find(state_matrix(:,StateMatCols.alive)==1 &...
                     state_matrix(:,StateMatCols.hiv_aware)~=1 &...
                     state_matrix(:,StateMatCols.vaccinated)==0 &...
                     state_matrix(:,StateMatCols.pox_aware)==0);

    elseif strcmp(p, 'vax13') % fixed vax size, all for Black
        % find those eligible for 1st dose
        % note that here we don't consider people aware of their pox status
        vax1_eligible_Ind = find(state_matrix(:,StateMatCols.alive)==1 &...
                     state_matrix(:,StateMatCols.race)==0 &...
                     state_matrix(:,StateMatCols.vaccinated)==0 &...
                     state_matrix(:,StateMatCols.pox_aware)==0);
        vax1_other_eligible_Ind = find(state_matrix(:,StateMatCols.alive)==1 &...
                     state_matrix(:,StateMatCols.race)~=0 &...
                     state_matrix(:,StateMatCols.vaccinated)==0 &...
                     state_matrix(:,StateMatCols.pox_aware)==0);

    elseif strcmp(p, 'vax14') % fixed vax size, all for Hispanic
        % find those eligible for 1st dose
        % note that here we don't consider people aware of their pox status
        vax1_eligible_Ind = find(state_matrix(:,StateMatCols.alive)==1 &...
                     state_matrix(:,StateMatCols.race)==1 &...
                     state_matrix(:,StateMatCols.vaccinated)==0 &...
                     state_matrix(:,StateMatCols.pox_aware)==0);
        vax1_other_eligible_Ind = find(state_matrix(:,StateMatCols.alive)==1 &...
                     state_matrix(:,StateMatCols.race)~=1 &...
                     state_matrix(:,StateMatCols.vaccinated)==0 &...
                     state_matrix(:,StateMatCols.pox_aware)==0);

    elseif strcmp(p, 'vax15') % fixed vax size, all for White
        % find those eligible for 1st dose
        % note that here we don't consider people aware of their pox status
        vax1_eligible_Ind = find(state_matrix(:,StateMatCols.alive)==1 &...
                     state_matrix(:,StateMatCols.race)==2 &...
                     state_matrix(:,StateMatCols.vaccinated)==0 &...
                     state_matrix(:,StateMatCols.pox_aware)==0);
        vax1_other_eligible_Ind = find(state_matrix(:,StateMatCols.alive)==1 &...
                     state_matrix(:,StateMatCols.race)~=2 &...
                     state_matrix(:,StateMatCols.vaccinated)==0 &...
                     state_matrix(:,StateMatCols.pox_aware)==0);
    end

    % run the following block only if one of above vax policies selected
    if sum(strcmp(p, {'vax12', 'vax13', 'vax14', 'vax15'}))>0 
        
        % prioritize people asypm and aware their pox status to get vax1
        vax1_pri_size = round(length(vax1_pri)*0.95);    
        % if the pri size is greater than the fixed size, then set it
        % equals to the fixed size since it should not exceed it
        if vax1_pri_size > vax1_size(t)  
            vax1_pri_size = vax1_size(t);
        end
        vax1_pri_ind = randsample(vax1_pri, vax1_pri_size);
        
        % update fixed vax1 size at time t for others
        vax1_size(t) = vax1_size(t) - vax1_pri_size;
        
        % randomly select people to be vaccinated (1st dose)
        % vax size is the same as in sq
        size_diff = vax1_size(t) - length(vax1_eligible_Ind);    
        % if more fixed vax1 left, equally distribute them to other groups
        if size_diff > 0 
            vax1_Ind = [vax1_pri_ind;...
                        vax1_eligible_Ind;...
                        randsample(vax1_other_eligible_Ind, size_diff)]; 
        else
            vax1_Ind = [vax1_pri_ind;...
                        randsample(vax1_eligible_Ind, vax1_size(t))];                    
        end

        state_matrix(vax1_Ind, StateMatCols.vaccinated) = 1;
        % find those eligible for 2nd dose
        vax2_eligible_Ind = find(state_matrix(:,StateMatCols.alive)==1 &...
                     state_matrix(:,StateMatCols.vax_wk)>=4 &...
                     state_matrix(:,StateMatCols.vax_wk)<=6 &...
                     state_matrix(:,StateMatCols.vaccinated)==1 &...
                     state_matrix(:,StateMatCols.pox_aware)==0);
        % randomly select people to be vaccinated (12nd dose)
        % vax size is the same as in sq
        vax2_Ind = randsample(vax2_eligible_Ind, vax2_size(t));            
        state_matrix(vax2_Ind, StateMatCols.vaccinated) = 2;

        % tally for all the demographic groups summed for this week
        vac1Tally = vac1Tally + length(vax1_Ind);
        vac2Tally = vac2Tally + length(vax2_Ind);
        vacTally = vac1Tally+ vac2Tally;
    end
    
    %%%%%%%%%%%%%%%%%%%% For each demographic group %%%%%%%%%%%%%%%%%%%%%%
    for demog_group_def = demog_table.'
        
        % Get state matrix row indices for people in this demographic group
        demog_group_idx = find_demog_rows2(state_matrix, StateMatCols, demog_group_def, DemogTblCols);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%% 4. aware %%%%%%%%%%%%%%%%%%%%%%%%%%        
        [state_matrix, awareTally1tmp] = transition5('pox_aware', aware1_transition_path, ...
                                    state_matrix, StateMatCols, demog_group_idx, ...
                                    future_state_param_names, 1);
        [state_matrix, awareTally2tmp] = transition5('pox_aware', aware2_transition_path, ...
                                    state_matrix, StateMatCols, demog_group_idx, ...
                                    future_state_param_names, 1);
               
        % tally for all the demographic groups summed for this week
        awareTally1 = awareTally1 + awareTally1tmp;     
        awareTally2 = awareTally2 + awareTally2tmp;  

        %%%%%%%%%%%%%%%%%%%%%%%% 5. vaccination %%%%%%%%%%%%%%%%%%%%%%%%
        % mpox vaccine requires 2 doses, vaccinated = 1 or 2        
        % run this block only if one of the above vax policies (12-15) are not selected
        if ~strcmp(p, {'vax12', 'vax13', 'vax14', 'vax15'})  % regular situation                
            % those sus or asymptomatic can be vaccinated (1st dose)
            % the aware asym and unaware asym have different prob of vaccinated
            [state_matrix, vac1Tallytmp] = transition5('vaccinated', vac1_transition_path, ...
                                        state_matrix, StateMatCols, demog_group_idx, ...
                                        future_state_param_names, 1);

            % people get 1st dose will get their 2nd dose on the 4th week       
            [state_matrix, vac2Tallytmp] = transition5('vaccinated', vac2_transition_path, ...
                                        state_matrix, StateMatCols, demog_group_idx, ...
                                        future_state_param_names, 1);

            % tally for all the demographic groups summed for this week
            vac1Tally = vac1Tally + vac1Tallytmp;
            vac2Tally = vac2Tally + vac2Tallytmp;
            vacTally = vac1Tally+ vac2Tally;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%% 6. pox status %%%%%%%%%%%%%%%%%%%%%%%%
        % transit from asymptomatic to recover if vaccinated within a week
        % only people infected and vaccinated this week have chances to recover
        % for these people, vaccinated=1 and vax_wk=0
        [state_matrix, asym2recoverTallytmp] = transition5('pox_status', asym2recover_transition_path, ...
                                state_matrix, StateMatCols, demog_group_idx, ...
                                future_state_param_names, 1);
                            
        asym2recoverTally = asym2recoverTally + asym2recoverTallytmp;
        
        % find indecies in the state matrix for those asym and sym
        asymInd = find_indices(state_matrix, 1:size(state_matrix,1), StateMatCols.pox_status, '=', 1);
        symInd = find_indices(state_matrix, 1:size(state_matrix,1), StateMatCols.pox_status, '=', 2);

        % create state matrix based on pox status
        asym_state_matrix = state_matrix(asymInd,:);
        sym_state_matrix = state_matrix(symInd,:);

        % create demographic group index
        asymInd_demog_idx = find_demog_rows2(asym_state_matrix, StateMatCols, demog_group_def, DemogTblCols);
        symInd_demog_idx = find_demog_rows2(sym_state_matrix, StateMatCols, demog_group_def, DemogTblCols);
             
        % transit from asymptomatic to symptomatic
        [asym_state_matrix, asym2symTallytmp] = transition5('pox_status', asym2sym_transition_path, ...
                                asym_state_matrix, StateMatCols, asymInd_demog_idx, ...
                                future_state_param_names, 1);
        % redefine the state matrix based on the index of the new matrix
        state_matrix(asymInd, StateMatCols.pox_status) = asym_state_matrix(:, StateMatCols.pox_status);      

        % transit from symptomatic to recover (with / without tretament)
        [sym_state_matrix, sym2recoverTallytmp] = transition5('pox_status', sym2recover_transition_path, ...
                                sym_state_matrix, StateMatCols, symInd_demog_idx, ...
                                future_state_param_names, 1);
        % redefine the state matrix based on the index of the new matrix
        state_matrix(symInd, StateMatCols.pox_status) = sym_state_matrix(:, StateMatCols.pox_status);
        
        % tally for all the demographic groups summed for this week
        asym2symTally = asym2symTally + asym2symTallytmp; 
        sym2recoverTally = sym2recoverTally + sym2recoverTallytmp;
   
        %%%%%%%%%%%%%%%%%%%%%%%%%% 7. treatment %%%%%%%%%%%%%%%%%%%%%%%%%%
        % people symptomatic and aware have the chance to start treatment
        % treatment and isolation status are independent
        
        % find indecies in the state matrix for those off treatment
        offtrt_aware_Ind = find(state_matrix(:,StateMatCols.alive)==1 &...
                         state_matrix(:,StateMatCols.treatment)==0 &...
                         state_matrix(:,StateMatCols.pox_aware)==1);
      
        % create state matrix based on treatment status
        offtrt_state_matrix = state_matrix(offtrt_aware_Ind,:);

        % create demographic group index
        offtrt_demog_idx = find_demog_rows2(offtrt_state_matrix, StateMatCols, demog_group_def, DemogTblCols);

        % only run the transition function if there're people off trt
        if isempty(offtrt_demog_idx)
            ontrtTallytmp = 0;
        else
            % for those off treatment, they might change status to on trt
            [offtrt_state_matrix, ontrtTallytmp] = transition5('treatment', ontrt_transition_path, ...
                                    offtrt_state_matrix, StateMatCols, offtrt_demog_idx, ...
                                    future_state_param_names, 1);
            % redefine the state matrix based on the index of the new matrix
            state_matrix(offtrt_aware_Ind, StateMatCols.treatment) = offtrt_state_matrix(:, StateMatCols.treatment);      
        end        
      
        % tally for all the demographic groups summed for this week
        ontrtTally = ontrtTally + ontrtTallytmp;     
        
        %%%%%%%%%%%%%%%%%%%%%%%%%% 8. isolation %%%%%%%%%%%%%%%%%%%%%%%%%%
        % people that are aware of their status have chance to isolate 
        % isolation continues until recovered
        
        % find indecies in the state matrix for those not isolated
        offiso_aware_Ind = find(state_matrix(:,StateMatCols.alive)==1 &...
                         state_matrix(:,StateMatCols.isolation)==0 &...
                         state_matrix(:,StateMatCols.pox_aware)==1);

        % create state matrix based on isolation status
        offiso_state_matrix = state_matrix(offiso_aware_Ind,:);

        % create demographic group index
        offiso_demog_idx = find_demog_rows2(offiso_state_matrix, StateMatCols, demog_group_def, DemogTblCols);

        % only run the transition function if this group is not empty
        if isempty(offiso_demog_idx)
            onhospTallytmp = 0;
        else
            % for those off iso, they might change status to on iso
            [offiso_state_matrix, onhospTallytmp] = transition5('isolation', iso_transition_path, ...
                                    offiso_state_matrix, StateMatCols, offiso_demog_idx, ...
                                    future_state_param_names, 1);
            % redefine the state matrix based on the index of the new matrix
            state_matrix(offiso_aware_Ind, StateMatCols.isolation) = offiso_state_matrix(:, StateMatCols.isolation);      
        end        

        % tally for all the demographic groups summed for this week
        onisoTally = onisoTally + onhospTallytmp;   

    end   
 
    % people newly to aware (including both asym and sym)
    new_aware = unique(find(sm_before_aware(:,StateMatCols.pox_aware)~=state_matrix(:, StateMatCols.pox_aware)));    
    
    % people who is hiv+ and aware of their mpox status
    new_aware_hiv = new_aware(state_matrix(new_aware, StateMatCols.hiv_status)~=0);
    aware_hiv = length(new_aware_hiv);
    
    aware_hiv_b = sum(state_matrix(new_aware_hiv, StateMatCols.race)==0);
    aware_hiv_h = sum(state_matrix(new_aware_hiv, StateMatCols.race)==1);
    aware_hiv_w = sum(state_matrix(new_aware_hiv, StateMatCols.race)==2);
    aware_hiv_ls40 = sum(state_matrix(new_aware_hiv, StateMatCols.age)<40);
    aware_hiv_ge40 = sum(state_matrix(new_aware_hiv, StateMatCols.age)>=40);
    
    % number of people in each race group turn to aware
    aware_b = sum(state_matrix(new_aware, StateMatCols.race)==0);
    aware_h = sum(state_matrix(new_aware, StateMatCols.race)==1);
    aware_w = sum(state_matrix(new_aware, StateMatCols.race)==2);
    
    % number of people in each age group turn to aware
    aware_a1 = sum(state_matrix(new_aware, StateMatCols.age_bucket)==1);
    aware_a2 = sum(state_matrix(new_aware, StateMatCols.age_bucket)==2);
    aware_a3 = sum(state_matrix(new_aware, StateMatCols.age_bucket)==3);
    aware_a4 = sum(state_matrix(new_aware, StateMatCols.age_bucket)==4);
    aware_a5 = sum(state_matrix(new_aware, StateMatCols.age_bucket)==5);
    
    % people newly to vax1
    new_vax1 = unique(find(sm_before_vax1(:,StateMatCols.vaccinated)==0 & state_matrix(:, StateMatCols.vaccinated)==1));
    new_vax2 = unique(find(sm_before_vax1(:,StateMatCols.vaccinated)==1 & state_matrix(:, StateMatCols.vaccinated)==2));
    
    % number of plwh get vax 1
    vax1_hiv = sum(state_matrix(new_vax1, StateMatCols.hiv_status)~=0);
    
    % number of people in each race group get vax 1
    vax1_b = sum(state_matrix(new_vax1, StateMatCols.race)==0);
    vax1_h = sum(state_matrix(new_vax1, StateMatCols.race)==1);
    vax1_w = sum(state_matrix(new_vax1, StateMatCols.race)==2);
    vax2_b = sum(state_matrix(new_vax2, StateMatCols.race)==0);
    vax2_h = sum(state_matrix(new_vax2, StateMatCols.race)==1);
    vax2_w = sum(state_matrix(new_vax2, StateMatCols.race)==2);
    
    vax_b = vax1_b + vax2_b;
    vax_h = vax1_h + vax2_h;
    vax_w = vax1_w + vax2_w;
    
    % number of people in each age group get vax 1
    vax1_a1 = sum(state_matrix(new_vax1, StateMatCols.age_bucket)==1);
    vax1_a2 = sum(state_matrix(new_vax1, StateMatCols.age_bucket)==2);
    vax1_a3 = sum(state_matrix(new_vax1, StateMatCols.age_bucket)==3);
    vax1_a4 = sum(state_matrix(new_vax1, StateMatCols.age_bucket)==4);
    vax1_a5 = sum(state_matrix(new_vax1, StateMatCols.age_bucket)==5);  
    vax2_a1 = sum(state_matrix(new_vax2, StateMatCols.age_bucket)==1);
    vax2_a2 = sum(state_matrix(new_vax2, StateMatCols.age_bucket)==2);
    vax2_a3 = sum(state_matrix(new_vax2, StateMatCols.age_bucket)==3);
    vax2_a4 = sum(state_matrix(new_vax2, StateMatCols.age_bucket)==4);
    vax2_a5 = sum(state_matrix(new_vax2, StateMatCols.age_bucket)==5);
    
    vax_a1 = vax1_a1 + vax2_a1;
    vax_a2 = vax1_a2 + vax2_a2;
    vax_a3 = vax1_a3 + vax2_a3;
    vax_a4 = vax1_a4 + vax2_a4;
    vax_a5 = vax1_a5 + vax2_a5;
       
    %%%%%%%%%%%%%%%%%%%%%%%%%% 8. AGING %%%%%%%%%%%%%%%%%%%%%%%%%%
    % if 52 weeks, add a year-age, and clear up the weekly age
    agw_wk_is52_id = find(state_matrix(:, StateMatCols.age_wk)==52);
    age_bucket_id = find(ismember(state_matrix(agw_wk_is52_id, StateMatCols.age), [29, 39, 49, 59]));
    
    if ~isempty(agw_wk_is52_id)
        state_matrix(agw_wk_is52_id, StateMatCols.age_wk) = 0;
        state_matrix(agw_wk_is52_id, StateMatCols.age) = state_matrix(agw_wk_is52_id, StateMatCols.age)+1;
        %%%% update aware_age_bucket, vax_age_bucket as needed %%%%%%%
        state_matrix(age_bucket_id, StateMatCols.age_bucket) = state_matrix(age_bucket_id, StateMatCols.age_bucket)+1;
    end
    % Add a week to alive MSM's weekly age
    state_matrix(find(alive==1), StateMatCols.age_wk) = state_matrix(find(alive==1), StateMatCols.age_wk) + 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% 9. Death %%%%%%%%%%%%%%%%%%%%%%%%%%
    % do this at end, want deaths happen after everything
    % For each demographic group
    for demog_group_def = demog_table.'      
        
        % create demographic group index
        demog_group_idx = find_demog_rows2(state_matrix, StateMatCols, demog_group_def, DemogTblCols);
 
        % If there's no one in this demographic group, continue
        if isempty(demog_group_idx)
            continue
        end 

        [state_matrix, deathTallytmp] = transition5('alive', deathNatural_transition_path, ...
                                        state_matrix, StateMatCols, demog_group_idx, ...
                                        future_state_param_names, 1);

        % tally for all the demographic groups summed for this week
        deathTally = deathTally + deathTallytmp; 
    end

    % Parameters from the updated state matrix
    alive = state_matrix(:,StateMatCols.alive);
    pox_status = state_matrix(:, StateMatCols.pox_status);
    pox_aware = state_matrix(:, StateMatCols.pox_aware);
    vaccinated = state_matrix(:, StateMatCols.vaccinated);
    treatment = state_matrix(:, StateMatCols.treatment);
    isolation = state_matrix(:, StateMatCols.isolation);
    vax_wk = state_matrix(:, StateMatCols.vax_wk);
    vax2_wk = state_matrix(:, StateMatCols.vax2_wk);
    hiv_status = state_matrix(:, StateMatCols.hiv_status);
    hiv_aware = state_matrix(:, StateMatCols.hiv_aware);
    hiv_age_bucket = state_matrix(:, StateMatCols.hiv_age_bucket);
    race = state_matrix(:, StateMatCols.race);
    ve = state_matrix(:, StateMatCols.ve);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% HIV INFECTION %%%%%%%%%%%%%%%%%%%%%%%%%%
    % assign people from HIV- to HIV+ every week
    %           15-29   30-49   50-64   65+
    % Black      3      3       2       1
    % Hispanic   6      7       5       2
    % White      2      2       2       1
    hiv_b_a1 = randsample(find(alive==1 & hiv_status==0 & race==0 & hiv_age_bucket==1), 3);
    hiv_b_a2 = randsample(find(alive==1 & hiv_status==0 & race==0 & hiv_age_bucket==2), 3);
    hiv_b_a3 = randsample(find(alive==1 & hiv_status==0 & race==0 & hiv_age_bucket==3), 2);
    hiv_b_a4 = randsample(find(alive==1 & hiv_status==0 & race==0 & hiv_age_bucket==4), 1);
    hiv_h_a1 = randsample(find(alive==1 & hiv_status==0 & race==1 & hiv_age_bucket==1), 6);
    hiv_h_a2 = randsample(find(alive==1 & hiv_status==0 & race==1 & hiv_age_bucket==2), 7);
    hiv_h_a3 = randsample(find(alive==1 & hiv_status==0 & race==1 & hiv_age_bucket==3), 5);
    hiv_h_a4 = randsample(find(alive==1 & hiv_status==0 & race==1 & hiv_age_bucket==4), 2);
    hiv_w_a1 = randsample(find(alive==1 & hiv_status==0 & race==2 & hiv_age_bucket==1), 2);
    hiv_w_a2 = randsample(find(alive==1 & hiv_status==0 & race==2 & hiv_age_bucket==2), 2);
    hiv_w_a3 = randsample(find(alive==1 & hiv_status==0 & race==2 & hiv_age_bucket==3), 2);
    hiv_w_a4 = randsample(find(alive==1 & hiv_status==0 & race==2 & hiv_age_bucket==4), 1);
    
    hiv_id = [hiv_b_a1; hiv_b_a2; hiv_b_a3; hiv_b_a4;...
        hiv_h_a1; hiv_h_a2; hiv_h_a3; hiv_h_a4;...
        hiv_w_a1; hiv_w_a2; hiv_w_a3; hiv_w_a4];
    
    state_matrix(hiv_id, StateMatCols.hiv_status) = 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% HIV Aware %%%%%%%%%%%%%%%%%%%%%%%%%%
    % assign people from HIV- to HIV+ every week
    %           15-29   30-49   50-64   65+
    % Black      2      3       2       0
    % Hispanic   5      7       4       1
    % White      2      3       2       0
    hiv_aware_b_a1 = randsample(find(alive==1 & hiv_status~=0 & hiv_aware==0 & race==0 & hiv_age_bucket==1), 2);
    hiv_aware_b_a2 = randsample(find(alive==1 & hiv_status~=0 & hiv_aware==0 & race==0 & hiv_age_bucket==2), 3);
    hiv_aware_b_a3 = randsample(find(alive==1 & hiv_status~=0 & hiv_aware==0 & race==0 & hiv_age_bucket==3), 2);
    hiv_aware_h_a1 = randsample(find(alive==1 & hiv_status~=0 & hiv_aware==0 & race==1 & hiv_age_bucket==1), 5);
    hiv_aware_h_a2 = randsample(find(alive==1 & hiv_status~=0 & hiv_aware==0 & race==1 & hiv_age_bucket==2), 7);
    hiv_aware_h_a3 = randsample(find(alive==1 & hiv_status~=0 & hiv_aware==0 & race==1 & hiv_age_bucket==3), 4);
    hiv_aware_h_a4 = randsample(find(alive==1 & hiv_status~=0 & hiv_aware==0 & race==1 & hiv_age_bucket==4), 1);
    hiv_aware_w_a1 = randsample(find(alive==1 & hiv_status~=0 & hiv_aware==0 & race==2 & hiv_age_bucket==1), 2);
    hiv_aware_w_a2 = randsample(find(alive==1 & hiv_status~=0 & hiv_aware==0 & race==2 & hiv_age_bucket==2), 3);
    hiv_aware_w_a3 = randsample(find(alive==1 & hiv_status~=0 & hiv_aware==0 & race==2 & hiv_age_bucket==3), 2);
    
    hiv_aware_id = [hiv_aware_b_a1; hiv_aware_b_a2; hiv_aware_b_a3;...
        hiv_aware_h_a1; hiv_aware_h_a2; hiv_aware_h_a3; hiv_aware_h_a4;...
        hiv_aware_w_a1; hiv_aware_w_a2; hiv_aware_w_a3];
    
    state_matrix(hiv_aware_id, StateMatCols.hiv_aware) = 1;
    
    %%%%%%%%%%%%%%% turn off treatment if recovered %%%%%%%%%%%%%%%%%%%%%%
    % find those recovered but trt is still on
    ontrt_rec_id = find(alive==1 & pox_status==3 & treatment==1);
    state_matrix(ontrt_rec_id, StateMatCols.treatment) = 0;

    %%%%%%%%%%%%%%% turn off aware if recovered %%%%%%%%%%%%%%%%%%%%%%
    % find those recovered but aware is still on
    aware_rec_id = find(alive==1 & pox_status==3 & pox_aware==1);
    state_matrix(aware_rec_id, StateMatCols.pox_aware) = 0;
    
    %%%%%%%%%%%%%%% turn off isolation if recovered %%%%%%
    % find those recovered but isolation is still on
    iso_rec_hos_id = find(alive==1 & pox_status==3 & isolation==1);
    state_matrix(iso_rec_hos_id, StateMatCols.isolation) = 0;
    
    % Parameters from the updated state matrix  
    treatment = state_matrix(:, StateMatCols.treatment);
    pox_aware = state_matrix(:, StateMatCols.pox_aware);
    isolation = state_matrix(:, StateMatCols.isolation);
    hiv_aware = state_matrix(:, StateMatCols.hiv_aware);

    %%%%%%%%%%%%%%% add a week to infectious_wks %%%%%%
    % find those mpox status is 2
    infectious_id = find(alive==1 & pox_status==2);
    state_matrix(infectious_id, StateMatCols.infectious_wks) = state_matrix(infectious_id, StateMatCols.infectious_wks)+1;

    %%%%%%%%%%%%% set infectious_wks=0 if not pox_status==2 %%%%%
    % find those mpox status is 2
    not_infectious_id = [find(alive==1 & pox_status~=2); find(alive==0)];
    state_matrix(not_infectious_id, StateMatCols.infectious_wks) = 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%% SAVE STATE MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%
    matrix_name = strcat(state_matrices_path, int2str(t));
%     csvwrite_with_headers([matrix_nsame, '.csv'], state_matrix, fieldnames(StateMatCols))
    save(matrix_name, 'state_matrix');

    infected_per_infectious = infectionTally / sum(alive_beginning==1 & pox_status_beginning==2);
    infectious_wks = state_matrix(:, StateMatCols.infectious_wks);
    r_t = infected_per_infectious * mean(infectious_wks(pox_status==2 & alive==1)); 

    % a vector for all tally
    tally(t+1,:) = [t, birthTally, infectionTally, infected_per_infectious,...
            r_t, max(infectious_wks),...
            infect_hiv, infect_vax1, infect_vax2,...
            infect_b, infect_h, infect_w,...
            infect_a1, infect_a2, infect_a3, infect_a4, infect_a5,...
            vac1Tally, vac2Tally, vacTally, vax1_hiv,...
            awareTally1, awareTally2, awareTally1+awareTally2, aware_hiv,...
            aware_hiv_b, aware_hiv_h, aware_hiv_w,...
            aware_hiv_ls40, aware_hiv_ge40,...
            aware_b, aware_h, aware_w,...
            aware_a1, aware_a2, aware_a3, aware_a4, aware_a5,...
            vax1_b, vax1_h, vax1_w,...
            vax1_a1, vax1_a2, vax1_a3, vax1_a4, vax1_a5,...
            vax2_b, vax2_h, vax2_w,...
            vax2_a1, vax2_a2, vax2_a3, vax2_a4, vax2_a5,...
            vax_b, vax_h, vax_w,...
            vax_a1, vax_a2, vax_a3, vax_a4, vax_a5,...
            asym2symTally, sym2recoverTally, asym2recoverTally,...
            ontrtTally, deathTally, onisoTally,...
            sum(alive),...
            sum(alive & pox_status==0), sum(alive & pox_status==1),...
            sum(alive & pox_status==2), sum(alive & pox_status==3),...
            sum(alive & pox_aware==1),...
            sum(alive & pox_aware==1 & pox_status==1),...
            sum(alive & pox_aware==1 & pox_status==2),...
            sum(alive & vaccinated==1), sum(alive & vaccinated==2),...
            sum(alive & vaccinated==1 & pox_status==0),...
            sum(alive & vaccinated==1 & pox_status==1),...
            sum(alive & vaccinated==1 & pox_status==2),...
            sum(alive & vaccinated==1 & pox_status==3),...
            sum(alive & vaccinated==2 & pox_status==0),...
            sum(alive & vaccinated==2 & pox_status==1),...
            sum(alive & vaccinated==2 & pox_status==2),...
            sum(alive & vaccinated==2 & pox_status==3),...
            sum(alive & vaccinated==1 & race==0), sum(alive & vaccinated==1 & race==1), sum(alive & vaccinated==1 & race==2),...
            sum(alive & vaccinated==2 & race==0), sum(alive & vaccinated==2 & race==1), sum(alive & vaccinated==2 & race==2),...
            sum(alive & vaccinated==1 & hiv_status~=0), sum(alive & vaccinated==2 & hiv_status~=0),...
            sum(alive & treatment==1), sum(alive & isolation==1), sum(alive & hiv_aware==1), sum(alive & hiv_aware==1 & pox_aware==1)];

    %%%%%%%%%%%%%%%%%%%%%%%%%% timer %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    completion = t/T*100;
    display(strcat(num2str(completion),'%'))

end

%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE Tally %%%%%%%%%%%%%%%%%%%%%%%%%%
matrix_name = [state_matrices_path,'Tally_',char(testVersion)];
tally_tbl = array2table(tally);
tally_tbl.Properties.VariableNames(1:end) = tallyHeaders;
writetable(tally_tbl, [matrix_name, '.csv'])






