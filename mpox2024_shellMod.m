% Main simulation script for mpox transmission and vaccination dynamics
% Simulates mpox spread in a population with vaccination, awareness, treatment, and isolation interventions

% Load simulation parameters
mpox2024_parameters;

% Initialize paths and directories
if exist('local_sim_dataDir', 'var')
    state_matrices_path = convertStringsToChars(local_sim_dataDir);
    S = local_S;
    WANING_VE_MODE = local_WANING_VE_MODE;
    ENABLE_SENSITIVITY = local_ENABLE_SENSITIVITY;
else
    state_matrices_path = convertStringsToChars(sim_dataDir);
end

% Load initial population & mixing matrix data
[state_matrix, StateMatCols] = read_table(init_pop_file);
[demog_table, DemogTblCols] = create_demog_groups(demog_var_def_file);
[mixing_table, MixingTblCols] = create_demog_groups(mixing_mat_def_file);
mixing_matrix = importdata(mixing_mat_file);
numPartners = importdata(num_partners_file);

% Initialize tracking variables
finalTransitionTally = [];
workingDir_home = pwd;

% Clear existing state matrices
cd(sim_dataDir);
delete *.mat
cd(workingDir_home)

% Save initial population state
matrix_name = strcat(state_matrices_path, int2str(0));
save(matrix_name, 'state_matrix');

% Define output tracking headers
tallyHeaders = {
    'Week', 'New Births', 'New Infections', 'infected_per_infectious',...
    'r_t', 'max(infectious_wks)',...
    'new infect|hiv', 'new infect|vax1', 'new infect|vax2',...
    'new infect b', 'new infect h', 'new infect w',...
    'new infect a1', 'new infect a2', 'new infect a3', 'new infect a4', 'new infect a5',...
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
    'on trt', 'on iso', 'hiv aware', 'pox + hiv aware'
};

% Initialize tally matrix
tally = zeros(T+1, length(tallyHeaders)); % +1 cuz include the initial wk as 0

% Initialize state tracking variables
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

% Calculate initial state tallies
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

%% %%%%%%%%%%%%%%%%%%% Scenario Settings %%%%%%%%%%%%%%%%%%%%%%%%
% randomly select weeks to assign low FoI 
rng("shuffle");
if S==17
    selected_wks = randsample(1:T, 22);
elseif S==18
    selected_wks = randsample(1:T, 13);
elseif S==19
    selected_wks = randsample(1:T, 9);
end

[IMPORTATION_TYPE, iso_transition_path, foi, const_import_num] = set_scenario_parameters(S);

% weekly input schedule: the week and its corresponding number of imports
if IMPORTATION_TYPE == 2
    import_wk = [1 3 5 7 8 12 13 15 16 20 21 22 23 24 26 27 28 30 31 32 33 35 37 39 41 45 46 50 52 53 54 57 59 60 61 63 64 65 67 68 69 70];
    import_num = [1 1 2 2 1 2 2 1 1 3 1 2 1 2 2 1 3 1 2 4 2 2 1 1 1 1 1 1 1 1 2 3 1 1 2 1 1 1 1 1 1 1];
end

% Set transition path for sensitivity analysis
if ENABLE_SENSITIVITY
    asym2sym_transition_path = asym2sym_transition_path_9;
end
% for sensitivity analysis
% import_num = import_num * 2;

%% %%%%%%%%%%%%%%% Main simulation loop for each week %%%%%%%%%%%%%%%%%%%%%%%
for t = 1:T
    % Initialize weekly tallies
    [infectionTally, vac1Tally, vac2Tally, awareTally1, awareTally2,...
     asym2symTally, sym2recoverTally, asym2recoverTally, ontrtTally, deathTally,...
     onisoTally] = deal(0);

    switch S
        case 11 % S11: week 4 – 16 (Apr-Jun) and week 26 – 34 (Sept-Oct), FOI = 0.7 
            if t>=4 && t<=16
                foi = 0.7;
            elseif t>=26 && t<=34
                foi = 0.7; 
            end
        case 15 % S15: For week 4 – 16 (Apr-Jun), FOI = 0.7
            if t>=4 && t<=16
                foi = 0.7;
            end
        case 16 % S16: For week 26 – 34 (Sept-Oct), FOI = 0.7
            if t>=26 && t<=34
                foi = 0.7; 
            end
        case 17 % S17: Randomly select 22 weeks with FoI = 0.7
            selected_wks = randsample(1:T, 22);
            if ismember(t, selected_wks)
                foi = 0.7;
            end
        case 18 % S18: Randomly select 13 weeks with FoI = 0.7
            selected_wks = randsample(1:T, 13);
            if ismember(t, selected_wks)
                foi = 0.7;
            end
        case 19 % S19: Randomly select 9 weeks with FoI = 0.7
            selected_wks = randsample(1:T, 9);
            if ismember(t, selected_wks)
                foi = 0.7;
            end
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%% update ve  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Vaccine effectiveness (VE) follows a piecewise function based on:
    % 1. HIV status (immunocompetent vs immunocompromised)
    % 2. Vaccination status (1st vs 2nd dose)
    % 3. Time since vaccination
    % 4. Waning mode (0: no waning, 1: wane to 0, 2: wane to 50%, 3: wane to 75%, 4: wane to 25%)
    % Resource: https://www.nejm.org/doi/full/10.1056/NEJMoa1817307
    state_matrix = update_ve(state_matrix, StateMatCols, WANING_VE_MODE);
    
    %% %%%%%%%%%%% 1. Update vax_wk and vax2_wk if vaccinated and alive %%%%%%%%%%%%%%%
    % for people die at the end of the week vax_wk is updated here before they die
    % for new birth this week, their vax_wk is 0
    vax_wk_id = find(alive & vaccinated ~= 0);
    state_matrix(vax_wk_id, StateMatCols.vax_wk) = state_matrix(vax_wk_id, StateMatCols.vax_wk)+1;
    
    vax2_wk_id = find(alive & vaccinated == 2);
    state_matrix(vax2_wk_id, StateMatCols.vax2_wk) = state_matrix(vax2_wk_id, StateMatCols.vax2_wk)+1;
 
    %% %%%%%%%%%%%%%% 2024 update: importation %%%%%%%%%%%%%%%%%%%%%
    if IMPORTATION_TYPE == 1
        % converting X people from sus to symptomatic every week starting from Jul 30 2023
        infect_id = find(alive==1 & pox_status==0);
        infect_id = randsample(infect_id, const_import_num);
        state_matrix(infect_id, StateMatCols.pox_status) = 2;
        state_matrix(infect_id, StateMatCols.pox_aware) = 1;
    elseif IMPORTATION_TYPE == 2
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
    end
    
    % Store initial state for tracking
    alive_beginning = state_matrix(:, StateMatCols.alive);
    pox_status_beginning = state_matrix(:, StateMatCols.pox_status);
    infectious_wks_beginning = state_matrix(:, StateMatCols.infectious_wks);
    
    %% %%%%%%%%%%%%%%%%%%%%%%%% 2. New births %%%%%%%%%%%%%%%%%%%%%%%%
    % Some 14 year-old enter the model
    [state_matrix, birthTally] = birth2(state_matrix, StateMatCols, age_def, race_def, race_prop, inflow);
    
    % update columns from the state matrix
    alive = state_matrix(:, StateMatCols.alive);
    
    %% %%%%%%%%%%%%%%%%%%%%%%% 3. Mpox INFECTION %%%%%%%%%%%%%%%%%%%%%%%%%
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
                                    numPartners, foi, hivInfect,...
                                    youngInfect, midInfect,...
                                    oldInfect, raceInfect, ...
                                    isolation_adherence);

        % tally for all the demographic groups summed for this week
        infectionTally = infectionTally + infectionTallytmp;
    end
    
    % Track new infections by demographic
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
    
    % calculate r_t
    infected_per_infectious = infectionTally / sum(alive_beginning==1 & pox_status_beginning==2);
    avg_infectious_wks = mean(infectious_wks_beginning(alive_beginning==1 & pox_status_beginning==2));
    r_t = infected_per_infectious * avg_infectious_wks; 

    %% %%%%%%%%%%%%%%%%%% vaccination under policies %%%%%%%%%%%%%%%%%%%%%% 
    % 95% of people aware of their pox status and asymp will be prioritized getting vax1
    vax1_pri = find(state_matrix(:,StateMatCols.alive)==1 &...
                 state_matrix(:,StateMatCols.pox_status)==1 &...
                 state_matrix(:,StateMatCols.vaccinated)==0 &...
                 state_matrix(:,StateMatCols.pox_aware)==1);

    %% %%%%%%%%%%%%%%%%%% For each demographic group %%%%%%%%%%%%%%%%%%%%%%
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

        %% %%%%%%%%%%%%%%%%%%%%%% 5. vaccination %%%%%%%%%%%%%%%%%%%%%%%%
        % mpox vaccine requires 2 doses, vaccinated = 1 or 2        
        % run this block only if one of the above vax policies (12-15) are not selected            
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
        
        %% %%%%%%%%%%%%%%%%%%%%%% 6. pox status %%%%%%%%%%%%%%%%%%%%%%%%
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
   
        %% %%%%%%%%%%%%%%%%%%%%%%%% 7. treatment %%%%%%%%%%%%%%%%%%%%%%%%%%
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
        
        %% %%%%%%%%%%%%%%%%%%%%%%%% 8. isolation %%%%%%%%%%%%%%%%%%%%%%%%%%
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
       
    %% %%%%%%%%%%%%%%%%%%%%%%%% 8. AGING %%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    %% %%%%%%%%%%%%%%%%%%%%%%%% 9. Death %%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    %% %%%%%%%%%%%%%%%%%%%%%%%% HIV Dynamics %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Count for HIV infection
    % Format: [Black, Hispanic, White] for each age group
    hiv_infection_rates = [
        3, 6, 2;  % 15-29
        3, 7, 2;  % 30-49
        2, 5, 2;  % 50-64
        1, 2, 1   % 65+
    ];
    
    % Count for HIV awareness
    % Format: [Black, Hispanic, White] for each age group
    hiv_awareness_rates = [
        2, 5, 2;  % 15-29
        3, 7, 3;  % 30-49
        2, 4, 2;  % 50-64
        0, 1, 0   % 65+
    ];
    
    % Initialize arrays to store IDs
    hiv_infection_ids = [];
    hiv_awareness_ids = [];
    
    % Process HIV infections
    for age_group = 1:4
        for race = 0:2
            % Get the number of new infections for this demographic
            num_infections = hiv_infection_rates(age_group, race + 1);
            if num_infections > 0
                % Find eligible individuals
                eligible_ids = find(alive==1 & hiv_status==0 & race==race & hiv_age_bucket==age_group);
                if ~isempty(eligible_ids)
                    % Randomly select individuals to infect
                    selected_ids = randsample(eligible_ids, min(num_infections, length(eligible_ids)));
                    hiv_infection_ids = [hiv_infection_ids; selected_ids];
                end
            end
        end
    end
    
    % Process HIV awareness
    for age_group = 1:4
        for race = 0:2
            % Get the number of new awareness cases for this demographic
            num_aware = hiv_awareness_rates(age_group, race + 1);
            if num_aware > 0
                % Find eligible individuals
                eligible_ids = find(alive==1 & hiv_status~=0 & hiv_aware==0 & race==race & hiv_age_bucket==age_group);
                if ~isempty(eligible_ids)
                    % Randomly select individuals to make aware
                    selected_ids = randsample(eligible_ids, min(num_aware, length(eligible_ids)));
                    hiv_awareness_ids = [hiv_awareness_ids; selected_ids];
                end
            end
        end
    end
    
    % Apply HIV infections and awareness
    state_matrix(hiv_infection_ids, StateMatCols.hiv_status) = 1;
    state_matrix(hiv_awareness_ids, StateMatCols.hiv_aware) = 1;
    
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

    % Update weekly tallies
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

%% %%%%%%%%%%%%%%%%%%%%%%%% SAVE Tally %%%%%%%%%%%%%%%%%%%%%%%%%%
matrix_name = [state_matrices_path,'Tally_',char(testVersion)];
tally_tbl = array2table(tally);
tally_tbl.Properties.VariableNames(1:end) = tallyHeaders;
writetable(tally_tbl, [matrix_name, '.csv'])

% Create memo file with scenario information
create_scenario_memo(S, testVerDir, testVersion, NUM_ITERATIONS, iso_transition_path, foi);

%% %%%%%%%%%%%%%%%%%%%%%%%% Helper Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function state_matrix = update_ve(state_matrix, StateMatCols, WANING_VE_MODE)
    % UPDATE_VE Updates vaccine effectiveness based on various factors
    %   Updates vaccine effectiveness (VE) based on:
    %   1. HIV status (immunocompetent vs immunocompromised)
    %   2. Vaccination status (1st vs 2nd dose)
    %   3. Time since vaccination
    %   4. Waning mode (0: no waning, 1: wane to 0, 2: wane to 50%, 3: wane to 75%, 4: wane to 25%)
    %
    %   Inputs:
    %   - state_matrix: Current state matrix
    %   - StateMatCols: Structure containing column indices
    %   - WANING_VE_MODE: Mode for VE waning (0, 1, 2, 3, or 4)
    %
    %   Output:
    %   - state_matrix: Updated state matrix with new VE values

    % Extract relevant columns
    alive = state_matrix(:, StateMatCols.alive);
    hiv_status = state_matrix(:, StateMatCols.hiv_status);
    vaccinated = state_matrix(:, StateMatCols.vaccinated);
    vax_wk = state_matrix(:, StateMatCols.vax_wk);

    % Constants for VE calculations
    vac1_plwh = 0.51;  % One dose VE for PLWH
    vac2_plwh = 0.702; % Two doses VE for PLWH
    vac1_normal = 0.721; % One dose VE for immunocompetent
    vac2_normal = 0.878; % Two doses VE for immunocompetent

    if WANING_VE_MODE == 0  % No waning VE
        % Set constant VE for PLWH and non-PLWH based on vaccination status      
        plwh_vax1_id = hiv_status~=0 & vaccinated==1 & vax_wk>=2 & alive==1; 
        plwh_vax2_id = hiv_status~=0 & vaccinated==2 & vax_wk>=6  & alive==1;
        state_matrix(plwh_vax1_id, StateMatCols.ve) = vac1_plwh;
        state_matrix(plwh_vax2_id, StateMatCols.ve) = vac2_plwh;
        
        normal_vax1_id = hiv_status==0 & vaccinated==1 & vax_wk>=2 & alive==1;
        normal_vax2_id = hiv_status==0 & vaccinated==2 & vax_wk>=6 & alive==1;
        state_matrix(normal_vax1_id, StateMatCols.ve) = vac1_normal;        
        state_matrix(normal_vax2_id, StateMatCols.ve) = vac2_normal;
    else
        % Define target VE ratios based on waning_ve
        if WANING_VE_MODE == 1
            target_ratio = 0;      % wanes to 0%
        elseif WANING_VE_MODE == 2
            target_ratio = 0.5;    % wanes to 50%
        elseif WANING_VE_MODE == 3
            target_ratio = 0.75;   % wanes to 75%
        elseif WANING_VE_MODE == 4
            target_ratio = 0.25;   % wanes to 25%
        end
        
        % Calculate target VE values
        target_ve_normal = vac1_normal * target_ratio;
        target_ve_plwh = vac1_plwh * target_ratio;
        target_ve2_normal = vac2_normal * target_ratio;
        target_ve2_plwh = vac2_plwh * target_ratio;
        
        % Define indices for different time periods
        % Single dose development (weeks 0-4)
        idx1 = hiv_status==0 & vaccinated==1 & vax_wk>=0 & vax_wk<=2 & alive==1;  % week 0-2
        idx2 = hiv_status==0 & vaccinated==1 & vax_wk>2 & vax_wk<=4 & alive==1;   % week 2-4
        idx1_plwh = hiv_status~=0 & vaccinated==1 & vax_wk>=0 & vax_wk<=2 & alive==1;
        idx2_plwh = hiv_status~=0 & vaccinated==1 & vax_wk>2 & vax_wk<=4 & alive==1;
        
        % Single dose waning (weeks 4-56)
        idx3 = hiv_status==0 & vaccinated==1 & vax_wk>4 & vax_wk<=56 & alive==1;
        idx3_plwh = hiv_status~=0 & vaccinated==1 & vax_wk>4 & vax_wk<=56 & alive==1;
        
        % Two dose development (weeks 4-6)
        idx4 = hiv_status==0 & vaccinated==2 & vax_wk>=4 & vax_wk<=6 & alive==1;
        idx4_plwh = hiv_status~=0 & vaccinated==2 & vax_wk>=4 & vax_wk<=6 & alive==1;
        
        % Two dose waning (weeks 6-58)
        idx5 = hiv_status==0 & vaccinated==2 & vax_wk>6 & vax_wk<=58 & alive==1;
        idx6 = hiv_status==0 & vaccinated==2 & vax_wk>58 & alive==1;
        idx5_plwh = hiv_status~=0 & vaccinated==2 & vax_wk>6 & vax_wk<=58 & alive==1;
        idx6_plwh = hiv_status~=0 & vaccinated==2 & vax_wk>58 & alive==1;
        
        % Calculate VE development rates
        normal_ve_rate = vac1_normal / 2;  % rate per week to reach vac1_normal in 2 weeks
        plwh_ve_rate = vac1_plwh / 2;
        normal_ve2_rate = (vac2_normal - vac1_normal) / 2;  % rate per week to reach vac2_normal in 2 weeks
        plwh_ve2_rate = (vac2_plwh - vac1_plwh) / 2;
        
        % Calculate VE waning rates
        normal_ve_drop_rate = (vac1_normal - target_ve_normal) / 52;  % rate per week to reach target VE
        plwh_ve_drop_rate = (vac1_plwh - target_ve_plwh) / 52;
        normal_ve2_drop_rate = (vac2_normal - target_ve2_normal) / 52;
        plwh_ve2_drop_rate = (vac2_plwh - target_ve2_plwh) / 52;
        
        % Apply VE development for single dose
        state_matrix(idx1, StateMatCols.ve) = normal_ve_rate * vax_wk(idx1);
        state_matrix(idx2, StateMatCols.ve) = vac1_normal;
        state_matrix(idx1_plwh, StateMatCols.ve) = plwh_ve_rate * vax_wk(idx1_plwh);
        state_matrix(idx2_plwh, StateMatCols.ve) = vac1_plwh;
        
        % Apply VE waning for single dose
        state_matrix(idx3, StateMatCols.ve) = max(round(-normal_ve_drop_rate * (vax_wk(idx3) - 4) + vac1_normal, 3), 0);
        state_matrix(idx3_plwh, StateMatCols.ve) = max(round(-plwh_ve_drop_rate * (vax_wk(idx3_plwh) - 4) + vac1_plwh, 3), 0);
        
        % Apply VE development for two doses
        state_matrix(idx4, StateMatCols.ve) = max(round(normal_ve2_rate * (vax_wk(idx4) - 4) + vac1_normal, 3), 0);
        state_matrix(idx4_plwh, StateMatCols.ve) = max(round(plwh_ve2_rate * (vax_wk(idx4_plwh) - 4) + vac1_plwh, 3), 0);
        
        % Apply VE waning for two doses
        state_matrix(idx5, StateMatCols.ve) = max(round(-normal_ve2_drop_rate * (vax_wk(idx5) - 6) + vac2_normal, 3), 0);
        state_matrix(idx6, StateMatCols.ve) = target_ve2_normal;
        state_matrix(idx5_plwh, StateMatCols.ve) = max(round(-plwh_ve2_drop_rate * (vax_wk(idx5_plwh) - 6) + vac2_plwh, 3), 0);
        state_matrix(idx6_plwh, StateMatCols.ve) = target_ve2_plwh;
    end
end

function [IMPORTATION_TYPE, iso_transition_path, foi, const_import_num] = set_scenario_parameters(S)
    % SET_SCENARIO_PARAMETERS Sets parameters based on scenario number
    % Inputs:
    %   S - Scenario number
    %   t - Current time step
    % Outputs:
    %   IMPORTATION_TYPE - Type of importation (0: none, 1: constant, 2: scheduled)
    %   iso_transition_path - Path to isolation transition file
    %   foi - Force of infection value
    %   const_import_num - Number of constant imports (only used for IMPORTATION_TYPE=1)
    
    % Initialize default values
    iso_transition_path = '../input/isolationOn_0.2.csv';
    foi = 0.7;
    const_import_num = 0;
    
    % Scenario groups
    no_const_import_scenarios = [0, 8, 9, 10];
    phylo_import_scenarios = [1:7, 12:14, 20:22];
    
    if ismember(S, no_const_import_scenarios)
        % No importation or constant importation scenarios
        if S == 0
            % Baseline scenario
            IMPORTATION_TYPE = 0;
        else
            % Constant importation scenarios
            IMPORTATION_TYPE = 1;
            switch S
                case 8
                    const_import_num = 5;
                case 9
                    const_import_num = 10;
                case 10
                    const_import_num = 15;
            end
        end
        
    elseif ismember(S, phylo_import_scenarios)
        % Importation scenarios with varying isolation and FOI
        IMPORTATION_TYPE = 2;
        
        % Set isolation and FOI based on scenario
        switch S
            case 1
                iso_transition_path = '../input/isolationOn_0.5.csv';
                foi = 0.7;
            case 2
                iso_transition_path = '../input/isolationOn_0.5.csv';
                foi = 1.45;
            case 3
                iso_transition_path = '../input/isolationOn_0.5.csv';
                foi = 2.2;
            case 4
                iso_transition_path = '../input/isolationOn_0.4.csv';
                foi = 2.2;
            case 5
                iso_transition_path = '../input/isolationOn_0.3.csv';
                foi = 2.2;
            case 6
                iso_transition_path = '../input/isolationOn_0.2.csv';
                foi = 2.2;
            case 7
                iso_transition_path = '../input/isolationOn_0.csv';
                foi = 2.2;
            case 12
                iso_transition_path = '../input/isolationOn_0.2.csv';
                foi = 1.45;
            case 13
                iso_transition_path = '../input/isolationOn_0.2.csv';
                foi = 0.7;
            case 14
                iso_transition_path = '../input/isolationOn_0.8.csv';
                foi = 2.2;
            case 20
                iso_transition_path = '../input/isolationOn_0.2.csv';
                foi = 1.6;
            case 21
                iso_transition_path = '../input/isolationOn_0.2.csv';
                foi = 1.8;
            case 22
                iso_transition_path = '../input/isolationOn_0.2.csv';
                foi = 2.0;
        end
    else % S11, S15-S19
        IMPORTATION_TYPE = 2;
        foi = 2.2; % Default FOI
        iso_transition_path = '../input/isolationOn_0.2.csv';
    end
end

function create_scenario_memo(S, testVerDir, testVersion, NUM_ITERATIONS, iso_transition_path, foi)
    % Create a memo file with scenario information
    % Inputs:
    %   S - Scenario number
    %   testVerDir - Simulation output directory
    %   iso_transition_path - Path to isolation transition file
    %   foi - Force of infection value
    
    % Create memo file path
    memo_path = fullfile(testVerDir, 'memo.txt');
    
    % Open file for writing
    fid = fopen(memo_path, 'w');
    
    % Write scenario information
    fprintf(fid, 'Scenario Information\n');
    fprintf(fid, '===================\n\n');
    fprintf(fid, 'Test Version: %s\n', testVersion);
    fprintf(fid, 'Number of Iterations: %d\n', NUM_ITERATIONS);
    fprintf(fid, 'Scenario Number: S%d\n', S);
    fprintf(fid, 'Isolation Transition Path: %s\n', iso_transition_path);
    fprintf(fid, 'Force of Infection: %.2f\n\n', foi);
    
    % Add detailed scenario description based on scenario number
    fprintf(fid, 'Scenario Description:\n');
    
    if S == 0
        fprintf(fid, 'Baseline scenario with no importation\n');
        fprintf(fid, '- Standard isolation (0.2)\n');
        fprintf(fid, '- Constant FOI = 0.7\n');
    elseif S >= 8 && S <= 10
        fprintf(fid, 'Constant importation scenario\n');
        fprintf(fid, '- Standard isolation (0.2)\n');
        fprintf(fid, '- Constant FOI = 0.7\n');
        if S == 8
            fprintf(fid, '- 5 new cases imported per week\n');
        elseif S == 9
            fprintf(fid, '- 10 new cases imported per week\n');
        elseif S == 10
            fprintf(fid, '- 15 new cases imported per week\n');
        end
    elseif (S >= 3 && S <= 7) || S == 14
        fprintf(fid, 'Phylogenetic importation scenario with varying isolation levels\n');
        fprintf(fid, '- Fixed FOI = 2.2\n');
        if S == 3
            fprintf(fid, '- Isolation = 0.5\n');
        elseif S == 4
            fprintf(fid, '- Isolation = 0.4\n');
        elseif S == 5
            fprintf(fid, '- Isolation = 0.3\n');
        elseif S == 6
            fprintf(fid, '- Isolation = 0.2\n');
        elseif S == 7
            fprintf(fid, '- Isolation = 0.0\n');
        elseif S == 14
            fprintf(fid, '- Isolation = 0.8\n');
        end
    elseif (S >= 12 && S <= 13) || S == 20 || S == 21 || S == 22
        fprintf(fid, 'Phylogenetic importation scenario with varying FOI\n');
        fprintf(fid, '- Standard isolation (0.2)\n');
        if S == 12
            fprintf(fid, '- FOI = 1.45\n');
        elseif S == 13
            fprintf(fid, '- FOI = 0.7\n');
        elseif S == 20
            fprintf(fid, '- FOI = 1.6\n');
        elseif S == 21
            fprintf(fid, '- FOI = 1.8\n');
        elseif S == 22
            fprintf(fid, '- FOI = 2.0\n');
        end
    elseif S == 11
        fprintf(fid, 'Time-varying FOI scenario\n');
        fprintf(fid, '- Standard isolation (0.2)\n');
        fprintf(fid, '- FOI = 0.7 during weeks 4-16 and 26-34\n');
        fprintf(fid, '- FOI = 2.2 during other weeks\n');
    elseif S >= 15 && S <= 19
        fprintf(fid, 'Time-varying FOI scenario\n');
        fprintf(fid, '- Standard isolation (0.2)\n');
        if S == 15
            fprintf(fid, '- FOI = 0.7 during weeks 4-16\n');
            fprintf(fid, '- FOI = 2.2 during other weeks\n');
        elseif S == 16
            fprintf(fid, '- FOI = 0.7 during weeks 26-34\n');
            fprintf(fid, '- FOI = 2.2 during other weeks\n');
        elseif S == 17
            fprintf(fid, '- FOI = 0.7 for 22 randomly selected weeks\n');
            fprintf(fid, '- FOI = 2.2 during other weeks\n');
        elseif S == 18
            fprintf(fid, '- FOI = 0.7 for 13 randomly selected weeks\n');
            fprintf(fid, '- FOI = 2.2 during other weeks\n');
        elseif S == 19
            fprintf(fid, '- FOI = 0.7 for 9 randomly selected weeks\n');
            fprintf(fid, '- FOI = 2.2 during other weeks\n');
        end
    else
        fprintf(fid, 'Custom scenario\n');
    end
    
    % Close the file
    fclose(fid);
end


