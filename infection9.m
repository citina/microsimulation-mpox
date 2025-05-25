% 2023
% This function will take in a particular demographic group, a state
% matrix, and a mixing matrix and will transition people in this
% demographic group to infection based on a uniform random number draw

% infection requries 2 state matrix. One is the state matrix being updated with state progression. (the output)
% the 2nd state matrix is the one that represents the status at the beginning of the time window (year). 
    % this is necessary because if the state matrix reflecting progression is being used ot calculated probability of infection,
    % the chance of infection is higher for each sequential demographic group, because those who were only susceptible before have become infected.

% P(infection) can be calibrated by a grand multiplier
% (calibrationConstant), age (YoungForceOfInfectionScale), and race
% (RaceForceOfInfectionScale). It can also be altered by quarantine
% adherence or the effeciency of vaccination

function [updated_state_matrix, infectionTally] = infection9(reference_demog_group_def, DemogTblCols, ...
                                                             mixing_matrix, mixing_table, MixingTblCols, ...
                                                             state_matrix, ...
                                                             StateMatCols, numPartners, ...
                                                             calibrationConstant, ...
                                                             hivInfect,...
                                                             YoungForceOfInfectionScale, ...
                                                             MidForceOfInfectionScale,...
                                                             oldForceOfInfectionScale, ...
                                                             RaceForceOfInfectionScale,...
                                                             isolation_adherence)
    
    % make a copy of the state matrix such that a tally can be made of the number of changes that occur
    state_matrix_prior_inf = state_matrix;
    state_matrix_toUpdate = state_matrix;

    % The reference dographic group is formatted as [age min, age max, race].
    % we can thus check the 2nd number and if it is under 25, we know to divide by our risk due to condom use
    % and apply it to our calibration constant... 
    % this works because our demographic group cutoff is at 25

    if reference_demog_group_def(1) >= 15 && reference_demog_group_def(2) <= 24
        % if young, divide by the scale (becuase scale less than one)
        ForceOfInfection = calibrationConstant / YoungForceOfInfectionScale;
    elseif reference_demog_group_def(1) >= 35 && reference_demog_group_def(2) <= 44       
        ForceOfInfection = calibrationConstant / MidForceOfInfectionScale;
    elseif reference_demog_group_def(1) >= 45
        ForceOfInfection = calibrationConstant / oldForceOfInfectionScale;
    else
        % otherwise, leave calibration constant as is
        ForceOfInfection = calibrationConstant;
    end
    
    % the force of infection must then be scaled by the race characteristic.
    % the race characteristic is the third one of demographic reference (add one cuz demog group starts at 0)
    ForceOfInfection = ForceOfInfection * RaceForceOfInfectionScale(reference_demog_group_def(3) + 1);

    % Returns the row index based on the mixing_mat_def
    % demographic group table
    demog_row = find_demog_rows(mixing_table, MixingTblCols, reference_demog_group_def, DemogTblCols, MixingTblCols);
 
    % Based on the demographic row, extract the mixing probabilities from the
    % corresponding row number in the mixing matrix
    mixingProbabilites = mixing_matrix(demog_row,:);
       
    % Placeholder for accounting for all the infection probabilities over all
    % demographic group our reference group is mixing with
    total_no_infec_prob = 1;
    
    % calculate the product part in the formula
    % For each of all possible demographic groups
    for demog_group_idx = 1:size(mixingProbabilites, 2)

        % Get the demographic group definition we are mixing with
        mix_demog_group_def = mixing_table(demog_group_idx,:);
        
        % pull mixing probability for the particular instnace (who the susceptible person is mixing with)
        mixingProb_group = mixingProbabilites(demog_group_idx);
        
        % pull number of partners associated with demographic row (the demographic of suceptible person)
        numPart = numPartners(demog_row);

        % determine probability of no infection when mixing with this particular race
        % account for the number of partners in the particular racial group
        % calibration constant is applied to the I/N calculation (calf_infec_prob)
        % I_N represent I/N
        I_N = calc_infec_prob4(state_matrix_prior_inf, StateMatCols, mix_demog_group_def, ...
                         MixingTblCols, MixingTblCols, isolation_adherence);
        no_infec_prob_1par = (1 - ForceOfInfection*I_N);
        if no_infec_prob_1par < 0
            no_infec_prob_1par = 0;
        end
        no_infec_prob = no_infec_prob_1par ^ (numPart * mixingProb_group);
        
        % Multiply this to the total no infection probability over all stratifications may mix with
        % groups we are mixing with
        total_no_infec_prob = total_no_infec_prob * no_infec_prob;  
    end  

    % Find the people in the reference demographic group
    reference_demog_indices = find_demog_rows2(state_matrix_toUpdate, StateMatCols, ...
                                    reference_demog_group_def, DemogTblCols);
    
    alive = state_matrix_toUpdate(reference_demog_indices, StateMatCols.alive);
    pox_status = state_matrix_toUpdate(reference_demog_indices, StateMatCols.pox_status);
    ve = state_matrix_toUpdate(reference_demog_indices, StateMatCols.ve);
    
    % Find eligible people to infect in the reference demographic group
    sus_idx = find(alive & pox_status==0);
    
    % some vaccinated people are removed from susceptible based on their ve.
    % ve is calculated for each individual based on their vaccinated dosage
    % time after vaccination, and HIV status (see "update ve" section)
    % the higher the ve, the more likely the person is removed from sus
    unif_prob = rand([1, length(sus_idx)]);
    keep_sus = ve(sus_idx)' < unif_prob;  
    keep_sus_idx = sus_idx(keep_sus);

    % The final suspectible indecies
    susceptible = keep_sus_idx;
    
    % find those HIV aware and unaware mpox-susceptible
    sus_hiv_aware = susceptible(state_matrix_toUpdate(susceptible,StateMatCols.hiv_aware)==1);
    sus_hiv_unaware = susceptible(state_matrix_toUpdate(susceptible,StateMatCols.hiv_aware)==0); %either unaware or HIV-
    
    % For aware HIV+, apply a relative risk of acquiring Mpox
    probInfection_hivneg = 1 - total_no_infec_prob;
    probInfection_hiv = min(probInfection_hivneg * hivInfect, 1);
    
    % people are draw to infected based on their P(infect)
    unif_prob_hiv = rand([1, length(sus_hiv_aware)]);
    unif_prob_hivneg = rand([1, length(sus_hiv_unaware)]);
    to_infect_hiv = probInfection_hiv > unif_prob_hiv;
    to_infect_hivneg = probInfection_hivneg > unif_prob_hivneg; %either unaware or HIV-
                        
    % Get indices to infect based on people who are eligible
    state_mat_infec_rows_hiv = sus_hiv_aware(to_infect_hiv);
    state_mat_infec_rows_hivneg = sus_hiv_unaware(to_infect_hivneg); %either unaware or HIV-
    
    % Infect chosen people
    state_matrix_toUpdate([state_mat_infec_rows_hiv; state_mat_infec_rows_hivneg], StateMatCols.pox_status) = 1;
 
    % we assume all individuals have same average number of partners
    % prob to infect when have multiple partners is binomial with p = probability infection by 1 person
    % we want 1 - p(infection by nobody) to show probability that infection by at least 1 person
    % this could need to be swithced to a distribution such that each person has dif number of partners
 
    % Return the updated state matrix
    updated_state_matrix = state_matrix_toUpdate;

    % return infection infection infectionTally
    infectionTally = sum(state_matrix_prior_inf(:, StateMatCols.pox_status) ~= state_matrix_toUpdate(:, StateMatCols.pox_status));
  
end


