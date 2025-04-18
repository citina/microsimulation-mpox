% This function will take in a state matrix and calculate the infection
% probabilities for each demographic group using the formula 
% P(infection) = (I/N) * S

function infection_probability = calc_infec_prob4(state_matrix, StateMatCols, ...
    mix_demog_group_def, mix_demog_col_struct, common_col_struct, isolation_adherence)

    % Find the people in the demographic group
    demog_mat_indices = find_demog_rows(state_matrix, StateMatCols, mix_demog_group_def, mix_demog_col_struct, common_col_struct);
       
    % If there are no people in the demographic group, return an infection
    % probability of 0
    if size(demog_mat_indices, 2) == 0
       infection_probability = 0;
       return;
    end
    
    % find total number of people in the demographic group (N) who are alive
    alive = find_indices(state_matrix, demog_mat_indices, StateMatCols.alive, '=', 1);
    N = size(alive, 2);
   
    % number of peole who are symptomatic infected (pox_status=2) and not
    % hospitalized, not adherent isolated, or not isolated
    I_all = find_indices(state_matrix, alive, StateMatCols.pox_status, '=', 2);
    I_niso = find_indices(state_matrix, I_all, StateMatCols.isolation, '=', 0);
    I_iso = find_indices(state_matrix, I_all, StateMatCols.isolation, '=', 1);
    I = size(I_niso, 2) + size(I_iso, 2) * (1-isolation_adherence);
    
    % if nobody alive in demographic group, add 0
    if N == 0
       infection_probability = 0;
       return;
    end
    infection_probability = I/N;

end